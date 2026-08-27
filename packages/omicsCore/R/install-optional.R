# Optional-Suggests groups. Each group lists CRAN- and Bioconductor-only
# packages that the corresponding analysis backend needs at runtime. Heavy
# Bioconductor packages live here (not in Imports) so that a fresh
# `install.packages("omicsCore")` stays slim on restricted environments.
#
# Note: `limma`, `clusterProfiler`, `msigdbr`, and `qs2` are now in Imports
# (auto-installed with omicsCore) — they are required by the default
# differential-expression, enrichment, and persistence paths.
OPTIONAL_GROUPS <- list(
  rnaseq      = c("DESeq2", "edgeR", "tximport", "GenomicFeatures"),
  proteomics  = c("missForest", "pcaMethods", "vsn"),
  enrichment  = c("fgsea", "GSVA"),
  imputation  = c("impute", "missForest", "pcaMethods"),
  viz         = c("ComplexHeatmap", "circlize", "ggrepel", "patchwork", "ggpubr"),
  persistence = character(0)
)

# Packages we expect to find on Bioconductor (rather than CRAN). Used to
# pick the right installer when `pak` is not available.
BIOC_PACKAGES <- c(
  "DESeq2", "edgeR", "tximport", "GenomicFeatures",
  "fgsea", "GSVA", "ComplexHeatmap", "pcaMethods", "vsn",
  "impute"
)

#' Resolve a Suggests group into its package list
#'
#' @param group One of `"rnaseq"`, `"proteomics"`, `"enrichment"`,
#'   `"imputation"`, `"viz"`, `"persistence"`, `"all"`.
#'
#' @return Character vector of package names.
#' @keywords internal
resolve_install_group <- function(group) {
  group <- match.arg(
    group,
    choices = c(names(OPTIONAL_GROUPS), "all"),
    several.ok = FALSE
  )
  if (identical(group, "all")) {
    return(unique(unlist(OPTIONAL_GROUPS, use.names = FALSE)))
  }
  OPTIONAL_GROUPS[[group]]
}

#' Install optional dependency groups
#'
#' Installs the on-demand `Suggests` packages required by a given backend
#' group. Uses `pak::pkg_install()` if available, otherwise falls back to
#' `BiocManager::install()` (for Bioconductor packages) and
#' `utils::install.packages()` (for CRAN packages).
#'
#' This function exists so that a fresh `install.packages("omicsCore")` can
#' stay small on restricted environments where Docker / system installs are
#' not available. Users opt in to heavy backends only when they need them.
#'
#' @param group One of `"rnaseq"`, `"proteomics"`, `"enrichment"`,
#'   `"imputation"`, `"viz"`, `"persistence"`, or `"all"`.
#' @param ask If `TRUE` (default in interactive sessions), prompt before
#'   installing.
#' @param upgrade If `TRUE`, allow upgrading already-installed packages.
#'   Defaults to `FALSE`.
#'
#' @return Invisibly returns a character vector of the packages that were
#'   targeted for installation.
#' @export
#' @family install
#' @examples
#' \dontrun{
#'   install_optional("rnaseq")
#'   install_optional("all", ask = FALSE)
#' }
install_optional <- function(
  group = c("rnaseq", "proteomics", "enrichment", "imputation", "viz",
            "persistence", "all"),
  ask = interactive(),
  upgrade = FALSE
) {
  group <- match.arg(group)
  pkgs <- resolve_install_group(group)

  missing_pkgs <- pkgs[!vapply(pkgs, is_installed, logical(1))]
  if (length(missing_pkgs) == 0L) {
    message("All packages in group '", group, "' are already installed.")
    return(invisible(pkgs))
  }

  if (isTRUE(ask)) {
    msg <- paste0(
      "About to install ", length(missing_pkgs),
      " package(s) for group '", group, "':\n  ",
      paste(missing_pkgs, collapse = ", "),
      "\nProceed?"
    )
    ans <- utils::menu(c("Yes", "No"), title = msg)
    if (!identical(ans, 1L)) {
      message("Skipped installation.")
      return(invisible(pkgs))
    }
  }

  if (is_installed("pak")) {
    install_via_pak(missing_pkgs, upgrade = upgrade)
  } else {
    install_via_fallback(missing_pkgs, upgrade = upgrade)
  }

  invisible(pkgs)
}

#' Report which optional dependency groups are installed
#'
#' Returns a `data.frame` with one row per requested group listing total
#' package count, how many are installed, and the names of any missing
#' packages. Useful as a pre-flight check in the Shiny app.
#'
#' @param features Character vector of group names to check.
#'
#' @return A `data.frame` with columns `group`, `n_total`, `n_installed`,
#'   `is_ready`, `missing`.
#' @export
#' @family install
check_install <- function(
  features = c("rnaseq", "proteomics", "enrichment", "imputation", "viz",
               "persistence")
) {
  rows <- lapply(features, function(g) {
    pkgs <- resolve_install_group(g)
    installed <- vapply(pkgs, is_installed, logical(1))
    data.frame(
      group = g,
      n_total = length(pkgs),
      n_installed = sum(installed),
      is_ready = all(installed),
      missing = paste(pkgs[!installed], collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---- internal helpers --------------------------------------------------

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

install_via_pak <- function(pkgs, upgrade = FALSE) {
  bioc_targets <- intersect(pkgs, BIOC_PACKAGES)
  cran_targets <- setdiff(pkgs, BIOC_PACKAGES)

  targets <- c(
    if (length(bioc_targets) > 0L) paste0("bioc::", bioc_targets) else NULL,
    cran_targets
  )
  pak <- asNamespace("pak")
  pak$pkg_install(targets, upgrade = upgrade, ask = FALSE)
}

install_via_fallback <- function(pkgs, upgrade = FALSE) {
  bioc_targets <- intersect(pkgs, BIOC_PACKAGES)
  cran_targets <- setdiff(pkgs, BIOC_PACKAGES)

  if (length(bioc_targets) > 0L) {
    if (!is_installed("BiocManager")) {
      utils::install.packages("BiocManager")
    }
    bm <- asNamespace("BiocManager")
    bm$install(bioc_targets, update = upgrade, ask = FALSE)
  }
  if (length(cran_targets) > 0L) {
    utils::install.packages(cran_targets)
  }
  invisible(TRUE)
}
