#!/usr/bin/env Rscript
#
# Check that the Dockerfile's pinned repositories can actually satisfy
# every version requirement in the dependency closure, before spending
# 45 minutes finding out that they cannot.
#
#   Rscript deploy/scripts/check_pins.R
#   Rscript deploy/scripts/check_pins.R 2025-01-15 2025-03-01   # compare
#
# This exists because "every package resolves" and "every package's
# requirements are satisfiable" are different questions, and only the
# second one predicts whether the build works. A snapshot two days after
# Bioconductor 3.20 shipped had all 228 packages present and still failed
# at BiocParallel, which wants BH >= 1.87.0 -- released later. The
# Bioconductor mirror serves the *final* state of the 3.20 branch, so the
# CRAN snapshot has to come from the end of that release's life, not its
# beginning.
#
# Run it whenever CRAN_SNAPSHOT, BIOC_MIRROR, the Bioconductor release,
# or either package list changes.

BIOC_RELEASE <- "3.20"
BIOC_MIRROR  <- Sys.getenv("BIOC_MIRROR",
                           "https://mirrors.westlake.edu.cn/bioconductor")

# Must match the two lists in deploy/docker/Dockerfile.
CRAN_PKGS <- c(
  "shiny", "bslib", "bsicons", "sass", "htmltools", "ggplot2", "plotly", "DT",
  "shinyWidgets", "shinyjs", "promises", "future",
  "tibble", "dplyr", "stringr", "readxl", "openxlsx", "jsonlite", "scales",
  "rlang", "msigdbr", "qs2",
  "ggrepel", "patchwork", "circlize", "here", "knitr", "rmarkdown",
  "missForest", "ActivePathways", "pkgload"
)
BIOC_PKGS <- c(
  "limma", "clusterProfiler", "DESeq2", "edgeR", "S4Vectors", "enrichplot",
  "fgsea", "GSVA", "ComplexHeatmap", "impute", "pcaMethods", "vsn"
)

# Versions the image depends on for reasons that are not obvious from the
# version number, so a snapshot that silently moves them is worth naming.
WATCH <- c(
  # ggtree calls getFromNamespace("check_linewidth", "ggplot2"); 4.0
  # removed it.
  ggplot2 = "< 4",
  # < 10 carries its data in the package. 10+ moves it to msigdbdf, from
  # a repository this image does not have, and the prewarm layer then
  # needs the network.
  msigdbr = "< 10"
)

base_pkgs <- rownames(installed.packages(priority = "base"))

repos_for <- function(snapshot) {
  bioc <- paste0(BIOC_MIRROR, "/packages/", BIOC_RELEASE)
  c(BioCsoft = paste0(bioc, "/bioc"),
    BioCann  = paste0(bioc, "/data/annotation"),
    BioCexp  = paste0(bioc, "/data/experiment"),
    CRAN     = paste0("https://p3m.dev/cran/__linux__/noble/", snapshot))
}

# Parse one "pkg (>= 1.2.3)" entry. Returns NULL for an unversioned
# dependency, which is always satisfiable.
parse_req <- function(item) {
  m <- regmatches(item, regexec(
    "^([A-Za-z0-9._]+)\\s*\\(\\s*([><=]+)\\s*([0-9.-]+)\\s*\\)$", item))[[1]]
  if (!length(m)) return(NULL)
  list(pkg = m[2], op = m[3], version = m[4])
}

satisfied <- function(have, op, want) {
  have <- package_version(have); want <- package_version(want)
  switch(op, ">=" = have >= want, ">" = have > want,
         "==" = have == want, "<=" = have <= want, "<" = have < want, TRUE)
}

check_snapshot <- function(snapshot) {
  ap <- tryCatch(available.packages(repos = repos_for(snapshot)),
                 error = function(e) NULL)
  if (is.null(ap)) {
    return(list(snapshot = snapshot, reachable = FALSE))
  }

  wanted <- c(CRAN_PKGS, BIOC_PKGS)
  absent <- setdiff(wanted, rownames(ap))
  deps <- tools::package_dependencies(
    wanted, db = ap, recursive = TRUE,
    which = c("Depends", "Imports", "LinkingTo"))
  closure <- setdiff(unique(c(wanted, unlist(deps))), base_pkgs)
  closure <- intersect(closure, rownames(ap))

  conflicts <- character(0)
  for (p in closure) {
    for (field in c("Depends", "Imports", "LinkingTo")) {
      spec <- ap[p, field]
      if (is.na(spec)) next
      for (item in trimws(strsplit(spec, ",")[[1]])) {
        req <- parse_req(item)
        if (is.null(req) || req$pkg == "R" || req$pkg %in% base_pkgs) next
        if (!req$pkg %in% rownames(ap)) {
          conflicts <- c(conflicts,
                         sprintf("%s needs %s, which is not in any repository",
                                 p, req$pkg))
          next
        }
        have <- ap[req$pkg, "Version"]
        if (!isTRUE(satisfied(have, req$op, req$version))) {
          conflicts <- c(conflicts,
                         sprintf("%s needs %s %s %s, but the snapshot has %s",
                                 p, req$pkg, req$op, req$version, have))
        }
      }
    }
  }

  watched <- vapply(names(WATCH), function(p) {
    if (p %in% rownames(ap)) ap[p, "Version"] else NA_character_
  }, character(1))

  list(snapshot = snapshot, reachable = TRUE, n = length(closure),
       absent = absent, conflicts = unique(conflicts), watched = watched)
}

report <- function(r) {
  if (!isTRUE(r$reachable)) {
    cat(sprintf("%s  repositories unreachable\n", r$snapshot))
    return(FALSE)
  }
  cat(sprintf("\n%s  closure %d packages\n", r$snapshot, r$n))
  for (p in names(r$watched)) {
    cat(sprintf("  %-10s %-10s  (image expects %s)\n",
                p, r$watched[[p]], WATCH[[p]]))
  }
  ok <- TRUE
  if (length(r$absent)) {
    ok <- FALSE
    cat("  MISSING from every repository:\n")
    cat(paste0("    ", r$absent, collapse = "\n"), "\n")
  }
  if (length(r$conflicts)) {
    ok <- FALSE
    cat(sprintf("  %d unsatisfiable version requirement(s):\n",
                length(r$conflicts)))
    cat(paste0("    ", r$conflicts, collapse = "\n"), "\n")
  } else {
    cat("  every version requirement is satisfiable\n")
  }
  ok
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  dockerfile <- file.path(dirname(dirname(
    normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(),
                                           value = TRUE)[1])))),
    "docker", "Dockerfile")
  args <- if (file.exists(dockerfile)) {
    sub(".*CRAN_SNAPSHOT=", "",
        grep("^ARG CRAN_SNAPSHOT=", readLines(dockerfile), value = TRUE)[1])
  } else {
    stop("pass a snapshot date, e.g. 2025-03-01")
  }
  cat("Snapshot taken from the Dockerfile:", args, "\n")
}

options(timeout = 300)
results <- lapply(args, function(d) report(check_snapshot(d)))
cat("\n")
if (!all(unlist(results))) {
  quit(status = 1)
}
