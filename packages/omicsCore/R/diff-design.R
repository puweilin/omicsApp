# Model design construction.
#
# Every column name and every factor level in `meta_df` comes from a file the
# user uploaded. Pasting either into a formula string hands the R parser
# something the uploader controls, and `as.formula()` + `model.matrix()` then
# evaluate it: a metadata column named `file.create("...")` runs before the
# model is ever fitted.
#
# So nothing user-supplied reaches the parser. Terms become `t001, t002, ...`
# and factor levels become `L001, L002, ...`, both generated here; the design
# matrix's column names are built only from those. The user's names survive in
# `term_labels()` / `level_labels()` for display, and the mapping is what turns
# a requested contrast back into a column position.
#
# Contrasts are numeric vectors indexed by column position. `makeContrasts()`
# parses its `contrasts=` string, which is the same exposure by another route,
# and matching a level by name breaks as soon as a legitimate group is called
# what an internal identifier is called.
#
# `model.matrix(~ ., data = df)` is safe and stays: `.` expands through symbols
# rather than through the parser, so a hostile column name in `df` is inert.
# It is the paste-a-string-then-parse pattern that has to go.

term_id <- function(i) sprintf("t%03d", i)
level_id <- function(i) sprintf("L%03d", i)

#' Rename metadata columns and factor levels to generated identifiers
#'
#' @param meta A data frame of sample metadata.
#' @param terms Character vector of column names to take from `meta`.
#'
#' @return List with `df` (renamed data frame, safe for `~ .`), `ids`,
#'   `terms` (id -> user column name) and `levels` (id -> user level labels).
#' @keywords internal
#' @noRd
safe_frame <- function(meta, terms) {
  terms <- as.character(terms %||% character(0))
  missing <- setdiff(terms, colnames(meta))
  if (length(missing) > 0L) {
    stop("Columns not found in metadata: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  ids <- if (length(terms)) term_id(seq_along(terms)) else character(0)
  df <- meta[, terms, drop = FALSE]
  names(df) <- ids

  level_labels <- list()
  for (i in seq_along(ids)) {
    col <- df[[i]]
    if (is.factor(col) || is.character(col)) {
      # An existing factor keeps its level order -- run_diff()'s backends
      # relevel so that control comes first, and re-factoring here would
      # silently drop that.
      f <- if (is.factor(col)) col else factor(col)
      level_labels[[ids[[i]]]] <- levels(f)
      levels(f) <- level_id(seq_along(levels(f)))
      df[[i]] <- f
    }
  }

  list(df = df, ids = ids,
       terms = stats::setNames(as.list(terms), ids),
       levels = level_labels)
}

#' Build a design matrix without exposing user strings to the parser
#'
#' @param meta Sample metadata.
#' @param terms Model terms, in order; the first is normally the primary one.
#' A formula that does not carry the frame it was built in
#'
#' `stats::as.formula()` on a string takes `env = parent.frame()`, and a
#' formula holds its environment for as long as anything holds the formula.
#' That is fine for a formula consumed and dropped, and a leak for one that
#' is stored: DESeq2 keeps the design inside the `DESeqDataSet`, so the
#' object returned by `run_deseq2_group()` referred back to that function's
#' frame -- which by the time it returned held the subsetted count matrix,
#' the standardized results, and the `DESeqDataSet` itself. A cycle.
#'
#' Measured at 63k features x 12 samples: 39 MB of actual data, 131 MB once
#' the environment is followed, of which the formula accounted for 98 MB.
#' `object.size()` does not follow environments and so reported the first
#' number; `future` does, and refused to send the bundle to a worker.
#'
#' The design terms are columns of the data frame passed alongside, so
#' nothing is ever looked up in the formula's environment. `baseenv()`
#' rather than `emptyenv()`, because code that calls `environment(f)` on a
#' formula expects something it can search.
#'
#' @param text The formula, as a string.
#' @keywords internal
#' @noRd
detached_formula <- function(text) {
  stats::as.formula(text, env = baseenv())
}

#' @param primary Which of `terms` carries the contrast, if any.
#' @param intercept Include an intercept column.
#'
#' @return A model matrix carrying `term_labels`, `level_labels` and
#'   `primary_id` attributes.
#' @keywords internal
#' @noRd
build_design <- function(meta, terms, primary = NULL, intercept = TRUE) {
  sf <- safe_frame(meta, terms)
  rhs <- paste(c(if (!intercept) "0", sf$ids), collapse = " + ")
  if (!nzchar(rhs)) rhs <- "1"
  mm <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data = sf$df)

  # The identifiers are generated, so a collision here means the generator is
  # wrong rather than that the data is unusual. Fail loudly: a duplicated
  # column name is how a contrast ends up pointing at a covariate.
  if (anyDuplicated(colnames(mm)) > 0L) {
    stop("Internal design column names collided; this is a bug in build_design().",
         call. = FALSE)
  }

  structure(mm,
    term_labels  = sf$terms,
    level_labels = sf$levels,
    primary_id   = if (is.null(primary)) NULL else sf$ids[[match(primary, terms)]]
  )
}

#' The design column holding one level of one term
#'
#' @param design A `build_design()` matrix.
#' @param id Term identifier.
#' @param label The user-facing level label.
#'
#' @return Column index, or `NA_integer_` when the level is absorbed into the
#'   intercept.
#' @keywords internal
#' @noRd
design_column <- function(design, id, label) {
  labels <- attr(design, "level_labels")[[id]]
  i <- match(label, labels)
  if (is.na(i)) {
    stop("Level not present in the model term: ", label, call. = FALSE)
  }
  match(paste0(id, level_id(i)), colnames(design))
}

#' Numeric contrast for `case - control` on the primary term
#'
#' @param design A `build_design()` matrix built with `intercept = FALSE`.
#' @param control_group,case_group User-facing level labels.
#'
#' @return A one-column contrast matrix suitable for `limma::contrasts.fit()`.
#' @keywords internal
#' @noRd
design_contrast <- function(design, control_group, case_group) {
  id <- attr(design, "primary_id")
  if (is.null(id)) {
    stop("This design has no primary term to contrast.", call. = FALSE)
  }
  k <- c(design_column(design, id, control_group),
         design_column(design, id, case_group))
  if (anyNA(k)) {
    stop("Contrast columns are missing from the design; the two groups are ",
         "not separately estimable.", call. = FALSE)
  }
  v <- numeric(ncol(design))
  v[[k[[1L]]]] <- -1
  v[[k[[2L]]]] <- 1
  matrix(v, ncol = 1L,
         dimnames = list(colnames(design),
                         paste(case_group, "-", control_group)))
}

#' The generated identifier standing for one user-facing level
#'
#' DESeq2 addresses a contrast as `c(variable, level, level)`, so it needs the
#' renamed level on its own rather than the design column name.
#'
#' @keywords internal
#' @noRd
safe_level <- function(sf, id, label) {
  i <- match(label, sf$levels[[id]])
  if (is.na(i)) {
    stop("Level not present in the model term: ", label, call. = FALSE)
  }
  level_id(i)
}

#' Coefficient name for one level of the primary term
#'
#' For backends that address a coefficient by name (`lm`, `edgeR`) rather than
#' by contrast matrix. The name is generated, so it cannot collide with a
#' covariate the way `paste0(group_col, case_group)` could.
#'
#' @keywords internal
#' @noRd
design_coef_name <- function(design_or_frame, id, label) {
  labels <- if (is.matrix(design_or_frame)) {
    attr(design_or_frame, "level_labels")[[id]]
  } else {
    design_or_frame$levels[[id]]
  }
  i <- match(label, labels)
  if (is.na(i)) {
    stop("Level not present in the model term: ", label, call. = FALSE)
  }
  paste0(id, level_id(i))
}
