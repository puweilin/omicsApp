# Contrasts beyond "each group against the control".
#
# A contrast is a set of weights over the levels of the group column that
# sum to zero: `B - A` is (A = -1, B = 1), `(B + C)/2 - A` is
# (A = -1, B = 0.5, C = 0.5). Every backend that fits one model over all
# the groups (limma, DESeq2, edgeR) can read any such contrast off that
# fit; the per-pair backends (t-test, lm) can only do simple pairs.
#
# Contrasts are written the way limma's makeContrasts() takes them, with
# one difference: a level that is not a valid R name is written in
# backticks (`` `Drug A` - Control ``), and the weights are worked out
# here rather than by makeContrasts() -- which is what lets "24h" or
# "KO-1" be a group at all.

#' Every pairwise comparison between groups
#'
#' @param levels Group labels, in the order they should be compared. Each
#'   later level is compared with each earlier one, so with a control
#'   listed first every comparison reads "treatment vs control".
#' @return A character vector of contrasts for `run_diff(contrasts = )`,
#'   e.g. `c("B - A", "C - A", "C - B")`.
#' @export
#' @family diff
#' @examples
#' pairwise_contrasts(c("Control", "TreatA", "TreatB"))
pairwise_contrasts <- function(levels) {
  assert_labels(levels, "levels")
  levels <- as.character(levels)
  if (length(levels) < 2L) {
    stop("`levels` needs at least two groups to compare.", call. = FALSE)
  }
  out <- character(0)
  for (j in seq_along(levels)[-1L]) {
    for (i in seq_len(j - 1L)) {
      out <- c(out, paste(contrast_level_token(levels[[j]]), "-",
                          contrast_level_token(levels[[i]])))
    }
  }
  out
}

# A level as it has to be written inside a contrast expression.
contrast_level_token <- function(x) {
  x <- as.character(x)
  if (identical(make.names(x), x)) x else paste0("`", gsub("`", "\\\\`", x), "`")
}

#' Parse contrast expressions into weights over the group levels
#'
#' @param contrasts Character vector of contrast expressions, or the
#'   single string `"pairwise"` for every pair of `levels`.
#' @param levels The levels present in the group column.
#' @return A list, one element per contrast, each a list with `label`
#'   (the `comparison` value the results carry), `spec` (the expression
#'   as written), `weights` (named numeric over `levels`, zero-sum),
#'   and `case` / `control` when the contrast is a simple pair.
#' @keywords internal
#' @noRd
parse_diff_contrasts <- function(contrasts, levels) {
  levels <- as.character(levels)
  if (identical(contrasts, "pairwise")) contrasts <- pairwise_contrasts(levels)
  if (!is.character(contrasts) || !length(contrasts) || anyNA(contrasts) ||
      !all(nzchar(trimws(contrasts)))) {
    stop("`contrasts` must be a character vector of contrast expressions ",
         "such as \"B - A\", or \"pairwise\".", call. = FALSE)
  }
  out <- lapply(contrasts, function(spec) {
    w <- contrast_weights(spec, levels)
    nz <- w[abs(w) > 1e-12]
    pair <- length(nz) == 2L && isTRUE(all.equal(sort(unname(nz)), c(-1, 1)))
    case <- if (pair) names(nz)[nz > 0] else NULL
    control <- if (pair) names(nz)[nz < 0] else NULL
    list(
      label = if (pair) paste0(case, "_vs_", control) else normalize_contrast_label(spec),
      spec = trimws(spec),
      weights = w,
      case = case,
      control = control
    )
  })
  labels <- vapply(out, `[[`, character(1), "label")
  if (anyDuplicated(labels)) {
    stop("The same comparison appears twice in `contrasts`: ",
         paste(unique(labels[duplicated(labels)]), collapse = ", "), ".",
         call. = FALSE)
  }
  out
}

# The pairs "each case against the control" as contrast specs.
case_control_contrasts <- function(control_group, case_group) {
  lapply(as.character(case_group), function(cg) {
    lv <- c(as.character(control_group), cg)
    w <- stats::setNames(c(-1, 1), lv)
    list(label = paste0(cg, "_vs_", control_group),
         spec = paste(contrast_level_token(cg), "-",
                      contrast_level_token(control_group)),
         weights = w, case = cg, control = as.character(control_group))
  })
}

# The levels a set of contrasts actually uses, in the order given.
contrast_levels <- function(specs, order = NULL) {
  used <- unique(unlist(lapply(specs, function(s) names(s$weights)[abs(s$weights) > 1e-12])))
  if (!is.null(order)) used <- c(intersect(as.character(order), used), setdiff(used, order))
  used
}

# Weights of one expression over `levels`: a small evaluator over the
# parsed expression that allows only what a contrast is made of -- group
# names, numbers, + - * / and brackets -- so a contrast can be checked
# for being a contrast, and so nothing in it is ever run as R code.
contrast_weights <- function(spec, levels) {
  expr <- tryCatch(str2lang(spec), error = function(e) {
    stop("Could not read the contrast \"", spec, "\": ", conditionMessage(e),
         ". Write group names that are not plain words in backticks, ",
         "e.g. `Drug A` - Control.", call. = FALSE)
  })
  zero <- stats::setNames(rep(0, length(levels)), levels)
  is_const <- function(x) is.numeric(x) && is.null(names(x))
  ev <- function(e) {
    if (is.numeric(e) && length(e) == 1L) return(as.numeric(e))
    if (is.name(e)) {
      nm <- as.character(e)
      if (!nm %in% levels) {
        stop("\"", nm, "\" in the contrast \"", spec, "\" is not a group. ",
             "Groups: ", paste(levels, collapse = ", "), ".", call. = FALSE)
      }
      v <- zero
      v[[nm]] <- 1
      return(v)
    }
    if (is.call(e)) {
      op <- as.character(e[[1L]])
      args <- lapply(as.list(e)[-1L], ev)
      if (op == "(" && length(args) == 1L) return(args[[1L]])
      if (op %in% c("+", "-") && length(args) == 1L) {
        return(if (op == "-") -args[[1L]] else args[[1L]])
      }
      if (op %in% c("+", "-") && length(args) == 2L) {
        a <- args[[1L]]; b <- args[[2L]]
        if (is_const(a) != is_const(b)) {
          stop("The contrast \"", spec, "\" adds a number to a group; a ",
               "contrast is a weighted sum of groups only.", call. = FALSE)
        }
        return(if (op == "+") a + b else a - b)
      }
      if (op == "*" && length(args) == 2L) {
        if (!is_const(args[[1L]]) && !is_const(args[[2L]])) {
          stop("The contrast \"", spec, "\" multiplies two groups.", call. = FALSE)
        }
        return(args[[1L]] * args[[2L]])
      }
      if (op == "/" && length(args) == 2L) {
        if (!is_const(args[[2L]]) || args[[2L]] == 0) {
          stop("The contrast \"", spec, "\" divides by something other than ",
               "a non-zero number.", call. = FALSE)
        }
        return(args[[1L]] / args[[2L]])
      }
    }
    stop("The contrast \"", spec, "\" uses something other than group names, ",
         "numbers, + - * / and brackets.", call. = FALSE)
  }
  w <- ev(expr)
  if (is_const(w)) {
    stop("The contrast \"", spec, "\" names no group.", call. = FALSE)
  }
  if (abs(sum(w)) > 1e-8) {
    stop("The weights of \"", spec, "\" sum to ", format(sum(w), digits = 3),
         ", not 0, so it is not a comparison between groups ",
         "(e.g. \"(A + B)/2 - C\", not \"A + B - C\").", call. = FALSE)
  }
  if (all(abs(w) < 1e-12)) {
    stop("The contrast \"", spec, "\" cancels out to nothing.", call. = FALSE)
  }
  w
}

normalize_contrast_label <- function(spec) {
  gsub("\\s+", " ", trimws(gsub("`", "", spec)))
}

# Weights as a (levels x contrasts) matrix in a given level order.
contrast_matrix_from_specs <- function(specs, levels) {
  m <- vapply(specs, function(s) {
    w <- stats::setNames(rep(0, length(levels)), levels)
    # The parsed weights may name every level of the column, with zeros
    # for the groups a contrast leaves out; only the fitted levels have a
    # row here. Assigning all of them grew the vector, and any contrast
    # not touching every group ("TreatA - Control" beside a TreatB)
    # failed with "values must be length 2".
    sw <- s$weights[abs(s$weights) > 0]
    missing <- setdiff(names(sw), levels)
    if (length(missing)) {
      stop("The contrast \"", s$label, "\" uses group(s) that are not in the fit: ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    w[names(sw)] <- sw
    w
  }, numeric(length(levels)))
  m <- matrix(m, nrow = length(levels),
              dimnames = list(levels, vapply(specs, `[[`, character(1), "label")))
  m
}

#' Check contrasts against a set of groups
#'
#' What [run_diff()] would call each contrast, or an error that says why
#' a contrast is not one -- for checking a contrast against a layer before
#' spending a model fit on it.
#'
#' @param contrasts Contrast expressions, as for `run_diff(contrasts = )`.
#' @param levels The groups available.
#' @return The `comparison` labels, one per contrast.
#' @export
#' @family diff
#' @examples
#' contrast_labels(c("B - A", "(B + C)/2 - A"), c("A", "B", "C"))
contrast_labels <- function(contrasts, levels) {
  assert_labels(levels, "levels")
  vapply(parse_diff_contrasts(contrasts, as.character(levels)),
         `[[`, character(1), "label")
}
