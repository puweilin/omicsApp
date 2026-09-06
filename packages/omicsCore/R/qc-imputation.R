# Imputation, using the method set the proteomics field already agrees on.
#
# These are DEP's methods, because DEP is what a proteomics analyst has
# read about and because divergence here is expensive: a result labelled
# "half_min" cannot be compared with anything in the literature, and an
# analyst who wants MinProb has no way to ask for it.
#
# DEP::impute() offers
#   bpca, knn, QRILC, MLE, MinDet, MinProb, man, min, zero, mixed, nbavg
# and delegates to MSnbase, which delegates to imputeLCMD for everything
# left-censored. We call imputeLCMD directly: identical arithmetic, none
# of the mzR/ProtGenerics stack, and two of its four dependencies
# (pcaMethods, impute) were already here.
#
# ---- why the default is MinProb and not "none" ----------------------
#
# Leaving NA is not neutral. limma drops a feature it cannot fit and PCA
# cannot scale a column it cannot centre, so "none" is itself a decision
# -- complete-case analysis -- taken silently and usually not the one
# the analyst would choose.
#
# Missingness in DIA/DDA is mostly MNAR: a protein is absent from the
# table because it fell below the detection limit, not at random. That
# is why the MAR methods (knn, MLE, bpca) are the wrong default here
# even though they are the familiar ones -- they infer the missing value
# from samples where the protein *was* seen, which is exactly the
# population it is not from.
#
# MinProb draws from a narrow distribution near the observed minimum,
# which is the left-censored story stated as arithmetic. DEP's own two
# vignettes use `fun = "MinProb", q = 0.01` as their worked example.
#
# Counts are a different world and keep "none": a zero count is an
# observation, and imputing it feeds a negative-binomial model numbers
# it never saw. resolve_impute_method() encodes that split.

#' Imputation methods, grouped by what they assume
#'
#' `IMPUTE_METHOD_ASSUMPTION` maps each method to the missingness it
#' assumes; `IMPUTE_METHODS` is its names, in the order the app lists
#' them, and is what [impute_matrix()] and [run_qc()] accept.
#'
#' @format A named character vector: names are the method, values are
#'   `"MNAR"`, `"MAR"`, `"either"` or `"none"`.
#' @export
#' @family qc
IMPUTE_METHOD_ASSUMPTION <- c(
  none     = "none",
  MinProb  = "MNAR",
  MinDet   = "MNAR",
  QRILC    = "MNAR",
  min      = "MNAR",
  zero     = "MNAR",
  knn      = "MAR",
  MLE      = "MAR",
  bpca     = "MAR",
  mixed    = "either",
  man      = "either"
)

#' @rdname IMPUTE_METHOD_ASSUMPTION
#' @export
IMPUTE_METHODS <- names(IMPUTE_METHOD_ASSUMPTION)

#' Which imputation a layer gets when the caller does not say
#'
#' Proteomics gets `MinProb`; anything counted gets `none`. Mirrors how
#' `run_qc()` already resolves `outlier_method` per omics type, and for
#' the same reason: the right answer differs by measurement, and a single
#' global default is wrong for one of them whichever it is.
#'
#' @param omics_type Modality string.
#' @return A single method name.
#' @export
#' @family qc
resolve_impute_method <- function(omics_type) {
  switch(omics_type %||% "",
    proteomics = "MinProb",
    rnaseq     = "none",
    "none"
  )
}

#' Impute missing values in an expression matrix
#'
#' The methods are DEP's, so a choice made here means the same thing it
#' means in the proteomics literature.
#'
#' * `"none"` — leave `NA`. Not neutral: downstream this becomes
#'   complete-case analysis.
#' * `"MinProb"` — random draws from a narrow gaussian near the observed
#'   minimum (MNAR). The default for proteomics.
#' * `"MinDet"` — a low quantile of the observed values, deterministic (MNAR).
#' * `"QRILC"` — quantile regression for left-censored data (MNAR).
#' * `"min"` — the observed minimum of that feature (MNAR).
#' * `"zero"` — zero (MNAR). Offered for parity; it distorts variance.
#' * `"knn"` — k-nearest neighbours (MAR).
#' * `"MLE"` — maximum likelihood (MAR).
#' * `"bpca"` — Bayesian PCA (MAR).
#' * `"mixed"` — MAR or MNAR per feature, chosen by a test.
#' * `"man"` — manual shift/scale, as in DEP's `shift`/`scale` arguments.
#'
#' MNAR methods assume log-scale data, which is the scale DEP's workflow
#' imputes on. On linear intensities they will place imputed values on a
#' scale the observed values are not on.
#'
#' @param mat A numeric matrix, features in rows.
#' @param method One of [IMPUTE_METHODS].
#' @param ... Forwarded to the backend (`q` for MinProb/MinDet, `k` for
#'   knn, `shift`/`scale` for man).
#'
#' @return A numeric matrix with the same dimensions and names as `mat`.
#' @export
#' @family qc
impute_matrix <- function(mat, method = IMPUTE_METHODS, ...) {
  method <- match.arg(method, IMPUTE_METHODS)
  mat <- as.matrix(mat)
  if (method == "none" || !anyNA(mat)) {
    return(mat)
  }

  out <- switch(method,
    min     = impute_row_min(mat),
    zero    = impute_constant(mat, 0),
    man     = impute_manual(mat, ...),
    MinProb = lcmd_call("impute.MinProb", mat, ...),
    MinDet  = lcmd_call("impute.MinDet", mat, ...),
    QRILC   = lcmd_qrilc(mat, ...),
    knn     = lcmd_call("impute.wrapper.KNN", mat, K = 10, ...),
    MLE     = lcmd_call("impute.wrapper.MLE", mat, ...),
    bpca    = impute_bpca(mat, ...),
    mixed   = impute_mixed(mat, ...)
  )
  dimnames(out) <- dimnames(mat)
  out
}

# ---- backends ----------------------------------------------------------

# imputeLCMD works on features-in-rows matrices, the same orientation as
# ours, so nothing is transposed here. Two of its dependencies
# (pcaMethods, impute) were already required by this package.
lcmd_require <- function(fn) {
  if (!is_installed("imputeLCMD")) {
    stop(
      "Package 'imputeLCMD' is required for this imputation method. ",
      "Install it with: omicsCore::install_optional('imputation') ",
      "or install.packages('imputeLCMD').",
      call. = FALSE
    )
  }
  getExportedValue("imputeLCMD", fn)
}

# Every draw-based method fixes its seed and puts the caller's stream
# back. Without it the same matrix imputes differently on each run and
# every downstream p-value moves with it -- so a script from
# export_script() would disagree with the report it came from.
with_fixed_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

# capture.output(), because impute.MinProb() prints an intermediate
# value. Harmless at a console and not harmless inside a Shiny render,
# where it lands in the log for every QC re-run.
lcmd_call <- function(fn, mat, ..., seed = 1234L) {
  f <- lcmd_require(fn)
  with_fixed_seed(seed, {
    out <- NULL
    utils::capture.output(out <- f(mat, ...))
    as.matrix(out)
  })
}

# QRILC returns a list; the imputed matrix is the first element.
lcmd_qrilc <- function(mat, ..., seed = 1234L) {
  f <- lcmd_require("impute.QRILC")
  with_fixed_seed(seed, {
    out <- NULL
    utils::capture.output(out <- f(mat, ...))
    as.matrix(out[[1L]])
  })
}

# DEP's "mixed": test each feature for MNAR and route it accordingly.
#
# imputeLCMD spells its MAR methods in capitals and dispatches on them
# with switch(), which has no default -- so a lowercase "knn" matches
# nothing, the branch never assigns, and the failure surfaces as
# `object 'dataSet.MCAR.imputed' not found` several frames away.
impute_mixed <- function(mat, mar = "KNN", mnar = "MinProb", ...) {
  sel <- lcmd_require("model.Selector")(mat)
  f <- lcmd_require("impute.MAR.MNAR")
  with_fixed_seed(1234L, {
    out <- NULL
    utils::capture.output(
      out <- f(mat, sel, method.MAR = mar, method.MNAR = mnar))
    as.matrix(out)
  })
}

impute_row_min <- function(mat) {
  # A feature with no observed value has no minimum. min() answers Inf
  # with a warning; the Inf was caught below and the warning was not,
  # and it was the only warning this file's tests emitted.
  row_min <- suppressWarnings(apply(mat, 1L, min, na.rm = TRUE))
  row_min[!is.finite(row_min)] <- 0
  idx <- which(is.na(mat), arr.ind = TRUE)
  mat[idx] <- row_min[idx[, 1L]]
  mat
}

impute_constant <- function(mat, value) {
  mat[is.na(mat)] <- value
  mat
}

# DEP's "man": draw from a gaussian shifted down from each sample's own
# distribution, with DEP's defaults for shift and scale.
impute_manual <- function(mat, shift = 1.8, scale = 0.3, seed = 1234L, ...) {
  with_fixed_seed(seed, {
    for (j in seq_len(ncol(mat))) {
      na_j <- is.na(mat[, j])
      if (!any(na_j)) next
      obs <- mat[!na_j, j]
      if (length(obs) < 2L) {
        mat[na_j, j] <- if (length(obs) == 1L) obs else 0
        next
      }
      mat[na_j, j] <- stats::rnorm(sum(na_j),
                                   mean = mean(obs) - shift * stats::sd(obs),
                                   sd   = scale * stats::sd(obs))
    }
    mat
  })
}

impute_bpca <- function(mat, n_pcs = 3, ...) {
  if (!is_installed("pcaMethods")) {
    stop(
      "Package 'pcaMethods' is required for method = 'bpca'. ",
      "Install it with: omicsCore::install_optional('imputation').",
      call. = FALSE
    )
  }
  pc <- pcaMethods::pca(t(mat), method = "bpca", nPcs = n_pcs, ...)
  t(pcaMethods::completeObs(pc))
}
