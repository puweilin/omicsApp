# Per-feature statistics for every feature at once.
#
# The lm, t-test and limma-continuous backends used to fit one model per
# feature in an R loop -- an lm() + summary() or a t.test() for each of
# tens of thousands of rows, repeated for every contrast. These compute
# the same quantities with one QR decomposition per pattern of missing
# values (or row-wise sums for the t-test). Each reproduces the base-R
# call it replaces -- lm() + summary.lm(), t.test(), cor.test() -- to
# all.equal tolerance, including their edge cases (an all-missing row, a
# factor left with one level, a constant row), which now give NA rather
# than stopping the run. Measured on 8,000 proteins x 60 samples: lm with
# three contrasts 33 s -> 3.7 s; a pairwise t-test 6.4 s -> 0.36 s; the
# limma continuous fit 8.9 s -> 0.9 s.

# Which columns of the design belong to which level of each factor, so a
# subset of samples gets the design model.frame(drop.unused.levels = TRUE)
# would build for it: absent levels lose their column, and when the
# reference level is absent the first present level becomes the
# reference. Without this a confounded subset is aliased differently
# from lm() and gives another estimate.
make_fac_map <- function(X, fac) {
  if (is.null(fac) || !ncol(fac)) return(list())
  lapply(names(fac), function(v) {
    f <- fac[[v]]
    levs <- if (is.factor(f)) levels(f) else if (is.logical(f)) c("FALSE", "TRUE")
            else sort(unique(as.character(f)))
    cols <- match(paste0(v, levs[-1L]), colnames(X))
    if (anyNA(cols)) stop("cannot map factor columns for ", v)
    list(var = v, levels = levs, cols = cols)
  })
}

# summary(lm(y ~ <design>)) for each row of Y.
#
# Y        features x samples (may contain NA; lm drops those samples)
# X        full model matrix, samples x k (no NA)
# coef_col name of the coefficient to report (or NULL)
# fac      data.frame of the factor / character design variables; when one
#          has fewer than 2 distinct values among a feature's observed
#          samples lm() errors, and the feature gets NA.
# build_X  optional function(obs) -> model matrix for that subset (for
#          data-dependent bases such as splines::ns()); NULL means X[obs, ].
#
# Rows sharing a pattern of missing values share one QR. qr() is LINPACK
# dqrdc2 with tol 1e-7, as in lm.fit(), so aliasing decisions match.
lm_rows <- function(Y, X, coef_col = NULL, fac = NULL, build_X = NULL,
                    want_resid = FALSE, fac_map = NULL) {
  p <- nrow(Y)
  beta <- tval <- pval <- adjr2 <- rep(NA_real_, p)
  resid <- if (want_resid) matrix(NA_real_, p, ncol(Y)) else NULL
  na <- is.na(Y)
  has_na <- rowSums(na) > 0L
  groups <- list()
  if (any(!has_na)) groups[[1L]] <- which(!has_na)
  if (any(has_na)) {
    idx <- which(has_na)
    key <- apply(na[idx, , drop = FALSE], 1L, function(r) paste(which(r), collapse = ","))
    groups <- c(groups, unname(split(idx, key)))
  }
  coef_name <- coef_col
  if (is.null(fac_map) && is.null(build_X)) fac_map <- make_fac_map(X, fac)
  for (rows in groups) {
    obs <- !na[rows[1L], ]
    n_obs <- sum(obs)
    if (n_obs == 0L) next                         # lm: 0 (non-NA) cases
    if (!is.null(fac) && ncol(fac)) {
      lv_ok <- vapply(fac, function(v) length(unique(v[obs])) >= 2L, logical(1))
      if (!all(lv_ok)) next                       # lm: contrasts error
    }
    if (is.null(build_X)) {
      drop <- integer(0)
      for (fm in fac_map) {
        present <- fm$levels[fm$levels %in% as.character(fac[[fm$var]][obs])]
        absent <- setdiff(fm$levels[-1L], present)
        drop <- c(drop, fm$cols[match(absent, fm$levels[-1L])])
        if (!fm$levels[1L] %in% present) drop <- c(drop, fm$cols[match(present[1L], fm$levels[-1L])])
      }
      Xg <- if (length(drop)) X[obs, -drop, drop = FALSE] else X[obs, , drop = FALSE]
    } else Xg <- build_X(obs)
    if (is.null(Xg)) next
    q <- qr(Xg)
    Yg <- t(Y[rows, obs, drop = FALSE])
    rk <- q$rank
    rdf <- n_obs - rk
    res <- qr.resid(q, Yg)
    fit <- Yg - res
    rss <- colSums(res^2)
    fbar <- colMeans(fit)
    mss <- colSums((fit - rep(fbar, each = nrow(fit)))^2)
    r2 <- mss / (mss + rss)
    # summary.lm: a fit where only the intercept is estimable reports
    # r.squared = adj.r.squared = 0 rather than the formula's value.
    adjr2[rows] <- if (rk == 1L) 0 else 1 - (1 - r2) * ((n_obs - 1L) / rdf)
    if (want_resid) resid[rows, obs] <- t(res)
    coef_j <- if (is.null(coef_name)) NA_integer_ else match(coef_name, colnames(Xg))
    if (!is.na(coef_j)) {
      pos <- match(coef_j, q$pivot)
      if (pos <= rk) {
        cf <- qr.coef(q, Yg)
        b <- if (is.matrix(cf)) cf[coef_j, ] else cf[coef_j]
        Rinv <- chol2inv(q$qr[seq_len(rk), seq_len(rk), drop = FALSE])
        se <- sqrt(Rinv[pos, pos] * rss / rdf)
        tt <- b / se
        beta[rows] <- b
        tval[rows] <- tt
        pval[rows] <- 2 * stats::pt(abs(tt), rdf, lower.tail = FALSE)
      }
    }
  }
  list(beta = beta, t = tval, p = pval, adj_r2 = adjr2, resid = resid)
}

# cor.test(y, x, method = "spearman")$estimate on the complete pairs of
# each row; a non-finite value becomes NA, as cor.test gives.
spearman_rows <- function(Y, x) {
  out <- rep(NA_real_, nrow(Y))
  na <- is.na(Y) | matrix(is.na(x), nrow(Y), length(x), byrow = TRUE)
  full <- rowSums(na) == 0L
  if (any(full)) {
    R <- t(apply(Y[full, , drop = FALSE], 1L, rank))
    if (sum(full) == 1L) R <- matrix(R, nrow = 1L)
    rx <- rank(x)
    rc <- R - rowMeans(R)
    xc <- rx - mean(rx)
    out[full] <- (rc %*% xc)[, 1] / sqrt(rowSums(rc^2) * sum(xc^2))
  }
  for (i in which(!full)) {
    ok <- !na[i, ]
    if (sum(ok) < 2L) next
    out[i] <- suppressWarnings(stats::cor(rank(Y[i, ok]), rank(x[ok])))
  }
  out[!is.finite(out)] <- NA_real_
  out
}

design_factors <- function(df) {
  keep <- vapply(df, function(v) is.factor(v) || is.character(v) || is.logical(v), logical(1))
  df[, keep, drop = FALSE]
}

# t.test(case, ctrl, paired, var.equal) for every row at once, including
# t.test's refusals ("not enough observations", "data are essentially
# constant"), which give NA.
ttest_rows <- function(case_m, ctrl_m, paired = FALSE, var_equal = FALSE) {
  eps10 <- 10 * .Machine$double.eps
  rvar <- function(m, n, mu) {
    v <- rowSums((m - mu)^2, na.rm = TRUE) / (n - 1)
    v[n < 2] <- NA_real_
    v
  }
  if (paired) {
    d <- case_m - ctrl_m
    nx <- rowSums(!is.na(d))
    mx <- rowMeans(d, na.rm = TRUE)
    vx <- rvar(d, nx, mx)
    stderr <- sqrt(vx / nx)
    bad <- nx < 2 | (!is.na(stderr) & stderr < eps10 * abs(mx))
    t <- mx / stderr
    df <- nx - 1
  } else {
    nx <- rowSums(!is.na(case_m)); ny <- rowSums(!is.na(ctrl_m))
    mx <- rowMeans(case_m, na.rm = TRUE); my <- rowMeans(ctrl_m, na.rm = TRUE)
    vx <- rvar(case_m, nx, mx); vy <- rvar(ctrl_m, ny, my)
    if (var_equal) {
      df <- nx + ny - 2
      v <- ifelse(nx > 1, (nx - 1) * vx, 0) + ifelse(ny > 1, (ny - 1) * vy, 0)
      v <- v / df
      stderr <- sqrt(v * (1 / nx + 1 / ny))
      bad <- nx < 1 | ny < 1 | nx + ny < 3
    } else {
      sx <- sqrt(vx / nx); sy <- sqrt(vy / ny)
      stderr <- sqrt(sx^2 + sy^2)
      df <- stderr^4 / (sx^4 / (nx - 1) + sy^4 / (ny - 1))
      bad <- nx < 2 | ny < 2
    }
    bad <- bad | (!is.na(stderr) & stderr < eps10 * pmax(abs(mx), abs(my)))
    t <- (mx - my) / stderr
  }
  p <- 2 * stats::pt(-abs(t), df)
  t[bad] <- NA_real_; p[bad] <- NA_real_
  list(t = unname(t), p = unname(p))
}

# The adjusted R^2 and Spearman's rho reported beside the limma
# continuous fit.
limma_continuous_extras <- function(expr_mat, cont_vals, adj_df, method, df) {
  adjustment_terms <- names(adj_df)
  base <- data.frame(cont = cont_vals, adj_df, check.names = FALSE)
  if (method == "spline") {
    f <- if (length(adjustment_terms)) paste("~ splines::ns(cont, df =", df, ") +",
                                             paste(adjustment_terms, collapse = " + "))
         else paste("~ splines::ns(cont, df =", df, ")")
    build <- function(obs) stats::model.matrix(stats::as.formula(f), data = base[obs, , drop = FALSE])
    X <- build(rep(TRUE, nrow(base)))
  } else {
    f <- if (length(adjustment_terms)) paste("~ cont +", paste(adjustment_terms, collapse = " + ")) else "~ cont"
    X <- stats::model.matrix(stats::as.formula(f), data = base)
    build <- NULL
  }
  adj_r2 <- lm_rows(expr_mat, X, NULL, fac = design_factors(adj_df), build_X = build)$adj_r2
  names(adj_r2) <- rownames(expr_mat)
  if (length(adjustment_terms) == 0L) {
    rho <- spearman_rows(expr_mat, cont_vals)
  } else {
    Xa <- stats::model.matrix(stats::as.formula(paste("~", paste(adjustment_terms, collapse = " + "))), data = base)
    cont_res <- unname(stats::lm.fit(Xa, cont_vals)$residuals)
    # Residualising needs every sample, so the partial correlation is
    # given for the features observed in all of them and NA for the rest,
    # as run_lm_continuous() does. qr.resid() over the whole matrix
    # stopped the run at the first missing value ("NA/NaN/Inf in foreign
    # function call"): any proteomics layer with a gap, adjusted for a
    # covariate or paired, could not be analysed against a continuous
    # variable at all.
    full <- rowSums(is.na(expr_mat)) == 0L
    rho <- rep(NA_real_, nrow(expr_mat))
    if (any(full)) {
      yres <- t(qr.resid(qr(Xa), t(expr_mat[full, , drop = FALSE])))
      rho[full] <- spearman_rows(yres, cont_res)
    }
  }
  names(rho) <- rownames(expr_mat)
  list(adj_r2 = adj_r2, rho = rho)
}

# Pairwise-complete Pearson correlation between columns by matrix
# products: the same matrix as cor(use = "pairwise.complete.obs") to
# 1e-15, 40x faster on 8,000 x 300 with missing values (the sample
# connectivity check runs it on every QC change).
pairwise_cor <- function(X) {
  M <- !is.na(X)
  X0 <- X - rep(colMeans(X, na.rm = TRUE), each = nrow(X))   # centred: no cancellation
  X0[!M] <- 0
  Mf <- M * 1
  n <- crossprod(Mf)
  sx <- crossprod(X0, Mf)
  sxx <- crossprod(X0^2, Mf)
  sxy <- crossprod(X0)
  num <- sxy - sx * t(sx) / n
  vx <- sxx - sx^2 / n
  r <- num / sqrt(vx * t(vx))
  r[n < 2] <- NA_real_
  dimnames(r) <- list(colnames(X), colnames(X))
  r
}

# TMM factors, remembered by the content of the count matrix. Every QC
# change, the QC PCA, each limma / t-test / lm run and GSVA put a count
# layer on the log scale, and calcNormFactors() took 1.3 s at 30k x 60
# and ~15 s at 60k x 300 each time.
.tmm_cache <- new.env(parent = emptyenv())
cached_tmm <- function(m) {
  key <- rlang::hash(m)
  nf <- .tmm_cache[[key]]
  if (is.null(nf)) {
    nf <- edgeR::calcNormFactors(edgeR::DGEList(counts = m))$samples$norm.factors
    if (length(ls(.tmm_cache)) > 20L) rm(list = ls(.tmm_cache), envir = .tmm_cache)
    assign(key, nf, envir = .tmm_cache)
  }
  nf
}
