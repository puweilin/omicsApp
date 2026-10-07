# log2 normalisation with the samples aligned on their medians.

nm_input <- function() {
  set.seed(8)
  n_feat <- 300
  base <- stats::rnorm(n_feat, 20, 2)
  m <- matrix(stats::rnorm(n_feat * 6, base, 0.4), n_feat, 6,
              dimnames = list(paste0("P", seq_len(n_feat)), paste0("S", 1:6)))
  # Two samples loaded with more material: brighter across the board
  m[, 5:6] <- m[, 5:6] + 1.5
  m[sample(length(m), 120)] <- NA
  meta <- data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(m))
  feat <- data.frame(feature_id = rownames(m), row.names = rownames(m))
  omics_input(2^m, meta, feat, omics_type = "proteomics", assay_type = "raw_intensity")
}

test_that("center = 'median' puts every sample on the median of the sample medians", {
  inp <- nm_input()
  plain <- suppressMessages(normalize_omics(inp, method = "log2"))
  centred <- suppressMessages(normalize_omics(inp, method = "log2", center = "median"))

  plain_med <- apply(plain$expr_mat, 2, stats::median, na.rm = TRUE)
  med <- apply(centred$expr_mat, 2, stats::median, na.rm = TRUE)
  expect_equal(unname(med), rep(stats::median(plain_med), 6))
  # On the log2-intensity scale, not centred on zero
  expect_gt(min(med), 15)
  # Within a sample nothing moves relative to anything else
  d <- centred$expr_mat - plain$expr_mat
  expect_true(all(apply(d, 2, function(x) diff(range(x, na.rm = TRUE))) < 1e-12))
  expect_identical(is.na(centred$expr_mat), is.na(plain$expr_mat))
  expect_identical(dimnames(centred$expr_mat), dimnames(inp$expr_mat))
  expect_identical(centred$assay_type, "normalized_intensity")
  expect_identical(centred$normalization$center, "median")
})

test_that("the default is unchanged", {
  inp <- nm_input()
  a <- suppressMessages(normalize_omics(inp, method = "log2"))
  b <- suppressMessages(normalize_omics(inp, method = "log2", center = "none"))
  expect_identical(a$expr_mat, b$expr_mat)
  expect_identical(a$expr_mat, log2(inp$expr_mat + 1))
  expect_identical(a$normalization$center, "none")
})

test_that("median centring goes with log2 only, and says so", {
  inp <- nm_input()
  expect_error(normalize_omics(inp, method = "vsn", center = "median"),
               "goes with `method = \"log2\"`", fixed = TRUE)
  expect_error(normalize_omics(inp, method = "log2", center = "mean"),
               "should be one of")
})

test_that("the exported script repeats the centring and reproduces the values", {
  inp <- nm_input()
  centred <- suppressMessages(normalize_omics(inp, method = "log2", center = "median"))
  centred$source_path <- "upload.xlsx"
  proj <- omics_project("n", list(proteomics = centred))
  lines <- export_script(proj)
  start <- grep("normalize_omics(", lines, fixed = TRUE)
  expect_length(start, 1L)
  end <- start + which(lines[-seq_len(start)] == ")")[1L]
  call <- lines[start:end]
  expect_true(any(grepl('center *= "median"', call)))
  expect_silent(parse(text = lines))

  # Run the emitted call on the file's values: the same matrix comes out.
  env <- new.env(parent = asNamespace("omicsCore"))
  var <- sub(" <- normalize_omics\\($", "", lines[start])
  assign(var, inp, envir = env)
  suppressMessages(eval(parse(text = call), envir = env))
  expect_identical(get(var, envir = env)$expr_mat, centred$expr_mat)

  # An uncentred layer's script does not mention centring
  plain <- suppressMessages(normalize_omics(inp, method = "log2"))
  plain$source_path <- "upload.xlsx"
  lines2 <- export_script(omics_project("n", list(proteomics = plain)))
  expect_false(any(grepl("center", lines2, fixed = TRUE)))
})

test_that("the report says the samples were median-aligned", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  inp <- nm_input()
  centred <- suppressMessages(normalize_omics(inp, method = "log2", center = "median"))
  proj <- omics_project("n", list(proteomics = centred))
  out <- tempfile(fileext = ".html")
  on.exit(unlink(out), add = TRUE)
  suppressMessages(suppressWarnings(export_report(proj, out)))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  # Pandoc wraps lines, so words are separated by any whitespace.
  expect_match(html, "median\\s+aligned\\s+to\\s+the\\s+median\\s+of\\s+all\\s+samples")
})
