# A spreadsheet's "Total" or "Mean" column beside the samples used to be
# imported as one more sample. They are dropped with a note, as summary
# rows already were -- by name when the name can only be a summary, by
# value when the column equals a statistic of the others -- and a real
# sample whose name merely contains such a word is kept.

GENES20 <- c("TP53", "EGFR", "MYC", "KRAS", "BRCA1", "PTEN", "AKT1", "GAPDH",
             "ACTB", "VEGFA", "IL6", "TNF", "CD4", "CD8A", "STAT3", "JUN",
             "FOS", "SOX2", "NOTCH1", "ESR1")

sample_block <- function(samples, seed = 4) {
  set.seed(seed)
  matrix(round(2^stats::rnorm(length(GENES20) * length(samples), 16, 1)),
         length(GENES20), dimnames = list(GENES20, samples))
}

read_block <- function(df) {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(gene = rownames(df), df, check.names = FALSE), path,
                   row.names = FALSE)
  suppressWarnings(read_omics(path, omics_type = "proteomics", assay_type = "raw_intensity",
                              orientation = "features_in_rows"))
}

test_that("a summary name is enough on its own, alone or with a qualifier", {
  m <- sample_block(paste0("S", 1:4))
  set.seed(9)
  # Values that are no statistic of the samples: the name decides.
  odd <- function() round(stats::runif(nrow(m), 1, 100), 1)
  df <- cbind(m, Total = odd(), `Mean intensity` = odd(), `Row total` = odd(),
              `CV (%)` = odd(), StdDev = odd(), Average = odd())
  res <- read_block(df)
  expect_identical(colnames(res$input$expr_mat), paste0("S", 1:4))
  expect_true(paste0("Dropped 6 summary column(s): Total, Mean intensity, Row total, ",
                     "CV (%), StdDev, Average.") %in% res$report$warnings)
})

test_that("a column that equals a statistic of the samples is dropped by its values", {
  m <- sample_block(paste0("S", 1:5))
  df <- cbind(m,
              # Named like a sample, but it is every row's sum.
              Pooled = rowSums(m),
              # Bare words that can be names; the values say otherwise.
              N = ncol(m), Max = apply(m, 1L, max),
              # A mean rounded to two decimals still matches.
              Mean = round(rowMeans(m), 2))
  res <- read_block(df)
  expect_identical(colnames(res$input$expr_mat), paste0("S", 1:5))
  note <- grep("^Dropped 4 summary column", res$report$warnings, value = TRUE)
  expect_length(note, 1L)
  for (nm in c("Pooled", "N", "Max", "Mean")) expect_match(note, nm, fixed = TRUE)
})

test_that("a sample whose name only contains a summary word is kept", {
  samples <- c("Control_Mean1", "Total_RNA_1", "N", "T", "Max_2", "Sum3")
  m <- sample_block(samples)
  res <- read_block(m)
  expect_identical(colnames(res$input$expr_mat), samples)
  expect_false(any(grepl("summary column", res$report$warnings)))
})

test_that("summary columns and summary rows go together, each with its own note", {
  m <- sample_block(paste0("S", 1:4))
  m <- rbind(m, Total = colSums(m))
  df <- cbind(m, Total = rowSums(m))
  res <- read_block(df)
  expect_identical(dim(res$input$expr_mat), c(20L, 4L))
  expect_true("Dropped 1 summary column(s): Total." %in% res$report$warnings)
  expect_true("Dropped 1 summary row(s): Total." %in% res$report$warnings)
})

test_that("checking the values stays quick on a large matrix", {
  set.seed(1)
  big <- matrix(stats::rnorm(20000 * 200, 20), 20000,
                dimnames = list(NULL, paste0("S", 1:200)))
  t <- system.time(drop <- find_summary_columns(big))[["elapsed"]]
  expect_false(any(drop))
  expect_lt(t, 5)
})
