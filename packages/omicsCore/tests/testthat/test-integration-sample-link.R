# Two layers of one study rarely share sample ids, and should not:
# RD001-C is a cheek biopsy and RD001_Folli is a hair follicle. They are
# different samples of the same person. Pairing them wrong produces
# correlations between different people, which look like results.

sl_input <- function(ids, donor = NULL) {
  m <- matrix(seq_len(4 * length(ids)) * 1000, nrow = 4,
              dimnames = list(paste0("F", 1:4), ids))
  meta <- data.frame(sample_id = ids,
                     condition = rep(c("G1", "G2"), length.out = length(ids)),
                     row.names = ids, stringsAsFactors = FALSE)
  if (!is.null(donor)) meta$donor <- donor
  omics_input(m, meta,
              data.frame(feature_id = rownames(m), row.names = rownames(m)),
              omics_type = "proteomics", assay_type = "raw_intensity")
}

sl_project <- function(a_ids, b_ids, a_donor = NULL, b_donor = NULL) {
  omics_project("p", list(prot = sl_input(a_ids, a_donor),
                          rna  = sl_input(b_ids, b_donor)))
}

test_that("a stated donor column pairs the layers", {
  p <- sl_project(c("RD001-C", "RD002-C"), c("RD001_Folli", "RD002_Folli"),
                  c("RD001", "RD002"), c("RD001", "RD002"))
  res <- sample_pairing_preview(p, "prot", "rna")

  expect_identical(res$source, "donor")
  expect_identical(res$pairs$donor_id, c("RD001", "RD002"))
  expect_identical(res$pairs$a, c("RD001-C", "RD002-C"))
  expect_identical(res$pairs$b, c("RD001_Folli", "RD002_Folli"))
})

test_that("without a donor column the ids are used to suggest one", {
  # The recovery path. A user who imported without donor would otherwise
  # have to re-import and re-run everything to reach Integration, and
  # nothing else in the analysis reads donor at all.
  p <- sl_project(c("RD001-C", "RD002-C"), c("RD001_Folli", "RD002_Folli"))
  res <- sample_pairing_preview(p, "prot", "rna")

  expect_identical(res$source, "suggested")
  expect_identical(res$pairs$donor_id, c("RD001", "RD002"))
})

test_that("an ambiguous stem yields no suggestion at all", {
  # A-1 and A-2 both stem to A, which would pair two samples of one layer
  # to one donor. A suggestion worse than none: acting on it still
  # produces a result.
  p <- sl_project(c("A-1", "A-2"), c("A_x", "A_y"))
  res <- sample_pairing_preview(p, "prot", "rna")

  expect_identical(res$source, "none")
  expect_equal(nrow(res$pairs), 0L)
})

test_that("layers that genuinely share ids need no donor", {
  p <- sl_project(c("S1", "S2"), c("S1", "S2"))
  res <- sample_pairing_preview(p, "prot", "rna")

  expect_identical(res$source, "sample_id")
  expect_identical(res$pairs$a, res$pairs$b)
})

test_that("a stated donor beats a guess the ids would have supported", {
  # The ids would stem to RD001/RD002, but the study says otherwise.
  # Whoever wrote the metadata knows something the id format does not.
  p <- sl_project(c("RD001-C", "RD002-C"), c("RD001_Folli", "RD002_Folli"),
                  c("subjB", "subjA"), c("subjA", "subjB"))
  res <- sample_pairing_preview(p, "prot", "rna")

  expect_identical(res$source, "donor")
  paired <- res$pairs[res$pairs$donor_id == "subjA", ]
  expect_identical(paired$a, "RD002-C")
  expect_identical(paired$b, "RD001_Folli")
})

test_that("samples missing from one layer are dropped, not paired blindly", {
  p <- sl_project(c("RD001-C", "RD002-C", "RD003-C"),
                  c("RD001_Folli", "RD003_Folli"),
                  c("RD001", "RD002", "RD003"), c("RD001", "RD003"))
  res <- sample_pairing_preview(p, "prot", "rna")

  expect_identical(sort(res$pairs$donor_id), c("RD001", "RD003"))
})

test_that("donor_column finds the usual spellings and nothing else", {
  for (nm in c("donor", "Donor", "subject_id", "patient")) {
    df <- data.frame(a = 1, stringsAsFactors = FALSE)
    names(df) <- nm
    expect_identical(donor_column(df), nm, info = nm)
  }
  expect_null(donor_column(data.frame(condition = "G1")))
  expect_null(donor_column(NULL))
})

test_that("sample_stem cuts at the first separator", {
  expect_identical(sample_stem(c("RD001-C", "RD001_Folli", "RD001.x", "RD001")),
                   rep("RD001", 4L))
})

test_that("derive_sample_link needs the column on both sides", {
  one_sided <- sl_project(c("RD001-C"), c("RD001_Folli"),
                          a_donor = "RD001")
  expect_null(derive_sample_link(one_sided))
})
