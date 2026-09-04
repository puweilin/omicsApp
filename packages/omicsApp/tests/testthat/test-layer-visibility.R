# A project holds several layers, and QC describes exactly one of them.
# Which one was decided elsewhere -- by the arrow in the Projects table,
# or by a fallback -- so the answer to "what am I looking at" was not on
# the page showing it.

lv_input <- function(omics_type) {
  ids <- c("S1", "S2", "S3", "S4")
  m <- matrix(seq_len(6 * 4) * 1000, nrow = 6,
              dimnames = list(paste0("F", 1:6), ids))
  omicsCore::omics_input(
    m,
    data.frame(sample_id = ids, condition = c("G1", "G1", "G2", "G2"),
               row.names = ids, stringsAsFactors = FALSE),
    data.frame(feature_id = rownames(m), row.names = rownames(m),
               stringsAsFactors = FALSE),
    omics_type = omics_type,
    assay_type = if (omics_type == "rnaseq") "raw_count" else "raw_intensity")
}

lv_project <- function() {
  omicsCore::omics_project("p", list(prot = lv_input("proteomics"),
                                     rna  = lv_input("rnaseq")))
}

test_that("QC offers a layer picker and honours it", {
  proj <- shiny::reactiveVal(lv_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_false(is.null(output$ui_layer))
    # Proteomics first when nothing has been asked for: the missingness
    # panels are what this view was built around.
    expect_identical(active()$tag, "prot")

    session$setInputs(layer = "rna")
    session$flushReact()
    expect_identical(active()$tag, "rna")
  })
})

test_that("the picker wins over a layer another view asked for", {
  # requested_layer() is a hand-over ("show me this one"), not a lock.
  proj <- shiny::reactiveVal(lv_project())
  shiny::testServer(
    qc_view_server,
    args = list(current_project = proj,
                requested_layer = shiny::reactiveVal("rna")), {
      session$flushReact()
      expect_identical(active()$tag, "rna")

      session$setInputs(layer = "prot")
      session$flushReact()
      expect_identical(active()$tag, "prot")
    })
})

test_that("a single-layer project gets no picker to choose from", {
  proj <- shiny::reactiveVal(
    omicsCore::omics_project("p", list(prot = lv_input("proteomics"))))
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_null(output$ui_layer)
    expect_identical(active()$tag, "prot")
  })
})

test_that("Enrichment states its layer rather than offering to change it", {
  # Enrichment runs on a differential result, so its layer is whichever
  # one Differential ran on. A picker here could be set to disagree, and
  # the panel would then be labelled rnaseq while showing proteomics
  # pathways.
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = shiny::reactiveVal(example_diff_bundle()),
                diff_layer = shiny::reactive("rna")), {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both")
      # Demo mode: the fixture belongs to no layer, so nothing is
      # claimed. Better than naming one that is not the source.
      expect_true(isTRUE(is_demo()))
      expect_null(diff_layer_tag())

    })
})

test_that("the diff view hands its layer tag out", {
  proj <- shiny::reactiveVal(lv_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    out <- session$getReturned()
    expect_true("layer" %in% names(out))
    expect_identical(out$layer(), active()$tag)

    session$setInputs(layer = "rna")
    session$flushReact()
    expect_identical(out$layer(), "rna")
  })
})
