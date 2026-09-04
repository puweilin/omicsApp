# The pairing card exists so that "these two samples are the same person"
# is something the reader sees before reading anything computed on it --
# and so that a user who imported without a donor column can fix it here
# rather than by re-importing and re-running everything.

pair_input <- function(ids, donor = NULL) {
  m <- matrix(seq_len(4 * length(ids)) * 1000, nrow = 4,
              dimnames = list(paste0("F", 1:4), ids))
  meta <- data.frame(sample_id = ids,
                     condition = rep(c("G1", "G2"), length.out = length(ids)),
                     row.names = ids, stringsAsFactors = FALSE)
  if (!is.null(donor)) meta$donor <- donor
  omicsCore::omics_input(
    m, meta,
    data.frame(feature_id = rownames(m), row.names = rownames(m)),
    omics_type = "proteomics", assay_type = "raw_intensity")
}

# as.character() on a tag yields one element per node, and expect_match
# defaults to requiring all of them -- the leading <svg> never matches.
note_text <- function(x) paste(as.character(x), collapse = "")

pair_project <- function(a, b, a_donor = NULL, b_donor = NULL) {
  omicsCore::omics_project("p", list(prot = pair_input(a, a_donor),
                                     rna  = pair_input(b, b_donor)))
}

test_that("a guessed pairing offers to be accepted, and saves when it is", {
  proj <- shiny::reactiveVal(
    pair_project(c("RD001-C", "RD002-C"), c("RD001_Folli", "RD002_Folli")))

  shiny::testServer(integration_view_server,
                    args = list(current_project = proj), {
    expect_identical(pairing()$source, "suggested")
    expect_equal(nrow(pairing()$pairs), 2L)
    # Offered, not applied: a wrong pairing correlates different people
    # and still produces a plausible-looking result.
    expect_false(is.null(output$pairing_action))

    session$setInputs(accept_pairing = 1)

    link <- current_project()$sample_link
    expect_equal(nrow(link), 4L)
    expect_identical(sort(unique(link$donor_id)), c("RD001", "RD002"))
    # Saved on the project, so it is no longer a guess.
    expect_identical(pairing()$source, "linked")
    expect_null(output$pairing_action)
  })
})

test_that("a stated donor column needs no confirmation", {
  proj <- shiny::reactiveVal(
    pair_project(c("RD001-C", "RD002-C"), c("RD001_Folli", "RD002_Folli"),
                 c("RD001", "RD002"), c("RD001", "RD002")))

  shiny::testServer(integration_view_server,
                    args = list(current_project = proj), {
    expect_identical(pairing()$source, "donor")
    expect_null(output$pairing_action)
  })
})

test_that("no pairing says what to do instead of failing later", {
  proj <- shiny::reactiveVal(pair_project(c("A-1", "A-2"), c("B_x", "B_y")))

  shiny::testServer(integration_view_server,
                    args = list(current_project = proj), {
    expect_equal(nrow(pairing()$pairs), 0L)
    html <- note_text(output$pairing_note)
    expect_match(html, "donor", fixed = TRUE)
    # The reassurance matters as much as the instruction: without it the
    # user assumes a re-import.
    expect_match(html, "re-run", ignore.case = TRUE)
  })
})

test_that("one layer is not a pairing problem, it is a missing layer", {
  proj <- shiny::reactiveVal(
    omicsCore::omics_project("p", list(prot = pair_input(c("S1", "S2")))))

  shiny::testServer(integration_view_server,
                    args = list(current_project = proj), {
    expect_null(pairing())
    expect_match(note_text(output$pairing_note), "second layer")
  })
})
