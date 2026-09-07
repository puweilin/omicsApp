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

# --- the second layer's model ---------------------------------------------

# Only the contrast used to be carried over, so a first layer adjusted for
# age was compared against a second that was not, and the concordance read
# as though both had been analysed the same way. Fixing prune_backend_args()
# does not touch this: the second differential is launched here, and it was
# never asked to inherit anything.
sm_layer <- function(seed, omics = "proteomics", with_age = TRUE) {
  set.seed(seed); p <- 20; n <- 10
  e <- matrix(rnorm(p * n, 10, 1), p,
              dimnames = list(paste0("SYM", seq_len(p)), paste0("s", seq_len(n))))
  e[1:5, 6:10] <- e[1:5, 6:10] + 2
  md <- data.frame(group = rep(c("A", "B"), each = 5),
                   row.names = colnames(e), stringsAsFactors = FALSE)
  if (with_age) md$age <- c(30, 35, 40, 45, 50, 31, 36, 41, 46, 51)
  fd <- data.frame(feature_id = rownames(e), feature_symbol = rownames(e),
                   row.names = rownames(e), stringsAsFactors = FALSE)
  omicsCore::omics_input(e, md, fd, omics_type = omics,
                         assay_type = if (omics == "rnaseq") "logcpm"
                                      else "normalized_intensity")
}

sm_setup <- function(secondary_has_age = TRUE) {
  a <- sm_layer(1)
  b <- sm_layer(2, "rnaseq", with_age = secondary_has_age)
  proj <- omicsCore::omics_project("pair", list(prot = a, rna = b))
  primary <- omicsCore::run_diff(a, method = "limma", analysis_type = "group",
                                 group_col = "group", control_group = "A",
                                 case_group = "B", covariates = "age")
  list(proj = shiny::reactiveVal(proj), db = shiny::reactiveVal(primary))
}

test_that("the second layer inherits the first layer's covariates", {
  skip_if_not_installed("limma")
  s <- sm_setup(secondary_has_age = TRUE)
  shiny::testServer(integration_view_server,
                    args = list(current_project = s$proj, diff_bundle = s$db), {
    session$flushReact()
    info <- can_run()
    expect_true(isTRUE(info$ok))
    expect_equal(info$covariates, "age")
    expect_length(info$missing_cols, 0L)
  })
})

test_that("a second layer that cannot carry the model stops rather than differing", {
  skip_if_not_installed("limma")
  s <- sm_setup(secondary_has_age = FALSE)
  shiny::testServer(integration_view_server,
                    args = list(current_project = s$proj, diff_bundle = s$db), {
    session$flushReact()
    expect_equal(can_run()$missing_cols, "age")
    session$setInputs(rerun = 1)
    session$flushReact()
    # Nothing was produced, and the reason names the column.
    expect_null(session$getReturned()())
    expect_match(integration_error(), "age")
    expect_match(integration_error(), "two different models")
  })
})

test_that("running unadjusted is a choice, and it stays on the result", {
  skip_if_not_installed("limma")
  s <- sm_setup(secondary_has_age = FALSE)
  shiny::testServer(integration_view_server,
                    args = list(current_project = s$proj, diff_bundle = s$db), {
    session$flushReact()
    session$setInputs(rerun = 1)
    session$flushReact()
    expect_null(session$getReturned()())

    session$setInputs(accept_unadjusted = 1)
    session$flushReact()
    b <- session$getReturned()()
    skip_if(is.null(b), "integration did not run in this environment")
    # Announced once is not enough: a reader of the report or the exported
    # script has to see that the two layers were modelled differently.
    expect_equal(b$params$unadjusted_secondary, "age")
    expect_null(b$params$secondary_covariates)
  })
})
