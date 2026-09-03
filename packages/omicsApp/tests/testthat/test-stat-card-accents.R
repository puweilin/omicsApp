# stat_card()'s accents live in two places -- match.arg() in R and
# `.stat-card.accent-*` in styles.scss -- and nothing connected them.
# The QC stats card asked for "warn" whenever a sample was flagged as an
# outlier and got match.arg()'s error instead, which reached the user as
# a red Shiny overlay the moment they switched outlier method to PCA.
# The import summary asks for the same accent above 50% missingness.

accents_in_r <- function() {
  eval(formals(stat_card)$accent)
}

accents_in_css <- function() {
  path <- system.file("app", "www", "styles.scss", package = "omicsApp")
  if (!nzchar(path)) path <- "../../inst/app/www/styles.scss"
  skip_if_not(file.exists(path), "styles.scss not found")
  css <- readLines(path, warn = FALSE)
  m <- regmatches(css, regexpr("\\.stat-card\\.accent-[a-z]+", css))
  unique(sub("^\\.stat-card\\.accent-", "", m))
}

test_that("every accent R accepts has a rule to style it", {
  # Without this an accepted accent renders as an unstyled card: no
  # error, just a stat card that silently looks like the others.
  expect_setequal(setdiff(accents_in_r(), accents_in_css()), character(0))
})

test_that("every accent the stylesheet defines is one R will accept", {
  # And the other direction, so a rule cannot be added for a value
  # match.arg() then rejects.
  expect_setequal(setdiff(accents_in_css(), accents_in_r()), character(0))
})

test_that("each accent renders the class its stylesheet rule targets", {
  for (a in accents_in_r()) {
    html <- as.character(stat_card("L", 1, accent = a))
    expect_true(grepl(paste0("accent-", a), html, fixed = TRUE), info = a)
  }
})

test_that("warn is accepted, since two callers already depend on it", {
  expect_silent(stat_card("Samples passing", "10 / 12", accent = "warn"))
})

# ---- the reported failure ---------------------------------------------

test_that("the QC stats card survives a sample being flagged", {
  # The flag is injected rather than provoked: none of the three
  # methods flags anything on the demo at its default settings, which
  # is exactly why this survived to production. Driving run_qc() until
  # it flagged something would make the test about the outlier
  # detectors instead of about the card, and would go quiet again the
  # day their thresholds moved.
  shiny::testServer(qc_view_server, args = list(), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "pca")
    b <- last_bundle()
    expect_false(is.null(b))
    b$results$qc_summary$recommended_filters$remove_samples <- "S03"
    last_bundle(b)
    # Outputs are memoised per flush cycle, so without this the card
    # read back is the one computed before the injection.
    session$flushReact()
    html <- paste(unlist(output$stats), collapse = " ")
    expect_true(grepl("accent-warn", html, fixed = TRUE))
    expect_true(grepl("1 flagged", html, fixed = TRUE))
  })
})

test_that("an import that is mostly missing still renders its summary", {
  # mod_import_view.R asks for "warn" past 50% missing. Nothing in the
  # fixtures is that sparse, so this path had never been drawn.
  expect_silent(stat_card("Missing", "72%", trend = "of all cells",
                          accent = "warn"))
})
