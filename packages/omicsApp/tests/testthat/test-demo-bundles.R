# The demo views draw through the same `omicsCore::plot_*()` dispatchers
# as the live ones. That only holds while the fixtures keep satisfying
# the result schemas -- and a fixture that drifts out of schema fails at
# render time in front of a user, not here, unless something checks.

test_that("the demo diff fixture plots through the dispatcher", {
  b <- example_diff_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_s3_class(omicsCore::plot_volcano(b), "ggplot")
})

test_that("the demo enrichment fixture is a bundle that satisfies the schema", {
  b <- example_enrich_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_enrichment")
  # plot_enrichment() runs the result through check_enrich_result_schema(),
  # so this fails loudly if the fixture loses a required column.
  expect_s3_class(omicsCore::plot_enrichment(b, view = "dot", top_n = 12L),
                  "ggplot")
})

test_that("the demo integration fixture is a bundle that satisfies the schema", {
  b <- example_integration_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_integration")
  expect_s3_class(omicsCore::plot_integration(b, view = "dual_volcano"),
                  "ggplot")
  expect_s3_class(omicsCore::plot_integration(b, view = "effect_pair"),
                  "ggplot")
})

test_that("demo bundles are cached rather than rebuilt", {
  # Each demo view resolves its own bundle; rebuilding would reseed the
  # fixtures and leave two views showing different synthetic data, which
  # is the bug the fixture cache was introduced for.
  expect_identical(example_enrich_bundle(), example_enrich_bundle())
  expect_identical(example_integration_bundle(), example_integration_bundle())
})

test_that("the app palette is the one omicsCore draws with", {
  # Two palettes means a figure on screen and the same figure in an
  # exported report are tinted differently.
  expect_identical(omics_colors, omicsCore::omics_colors)
})
