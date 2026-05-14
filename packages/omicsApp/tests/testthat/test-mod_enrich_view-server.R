# testServer harness for the Enrichment view (slice 3E).
#
# Three scenarios:
#   1. NULL diff_bundle → demo fixture fills the dot/hits cards
#      and a notice nudges the user to run a diff.
#   2. Live diff_bundle but clusterProfiler not installed →
#      `enrich_error()` carries the install hint; demo fallback.
#   3. Live diff_bundle + clusterProfiler installed →
#      `enrich_bundle()` is an analysis_bundle from run_enrichment.

test_that("enrich view falls back to demo when diff_bundle is NULL", {
  diff_bundle <- shiny::reactiveVal(NULL)
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = diff_bundle),
    {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both")
      expect_true(isTRUE(is_demo()))
      expect_null(enrich_bundle())
      # table_data() should fall through to example_enrich_table.
      expect_true(nrow(table_data()) > 0L)
    }
  )
})

test_that("enrich view surfaces clusterProfiler install hint when needed", {
  skip_if(requireNamespace("clusterProfiler", quietly = TRUE),
          "clusterProfiler is installed; cannot exercise the gating notice.")
  diff_bundle <- shiny::reactiveVal(omicsApp:::example_diff_bundle())
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = diff_bundle),
    {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both", rerun = 1)
      expect_true(isTRUE(is_demo()))
      expect_match(enrich_error(), "clusterProfiler", fixed = TRUE)
    }
  )
})

test_that("enrich view runs run_enrichment when clusterProfiler is installed", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  diff_bundle <- shiny::reactiveVal(omicsApp:::example_diff_bundle())
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = diff_bundle),
    {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both", rerun = 1)
      b <- enrich_bundle()
      expect_s3_class(b, "analysis_bundle")
      expect_identical(b$analysis_name, "run_enrichment")
      expect_false(isTRUE(is_demo()))
    }
  )
})
