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

# ---- full-table download ---------------------------------------------
# The CSV is what gets opened in Excel a week later. It carries every
# pathway and every column -- including the gene lists, which are what a
# follow-up needs and the only thing the on-screen table cannot show.

test_that("the enrichment download writes the whole table, not the filtered view", {
  diff_bundle <- shiny::reactiveVal(NULL)
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = diff_bundle),
    {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both",
                        # Strict enough to empty the table on screen.
                        # The file must not follow it: a CSV silently
                        # truncated to one threshold is the kind of thing
                        # nobody notices until a number cannot be
                        # reproduced.
                        show_p = "adjusted", show_cutoff = 1e-12)

      full <- table_data()
      expect_true(nrow(full) > 0L)

      written <- utils::read.csv(output$download_table, check.names = FALSE)

      expect_equal(nrow(written), nrow(full))
      expect_equal(ncol(written), ncol(full))
      expect_true(all(c("pathway_name", "p_value", "adj_p_value") %in%
                        names(written)))
      # The gene lists are the reason to export at all.
      expect_true(any(c("overlap_features", "leading_features") %in%
                        names(written)))
    }
  )
})

test_that("the exported table is sorted so the file opens on the strongest result", {
  diff_bundle <- shiny::reactiveVal(NULL)
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = diff_bundle),
    {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both")
      written <- utils::read.csv(output$download_table, check.names = FALSE)
      by_db <- split(written$adj_p_value, written$database)
      for (v in by_db) {
        expect_false(is.unsorted(v[!is.na(v)]))
      }
    }
  )
})
