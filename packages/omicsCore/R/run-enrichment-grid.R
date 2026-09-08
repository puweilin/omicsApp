# Every enrichment a differential supports, in one pass.
#
# The Enrich view runs one combination at a time, which is what you want
# while you are looking at a result. It is not what you want when you are
# finished: getting ORA and GSEA across seven databases out of the
# application meant fourteen visits to the same view, each with a download.
#
# The grid runs them all and stacks the rows. They already carry
# `database`, `direction` and `result_type`, so a stacked table loses
# nothing -- it is the same rows with the combination named in columns
# rather than in a filename.

#' Run every enrichment a differential supports
#'
#' ORA and GSEA across the requested databases, returned as one bundle
#' whose `enrich_result_df` holds every row.
#'
#' GSEA is run once per database rather than once per direction. Its
#' `direction` argument filters the standardized output and does not
#' change the analysis -- it is a ranked-list method that reads every
#' feature -- so three runs would be the same computation filtered three
#' ways, and the rows say which way each one went. ORA is different: the
#' direction decides which genes go into the test, so `up`, `down` and
#' `both` are three analyses and all three are run.
#'
#' A combination that fails does not stop the others. Missing gene sets,
#' an organism a database does not cover, a differential with nothing
#' significant in one direction -- each is one empty row in the log and
#' the rest of the grid still runs. `results$grid_log` is that log: one
#' row per combination, with what happened.
#'
#' @param diff_bundle An [`analysis_bundle`][is_analysis_bundle()] from
#'   [run_diff()].
#' @param types Which tests to run.
#' @param databases Which databases, defaulting to all supported.
#' @param ora_directions Directions to run ORA in.
#' @param on_progress Optional `function(done, total, label)`, called
#'   before each combination. The application passes a progress bar; the
#'   grid takes tens of seconds and a caller that cannot say where it has
#'   got to looks like one that has stopped.
#' @param ... Passed to [run_enrichment()].
#'
#' @return An `analysis_bundle` with `results$enrich_result_df` (every
#'   row, stacked) and `results$grid_log`.
#' @export
#' @family enrichment
#' @examples
#' \dontrun{
#'   grid <- run_enrichment_grid(diff_bundle)
#'   table(grid$results$enrich_result_df$database)
#' }
run_enrichment_grid <- function(
  diff_bundle,
  types = c("ora", "gsea"),
  databases = SUPPORTED_ENRICH_DATABASES,
  ora_directions = c("both", "up", "down"),
  on_progress = NULL,
  ...
) {
  if (!is_analysis_bundle(diff_bundle) ||
      !identical(diff_bundle$analysis_name, "run_diff")) {
    stop("`diff_bundle` must be an analysis_bundle from run_diff().")
  }
  types <- match.arg(types, c("ora", "gsea"), several.ok = TRUE)
  databases <- match.arg(databases, SUPPORTED_ENRICH_DATABASES,
                         several.ok = TRUE)
  ora_directions <- match.arg(ora_directions, c("both", "up", "down"),
                              several.ok = TRUE)

  plan <- list()
  for (db in databases) {
    if ("ora" %in% types) {
      for (dir_ in ora_directions) {
        plan[[length(plan) + 1L]] <- list(type = "ora", database = db,
                                          direction = dir_)
      }
    }
    # Once, at "both". See above.
    if ("gsea" %in% types) {
      plan[[length(plan) + 1L]] <- list(type = "gsea", database = db,
                                        direction = "both")
    }
  }

  frames <- list()
  log <- vector("list", length(plan))
  for (i in seq_along(plan)) {
    step <- plan[[i]]
    label <- sprintf("%s / %s / %s", toupper(step$type), step$database,
                     step$direction)
    if (is.function(on_progress)) on_progress(i - 1L, length(plan), label)

    t0 <- Sys.time()
    # Not suppressWarnings(). The backends answer a failure with a warning
    # and NULL rather than an error, so suppressing warnings here made a
    # database that had *broken* indistinguishable from one that found
    # nothing: `status = "empty"`, no message, and a reader concluding
    # there is nothing there. The classed warning is collected; ordinary
    # ones are still quietened, because a grid of 28 runs otherwise
    # reports every DESeq2 fit note 28 times.
    backend_failure <- NULL
    res <- tryCatch(
      withCallingHandlers(
        run_enrichment(diff_bundle, type = step$type,
                       database = step$database,
                       direction = step$direction, ...),
        omicsCore_enrich_backend_failure = function(w) {
          backend_failure <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        },
        warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) e)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (!is.null(backend_failure) && !inherits(res, "error")) {
      res <- simpleError(backend_failure)
    }
    if (inherits(res, "error")) {
      log[[i]] <- data.frame(
        type = step$type, database = step$database,
        direction = step$direction, rows = 0L, seconds = round(secs, 2),
        status = "failed", message = conditionMessage(res),
        stringsAsFactors = FALSE)
      next
    }
    df <- res$results$enrich_result_df
    n <- if (is.null(df)) 0L else nrow(df)
    if (n > 0L) frames[[length(frames) + 1L]] <- df
    log[[i]] <- data.frame(
      type = step$type, database = step$database,
      direction = step$direction, rows = n, seconds = round(secs, 2),
      status = if (n > 0L) "ok" else "empty", message = NA_character_,
      stringsAsFactors = FALSE)
  }
  if (is.function(on_progress)) {
    on_progress(length(plan), length(plan), "done")
  }

  grid_log <- do.call(rbind, log)
  # rbind on frames from different databases: the standardized schema is
  # the same for all of them, which is what makes stacking legitimate
  # rather than a coincidence that holds until a backend adds a column.
  stacked <- if (length(frames) == 0L) {
    empty_enrich_result_df()
  } else {
    cols <- Reduce(intersect, lapply(frames, names))
    do.call(rbind, lapply(frames, function(d) d[, cols, drop = FALSE]))
  }

  new_analysis_bundle(
    analysis_name = "run_enrichment_grid",
    input_info = diff_bundle$input_info,
    # The parameters it actually ran with, `...` included. Recording only
    # the grid's own arguments left the thresholds -- the numbers that
    # decide which genes go into ORA -- written down nowhere, so a batch
    # run could not be told from any other batch run, and the script
    # export had nothing to replay.
    params = c(
      list(
        types = types,
        databases = databases,
        ora_directions = ora_directions,
        comparison = diff_bundle$params$comparison,
        n_combinations = length(plan)
      ),
      list(...)
    ),
    results = list(
      enrich_result_df = stacked,
      grid_log = grid_log
    ),
    warnings = {
      failed <- grid_log[grid_log$status == "failed", , drop = FALSE]
      if (nrow(failed) == 0L) character(0) else
        sprintf("%s / %s / %s: %s", failed$type, failed$database,
                failed$direction, failed$message)
    }
  )
}

# An enrichment table with the right columns and no rows, for the case
# where nothing came back at all. Returning NULL instead would make every
# caller test for it.
empty_enrich_result_df <- function() {
  cols <- ENRICH_RESULT_REQUIRED_COLS
  out <- as.data.frame(
    stats::setNames(replicate(length(cols), character(0), simplify = FALSE),
                    cols),
    stringsAsFactors = FALSE)
  out
}
