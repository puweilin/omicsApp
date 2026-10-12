# Committing an import: the layer name, the replace check, normalization
# and archiving the upload.
#
# Split out of mod_import_view.R and called from inside its module
# server, so the ids stay "import-confirm", "import-layer_name". The
# confirmed input it writes is the module's, and the module returns it.
import_commit_server <- function(input, output, session, current_project, parsed, parse_ok,
                                 normalizable, apply_design, confirmed_input) {
  ns <- session$ns

  # ---- confirm, guarded when it would replace live data -------------
  # Committing an input whose omics_type already exists in the project
  # replaces that layer, and every analysis computed on it becomes
  # meaningless. Three cases, only the last of which is destructive:
  #
  #   no such layer yet     -> commit straight away
  #   same file re-selected -> nothing changed, keep the results
  #   different data        -> ask first, then commit and clear
  #
  # The middle case is why the fingerprint is worth carrying: a
  # mis-click that re-picks the same file must not cost the user an
  # afternoon of analysis.
  layer_being_replaced <- function(cand) {
    proj <- current_project()
    if (is.null(proj) || is.null(cand)) return(NULL)
    tag <- target_tag(cand)
    if (!tag %in% names(proj$experiments)) return(NULL)
    proj$experiments[[tag]]
  }

  # The name the layer goes in under. A project can hold several layers
  # of one omics type (two proteomics batches), so the name is the
  # user's: it defaults to the omics type, and a name already in the
  # project means "replace that layer".
  target_tag <- function(cand) {
    clean_layer_name(input$layer_name, cand$omics_type %||% "experiment")
  }
  output$layer_name_ui <- shiny::renderUI({
    shiny::req(parse_ok())
    default <- parsed()$input$omics_type %||% "experiment"
    htmltools::tags$div(
      class = "inline-control", style = "max-width:200px",
      shiny::textInput(ns("layer_name"), label = "Layer name",
                       value = shiny::isolate(input$layer_name) %||% default,
                       placeholder = default)
    )
  })
  # The omics type is the default name, and changes with the type
  # control until the user writes a name of their own.
  shiny::observeEvent(parsed()$input$omics_type, {
    typ <- parsed()$input$omics_type
    cur <- input$layer_name
    if (is.null(cur) || !nzchar(trimws(cur)) || cur %in% c("proteomics", "rnaseq")) {
      shiny::updateTextInput(session, "layer_name", value = typ)
    }
  }, ignoreInit = TRUE)

  # Archive the upload alongside the parsed input. Only on commit:
  # parsing happens on every file pick and radio change, most of which
  # the user never confirms. Archiving failing must not stop the
  # import, so its outcome is surfaced but not acted on.
  # Which normalization the current controls would apply. "none" when the
  # layer is not normalizable or the box is unticked -- both have to be
  # distinguishable in the fingerprint from an actual method.
  pending_normalize_method <- function() {
    if (normalizable() && isTRUE(input$normalize %||% TRUE)) {
      input$normalize_method %||% "vsn"
    } else {
      "none"
    }
  }

  # The identity the data would have once committed. Computed before the
  # replace check as well as inside commit(), so that re-picking the same
  # file with the same settings still reads as "nothing changed" rather than
  # as new data.
  stamp_fingerprint <- function(cand) {
    f <- input$file
    if (is.null(f)) return(cand)
    cand$source_fingerprint <- input_fingerprint(
      f$datapath, cand$omics_type, cand$assay_type,
      normalize = pending_normalize_method())
    # A different sample sheet is different data.
    sf <- input$sample_file
    if (!is.null(sf)) {
      cand$source_fingerprint <- paste0(cand$source_fingerprint, ":sheet=",
                                        unname(tools::md5sum(sf$datapath)))
    }
    cand
  }

  commit <- function(cand, tag = target_tag(cand)) {
    method <- pending_normalize_method()
    do_normalize <- !identical(method, "none")
    cand <- stamp_fingerprint(cand)
    f <- input$file

    if (do_normalize) {
      normalized <- tryCatch(
        suppressMessages(do.call(omicsCore::normalize_omics,
                                 c(list(cand), normalize_args(method)))),
        error = function(e) e
      )
      if (inherits(normalized, "error")) {
        # Importing raw and telling the user beats importing something the
        # backends will silently mistreat
        shiny::showNotification(
          paste0("Normalization failed, importing unnormalized: ",
                 conditionMessage(normalized)),
          type = "error", duration = 12
        )
      } else {
        # normalize_omics() keeps the pre-normalization matrix in raw_mat
        normalized$source_fingerprint <- cand$source_fingerprint
        cand <- normalized
        shiny::showNotification(
          sprintf("Normalized with %s; values are now '%s'.",
                  names(NORMALIZE_CHOICES)[NORMALIZE_CHOICES == method],
                  cand$assay_type),
          type = "message", duration = 6
        )
      }
    }

    cand <- apply_design(cand)

    # Quantification files are archived one by one, with the names
    # they were uploaded under: export_script() reads them back with
    # read_quant_files(), and the names are the sample names (the
    # archived copies are renamed to stay unique).
    is_quant <- !is.null(parsed()$report$suggested_input$quant_format)
    if (!is.null(f) && is_quant) {
      res <- lapply(seq_along(f$datapath), function(i) {
        store_raw_upload(f$datapath[[i]], f$name[[i]], unname(tools::md5sum(f$datapath[[i]])))
      })
      if (all(vapply(res, function(r) isTRUE(r$ok), logical(1)))) {
        cand$quant_source <- list(paths = vapply(res, `[[`, character(1), "path"),
                                  names = as.character(f$name))
      } else if (any(grepl("quota", vapply(res, function(r) r$message %||% "", character(1)),
                           fixed = TRUE))) {
        shiny::showNotification(res[[1L]]$message, type = "warning", duration = 8)
      }
    } else if (!is.null(f) && length(f$datapath) == 1L) {
      res <- store_raw_upload(f$datapath, f$name, cand$source_fingerprint)
      if (isTRUE(res$ok)) {
        # Recorded so `omicsCore::export_script()` can point its
        # read_omics() line at a file that exists, which is the
        # difference between a script that runs and one that only
        # documents what was run.
        cand$source_path <- res$path
      } else if (grepl("quota", res$message, fixed = TRUE)) {
        shiny::showNotification(res$message, type = "warning", duration = 8)
      }
    }
    sf <- input$sample_file
    if (!is.null(sf) && !is.null(cand$sample_sheet_path)) {
      res <- store_raw_upload(sf$datapath, sf$name,
                              unname(tools::md5sum(sf$datapath)))
      cand$sample_sheet_path <- if (isTRUE(res$ok)) res$path else NULL
    }
    cand$layer_tag <- tag
    confirmed_input(cand)
  }

  shiny::observeEvent(input$confirm, {
    shiny::req(parse_ok())
    cand <- stamp_fingerprint(parsed()$input)
    existing <- layer_being_replaced(cand)

    if (is.null(existing)) {
      commit(cand)
      return()
    }
    if (fingerprints_match(existing, cand)) {
      shiny::showNotification(
        "That is the file already loaded \u2014 nothing to re-import.",
        type = "message"
      )
      return()
    }
    tag <- target_tag(cand)
    shiny::showModal(
      replace_layer_modal(ns, tag, current_project(),
                          keep_both = next_free_tag(tag, names(current_project()$experiments)))
    )
  })

  shiny::observeEvent(input$confirm_replace, {
    shiny::removeModal()
    shiny::req(parse_ok())
    commit(parsed()$input)
  })
  # Both: the new data as a layer of its own beside the old one, whose
  # results stay.
  shiny::observeEvent(input$confirm_keep_both, {
    shiny::removeModal()
    shiny::req(parse_ok())
    cand <- parsed()$input
    tag <- next_free_tag(target_tag(cand), names(current_project()$experiments))
    shiny::updateTextInput(session, "layer_name", value = tag)
    commit(cand, tag = tag)
  })
  invisible()
}

# Human-readable names for the bundles about to be discarded, so the
# dialog names what is at stake rather than saying "analyses".
BUNDLE_LABELS <- c(
  qc          = "Quality control",
  diff        = "Differential analysis",
  enrich      = "Pathway enrichment",
  integration = "Multi-omics integration"
)

# A layer name as typed, made safe to file under: letters, digits and
# `_ - .` (a space becomes `_`), at most 40 characters, the omics type
# when nothing usable is left.
clean_layer_name <- function(x, default) {
  x <- trimws(as.character(x %||% ""))
  if (!length(x) || !nzchar(x[[1L]])) return(default)
  x <- gsub("[^[:alnum:]_.-]+", "_", x[[1L]], perl = TRUE)
  x <- gsub("^_+|_+$", "", substr(x, 1L, 40L))
  if (nzchar(x)) x else default
}

# `tag`, or `tag_2`, `tag_3`, ... -- the first name no layer has.
next_free_tag <- function(tag, taken) {
  if (!tag %in% taken) return(tag)
  i <- 2L
  while (paste0(tag, "_", i) %in% taken) i <- i + 1L
  paste0(tag, "_", i)
}

replace_layer_modal <- function(ns, tag, project, keep_both = NULL) {
  # The results on this layer, not every result in the project: another
  # layer's analyses are not touched by replacing this one.
  all_b <- project$bundles %||% list()
  bundles <- names(all_b)[vapply(all_b, function(b) {
    !omicsCore::is_analysis_bundle(b) ||
      tag %in% (b$params$experiments %||% character(0)) ||
      identical(omicsCore::bundle_layer(project, b), tag)
  }, logical(1))]
  losing <- if (length(bundles) == 0L) {
    htmltools::tags$p(
      class = "muted",
      "No analyses have been run on the current data yet."
    )
  } else {
    htmltools::tagList(
      htmltools::tags$p("These results were computed on the current data ",
                        "and will be cleared:"),
      htmltools::tags$ul(
        lapply(bundles, function(b) {
          htmltools::tags$li(unname(BUNDLE_LABELS[b]) %||% b)
        })
      )
    )
  }
  shiny::modalDialog(
    title = sprintf("Replace the %s layer?", tag),
    losing,
    htmltools::tags$p(
      class = "muted",
      style = "font-size:12px",
      "The uploaded file differs from the one currently loaded. Keeping ",
      "results computed on the previous data would misreport them as ",
      "belonging to the new data."
    ),
    easyClose = FALSE,
    footer = htmltools::tagList(
      shiny::modalButton("Cancel"),
      if (!is.null(keep_both)) {
        shiny::actionButton(ns("confirm_keep_both"),
                            sprintf("Keep both (add as '%s')", keep_both),
                            class = "btn btn-outline-primary")
      },
      shiny::actionButton(ns("confirm_replace"), "Replace and clear",
                          class = "btn btn-danger")
    )
  )
}
