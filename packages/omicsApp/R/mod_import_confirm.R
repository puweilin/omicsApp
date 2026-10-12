# The "What will be imported" card of the import view.
#
# Split out of mod_import_view.R and called from inside its module
# server, so the ids stay "import-confirm_shape", "import-design_group".
# The parse state and do_parse() are the module's; this returns
# apply_design(), which commit uses, and normalizable(), which says
# whether the normalization controls apply.
import_confirm_server <- function(input, output, session, parsed, parse_gen,
                                  orientation_override, role_overrides, do_parse,
                                  parse_ok) {
  ns <- session$ns

  # ---- confirmation card --------------------------------------------
  output$confirm_shape <- shiny::renderUI({
    cand <- parsed()
    if (is.null(cand)) return(NULL)
    htmltools::tagList(
      confirm_shape_ui(cand$input, cand$report),
      orientation_picker_ui(ns, cand$report)
    )
  })

  # Samples as columns or rows. Asked outright rather than left to the
  # role dropdowns: a transposed matrix imports cleanly and analyses
  # features as samples, and the only sign is a count the user has to
  # notice.
  shiny::observeEvent(input$orientation_pick, {
    cand <- parsed()
    shiny::req(cand)
    current <- cand$report$suggested_input$orientation
    if (identical(input$orientation_pick, current)) return()
    orientation_override(input$orientation_pick)
    do_parse()
  }, ignoreInit = TRUE)

  output$confirm_roles <- shiny::renderUI({
    cand <- parsed()
    if (is.null(cand)) return(NULL)
    confirm_roles_ui(ns, cand$report, gen = parse_gen())
  })

  # One observer per sheet row, created after the report is known. The
  # dropdowns are rendered by confirm_roles_ui(), so the number of them is
  # data-dependent; observers are registered once and read whatever exists.
  shiny::observe({
    cand <- parsed()
    shiny::req(cand, cand$report$sheets)
    sheets <- cand$report$sheets

    gen <- parse_gen()
    chosen <- vapply(seq_len(nrow(sheets)), function(i) {
      val <- input[[role_input_id(gen, i)]]
      if (is.null(val)) NA_character_ else val
    }, character(1))
    if (all(is.na(chosen))) return()

    current <- sheets$role
    changed <- !is.na(chosen) & chosen != current
    if (!any(changed)) return()

    # Carry every explicit choice, not just the changed one: re-parsing
    # rebuilds the table from the classifier, and an earlier override would
    # otherwise be undone by the next one.
    overrides <- stats::setNames(chosen[!is.na(chosen)],
                                 sheets$name[!is.na(chosen)])
    role_overrides(overrides)
    do_parse()
  })

  # Which column the gene symbol came from, or that none was found.
  #
  # The import matches it by heading, and a heading can be wrong or
  # absent. Getting it wrong relabels every feature; getting nothing
  # leaves feature_symbol as the accession, and enrichment then
  # returns empty rather than failing -- clusterProfiler answers "No
  # gene can be mapped" and hands back NULL, which arrives as a result
  # with nothing in it. Either way the user cannot tell from the
  # outcome, so it is said here, before anything is run.
  output$confirm_symbol_source <- shiny::renderUI({
    cand <- parsed()
    if (is.null(cand) || is.null(cand$input)) return(NULL)
    src <- attr(cand$input$feature_df, "symbol_column") %||% NA_character_
    syms <- cand$input$feature_df$feature_symbol
    share <- if (length(syms)) mean(is_gene_symbol(syms)) else 0
    # Judged on the values, not on whether a column was named like a
    # gene column: ids that are gene symbols, or Ensembl ids mapped to
    # symbols, were told enrichment "will return nothing".
    if (is.na(src) && share >= 0.5) {
      example <- utils::head(syms[is_gene_symbol(syms)], 3L)
      return(htmltools::tags$div(
        class = "muted", style = "font-size:12.5px;margin:4px 0 12px",
        htmltools::tags$strong("Gene symbol"), htmltools::HTML(" &middot; "),
        sprintf("%.0f%% of features carry a gene symbol ", 100 * share),
        htmltools::tags$span(class = "text-mono", paste(example, collapse = ", "))))
    }

    if (is.na(src)) {
      return(htmltools::tags$div(
        class = "muted",
        style = "font-size:12.5px;margin:4px 0 12px",
        htmltools::tags$strong("Gene symbol"),
        htmltools::HTML(" &middot; "),
        "no gene column found \u2014 feature IDs will be used instead. ",
        htmltools::tags$span(
          style = "color:var(--warn)",
          "Enrichment needs gene symbols, so it will return nothing."
        )
      ))
    }
    example <- utils::head(
      stats::na.omit(cand$input$feature_df$feature_symbol), 3L)
    htmltools::tags$div(
      class = "muted",
      style = "font-size:12.5px;margin:4px 0 12px",
      htmltools::tags$strong("Gene symbol"),
      htmltools::HTML(" &middot; "),
      "taken from ",
      htmltools::tags$span(class = "text-mono", sprintf("\u201c%s\u201d", src)),
      if (length(example)) {
        htmltools::tagList(
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(class = "text-mono",
                               paste(example, collapse = ", "))
        )
      }
    )
  })

  # Integer values with gene ids imported as proteomics: almost always
  # RNA-seq counts on the wrong radio button, which then go through vsn
  # and limma as intensities.
  output$omics_hint <- shiny::renderUI({
    cand <- parsed()
    if (is.null(cand) || is.null(cand$input)) return(NULL)
    if (!identical(input$omics_type %||% "proteomics", "proteomics")) return(NULL)
    m <- cand$input$expr_mat
    v <- m[!is.na(m)]
    v <- v[seq_len(min(length(v), 20000L))]
    ids <- utils::head(rownames(m), 200L)
    if (length(v) && all(v >= 0) && all(v == round(v)) && max(v) > 100 &&
        mean(grepl("^ENS[A-Z]*G[0-9]+", ids)) > 0.5) {
      notice(title = "This looks like RNA-seq counts",
             detail = "Whole-number values on Ensembl gene ids. If these are read counts, choose RNA-seq above so they are modelled as counts.",
             kind = "warn")
    }
  })

  output$confirm_matrix_preview <- shiny::renderTable({
    cand <- parsed()
    shiny::req(cand, cand$input)
    preview_matrix(cand$input$expr_mat)
  }, striped = TRUE, spacing = "xs", width = "100%", digits = 2)

  # ---- study design ---------------------------------------------------
  output$confirm_design <- shiny::renderUI({
    cand <- parsed()
    shiny::req(cand, cand$input)
    meta <- cand$input$meta_df
    cands <- grouping_candidates(meta)
    if (!length(cands)) {
      info_cols <- setdiff(names(meta), "sample_id")
      msg <- if (!length(info_cols)) {
        paste("This file has no sample information, so there are no groups to",
              "compare. Add a sample sheet above (one row per sample, a column",
              "naming its group), or put it in a second sheet of the workbook.")
      } else {
        counts <- lapply(meta[info_cols], function(x) table(as.character(x)))
        single <- unique(unlist(lapply(counts, function(t) names(t)[t < 2L])))
        paste0("No column splits the samples into groups of two or more",
               if (length(single)) sprintf(" (only one sample in: %s)",
                                           paste(utils::head(single, 5L), collapse = ", ")),
               ". A group needs at least two samples to be compared.")
      }
      return(htmltools::tags$div(
        class = "notice notice-warn", style = "font-size:12.5px;margin:6px 0 10px",
        msg))
    }
    sel <- shiny::isolate(input$design_group)
    if (is.null(sel) || !sel %in% c(cands, "")) sel <- cands[[1L]]
    htmltools::tags$div(
      class = "design-picker",
      htmltools::tags$h5("Study design",
                         info_tip(paste("Which column holds the groups, and which group",
                                        "is the control. QC colours by it and the",
                                        "Differential view compares against the control",
                                        "by default; both can still be changed there."))),
      htmltools::tags$div(
        class = "row-grid r-6-6",
        shiny::selectInput(ns("design_group"), "Group column",
                           choices = c(cands, "(none)" = ""), selected = sel),
        shiny::uiOutput(ns("design_reference_ui"))
      )
    )
  })

  output$design_reference_ui <- shiny::renderUI({
    cand <- parsed()
    gc <- input$design_group
    shiny::req(cand, cand$input, gc, nzchar(gc), gc %in% names(cand$input$meta_df))
    lv <- sort(unique(as.character(stats::na.omit(cand$input$meta_df[[gc]]))))
    shiny::selectInput(ns("design_reference"), "Control (reference) group",
                       choices = lv, selected = default_control_level(lv))
  })

  # The design as chosen, checked against the metadata it names.
  apply_design <- function(inp) {
    gc <- input$design_group
    if (is.null(gc) || !nzchar(gc) || !gc %in% names(inp$meta_df)) return(inp)
    ref <- input$design_reference
    lv <- unique(as.character(stats::na.omit(inp$meta_df[[gc]])))
    if (is.null(ref) || !ref %in% lv) ref <- default_control_level(sort(lv))
    tryCatch(omicsCore::set_study_design(inp, gc, ref), error = function(e) inp)
  }

  output$confirm_meta_preview <- shiny::renderTable({
    cand <- parsed()
    shiny::req(cand, cand$input)
    preview_metadata(cand$input$meta_df)
  }, striped = TRUE, spacing = "xs", width = "100%")

  # ---- normalization, applied on commit -----------------------------
  # This belongs to import rather than QC because it is what turns a file
  # into something the analysis backends can read: limma applies no
  # transform of its own, so an un-normalized layer means limma runs on raw
  # instrument output. The legacy framework normalized in its
  # data-input layer for the same reason; that layer is what did not survive
  # the port into this package.
  #
  # RNA-seq is deliberately excluded: DESeq2 and edgeR model raw counts
  # directly, and the t-test / lm backends log-transform "raw_count"
  # themselves.
  normalizable <- shiny::reactive({
    if (!parse_ok()) return(FALSE)
    chosen <- input$assay_type %||% parsed()$input$assay_type
    identical(parsed()$input$omics_type, "proteomics") &&
      !is.null(chosen) &&
      !chosen %in% omicsCore::LOG_SCALE_ASSAY_TYPES
  })

  output$normalize_controls <- shiny::renderUI({
    if (!parse_ok()) return(NULL)
    if (!normalizable()) {
      if (!identical(parsed()$input$omics_type, "proteomics")) return(NULL)
      return(htmltools::tags$div(
        class = "muted",
        style = "font-size:12px;margin-bottom:10px",
        "Already on a transformed scale \u2014 nothing to normalize."
      ))
    }
    htmltools::tagList(
      shiny::checkboxInput(
        ns("normalize"),
        label = "Normalize on import",
        value = TRUE
      ),
      shiny::conditionalPanel(
        condition = "input.normalize",
        ns = ns,
        shiny::selectInput(
          ns("normalize_method"),
          label = "Method",
          choices = NORMALIZE_CHOICES,
          selected = "vsn"
        )
      )
    )
  })

  list(apply_design = apply_design, normalizable = normalizable)
}

# ---- confirmation card -----------------------------------------------

# The schema card answers "what did the classifier decide". This one answers
# "what does that decision mean for my data", which is the question a user can
# actually check. A sheet assignment can be wrong at high confidence and still
# yield an omics_input that analyses cleanly -- metadata read as the matrix
# gives numbers, dimensions, and a full result table, all meaningless. Nothing
# downstream errors on it, so this is the last point where it is catchable.
import_confirm_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "What will be imported"),
      htmltools::tags$span(
        class = "card-sub",
        "check this against what you know about the file"
      )
    ),
    bslib::card_body(
      shiny::uiOutput(ns("confirm_shape")),
      shiny::uiOutput(ns("confirm_roles")),
      shiny::uiOutput(ns("confirm_symbol_source")),
      # The groups and the control, stated once here by the person who
      # knows the study, and read by QC, Differential and Integration.
      shiny::uiOutput(ns("confirm_design")),
      htmltools::tags$div(
        class = "row-grid r-6-6",
        htmltools::tags$div(
          htmltools::tags$h5("Expression matrix"),
          htmltools::tags$div(class = "muted",
                              style = "font-size:12px;margin-bottom:6px",
                              "first rows and columns, as parsed"),
          shiny::tableOutput(ns("confirm_matrix_preview"))
        ),
        htmltools::tags$div(
          htmltools::tags$h5("Sample metadata"),
          htmltools::tags$div(class = "muted",
                              style = "font-size:12px;margin-bottom:6px",
                              "columns available for grouping and covariates"),
          shiny::tableOutput(ns("confirm_meta_preview"))
        )
      )
    )
  )
}

# Numbers first, because "3 features x 240 samples" on a file the user knows
# has 240 features is the fastest way to catch a transposed matrix.
confirm_shape_ui <- function(input_obj, report) {
  if (is.null(input_obj)) {
    return(notice(
      title  = "Nothing to import yet",
      detail = "Upload a file, or correct the sheet roles above if the classifier could not find an expression matrix.",
      kind   = "warn"
    ))
  }
  mat <- input_obj$expr_mat
  n_missing <- sum(is.na(mat))
  orientation <- report$suggested_input$orientation %||% "features_in_rows"

  htmltools::tags$div(
    class = "stat-grid",
    style = "margin-bottom:16px",
    stat_card(
      label = "Features", value = format(nrow(mat), big.mark = ","),
      trend = "rows of the matrix", mono = TRUE
    ),
    stat_card(
      label = "Samples", value = format(ncol(mat), big.mark = ","),
      trend = "columns of the matrix", mono = TRUE
    ),
    stat_card(
      label = "Missing", value = sprintf("%.1f%%", 100 * n_missing / length(mat)),
      trend = sprintf("%s cells", format(n_missing, big.mark = ",")),
      accent = if (n_missing / length(mat) > 0.5) "warn" else "ok"
    ),
    stat_card(
      label = "Orientation",
      value = if (identical(orientation, "features_in_rows")) "features in rows"
              else "samples in rows",
      trend = if (identical(report$suggested_input$orientation_source, "user")) "as you set it"
              else sprintf("detected (confidence %.2f)",
                           report$suggested_input$orientation_confidence %||% NA_real_),
      accent = if ((report$suggested_input$orientation_confidence %||% 1) < 0.6) "warn" else "ok"
    )
  )
}

# The control that corrects the orientation, prominent when it was a guess.
orientation_picker_ui <- function(ns, report) {
  sug <- report$suggested_input
  if (is.null(sug$orientation) || identical(sug$matrix_sheet, "rds")) return(NULL)
  guessed <- (sug$orientation_confidence %||% 1) < 0.6 &&
    !identical(sug$orientation_source, "user")
  htmltools::tags$div(
    style = "margin:-6px 0 14px",
    if (guessed) {
      notice(title = "Check the orientation",
             detail = paste("The file did not make it clear whether samples are",
                            "columns or rows. If the Samples count above is really",
                            "the number of features, switch it here."),
             kind = "warn")
    },
    shiny::radioButtons(
      ns("orientation_pick"), label = "Samples are the matrix's",
      choices = c("columns (features in rows)" = "features_in_rows",
                  "rows (samples in rows, e.g. Olink NPX)" = "samples_in_rows"),
      selected = sug$orientation, inline = TRUE)
  )
}

# Which sheet became what, with a dropdown to say otherwise. Sheets the
# classifier could not place are worth showing too: an "unknown" sheet is
# often the metadata, and silently dropping it is how a grouping column goes
# missing later.
# Input id for one sheet's role dropdown. The generation is part of the id so
# a new upload gets fresh controls rather than inheriting the last file's.
role_input_id <- function(gen, i) paste0("role_", gen, "_", i)

confirm_roles_ui <- function(ns, report, gen = 0L) {
  sheets <- report$sheets
  if (is.null(sheets) || nrow(sheets) == 0L) return(NULL)

  choices <- c("expression matrix" = "matrix",
               "sample metadata" = "metadata",
               "feature annotation" = "feature_annot",
               "ignore" = "unknown")

  rows <- lapply(seq_len(nrow(sheets)), function(i) {
    nm <- sheets$name[i]
    low_conf <- !is.na(sheets$confidence[i]) && sheets$confidence[i] < 0.5
    htmltools::tags$div(
      style = "display:flex;align-items:center;gap:10px;margin-bottom:6px",
      htmltools::tags$code(style = "min-width:150px", nm),
      htmltools::tags$span(
        class = "muted", style = "font-size:12px;min-width:110px",
        sprintf("%s x %s", sheets$n_rows[i], sheets$n_cols[i])
      ),
      shiny::selectInput(
        ns(role_input_id(gen, i)),
        label = htmltools::tags$span(class = "visually-hidden-label", sprintf("Role of sheet %s", nm)),
        choices = choices,
        selected = sheets$role[i], width = "180px"
      ),
      if (low_conf) {
        pill("low confidence", kind = "warn")
      } else if (identical(sheets$notes[i], "role set by user")) {
        pill("you set this", kind = "ok")
      } else NULL
    )
  })

  htmltools::tagList(
    htmltools::tags$h5("Sheet roles"),
    htmltools::tags$div(class = "muted",
                        style = "font-size:12px;margin-bottom:8px",
                        "Changing a role re-reads the file with your assignment."),
    rows,
    htmltools::tags$hr(style = "margin:14px 0")
  )
}

# A corner of the matrix. Seeing the actual values is what catches a header
# row parsed as data, or an ID column read as a sample.
preview_matrix <- function(mat, n_row = 5L, n_col = 4L) {
  if (is.null(mat) || nrow(mat) == 0L) return(NULL)
  sub <- mat[seq_len(min(n_row, nrow(mat))),
             seq_len(min(n_col, ncol(mat))), drop = FALSE]
  df <- as.data.frame(round(sub, 2))
  df <- cbind(feature = rownames(sub), df)
  rownames(df) <- NULL
  df
}

preview_metadata <- function(meta, n_row = 5L) {
  if (is.null(meta) || nrow(meta) == 0L) return(NULL)
  df <- meta[seq_len(min(n_row, nrow(meta))), , drop = FALSE]
  cbind(sample = rownames(df), as.data.frame(df, stringsAsFactors = FALSE))
}
