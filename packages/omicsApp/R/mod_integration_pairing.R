# The sample-pairing and feature-matching cards of the integration view.
#
# Split out of mod_integration_view.R. These are plain functions called
# from inside integration_view_server()'s moduleServer(), not modules, so
# every input and output keeps its id ("integration-accept_pairing",
# "integration-link_file").

# Which sample of one layer is which sample of the other. Returns the
# pairing preview.
integration_pairing_server <- function(input, output, session, current_project, layers) {
  # ---- sample pairing ------------------------------------------------
  # Shown whether or not the current method needs it. Concordance
  # matches features, not samples, so it runs without a pairing -- but
  # a user who gets that far and then cannot correlate should find out
  # here rather than after the run, and should be able to fix it here
  # rather than by re-importing.
  pairing <- shiny::reactive({
    proj <- current_project()
    l <- layers()
    if (is.null(proj) || is.null(l)) return(NULL)
    omicsCore::sample_pairing_preview(proj, l$primary, l$partner)
  })

  output$pairing_table <- DT::renderDT({
    p <- pairing()
    shiny::req(p, nrow(p$pairs) > 0L)
    l <- layers()
    out <- p$pairs
    names(out) <- c("Donor", l$primary, l$partner)
    DT::datatable(out, rownames = FALSE, selection = "none",
                  options = list(pageLength = 8, dom = "tip"))
  }, server = TRUE)

  output$pairing_note <- shiny::renderUI({
    p <- pairing()
    if (is.null(p)) {
      return(notice("Import a second layer to pair samples across omics.",
                    kind = "info"))
    }
    n <- nrow(p$pairs)
    if (n == 0L) {
      return(notice(
        "No pairing found",
        paste0(
          "Sample-level integration needs to know which sample in each ",
          "layer came from the same person. Add a `donor` column to ",
          "each layer's sample metadata and re-import those layers, or ",
          "rename the samples so they share a leading id \u2014 RD001-C ",
          "and RD001_Folli both give RD001. Re-importing a layer clears ",
          "the results computed on it, so do this before the other ",
          "analyses if you can. Fold-change concordance does not need a pairing."),
        kind = "warn"))
    }
    n_donor <- length(unique(p$pairs$donor_id))
    dup <- n > n_donor
    extra <- if (dup) {
      sprintf(" %d donor(s) have more than one sample in a layer; correlation uses one sample per donor.",
              n - n_donor)
    } else ""
    switch(
      p$source,
      linked = notice(sprintf("%d pairs, saved with this project.%s", n, extra),
                      kind = if (dup) "warn" else "info"),
      donor = notice(sprintf("%d pairs, from the donor column in both layers.%s", n, extra),
                     kind = if (dup) "warn" else "info"),
      sample_id = notice(sprintf("%d pairs, from sample ids that match outright.", n),
                         kind = "info"),
      suggested = notice(
        sprintf("%d pairs, guessed from the sample ids \u2014 check them", n),
        paste0("Nothing states that these are the same person; the ids ",
               "merely share a leading part. Accept the pairing to keep ",
               "it with the project (sample-level correlation will not ",
               "use a guess), or state it properly with a donor column."),
        kind = "warn"),
      notice(sprintf("%d pairs.", n), kind = "info")
    )
  })

  # A guess becomes a decision only when someone says so, and then it is
  # saved with the project rather than re-derived every time.
  output$pairing_action <- shiny::renderUI({
    p <- pairing()
    if (is.null(p) || nrow(p$pairs) == 0L) return(NULL)
    if (!identical(p$source, "suggested")) return(NULL)
    shiny::actionButton(session$ns("accept_pairing"), "Accept pairing",
                        class = "btn btn-sm btn-primary")
  })

  shiny::observeEvent(input$accept_pairing, {
    proj <- current_project()
    p <- pairing()
    l <- layers()
    shiny::req(proj, p, l, nrow(p$pairs) > 0L)
    new_rows <- rbind(
      data.frame(tag = l$primary, sample_id = p$pairs$a,
                 donor_id = p$pairs$donor_id, stringsAsFactors = FALSE),
      data.frame(tag = l$partner, sample_id = p$pairs$b,
                 donor_id = p$pairs$donor_id, stringsAsFactors = FALSE)
    )
    # Added to what the project already states, not in place of it: a
    # third layer's saved pairing is not this pair's to throw away.
    old <- proj$sample_link
    if (!is.null(old) && nrow(old)) {
      clash <- paste(old$tag, old$sample_id) %in%
        paste(new_rows$tag, new_rows$sample_id)
      new_rows <- rbind(old[!clash, c("tag", "sample_id", "donor_id"), drop = FALSE],
                        new_rows)
    }
    rownames(new_rows) <- NULL
    proj$sample_link <- new_rows
    current_project(proj)
    shiny::showNotification(
      sprintf("Pairing saved: %d donors.", nrow(p$pairs)),
      type = "message")
  })

  pairing
}

# How the features of the two layers are matched, and the mapping-table
# upload. Returns the reactives the module keeps.
integration_feature_server <- function(input, output, session, current_project, layers) {
  # ---- feature matching ----------------------------------------------
  # How the features of the two layers will be matched: by gene symbol,
  # or by a mapping table the user uploaded (saved with the project as
  # its feature_link, and archived so the exported script reads it).
  feature_pairing <- shiny::reactive({
    proj <- current_project()
    l <- layers()
    if (is.null(proj) || is.null(l)) return(NULL)
    tryCatch(omicsCore::feature_pairing_preview(proj, l$primary, l$partner),
             error = function(e) list(error = conditionMessage(e)))
  })

  output$feature_note <- shiny::renderUI({
    fp <- feature_pairing()
    l <- layers()
    if (is.null(fp) || is.null(l)) {
      return(notice("Import a second layer to match features across omics.",
                    kind = "info"))
    }
    if (!is.null(fp$error)) {
      return(notice("The saved mapping table cannot be used", fp$error, kind = "warn"))
    }
    feature_pairing_notice(fp, l$primary, l$partner)
  })

  # An uploaded table, read but not yet in use: its columns, and the
  # guess of which holds which layer's identifiers.
  link_upload <- shiny::reactive({
    f <- input$link_file
    if (is.null(f)) return(NULL)
    tab <- tryCatch(omicsCore::read_feature_link(f$datapath, columns = NULL),
                    error = function(e) e)
    if (inherits(tab, "error")) return(list(error = conditionMessage(tab)))
    if (ncol(tab) < 2L) {
      return(list(error = "The table needs two columns: one for each layer's identifiers."))
    }
    proj <- current_project()
    l <- layers()
    guess <- if (!is.null(proj) && !is.null(l)) guess_link_columns(tab, proj, l$primary, l$partner)
    list(tab = tab, name = f$name, path = f$datapath, guess = guess %||% names(tab)[1:2])
  })

  output$link_columns <- shiny::renderUI({
    up <- link_upload()
    l <- layers()
    if (is.null(up) || is.null(l)) return(NULL)
    if (!is.null(up$error)) {
      return(notice("The table could not be read", up$error, kind = "warn"))
    }
    cols <- names(up$tab)
    htmltools::tagList(
      htmltools::tags$div(
        class = "row-grid r-6-6",
        shiny::selectInput(session$ns("link_col_a"),
                           sprintf("Column with the %s identifiers", l$primary),
                           choices = cols, selected = up$guess[[1L]]),
        shiny::selectInput(session$ns("link_col_b"),
                           sprintf("Column with the %s identifiers", l$partner),
                           choices = cols, selected = up$guess[[2L]])),
      shiny::uiOutput(session$ns("link_preview")),
      shiny::actionButton(session$ns("use_link"), "Use this table",
                          class = "btn btn-sm btn-primary"))
  })

  # The table as it would be used with the columns chosen.
  pending_link <- shiny::reactive({
    up <- link_upload()
    l <- layers()
    a <- input$link_col_a
    b <- input$link_col_b
    if (is.null(up) || !is.null(up$error) || is.null(l) || is.null(a) || is.null(b) ||
        !all(c(a, b) %in% names(up$tab))) return(NULL)
    if (identical(a, b)) return(list(error = "Choose a different column for each layer."))
    link <- stats::setNames(up$tab[c(a, b)], c(l$primary, l$partner))
    prev <- tryCatch(omicsCore::feature_pairing_preview(current_project(), l$primary,
                                                        l$partner, feature_link = link),
                     error = function(e) list(error = conditionMessage(e)))
    list(link = link, preview = prev)
  })

  output$link_preview <- shiny::renderUI({
    pl <- pending_link()
    l <- layers()
    if (is.null(pl)) return(NULL)
    if (!is.null(pl$error)) return(notice(pl$error, kind = "warn"))
    if (!is.null(pl$preview$error)) return(notice(pl$preview$error, kind = "warn"))
    htmltools::tags$div(
      class = "muted", style = "font-size:12.5px;margin:4px 0 8px",
      paste("With this table:", feature_pairing_sentence(pl$preview, l$primary, l$partner)))
  })

  shiny::observeEvent(input$use_link, {
    up <- link_upload()
    pl <- pending_link()
    proj <- current_project()
    l <- layers()
    shiny::req(up, pl, proj, l, is.null(pl$error), is.null(pl$preview$error))
    # Archived like a sample sheet, so the exported script can read the
    # same file; read back from the archived copy so the link remembers
    # where it lives. Archiving fails soft (quota): the link is then
    # used from the upload and the script asks for the file instead.
    res <- store_raw_upload(up$path, up$name, unname(tools::md5sum(up$path)))
    src <- if (isTRUE(res$ok)) res$path else up$path
    link <- tryCatch(
      omicsCore::read_feature_link(
        src, columns = stats::setNames(c(input$link_col_a, input$link_col_b),
                                       c(l$primary, l$partner))),
      error = function(e) e)
    if (inherits(link, "error")) {
      shiny::showNotification(conditionMessage(link), type = "error")
      return()
    }
    if (!isTRUE(res$ok)) {
      attr(link, "source") <- NULL
      if (grepl("quota", res$message %||% "", fixed = TRUE)) {
        shiny::showNotification(res$message, type = "warning", duration = 8)
      }
    }
    proj$feature_link <- link
    current_project(proj)
    shiny::showNotification(
      sprintf("Mapping table saved with the project: %s rows.",
              format(nrow(link), big.mark = ",")),
      type = "message")
  })

  output$link_action <- shiny::renderUI({
    fp <- feature_pairing()
    if (is.null(fp) || !identical(fp$source, "project")) return(NULL)
    shiny::actionButton(session$ns("drop_link"), "Match by gene symbol instead",
                        class = "btn btn-sm btn-ghost")
  })

  shiny::observeEvent(input$drop_link, {
    proj <- current_project()
    shiny::req(proj)
    proj$feature_link <- NULL
    current_project(proj)
    shiny::showNotification("Mapping table removed; features are matched by gene symbol.",
                            type = "message")
  })

  list(feature_pairing = feature_pairing, link_upload = link_upload,
       pending_link = pending_link)
}

integration_pairing_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Sample pairing"),
      htmltools::tags$span(class = "card-sub",
                           "which sample is which person"),
      shiny::uiOutput(ns("pairing_action"), inline = TRUE)
    ),
    bslib::card_body(
      shiny::uiOutput(ns("pairing_note")),
      # Folded: the note says how many pairs and where they came from,
      # which is what most readers need; the rows are one click away for
      # the reader who wants to check them.
      htmltools::tags$details(
        class = "pairing-details",
        htmltools::tags$summary("Show the pairs"),
        DT::DTOutput(ns("pairing_table"))
      )
    )
  )
}

integration_feature_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Feature matching"),
      htmltools::tags$span(class = "card-sub", "which protein is which gene"),
      shiny::uiOutput(ns("link_action"), inline = TRUE)
    ),
    bslib::card_body(
      shiny::uiOutput(ns("feature_note")),
      htmltools::tags$details(
        class = "pairing-details",
        htmltools::tags$summary("Match with your own table"),
        htmltools::tags$div(
          class = "muted", style = "font-size:12px;margin:6px 0 8px",
          paste("A table with two columns: the identifiers of one layer (for",
                "example UniProt accessions) and the matching identifiers of",
                "the other (for example gene symbols or Ensembl gene ids), one",
                "pair per row. Features are then matched by this table instead",
                "of by gene symbol, and features it does not list are left out.",
                "An accession without an isoform suffix (P04637) also covers its",
                "isoforms (P04637-2). The table is saved with the project and",
                "replaces any table saved before.")),
        shiny::fileInput(
          ns("link_file"), label = NULL, multiple = FALSE,
          accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
          placeholder = "mapping.csv: two columns"),
        shiny::uiOutput(ns("link_columns"))
      )
    )
  )
}

# How the two layers' features meet, in one sentence: how many of each
# were matched, and how many share their partner with another feature of
# their layer (isoforms of one gene), each such pair being compared on
# its own.
feature_pairing_sentence <- function(fp, tag_a, tag_b) {
  n <- function(x) format(x, big.mark = ",")
  out <- sprintf("%s of %s features of %s are matched with %s of %s features of %s (%s pairs).",
                 n(fp$n_paired_a), n(fp$n_features_a), tag_a,
                 n(fp$n_paired_b), n(fp$n_features_b), tag_b, n(fp$n_pairs))
  if (fp$n_a_sharing > 0L) {
    out <- paste(out, sprintf(
      "%s features of %s share a match in %s with another (%s of its features have more than one).",
      n(fp$n_a_sharing), tag_a, tag_b, n(fp$n_b_shared)))
  }
  if (fp$n_b_sharing > 0L) {
    out <- paste(out, sprintf(
      "%s features of %s share a match in %s with another (%s of its features have more than one).",
      n(fp$n_b_sharing), tag_b, tag_a, n(fp$n_a_shared)))
  }
  out
}

feature_pairing_notice <- function(fp, tag_a, tag_b) {
  how <- switch(fp$source,
                project = "Matched by the mapping table saved with this project.",
                supplied = "Matched by the mapping table.",
                "Matched by gene symbol.")
  if (fp$n_pairs == 0L) {
    return(notice(
      "No features matched",
      paste(how, "Check that both layers carry gene symbols from the same",
            "organism, or match them with your own table below."),
      kind = "warn"))
  }
  shared <- fp$n_a_sharing > 0L || fp$n_b_sharing > 0L
  notice(
    paste(how, feature_pairing_sentence(fp, tag_a, tag_b)),
    if (shared) paste("Every feature is kept: concordance and correlation",
                      "compare each matched pair on its own, and ActivePathways",
                      "scores a gene by its most significant feature, allowing",
                      "for how many it has."),
    kind = "info")
}

# The likeliest columns of an uploaded table for each layer: the pair
# that matches the most features. Only the first few columns are tried;
# a mapping table is two columns, perhaps with a description beside them.
guess_link_columns <- function(tab, proj, tag_a, tag_b) {
  cols <- utils::head(names(tab), 6L)
  best <- cols[1:2]
  score <- -1
  for (a in cols) for (b in cols) {
    if (identical(a, b)) next
    link <- stats::setNames(tab[c(a, b)], c(tag_a, tag_b))
    fp <- tryCatch(omicsCore::feature_pairing_preview(proj, tag_a, tag_b, feature_link = link),
                   error = function(e) NULL)
    if (is.null(fp)) next
    sc <- fp$n_paired_a + fp$n_paired_b
    if (sc > score) {
      score <- sc
      best <- c(a, b)
    }
  }
  best
}
