#' Internal helpers shared across UI modules
#'
#' Tag-list factories that wrap the CSS classes declared in
#' `inst/app/www/styles.scss`. Keeping these centralised means every
#' view module emits identical markup for the same widget — and lets
#' the SCSS file remain the single source of truth for visual rules.
#'
#' All helpers are package-internal (`@keywords internal`, `@noRd`).
#' They're meant to be called from inside `*_view_ui()` definitions,
#' not by end-users.
#'
#' @keywords internal
#' @name ui-helpers
#' @noRd
NULL

# ---- view scaffolding ------------------------------------------------

#' Empty styled card with title + subtitle
#'
#' Used by the slice-2A view stubs while real per-view content is
#' deferred to 2C-2F. Once a view ships real content, the call to
#' `view_placeholder()` is replaced inline; the helper itself stays
#' so other views can keep using it.
#'
#' @param title View name (shown in the card header).
#' @param subtitle Optional one-liner under the title.
#'
#' @return A `bslib::card` tag list.
#'
#' @keywords internal
#' @noRd
view_placeholder <- function(title, subtitle = NULL) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$div(class = "card-title", title),
      if (!is.null(subtitle)) htmltools::tags$div(class = "card-sub", subtitle)
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "placeholder-pad",
        htmltools::tags$h3(paste0(title, " view")),
        htmltools::tags$p(
          class = "muted",
          "Content will land in a later Phase 2 slice."
        )
      )
    )
  )
}

#' Page-level header row for a view
#'
#' Renders the `.page-head` block defined in styles.scss: title on the
#' left, optional subtitle directly below it, optional action buttons
#' pushed to the right.
#'
#' Reference markup: `mockup/index.html:201-204`.
#'
#' @param title Title text (rendered inside `<h1 class="page-title">`).
#' @param subtitle Optional secondary line (rendered inside
#'   `<div class="page-sub">`).
#' @param actions Optional `htmltools::tag` / `tagList` of right-aligned
#'   action buttons. NULL skips the actions slot entirely.
#'
#' @return A single `<div class="page-head">` tag.
#'
#' @keywords internal
#' @noRd
view_header <- function(title, subtitle = NULL, actions = NULL) {
  htmltools::tags$div(
    class = "page-head",
    htmltools::tags$div(
      htmltools::tags$h1(class = "page-title", title),
      if (!is.null(subtitle)) htmltools::tags$div(class = "page-sub", subtitle)
    ),
    if (!is.null(actions)) htmltools::tags$div(class = "page-actions", actions)
  )
}

# ---- atomic widgets --------------------------------------------------

#' Stat card (KPI tile)
#'
#' Emits `<div class="card stat-card accent-{accent}">` with a label,
#' a big value, and an optional trend line below. Used at the top of
#' the Project, QC, and Differential views.
#'
#' Reference markup: `mockup/index.html:223-234`.
#'
#' @param label Short caption (uppercased by CSS).
#' @param value Headline number / text. Strings are rendered as-is.
#' @param trend Optional trend line under the value. May be a string
#'   or a `htmltools::tag`.
#' @param accent Border-top color: one of `"brand"`, `"up"`, `"down"`,
#'   `"ok"`. Defaults to `"brand"`.
#' @param mono Logical. If `TRUE`, render the value in the monospace
#'   font (for IDs, hashes, ratios).
#'
#' @return A `<div class="card stat-card ...">` tag.
#'
#' @keywords internal
#' @noRd
#' @section Accents:
#' `"warn"` is here because two callers already asked for it and got
#' `match.arg()`'s error instead: the QC stats card whenever a sample is
#' flagged as an outlier, and the import summary whenever a matrix is
#' more than half missing. Both are warnings, and mapping them onto
#' `"up"` for its red would have spent the volcano's up-regulated colour
#' on something that is not a direction.
#'
#' Every value here needs a matching `.stat-card.accent-*` rule in
#' `inst/app/www/styles.scss`; the two are checked against each other in
#' `test-stat-card-accents.R`.
stat_card <- function(label, value, trend = NULL,
                      accent = c("brand", "up", "down", "ok", "warn"),
                      mono = FALSE) {
  accent <- match.arg(accent)
  value_class <- paste(c("stat-value", if (isTRUE(mono)) "mono"), collapse = " ")
  htmltools::tags$div(
    class = paste0("card stat-card accent-", accent),
    htmltools::tags$div(class = "stat-label", label),
    htmltools::tags$div(class = value_class, value),
    if (!is.null(trend)) htmltools::tags$div(class = "stat-trend", trend)
  )
}

#' Status pill
#'
#' `<span class="pill pill-{kind}">…</span>`. Used for direction
#' indicators (up/down/ns) in tables and for status badges
#' (ok/warn/info) in the import wizard.
#'
#' Reference markup: `mockup/index.html:255-261, 827`.
#'
#' @param text Pill content (string or inline tag).
#' @param kind One of `"info"`, `"up"`, `"down"`, `"ns"`, `"ok"`,
#'   `"warn"`. Defaults to `"info"`.
#'
#' @return A `<span>` tag.
#'
#' @keywords internal
#' @noRd
pill <- function(text, kind = c("info", "up", "down", "ns", "ok", "warn")) {
  kind <- match.arg(kind)
  htmltools::tags$span(
    class = paste0("pill pill-", kind),
    text
  )
}

#' Confidence bar
#'
#' Renders `<div class="conf-bar"><div class="fill"/></div>` with the
#' inner fill width set to `value * 100%` and its background colour
#' driven by an accent CSS variable. Used by the Import view (slice
#' 2C) to show the smart-parse confidence per inferred sheet role.
#'
#' Reference markup: `mockup/index.html:331-332, 640, 651`.
#'
#' @param value Numeric in `[0, 1]`; clamped if outside.
#' @param accent One of `"ok"`, `"warn"`, `"err"`, `"brand"`. Picks
#'   the CSS variable backing the fill colour.
#'
#' @return A `<div class="conf-bar">` tag.
#'
#' @keywords internal
#' @noRd
confidence_bar <- function(value, accent = c("ok", "warn", "err", "brand")) {
  accent <- match.arg(accent)
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    stop("`value` must be a numeric scalar.", call. = FALSE)
  }
  pct <- max(0, min(100, round(100 * value)))
  var_name <- switch(accent,
                     ok    = "--ok",
                     warn  = "--warn",
                     err   = "--err",
                     brand = "--brand-600")
  htmltools::tags$div(
    class = "conf-bar",
    htmltools::tags$div(
      class = "fill",
      style = sprintf("width:%d%%;background:var(%s)", pct, var_name)
    )
  )
}

# ---- import-wizard widgets ------------------------------------------

#' Single step in a `.steps` wizard strip
#'
#' Emits `<div class="step {state}">` containing a small numbered
#' badge and a (label, desc) pair. `state = "done"` swaps the number
#' for a check glyph so a strip of steps reads at a glance.
#'
#' Reference markup: `mockup/index.html:573-587`.
#'
#' @param index Step number (rendered as text when `state != "done"`).
#' @param label One-line step name.
#' @param desc Optional sub-line describing the step's current state.
#' @param state One of `"pending"`, `"active"`, `"done"`. The pending
#'   state adds no extra class — that matches the bare `<div class=
#'   "step">` in the mockup.
#'
#' @return A `<div class="step ...">` tag.
#'
#' @keywords internal
#' @noRd
step_item <- function(index, label, desc = NULL,
                      state = c("pending", "active", "done")) {
  state <- match.arg(state)
  classes <- c("step", if (state != "pending") state)
  ix_text <- if (state == "done") "\u2713" else as.character(index)
  htmltools::tags$div(
    class = paste(classes, collapse = " "),
    htmltools::tags$div(class = "ix", ix_text),
    htmltools::tags$div(
      htmltools::tags$div(class = "label", label),
      if (!is.null(desc)) htmltools::tags$div(class = "desc", desc)
    )
  )
}

#' One row inside the "Inferred schema" card on the Import view
#'
#' Lays out the five-column grid the mockup uses: index tag, title
#' and description, role pill, confidence bar + percent, action slot.
#' Reuses `pill()` for the role tag and `confidence_bar()` for the
#' fill so a single visual rule governs both widgets.
#'
#' Reference markup: `mockup/index.html:633-664`.
#'
#' @param ix Short index tag (e.g. `"S1"`).
#' @param title Bold one-liner describing the inferred role.
#' @param desc Optional secondary line with the parse evidence.
#' @param role Pill text (e.g. `"expression matrix"`).
#' @param confidence Numeric in `[0, 1]`; passed through to
#'   `confidence_bar()` and rendered as `XX%` next to the bar.
#' @param accent One of `"ok"`, `"warn"`, `"err"`, `"brand"` —
#'   forwarded to `confidence_bar()`.
#' @param role_kind One of `"info"`, `"ok"`, `"warn"`, `"up"`,
#'   `"down"`, `"ns"` — forwarded to `pill()`.
#' @param actions Optional tag (typically a button) rendered in the
#'   right-most slot. `NULL` leaves the cell empty.
#'
#' @return A `<div class="schema-row">` tag.
#'
#' @keywords internal
#' @noRd
schema_row <- function(ix, title, desc = NULL, role, confidence,
                       accent    = c("ok", "warn", "err", "brand"),
                       role_kind = c("info", "ok", "warn", "up", "down", "ns"),
                       actions   = NULL) {
  accent    <- match.arg(accent)
  role_kind <- match.arg(role_kind)
  if (!is.numeric(confidence) || length(confidence) != 1L || is.na(confidence)) {
    stop("`confidence` must be a numeric scalar.", call. = FALSE)
  }
  pct <- max(0, min(100, round(100 * confidence)))
  htmltools::tags$div(
    class = "schema-row",
    htmltools::tags$div(class = "ix", ix),
    htmltools::tags$div(
      htmltools::tags$div(class = "title", title),
      if (!is.null(desc)) htmltools::tags$div(class = "desc", desc)
    ),
    pill(role, kind = role_kind),
    htmltools::tags$div(
      class = "conf",
      confidence_bar(confidence, accent = accent),
      htmltools::tags$span(
        class = "text-mono",
        style = "font-size:11px",
        sprintf("%d%%", pct)
      )
    ),
    actions
  )
}

#' Inline notice / alert block
#'
#' Emits `<div class="notice notice-{kind}">` with a leading bsicon
#' and a (title, detail) body. Used in the Import view to flag
#' schema warnings and elsewhere for inline status messages.
#'
#' Reference markup: `mockup/index.html:666-674`.
#'
#' @param title Bold lead line.
#' @param detail Optional secondary line rendered in `.muted`.
#' @param kind One of `"info"`, `"warn"`, `"error"`. Picks the background
#'   colour and the leading icon.
#' @param technical Optional raw message (an R condition message, say).
#'   Folded away under "Technical details" so the notice leads with what
#'   happened in words a user can act on, and the message a developer
#'   needs is still one click away.
#'
#' @return A `<div class="notice ...">` tag.
#'
#' @keywords internal
#' @noRd
notice <- function(title, detail = NULL, kind = c("info", "warn", "error"),
                   technical = NULL) {
  kind <- match.arg(kind)
  icon_name <- switch(kind,
                      info  = "info-circle",
                      warn  = "exclamation-triangle",
                      error = "x-octagon")
  htmltools::tags$div(
    class = paste0("notice notice-", kind),
    role = if (kind == "error") "alert" else NULL,
    bsicons::bs_icon(icon_name, class = "icon"),
    htmltools::tags$div(
      htmltools::tags$strong(title),
      if (!is.null(detail)) htmltools::tags$div(class = "muted", detail),
      if (!is.null(technical) && length(technical) && nzchar(technical[[1L]])) {
        htmltools::tags$details(
          class = "notice-technical",
          htmltools::tags$summary("Technical details"),
          htmltools::tags$code(technical)
        )
      }
    )
  )
}

#' Single uploaded-file row in the Import view's file list
#'
#' Renders `<div class="file-row">` with a leading icon, a name +
#' optional sub-line, and a right-aligned size label.
#'
#' Reference markup: `mockup/index.html:595-610`.
#'
#' @param name File name (rendered in `.name`).
#' @param meta Optional sub-line (rendered in muted text).
#' @param size Optional right-aligned size string (rendered in
#'   monospace via `.meta`).
#' @param icon Leading icon tag. Defaults to a bsicon file glyph.
#'
#' @return A `<div class="file-row">` tag.
#'
#' @keywords internal
#' @noRd
file_row <- function(name, meta = NULL, size = NULL,
                     icon = bsicons::bs_icon("file-earmark-text", class = "icon")) {
  htmltools::tags$div(
    class = "file-row",
    icon,
    htmltools::tags$div(
      htmltools::tags$div(class = "name", name),
      if (!is.null(meta)) htmltools::tags$div(
        class = "muted",
        style = "font-size:11.5px",
        meta
      )
    ),
    if (!is.null(size)) htmltools::tags$div(class = "meta", size)
  )
}

# ---- diff-view sidebar widgets --------------------------------------

#' Single labelled group inside a `.param-stack` sidebar
#'
#' Wraps a section header (`<h4>`) and arbitrary form controls in
#' the `.param-group` block defined by `styles.scss`. Used by the
#' Diff view (slice 2D) to lay out the Method / Contrast /
#' Covariates / Thresholds groups.
#'
#' Reference markup: `mockup/index.html:761-810`.
#'
#' @param title Group label (uppercased by CSS).
#' @param ... Form controls (selectInput / sliderInput / switches).
#'
#' @return A `<div class="param-group">` tag.
#'
#' @keywords internal
#' @noRd
param_group <- function(title, ..., help = NULL) {
  htmltools::tags$div(
    class = "param-group",
    htmltools::tags$h4(title, if (!is.null(help)) info_tip(help)),
    ...
  )
}

#' A question-mark icon that explains a control on hover or focus
#'
#' For the statistical choices a user cannot be expected to know the
#' meaning of (what "auto" picks, what a covariate does, ORA vs GSEA).
#' Keyboard-reachable: the icon takes focus and the tooltip opens on it.
#'
#' @param text The explanation, one or two sentences.
#' @return A tag.
#' @keywords internal
#' @noRd
info_tip <- function(text) {
  bslib::tooltip(
    htmltools::tags$span(
      class = "info-tip", tabindex = "0", role = "img",
      `aria-label` = text,
      bsicons::bs_icon("question-circle")
    ),
    text,
    placement = "right"
  )
}

#' Inline legend chip
#'
#' Renders `<span><span class="swatch" .../>{label}</span>` for use
#' inside a `.legend` row. Centralising it keeps the inline-style
#' sprawl out of the view modules.
#'
#' Reference markup: `mockup/index.html:840-845`.
#'
#' @param label Legend text.
#' @param color CSS color or `var(--brand-600)`-style variable.
#'
#' @return A `<span>` tag.
#'
#' @keywords internal
#' @noRd
legend_swatch <- function(label, color) {
  htmltools::tags$span(
    htmltools::tags$span(
      class = "swatch",
      style = sprintf("background:%s", color)
    ),
    label
  )
}

# Static plots are drawn at 96 dpi rather than Shiny's 72: the same
# pixel size, with text a third larger. At 72 the axis text was about
# 7 px on a desktop and 5-6 px on a phone.
PLOT_RES <- 96

# Height for one row per comparison, each as tall as its label: a long
# name wraps over up to three lines (omicsCore::wrap_comparison()), and
# a fixed row height let wrapped labels run into each other.
label_rows_px <- function(comparisons, pad_px) {
  lines <- omicsCore::comparison_label_lines(comparisons)
  as.integer(sum(pad_px + 19L * lines))
}

# A run button that cannot be pressed again while its run is in flight,
# the pattern the Differential and Integration buttons already had.
# Enrichment's Re-run, "Enrich every comparison" and "Run global test"
# stayed live through an 8 s GSEA, and a second press queued a second
# run behind the first. `flag` is a reactiveVal the module keeps, so a
# button re-rendered by renderUI mid-run comes back disabled too.
set_button_busy <- function(id, busy, flag = NULL) {
  if (!is.null(flag)) flag(busy)
  tryCatch(if (busy) shinyjs::disable(id) else shinyjs::enable(id),
           error = function(e) NULL)
  invisible(busy)
}

disabled_if <- function(tag, busy) {
  if (isTRUE(busy)) htmltools::tagAppendAttributes(tag, disabled = NA) else tag
}

# On a phone a plot is about 280 px wide, and at 96 dpi its titles and
# legends no longer fit: the "Hits per comparison" legend and the PCA
# group names ran off the right edge. Below NARROW_PLOT_PX the text is
# set smaller and legends stack, rather than the figure being cut.
NARROW_PLOT_PX <- 420

fit_to_width <- function(output_id, p, session = shiny::getDefaultReactiveDomain()) {
  w <- tryCatch(session$clientData[[paste0("output_", session$ns(output_id), "_width")]],
                error = function(e) NULL)
  if (!is.numeric(w) || !length(w) || w >= NARROW_PLOT_PX) return(p)
  # Stacked only where a legend sits under the panel; one above it (up /
  # down) fits on a line.
  top <- identical(tryCatch(p$theme$legend.position, error = function(e) NULL), "top")
  shrink <- ggplot2::theme(text = ggplot2::element_text(size = 9))
  # A legend at the side took most of a 280 px panel: the enrichment dot
  # plot's points were squeezed into a sliver. Under the panel instead,
  # its entries stacked.
  # Legends run across, each on its own row: two upright legends side by
  # side (the enrichment colour bar and its size key) did not fit in the
  # width and were cut off. Group legends keep one entry per line (their
  # guide asks for one column), so long group names do not run off.
  if (!top) {
    shrink <- shrink + ggplot2::theme(legend.position = "bottom",
                                      legend.direction = "horizontal",
                                      legend.box = "vertical")
    # Titles above their keys, not beside them: "-log10(adjusted p)"
    # beside a colour bar was wider than the panel. (ggplot2 >= 3.5.)
    if ("legend.title.position" %in% names(ggplot2::get_element_tree())) {
      shrink <- shrink + ggplot2::theme(legend.title.position = "top")
    }
  }
  if (inherits(p, "patchwork")) p & shrink
  else if (inherits(p, c("gg", "ggplot"))) p + shrink
  else p
}
