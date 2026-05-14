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
stat_card <- function(label, value, trend = NULL,
                      accent = c("brand", "up", "down", "ok"),
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
