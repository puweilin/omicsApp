#' Internal helpers shared across UI modules
#'
#' These helpers will grow in slice 2B (`stat_card()`, `view_header()`,
#' `pill()`, `confidence_bar()`). For 2A we ship only the placeholder
#' used by every empty view stub.
#'
#' @keywords internal
#' @name ui-helpers
#' @noRd
NULL

#' Empty styled card with title + subtitle
#'
#' Used by slice-2A view stubs. Will be replaced per-view in 2C-2F.
#'
#' @param title View name (shown as card header).
#' @param subtitle One-liner describing what the view will do.
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
