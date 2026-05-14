#' Report view module (Slice 2A stub)
#'
#' Phase 2 slice 2A scaffold. Report builder lands in slice 2F.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
report_view_ui <- function(id) {
  ns <- shiny::NS(id)
  view_placeholder("Report", "Build a shareable PDF / HTML report from the current bundle.")
}

#' @rdname report_view_ui
#' @keywords internal
#' @noRd
report_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
  })
}
