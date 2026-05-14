#' Project view module (Slice 2A stub)
#'
#' Phase 2 slice 2A scaffold. The UI returns a placeholder card; the
#' server is empty. Real content lands in slice 2C.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
project_view_ui <- function(id) {
  ns <- shiny::NS(id)
  view_placeholder("Project", "Manage sessions, omics layers, and global parameters.")
}

#' @rdname project_view_ui
#' @keywords internal
#' @noRd
project_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # No reactives yet (slice 2A).
  })
}
