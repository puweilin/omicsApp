#' Import view module (Slice 2A stub)
#'
#' Phase 2 slice 2A scaffold. Real content (dropzone, smart-parse
#' schema review) lands in slice 2C.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
import_view_ui <- function(id) {
  ns <- shiny::NS(id)
  view_placeholder("Import", "Drop files, review the inferred schema, run the importer.")
}

#' @rdname import_view_ui
#' @keywords internal
#' @noRd
import_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
  })
}
