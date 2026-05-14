#' QC view module (Slice 2A stub)
#'
#' Phase 2 slice 2A scaffold. PCA, missingness, intensity QC land in
#' slice 2D.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
qc_view_ui <- function(id) {
  ns <- shiny::NS(id)
  view_placeholder("QC", "Sample-level diagnostics: PCA, missingness, intensity distribution.")
}

#' @rdname qc_view_ui
#' @keywords internal
#' @noRd
qc_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
  })
}
