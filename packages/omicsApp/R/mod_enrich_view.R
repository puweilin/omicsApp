#' Enrichment view module (Slice 2A stub)
#'
#' Phase 2 slice 2A scaffold. Dotplot + ORA / GSEA tables land in
#' slice 2E.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
enrich_view_ui <- function(id) {
  ns <- shiny::NS(id)
  view_placeholder("Enrichment", "ORA, GSEA, GSVA over the active differential bundle.")
}

#' @rdname enrich_view_ui
#' @keywords internal
#' @noRd
enrich_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
  })
}
