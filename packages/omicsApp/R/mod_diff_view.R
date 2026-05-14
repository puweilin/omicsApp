#' Differential analysis view module (Slice 2A stub)
#'
#' Phase 2 slice 2A scaffold. Volcano + top-hits + param controls land
#' in slice 2D.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
diff_view_ui <- function(id) {
  ns <- shiny::NS(id)
  view_placeholder("Differential", "Choose a contrast, run a method, inspect the volcano.")
}

#' @rdname diff_view_ui
#' @keywords internal
#' @noRd
diff_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
  })
}
