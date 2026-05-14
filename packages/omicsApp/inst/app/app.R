# Shiny entry point for omicsApp.
#
# `omicsApp::launch()` runs `shiny::runApp()` against this directory.
# Everything that defines the UI / server lives in the package's
# `R/` folder so it benefits from roxygen, R CMD check, and reuse
# from tests. Keep this file thin: build UI, attach shinyjs, hand
# both halves to `shinyApp()`.

shiny::shinyApp(
  ui = htmltools::tagList(
    shinyjs::useShinyjs(),
    omicsApp:::app_ui()
  ),
  server = omicsApp:::app_server
)
