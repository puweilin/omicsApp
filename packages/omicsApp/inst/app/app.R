# Shiny entry point for omicsApp.
#
# `omicsApp::launch()` runs `shiny::runApp()` against this directory.
# Everything that defines the UI / server lives in the package's
# `R/` folder so it benefits from roxygen, R CMD check, and reuse
# from tests. Keep this file thin: build UI, attach shinyjs, hand
# both halves to `shinyApp()`.

# Development hook. When OMICSAPP_DEV_ROOT names the monorepo's
# `packages/` directory, both packages are loaded from source instead of
# from the library. The shinytest2 harness sets it, because otherwise a
# browser test exercises whatever was last installed -- which was four
# months old the day this was added -- and its verdict says nothing
# about the code being changed.
dev_root <- Sys.getenv("OMICSAPP_DEV_ROOT", "")
if (nzchar(dev_root)) {
  pkgload::load_all(file.path(dev_root, "omicsCore"), quiet = TRUE)
  pkgload::load_all(file.path(dev_root, "omicsApp"), quiet = TRUE)
}

shiny::shinyApp(
  ui = htmltools::tagList(
    shinyjs::useShinyjs(),
    omicsApp:::app_ui()
  ),
  server = omicsApp:::app_server
)
