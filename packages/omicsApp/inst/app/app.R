# Placeholder Shiny app for Phase 0.
# Phase 2 will replace this with the modular bslib-based UI.

library(shiny)
library(bslib)

ui <- page_fluid(
  theme = bs_theme(version = 5, bootswatch = "minty"),
  card(
    card_header("omicsApp"),
    card_body(
      h3("Phase 0 skeleton"),
      p("The full multi-omics analysis interface will be implemented in Phase 2."),
      p("See ", tags$code("docs/export-manifest.md"), " for the planned API.")
    )
  )
)

server <- function(input, output, session) {
  # No reactives yet.
}

shinyApp(ui, server)
