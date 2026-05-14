#' Import view module
#'
#' Three-step wizard that walks the user from "drop a file" to
#' "build an `omics_input`". Slice 2C renders the **review** step
#' (step 2) literally — a snapshot of two fake uploaded files plus
#' the inferred-schema card. No upload handler is wired yet; the
#' dropzone is a visual placeholder and the "Edit" / "Re-detect"
#' buttons are inert. Real wiring lands in a later phase.
#'
#' Reference markup: `omicsApp/mockup/index.html:565-683`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
import_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    view_header(
      title    = "Import data",
      subtitle = "Auto-detect expression matrix, sample metadata, and feature annotation"
    ),
    import_steps_strip(),
    htmltools::tags$div(
      class = "row-grid r-4-8",
      import_files_card(),
      import_schema_card()
    )
  )
}

#' @rdname import_view_ui
#' @keywords internal
#' @noRd
import_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Slice 2C: static mock; no upload handler yet.
  })
}

# ---- internal helpers ------------------------------------------------

# Chevron between two steps. Inline SVG so we don't depend on a
# specific bsicon name for what is essentially a typographic glyph.
step_arrow <- function() {
  htmltools::tags$svg(
    class    = "icon step-arrow",
    viewBox  = "0 0 24 24",
    fill     = "none",
    stroke   = "currentColor",
    width    = "16",
    height   = "16",
    htmltools::tags$path(d = "M9 6l6 6-6 6")
  )
}

import_steps_strip <- function() {
  htmltools::tags$div(
    class = "steps",
    step_item(1L, "Upload", "2 files \u00B7 14.6 MB", state = "done"),
    step_arrow(),
    step_item(2L, "Review inferred schema", "Confidence 92%", state = "active"),
    step_arrow(),
    step_item(3L, "Confirm & import", "Creates omics_input objects",
              state = "pending")
  )
}

import_files_card <- function() {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Uploaded files")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "file-tree",
        file_row("CHISSS2_Cheek.xlsx",
                 meta = "3 sheets \u00B7 proteomics",
                 size = "2.4 MB"),
        file_row("CHISSS2_Cheek_RNA.xlsx",
                 meta = "2 sheets \u00B7 rna-seq",
                 size = "12.2 MB")
      ),
      htmltools::tags$div(class = "divider"),
      htmltools::tags$div(
        class = "dropzone",
        style = "padding:22px 16px",
        bsicons::bs_icon("cloud-upload", class = "icon icon-lg", size = "36px"),
        htmltools::tags$div(
          style = "font-weight:500;font-size:13px",
          "Drop more files"
        ),
        htmltools::tags$div(
          class = "muted",
          style = "font-size:12px",
          "Excel \u00B7 CSV \u00B7 TSV \u00B7 Salmon dir \u00B7 MaxQuant"
        )
      )
    )
  )
}

import_schema_card <- function() {
  # Inert static buttons: slice 2C doesn't wire any Re-detect / Edit
  # handlers, so we render plain HTML buttons instead of shiny
  # actionButtons (which would force a unique inputId for each row).
  ghost_btn <- function(label) {
    htmltools::tags$button(
      type  = "button",
      class = "btn btn-ghost",
      style = "padding:5px 10px;font-size:12px",
      label
    )
  }
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(
        class = "card-title",
        "Inferred schema \u00B7 ",
        htmltools::tags$span(class = "text-mono", "CHISSS2_Cheek.xlsx")
      ),
      htmltools::tags$span(class = "card-sub", "click any row to override"),
      htmltools::tags$div(
        class = "card-actions",
        ghost_btn("Re-detect")
      )
    ),
    bslib::card_body(
      schema_row(
        ix    = "S1",
        title = "Matrix \u00B7 features \u00D7 samples",
        desc  = "1,842 rows \u00B7 36 numeric columns \u00B7 orientation auto-detected",
        role  = "expression matrix",
        confidence = 0.96,
        actions    = ghost_btn("Edit")
      ),
      schema_row(
        ix    = "S2",
        title = "Pheno \u00B7 sample metadata",
        desc  = htmltools::tagList(
          "primary key ",
          htmltools::tags$span(class = "text-mono", "label"),
          " \u00B7 grouping vars: ",
          htmltools::tags$span(class = "text-mono", "Group, Age, Sex")
        ),
        role  = "sample meta",
        confidence = 0.91,
        actions    = ghost_btn("Edit")
      ),
      schema_row(
        ix    = "S3",
        title = "GeneMap \u00B7 feature annotation",
        desc  = "UniProt ID \u2192 Gene Symbol \u00B7 1,842 rows",
        role  = "feature annot",
        confidence = 0.88,
        actions    = ghost_btn("Edit")
      ),
      htmltools::tags$div(
        style = "margin-top:14px",
        notice(
          title  = "3 samples found in matrix but missing from metadata",
          detail = htmltools::tagList(
            htmltools::tags$span(class = "text-mono", "S023, S041, S099"),
            htmltools::HTML(" &middot; "),
            "will be excluded unless you add them to Pheno."
          ),
          kind   = "warn"
        )
      ),
      htmltools::tags$div(
        style = "display:flex;gap:8px;justify-content:flex-end;margin-top:18px",
        htmltools::tags$button(type = "button", class = "btn btn-ghost",
                               "Back"),
        htmltools::tags$button(type = "button", class = "btn btn-primary",
                               "Confirm & build omics_input")
      )
    )
  )
}
