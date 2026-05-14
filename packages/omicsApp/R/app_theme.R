#' Build the omicsApp Bootstrap 5 theme
#'
#' Returns a [bslib::bs_theme()] configured with the design tokens
#' shared with `mockup/styles.scss`. The companion custom-class rules
#' (stat cards, pills, schema rows, ...) are layered on via
#' [bslib::bs_add_rules()] reading `inst/app/www/styles.scss`.
#'
#' Keeping design tokens in `bs_theme()` arguments — rather than only
#' in the SCSS file — means Bootstrap's own variables
#' (`--bs-primary`, `--bs-body-bg`, ...) inherit our palette, so any
#' stock Bootstrap component picks up the brand colors without
#' per-component overrides.
#'
#' @return A `bs_theme` object suitable for the `theme` argument of
#'   bslib page constructors.
#'
#' @keywords internal
#' @noRd
app_theme <- function() {
  theme <- bslib::bs_theme(
    version = 5,
    bg = "#F5F6FA",
    fg = "#0F1320",
    primary   = "#2C3E99",
    secondary = "#14A085",
    success   = "#10B981",
    info      = "#3B82F6",
    warning   = "#F59E0B",
    danger    = "#EF4444",
    "body-bg"           = "#F5F6FA",
    "body-color"        = "#0F1320",
    "border-color"      = "#E5E7EB",
    "border-radius"     = "8px",
    "card-border-color" = "#E5E7EB",
    "card-cap-bg"       = "#FFFFFF",
    "font-family-base"      = "Inter, system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif",
    "font-family-monospace" = "'JetBrains Mono', 'SF Mono', Menlo, Consolas, monospace"
  )

  scss_path <- system.file("app", "www", "styles.scss", package = "omicsApp")
  if (nzchar(scss_path) && file.exists(scss_path)) {
    theme <- bslib::bs_add_rules(theme, sass::sass_file(scss_path))
  }
  theme
}
