# Tests for the shared UI tag-list factories. They guard against
# silent CSS-class drift: if a class name changes here but the SCSS
# file still references the old name, the rendered widget looks
# unstyled in production.

test_that("view_header() renders title, subtitle, and action slot", {
  html <- render_html(view_header(
    title    = "QC",
    subtitle = "Sample-level diagnostics",
    actions  = shiny::actionButton("act", "Run")
  ))
  expect_match(html, 'class="page-head"', fixed = TRUE)
  expect_match(html, 'class="page-title"', fixed = TRUE)
  expect_match(html, ">QC</h1>", fixed = TRUE)
  expect_match(html, 'class="page-sub"', fixed = TRUE)
  expect_match(html, "Sample-level diagnostics", fixed = TRUE)
  expect_match(html, 'class="page-actions"', fixed = TRUE)
  expect_match(html, 'id="act"', fixed = TRUE)
})

test_that("view_header() omits actions when not provided", {
  html <- render_html(view_header("Project"))
  expect_match(html, 'class="page-head"', fixed = TRUE)
  expect_false(grepl("page-actions", html, fixed = TRUE))
  expect_false(grepl("page-sub", html, fixed = TRUE))
})

test_that("stat_card() applies label, value, trend, and accent class", {
  html <- render_html(stat_card(
    label = "Experiments", value = 2,
    trend = "Proteomics + RNA-seq",
    accent = "brand"
  ))
  expect_match(html, "stat-card accent-brand", fixed = TRUE)
  expect_match(html, ">Experiments</div>", fixed = TRUE)
  expect_match(html, ">2</div>", fixed = TRUE)
  expect_match(html, ">Proteomics + RNA-seq</div>", fixed = TRUE)
})

test_that("stat_card() mono = TRUE switches the value to monospace class", {
  html <- render_html(stat_card("Project ID", "abc-123", mono = TRUE))
  expect_match(html, '"stat-value mono"', fixed = TRUE)
})

test_that("stat_card() rejects unknown accent values", {
  expect_error(stat_card("X", 1, accent = "unknown"), "should be one of")
})

test_that("pill() picks the right kind class", {
  for (k in c("info", "up", "down", "ns", "ok", "warn")) {
    html <- render_html(pill("x", kind = k))
    expect_match(html, paste0("pill pill-", k), fixed = TRUE)
  }
})

test_that("confidence_bar() clamps and renders the right CSS variable", {
  expect_match(render_html(confidence_bar(0.92, "ok")),
               "width:92%;background:var(--ok)", fixed = TRUE)
  expect_match(render_html(confidence_bar(-0.5, "err")),
               "width:0%;background:var(--err)", fixed = TRUE)
  expect_match(render_html(confidence_bar(1.5, "warn")),
               "width:100%;background:var(--warn)", fixed = TRUE)
  expect_match(render_html(confidence_bar(0.5, "brand")),
               "background:var(--brand-600)", fixed = TRUE)
})

test_that("confidence_bar() rejects non-numeric input", {
  expect_error(confidence_bar("hi"), "numeric scalar")
  expect_error(confidence_bar(NA_real_), "numeric scalar")
  expect_error(confidence_bar(c(0.1, 0.2)), "numeric scalar")
})

test_that("view_placeholder() still emits a card with title + body", {
  # Slice 2A's placeholder must keep working — the seven view stubs
  # still rely on it.
  html <- render_html(view_placeholder("Report", "Build a shareable report"))
  expect_match(html, "card-title", fixed = TRUE)
  expect_match(html, "card-sub", fixed = TRUE)
  expect_match(html, "Report view", fixed = TRUE)
  expect_match(html, "placeholder-pad", fixed = TRUE)
})
