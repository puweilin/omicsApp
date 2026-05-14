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

# ---- slice 2C helpers ------------------------------------------------

test_that("step_item() reflects state in the outer class", {
  expect_match(render_html(step_item(1, "Upload", state = "done")),
               'class="step done"', fixed = TRUE)
  expect_match(render_html(step_item(2, "Review", state = "active")),
               'class="step active"', fixed = TRUE)
  pending <- render_html(step_item(3, "Confirm", state = "pending"))
  expect_match(pending, 'class="step"', fixed = TRUE)
  # pending must NOT carry the active/done modifier
  expect_false(grepl("step active", pending, fixed = TRUE))
  expect_false(grepl("step done",   pending, fixed = TRUE))
})

test_that("step_item() done-state replaces the index with a check mark", {
  html_done <- render_html(step_item(1, "Upload", state = "done"))
  expect_match(html_done, "\u2713", fixed = TRUE)
  expect_false(grepl(">1<", html_done, fixed = TRUE))

  html_active <- render_html(step_item(2, "Review", state = "active"))
  expect_match(html_active, ">2<", fixed = TRUE)
})

test_that("step_item() omits desc when not provided", {
  html <- render_html(step_item(1, "Upload", state = "pending"))
  expect_match(html, 'class="label"', fixed = TRUE)
  expect_false(grepl('class="desc"', html, fixed = TRUE))
})

test_that("schema_row() renders the five-column grid", {
  html <- render_html(schema_row(
    ix = "S1", title = "Matrix",
    desc = "1842 rows",
    role = "expression matrix",
    confidence = 0.96,
    accent = "ok", role_kind = "info"
  ))
  expect_match(html, 'class="schema-row"',  fixed = TRUE)
  expect_match(html, 'class="ix">S1</div>', fixed = TRUE)
  expect_match(html, 'class="title">Matrix', fixed = TRUE)
  expect_match(html, 'class="desc">1842 rows', fixed = TRUE)
  expect_match(html, "pill pill-info", fixed = TRUE)
  expect_match(html, 'class="conf"',   fixed = TRUE)
  # conf-bar fill width + colour from confidence_bar()
  expect_match(html, "width:96%;background:var(--ok)", fixed = TRUE)
  # percentage label
  expect_match(html, ">96%</span>", fixed = TRUE)
})

test_that("schema_row() forwards role_kind and accent overrides", {
  html <- render_html(schema_row(
    ix = "S2", title = "Pheno", role = "sample meta",
    confidence = 0.4,
    accent = "warn", role_kind = "warn"
  ))
  expect_match(html, "pill pill-warn",  fixed = TRUE)
  expect_match(html, "background:var(--warn)", fixed = TRUE)
})

test_that("schema_row() renders the actions slot when provided", {
  edit <- htmltools::tags$button(type = "button", id = "rowEdit", "Edit")
  html <- render_html(schema_row(
    ix = "S1", title = "Matrix",
    role = "expression matrix",
    confidence = 0.9, actions = edit
  ))
  expect_match(html, 'id="rowEdit"', fixed = TRUE)
})

test_that("schema_row() rejects non-numeric confidence", {
  expect_error(
    schema_row("S1", "Matrix", role = "expression matrix",
               confidence = "high"),
    "numeric scalar"
  )
})

test_that("notice() emits the right modifier class and icon", {
  warn <- render_html(notice("Heads up", "details", kind = "warn"))
  expect_match(warn, "notice notice-warn", fixed = TRUE)
  # bsicons sets a class like "bi bi-exclamation-triangle"
  expect_match(warn, "bi-exclamation-triangle", fixed = TRUE)
  expect_match(warn, "<strong>Heads up</strong>", fixed = TRUE)
  expect_match(warn, 'class="muted">details', fixed = TRUE)

  info <- render_html(notice("FYI", kind = "info"))
  expect_match(info, "notice notice-info", fixed = TRUE)
  expect_match(info, "bi-info-circle", fixed = TRUE)
})

test_that("notice() omits detail when not provided", {
  html <- render_html(notice("Just the title", kind = "info"))
  expect_false(grepl('class="muted"', html, fixed = TRUE))
})

test_that("file_row() renders name, meta, size, and default icon", {
  html <- render_html(file_row("data.xlsx", meta = "3 sheets", size = "2.4 MB"))
  expect_match(html, 'class="file-row"', fixed = TRUE)
  expect_match(html, 'class="name">data.xlsx', fixed = TRUE)
  expect_match(html, "3 sheets", fixed = TRUE)
  expect_match(html, 'class="meta">2.4 MB', fixed = TRUE)
  # default file icon
  expect_match(html, "bi-file-earmark-text", fixed = TRUE)
})

test_that("file_row() omits optional slots when not provided", {
  html <- render_html(file_row("data.xlsx"))
  expect_match(html, 'class="name">data.xlsx', fixed = TRUE)
  expect_false(grepl('class="meta">', html, fixed = TRUE))
  # the muted sub-line is also gone
  expect_false(grepl("font-size:11.5px", html, fixed = TRUE))
})
