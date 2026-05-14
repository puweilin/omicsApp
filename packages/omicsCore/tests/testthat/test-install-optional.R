test_that("resolve_install_group expands 'all' to the union", {
  all_pkgs <- omicsCore:::resolve_install_group("all")
  expect_true(all(omicsCore:::OPTIONAL_GROUPS$rnaseq %in% all_pkgs))
  expect_true(all(omicsCore:::OPTIONAL_GROUPS$enrichment %in% all_pkgs))
  expect_false(anyDuplicated(all_pkgs) > 0)
})

test_that("resolve_install_group rejects unknown groups", {
  expect_error(omicsCore:::resolve_install_group("nope"))
})

test_that("check_install returns a row per group", {
  res <- check_install(features = c("rnaseq", "enrichment"))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_named(res, c("group", "n_total", "n_installed", "is_ready", "missing"))
  expect_true(all(res$n_total > 0))
})

test_that("install_optional() is a no-op when ask=FALSE and group is empty", {
  # If a user already has every package in a group, we return invisibly
  # without trying to install anything. We can't easily guarantee this in
  # a sandbox, so we only test the error path: an invalid group.
  expect_error(install_optional("nope"))
})

test_that("new_artifact_registry has the expected schema", {
  reg <- omicsCore:::new_artifact_registry()
  expect_s3_class(reg, "data.frame")
  expect_named(reg, c("artifact_type", "label", "path", "created_at"))
  expect_equal(nrow(reg), 0)
})

test_that("register_artifact appends a row", {
  reg <- omicsCore:::register_artifact(
    NULL, artifact_type = "plot", label = "volcano",
    path = "/tmp/volcano.pdf"
  )
  expect_equal(nrow(reg), 1)
  expect_equal(reg$artifact_type, "plot")
})

test_that("schema validators catch missing columns", {
  expect_error(omicsCore:::check_diff_result_schema(data.frame(x = 1)),
               "Missing required diff")
  expect_error(omicsCore:::check_enrich_result_schema(data.frame(x = 1)),
               "Missing required enrich")
})

test_that("new_analysis_bundle() returns a valid bundle", {
  b <- omicsCore:::new_analysis_bundle(
    analysis_name = "demo",
    results = list(diff_result_df = data.frame(x = 1))
  )
  expect_true(is_analysis_bundle(b))
  expect_equal(b$analysis_name, "demo")
  expect_s3_class(b$artifacts, "data.frame")
})
