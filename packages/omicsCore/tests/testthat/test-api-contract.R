# Every public function answers a wrong argument by name.
#
# Calling each exported function with one argument wrong at a time --
# NULL, NA, a negative number, a string, an empty vector, a list, a
# function, an object with the right class and nothing in it -- found
# 372 answers that named nothing the caller typed: "argument is of
# length zero", "$ operator is invalid for atomic vectors", "missing
# value where TRUE/FALSE needed". Those reach a user as the app's
# error notice and a script author as a stack trace into the engine.
# Now each names the argument, what it had to be, and what arrived.

skip_if_not_installed("limma")

# Messages R produces when code is handed something it never checked.
internal_error <- c(
  "subscript out of bounds", "is not subsettable", "argument is of length zero",
  "missing value where TRUE/FALSE needed", "non-numeric argument",
  "operator is invalid for atomic vectors", "object '[^']+' not found",
  "arguments imply differing number of rows", "undefined columns selected",
  "incorrect number of dimensions", "attempt to apply non-function",
  "non-character argument", "is missing, with no default", "cannot coerce",
  "unused argument", "no applicable method", "NA/NaN/Inf in", "invalid argument",
  "replacement has", "length of 'dimnames'", "invalid for factors",
  "invalid 'type'", "is.na\\(\\) applied", "invalid subscript", "need finite",
  "is not TRUE$"
)
is_internal <- function(msg) any(vapply(internal_error, grepl, logical(1), x = msg))

hostile_values <- list(
  `NULL` = NULL, `NA` = NA, `-1` = -1, `"x"` = "x", `""` = "",
  `character(0)` = character(0), `list()` = list(), `function` = function() NULL,
  `data.frame()` = data.frame(),
  `hollow omics_input` = structure(list(), class = "omics_input"),
  `hollow bundle` = structure(list(), class = "analysis_bundle"),
  `hollow project` = structure(list(), class = "omics_project")
)

# Calls each function with `base` and then with every formal replaced by
# every hostile value in turn. Each answer is a success or an error
# whose message is not one of R's own.
sweep_contract <- function(calls, values = hostile_values) {
  offences <- character(0)
  for (fn_name in names(calls)) {
    fn <- get(fn_name, envir = asNamespace("omicsCore"))
    base <- calls[[fn_name]]
    ok <- tryCatch({ suppressWarnings(suppressMessages(do.call(fn, base))); TRUE },
                   error = function(e) conditionMessage(e))
    if (!isTRUE(ok)) {
      offences <- c(offences, sprintf("%s: base call failed: %s", fn_name, ok))
      next
    }
    for (arg in setdiff(names(formals(fn)), "...")) {
      for (v in names(values)) {
        args <- base
        args[arg] <- list(values[[v]])
        msg <- tryCatch({ suppressWarnings(suppressMessages(do.call(fn, args))); NULL },
                        error = function(e) conditionMessage(e))
        if (!is.null(msg) && is_internal(msg)) {
          offences <- c(offences, sprintf("%s(%s = %s): %s", fn_name, arg, v,
                                          gsub("\n", " ", msg)))
        }
      }
    }
  }
  offences
}

prot <- realistic_input("proteomics")
prot$meta_df$age <- seq_len(ncol(prot$expr_mat)) * 3 + 20
rna <- realistic_input("rnaseq")
proj <- omics_project("p", experiments = list(proteomics = prot, rnaseq = rna))
db <- run_diff(prot, method = "limma", analysis_type = "group", group_col = "group",
               control_group = "G1", case_group = "G2")
qcb <- run_qc(prot)
proj$bundles <- list(diff = db)
sheet <- data.frame(feature_id = rownames(prot$expr_mat)[1:5], S01 = 1:5, S02 = 2:6)
raw_input <- omics_input(2^prot$expr_mat, prot$meta_df, prot$feature_df,
                         omics_type = "proteomics", assay_type = "raw_intensity")

light_calls <- list(
  omics_input = list(expr_mat = prot$expr_mat, meta_df = prot$meta_df,
                     feature_df = prot$feature_df, omics_type = "proteomics",
                     assay_type = "normalized_intensity"),
  validate_omics_input = list(x = prot),
  summarize_omics = list(x = prot),
  subset_omics = list(input = prot, samples = colnames(prot$expr_mat)[1:4],
                      features = rownames(prot$expr_mat)[1:10]),
  subset_omics_samples = list(omics_input = prot, sample_ids = colnames(prot$expr_mat)[1:4]),
  subset_omics_features = list(omics_input = prot, feature_ids = rownames(prot$expr_mat)[1:10]),
  select_complete_cases = list(omics_input = prot, feature_missing_cutoff = 0.5),
  drop_meta_na = list(omics_input = prot, cols = "group"),
  run_qc = list(input = prot),
  qc_missingness = list(input = prot),
  qc_outliers = list(input = prot),
  qc_depth = list(input = rna),
  qc_depth_outliers = list(depth_df = qc_depth(rna)),
  impute_matrix = list(mat = prot$expr_mat, method = "min"),
  winsorize_counts = list(count_mat = rna$expr_mat, k = 3),
  normalize_omics = list(input = raw_input, method = "log2"),
  infer_assay_type = list(expr_mat = prot$expr_mat, omics_type = "proteomics"),
  check_assay_scale = list(x = prot),
  applicable_diff_methods = list(input = prot, analysis_type = "group"),
  run_diff = list(input = prot, method = "ttest", analysis_type = "group",
                  group_col = "group", control_group = "G1", case_group = "G2"),
  run_diff_continuous = list(input = prot, method = "lm", continuous_col = "age"),
  filter_diff_results = list(result_df = db$results$diff_result_df, p_cutoff = 0.05),
  make_ranked_features = list(result_df = db$results$diff_result_df),
  plot_volcano = list(bundle = db),
  plot_ma = list(bundle = db),
  plot_pca = list(input = prot, color_by = "group"),
  plot_qc = list(bundle = qcb, view = "pca"),
  plot_heatmap = list(x = db, input = prot, n_top = 20),
  plot_feature_expression = list(input = prot, features = rownames(prot$expr_mat)[1:2],
                                 group_by = "group"),
  export_script = list(project = proj),
  omics_project = list(name = "p", experiments = list(proteomics = prot)),
  add_experiment = list(project = proj, name = "x", input = prot),
  remove_experiment = list(project = proj, name = "proteomics"),
  experiment_tags = list(project = proj),
  derive_sample_link = list(project = proj),
  suggest_sample_link = list(project = proj, tag_a = "proteomics", tag_b = "rnaseq"),
  sample_pairing_preview = list(project = proj, tag_a = "proteomics", tag_b = "rnaseq"),
  classify_sheet_role = list(df = sheet, name = "expression"),
  detect_orientation = list(df = sheet),
  detect_id_columns = list(df = sheet),
  new_analysis_bundle = list(analysis_name = "x"),
  new_import_report = list(),
  import_report_sheets = list(report = new_import_report()),
  import_report_warnings = list(report = new_import_report()),
  map_ensembl_symbols = list(ids = "ENSG00000141510"),
  resolve_impute_method = list(omics_type = "proteomics"),
  theme_omicsCore = list(),
  check_install = list(features = "proteomics"),
  is_omics_input = list(x = prot),
  is_analysis_bundle = list(x = db),
  is_omics_project = list(x = proj),
  is_import_report = list(x = new_import_report())
)

test_that("no public function answers a wrong argument with an internal error", {
  offences <- sweep_contract(light_calls)
  expect_identical(offences, character(0))
})

test_that("the file-writing functions hold the same line", {
  skip_if_not_installed("openxlsx")
  dir <- withr::local_tempdir()
  xlsx <- file.path(dir, "prot.xlsx")
  openxlsx::write.xlsx(list(
    expression = data.frame(feature_id = rownames(prot$expr_mat), prot$expr_mat,
                            check.names = FALSE),
    metadata = data.frame(sample_id = rownames(prot$meta_df), prot$meta_df),
    features = prot$feature_df), xlsx)
  omp <- file.path(dir, "p.omp")
  save_project(proj, omp)
  offences <- sweep_contract(list(
    read_omics = list(path = xlsx, omics_type = "proteomics",
                      assay_type = "normalized_intensity"),
    export_bundle = list(bundle = db, dir = file.path(dir, "bundle"), formats = "tsv"),
    save_project = list(project = proj, path = file.path(dir, "q.omp"), overwrite = TRUE),
    load_project = list(path = omp)
  ))
  expect_identical(offences, character(0))
})

# The engines that take seconds per call, and every argument of every
# function. Run with OMICSCORE_FUZZ_TESTS=1.
test_that("the heavy engines hold the same line", {
  skip_if(!nzchar(Sys.getenv("OMICSCORE_FUZZ_TESTS", "")),
          "set OMICSCORE_FUZZ_TESTS=1 to sweep the heavy engines")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("edgeR")
  dbe <- run_diff(rna, method = "edger", analysis_type = "group", group_col = "group",
                  control_group = "G1", case_group = "G2")
  enr <- run_enrichment(db, type = "ora", database = "hallmark")
  intb <- run_integration(proj, method = "concordance",
                          experiments = c("proteomics", "rnaseq"),
                          diff_bundles = list(proteomics = db, rnaseq = dbe))
  heavy <- list(
    run_enrichment = list(diff_bundle = db, type = "ora", database = "hallmark"),
    filter_enrich_results = list(enrich_df = enr$results$enrich_result_df, p_cutoff = 0.05),
    plot_enrichment = list(bundle = enr),
    run_integration = list(project = proj, method = "concordance",
                           experiments = c("proteomics", "rnaseq"),
                           diff_bundles = list(proteomics = db, rnaseq = dbe)),
    plot_integration = list(bundle = intb),
    list_gene_sets = list(database = "hallmark"),
    geneset_cache_status = list()
  )
  if (requireNamespace("GSVA", quietly = TRUE)) {
    gsv <- run_gsva(prot, database = "hallmark")
    heavy$run_gsva <- list(input = prot, database = "hallmark")
    heavy$plot_gsva_heatmap <- list(bundle = gsv, top_n = 10)
  }
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    heavy$export_report <- list(project = proj, path = tempfile(fileext = ".html"),
                                overwrite = TRUE)
  }
  expect_identical(sweep_contract(heavy), character(0))
})

# ---- what is refused, and in which words -------------------------------

test_that("a threshold that is not a probability is refused, not applied", {
  df <- db$results$diff_result_df
  expect_error(filter_diff_results(df, p_cutoff = NA),
               "`p_cutoff` must be a single number between 0 and 1, not NA.", fixed = TRUE)
  expect_error(filter_diff_results(df, p_cutoff = -1), "not -1", fixed = TRUE)
  expect_error(filter_diff_results(df, p_cutoff = "0.05"), "not '0.05'", fixed = TRUE)
  expect_error(filter_diff_results(df, effect_cutoff = -2), "`effect_cutoff` must be a single number of at least 0")
  expect_error(run_qc(prot, missing_threshold = "half"),
               "`missing_threshold` must be a single number between 0 and 1, not 'half'.", fixed = TRUE)
  expect_error(run_qc(prot, outlier_sd_threshold = list()), "not a list of length 0", fixed = TRUE)
  expect_error(qc_missingness(prot, feature_missing_cutoff = 2), "not 2", fixed = TRUE)
})

test_that("a count that is not a whole number is refused", {
  expect_error(plot_volcano(db, top_n = -1), "`top_n` must be a single whole number of at least 0, not -1.", fixed = TRUE)
  expect_error(plot_volcano(db, top_n = 2.5), "not 2.5", fixed = TRUE)
  expect_error(plot_heatmap(db, prot, n_top = 0), "at least 1", fixed = TRUE)
  expect_error(detect_id_columns(sheet, max_check = 0), "`max_check`", fixed = TRUE)
})

test_that("a name that is not a string is refused", {
  expect_error(run_diff(prot, method = "ttest", group_col = list("group"),
                        control_group = "G1", case_group = "G2"),
               "`group_col` must be a single non-empty string, not a list of length 1.", fixed = TRUE)
  expect_error(run_diff(prot, method = "ttest", group_col = "group",
                        control_group = "G1", case_group = "G2", covariates = function() NULL),
               "`covariates` must be a character vector of one or more non-empty names, not a function.", fixed = TRUE)
  expect_error(remove_experiment(proj, NULL), "`name` must be a single non-empty string, not NULL.", fixed = TRUE)
  expect_error(drop_meta_na(prot, data.frame()), "`cols`", fixed = TRUE)
  expect_error(make_ranked_features(db$results$diff_result_df, rank_col = character(0)),
               "not an empty character vector", fixed = TRUE)
})

test_that("a level may be a number, because a metadata column may hold one", {
  i <- prot
  i$meta_df$arm <- ifelse(i$meta_df$group == "G1", 1L, 2L)
  b <- run_diff(i, method = "ttest", group_col = "arm", control_group = 1, case_group = 2)
  expect_s3_class(b, "analysis_bundle")
})

test_that("an empty covariate selection means none", {
  b <- run_diff(prot, method = "limma", group_col = "group", control_group = "G1",
                case_group = "G2", covariates = character(0))
  expect_null(b$params$covariates)
})

test_that("a flag that is not TRUE or FALSE is refused", {
  expect_error(export_script(proj, include_plots = "yes"), "`include_plots` must be TRUE or FALSE, not 'yes'.", fixed = TRUE)
  expect_error(plot_heatmap(db, prot, cluster_rows = NA), "`cluster_rows` must be TRUE or FALSE, not NA.", fixed = TRUE)
})

test_that("an object of the wrong kind is described, not dereferenced", {
  expect_error(qc_depth(1), "omics_input")
  expect_error(qc_depth_outliers(1), "`depth_df` must be a data.frame, not 1.", fixed = TRUE)
  expect_error(winsorize_counts("x"), "`count_mat` must be a numeric matrix, not 'x'.", fixed = TRUE)
  expect_error(derive_sample_link(list()), "`project` must be an `omics_project`, not a list of length 0.", fixed = TRUE)
  expect_error(plot_heatmap(NA), "`x` must be an `omics_input` or an analysis_bundle from run_diff(), not NA.", fixed = TRUE)
  expect_error(import_report_sheets(db), "`report` must be an `ImportReport`, not an object of class analysis_bundle.", fixed = TRUE)
})

test_that("gene-set sizes must be whole, positive and ordered", {
  skip_if_not_installed("clusterProfiler")
  expect_error(run_enrichment(db, min_size = 0), "`min_size` must be a single whole number of at least 1")
  expect_error(run_enrichment(db, min_size = 100, max_size = 10), "`min_size` must not exceed `max_size`.", fixed = TRUE)
  expect_error(run_enrichment(db, p_adjust_method = "bh"), "`p_adjust_method` must be one of")
})

test_that("describe_value names what arrived", {
  expect_identical(describe_value(NULL), "NULL")
  expect_identical(describe_value(NA), "NA")
  expect_identical(describe_value("x"), "'x'")
  expect_identical(describe_value(2.5), "2.5")
  expect_identical(describe_value(character(0)), "an empty character vector")
  expect_identical(describe_value(1:3), "a integer vector of length 3")
  expect_identical(describe_value(list(1)), "a list of length 1")
  expect_identical(describe_value(data.frame(a = 1)), "a data.frame with 1 row(s)")
  expect_identical(describe_value(matrix(1:4, 2)), "a 2 x 2 matrix")
  expect_identical(describe_value(mean), "a function")
  expect_identical(describe_value(db), "an object of class analysis_bundle")
})
