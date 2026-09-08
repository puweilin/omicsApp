suppressPackageStartupMessages(pkgload::load_all('packages/omicsCore', quiet = TRUE))
helpers <- new.env(parent = baseenv())
sys.source('packages/omicsApp/R/async_helpers.R', envir = helpers)
set.seed(708)
n_genes <- 1200L
n_samples <- 200L
groups <- rep(c('ctrl', 'case'), each = n_samples / 2)
base_mu <- exp(runif(n_genes, log(5), log(500)))
mu <- matrix(rep(base_mu, n_samples), nrow = n_genes)
mu[1:80, groups == 'case'] <- mu[1:80, groups == 'case'] * 2
counts <- matrix(rnbinom(length(mu), mu = as.vector(mu), size = 10),
                 nrow = n_genes,
                 dimnames = list(paste0('g', seq_len(n_genes)),
                                 paste0('s', seq_len(n_samples))))
metadata <- data.frame(group = groups, row.names = colnames(counts))
features <- data.frame(feature_id = rownames(counts),
                       feature_symbol = paste0('S', seq_len(n_genes)),
                       row.names = rownames(counts))
x <- omicsCore::omics_input(counts, metadata, features,
                            omics_type = 'rnaseq', assay_type = 'raw_count')
args <- list(input = x, method = 'deseq2', analysis_type = 'group',
              group_col = 'group', control_group = 'ctrl', case_group = 'case',
              keep_model = TRUE)
cat('SYNTHETIC DATA; local runtime, not server benchmark\n')
cat('genes=', n_genes, '; samples=', n_samples, '\n', sep = '')
cat(R.version.string, '\n')
cat('DESeq2=', as.character(packageVersion('DESeq2')),
    '; future=', as.character(packageVersion('future')), '\n', sep = '')
cat('Fitting detached design ...\n')
t_detached <- system.time(b_detached <- do.call(omicsCore::run_diff, args))
cat('Fitting original attached design ...\n')
t_attached <- system.time(b_attached <- testthat::with_mocked_bindings(
  do.call(omicsCore::run_diff, args),
  detached_formula = function(text) stats::as.formula(text, env = parent.frame()),
  .package = 'omicsCore'))
stopifnot(identical(b_detached$results$diff_result_df,
                    b_attached$results$diff_result_df))
stopifnot(identical(b_detached$results$diff_raw_df,
                    b_attached$results$diff_raw_df))
stopifnot(identical(environment(DESeq2::design(b_detached$results$diff_object)), baseenv()))
stopifnot(!identical(environment(DESeq2::design(b_attached$results$diff_object)), baseenv()))
cat('Standardized and raw differential tables: IDENTICAL\n')
args_default <- args
args_default$keep_model <- FALSE
b_default <- do.call(omicsCore::run_diff, args_default)
stopifnot(is.null(b_default$results$diff_object))
stopifnot(identical(b_default$results$diff_result_df, b_detached$results$diff_result_df))
stopifnot(identical(b_default$results$diff_raw_df, b_detached$results$diff_raw_df))
cat('Actual default keep_model=FALSE: model absent, both tables IDENTICAL\n')
cat('Fit elapsed: detached=', t_detached[['elapsed']],
    '; attached=', t_attached[['elapsed']], '\n', sep = '')
b_no_model <- b_attached
b_no_model$results['diff_object'] <- list(NULL)
b_light <- b_no_model
b_light$results <- list(diff_result_df = b_no_model$results$diff_result_df)
b_light$params <- b_no_model$params['comparison']
make_func <- function(bundle) helpers$detached_call(function() {
  if (identical(type, 'ora')) {
    omicsCore::run_enrichment(diff_bundle = bundle, type = type,
      database = db_arg, direction = dir_, p_cutoff = thr$p_cutoff,
      p_preference = thr$p_preference, effect_cutoff = thr$effect_cutoff,
      output_p_cutoff = 1)
  } else {
    omicsCore::run_enrichment(diff_bundle = bundle, type = type,
      database = db_arg, direction = dir_, output_p_cutoff = 1)
  }
}, bundle = bundle, type = 'ora', db_arg = 'hallmark', dir_ = 'both',
   thr = list(p_cutoff = 0.05, p_preference = 'adjusted', effect_cutoff = 1))
measure <- function(bundle, label) {
  func <- make_func(bundle)
  e <- list2env(list(func = func), parent = baseenv())
  gp <- future::getGlobalsAndPackages(
    quote({ tryCatch(func(), error = function(e) e) }), envir = e, maxSize = 8 * 1024^3)
  sz <- vapply(gp$globals, future:::objectSize, numeric(1))
  cat(label, 'globals:', paste(names(sz), collapse = ', '), '\n')
  data.frame(case = label, bundle_base_MiB = as.numeric(object.size(bundle)) / 1024^2,
             bundle_future_MiB = future:::objectSize(bundle) / 1024^2,
             closure_serialized_MiB = length(serialize(func, NULL)) / 1024^2,
             globals_total_MiB = attr(gp$globals, 'total_size') / 1024^2,
             func_MiB = unname(sz['func']) / 1024^2,
             exported_bundle_MiB = unname(sz['bundle']) / 1024^2)
}
out <- do.call(rbind, list(measure(b_attached, 'original_full_model'),
                            measure(b_detached, 'detached_full_model'),
                            measure(b_no_model, 'no_model'),
                            measure(b_light, 'only_enrichment_inputs')))
print(out, row.names = FALSE)
write.csv(out, 'payload-measurements.csv', row.names = FALSE)
cat('All probes complete.\n')
