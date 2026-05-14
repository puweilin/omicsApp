# Tiny synthetic xlsx workbook for Import-view tests.
#
# Produces a 3-sheet workbook (expression / metadata / feature_annot)
# that `omicsCore::read_omics()` parses without warnings under
# `omics_type = "proteomics"`. Used by `test-mod_import_view-server.R`
# and reusable from any future view test that needs a real `omics_input`.

write_tiny_omics_xlsx <- function(path,
                                  n_features = 5L,
                                  n_samples  = 6L,
                                  seed       = 2031L) {
  set.seed(seed)
  feat_ids <- sprintf("P%04d", seq_len(n_features))
  samp_ids <- sprintf("S%02d", seq_len(n_samples))

  expr <- matrix(
    rnorm(n_features * n_samples, mean = 18, sd = 1.2),
    nrow = n_features,
    dimnames = list(feat_ids, samp_ids)
  )
  expr_df <- data.frame(
    feature_id = feat_ids,
    expr,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  meta_df <- data.frame(
    sample_id = samp_ids,
    group     = rep(c("G1", "G2"), length.out = n_samples),
    age       = seq.int(30L, 30L + n_samples - 1L),
    stringsAsFactors = FALSE
  )
  feat_df <- data.frame(
    feature_id     = feat_ids,
    feature_symbol = paste0("SYM", seq_len(n_features)),
    description    = paste("desc for", feat_ids),
    stringsAsFactors = FALSE
  )
  openxlsx::write.xlsx(
    list(expression = expr_df,
         metadata   = meta_df,
         feature_annot = feat_df),
    file = path,
    overwrite = TRUE
  )
  invisible(path)
}
