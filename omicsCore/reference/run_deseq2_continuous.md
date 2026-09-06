# DESeq2 continuous-variable differential test

DESeq2 continuous-variable differential test

## Usage

``` r
run_deseq2_continuous(
  input,
  continuous_col,
  covariates = NULL,
  paired_col = NULL
)
```

## Arguments

- input:

  A validated RNA-seq `omics_input` built from raw counts.

- continuous_col:

  Continuous metadata column.

- covariates:

  Optional covariate column names.

- paired_col:

  Optional pairing column.

## Value

List with `results_raw`, `results_std`, `model_object` (`DESeqDataSet`),
and `analysis_info`.
