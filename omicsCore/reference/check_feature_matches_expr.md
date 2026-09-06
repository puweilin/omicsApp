# Check that feature metadata matches an expression matrix

Check that feature metadata matches an expression matrix

## Usage

``` r
check_feature_matches_expr(expr_mat, feature_df, feature_id_col = "feature_id")
```

## Arguments

- expr_mat:

  Expression matrix.

- feature_df:

  Feature metadata.

- feature_id_col:

  Feature identifier column in `feature_df`.

## Value

Invisibly returns `TRUE` on success.
