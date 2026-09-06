# Plot colour palette

Named colours used by every `plot_*()` function, mirroring the design
tokens in the Shiny front end so a figure on screen and the same figure
in an exported report are tinted identically.

## Usage

``` r
omics_colors
```

## Format

A named list of hex colour strings.

- up, down, ns:

  Direction of change, and features that reach no threshold.

- fg_dark, border:

  Axis text and rule colours.

- scale_low, scale_high:

  Endpoints of continuous colour scales, low to high significance.

- conc_up_up, conc_down_down, conc_up_down, conc_down_up:

  The four concordance quadrants of a two-omics comparison.

- shared, unique\_:

  Pathways found in both layers versus one.

## Examples

``` r
omics_colors$up
#> [1] "#C0392B"
```
