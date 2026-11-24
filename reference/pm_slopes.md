# Slopes of a Covariate by Group

The slope (rate of change) of the fitted response in a numeric covariate
within each level of a model term, with standard errors, degrees of
freedom, confidence limits and tests, and pairwise comparisons of the
slopes, as
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html)
and SAS `ESTIMATE` statements on slope terms.

## Usage

``` r
pm_slopes(model, covariate, modelterm, at = NULL, conf_level = 0.95, adj = "none",
          df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment",
            "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
          vcov_adjust = TRUE,
          weights = c("equal", "proportional", "outer", "cells", "flat"),
          vcov = NULL, ...)
```

## Arguments

- model:

  A fitted model, as for
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- covariate:

  Name (in "quotes") of the numeric covariate.

- modelterm:

  Name (in "quotes") of the factor term whose levels get a slope each,
  e.g. `"Variety"`.

- at:

  The value of the covariate where the slopes are taken (matters only
  when the covariate enters the model non-linearly); the default is its
  mean.

- conf_level:

  Confidence level of the limits, the default is 0.95 (`level = 0.05`
  before 2.0).

- adj:

  Adjustment of the p-values of the pairwise comparisons of slopes, as
  in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- df_method, vcov_adjust, weights, vcov:

  As in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- ...:

  Not used: `level` (before 2.0) stops with its 2.0 spelling.

## Value

An object of class `"pm_slopes"`, a list in the layout of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
(`tidy()` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
its tables, `what = "slopes"` or `"pairs"`):

- slopes:

  one row per level of `modelterm`: the factor columns, `estimate` (the
  slope), `std.error`, `df`, `statistic` (t), `p.value` (of slope = 0),
  `conf.low` and `conf.high`.

- pairs:

  one row per pair of levels: `level1`, `level2`, `estimate` (slope of
  `level1` minus slope of `level2`), `std.error`, `df`, `statistic`,
  `p.value`, `p.adj`, `conf.low`, `conf.high` and, with the adjustments
  that have them (tukey, bonferroni, sidak, scheffe, mvt),
  `adj.conf.low`, `adj.conf.high`: simultaneous limits that leave out 0
  exactly when `p.adj` is below `1 - conf_level`.

- settings:

  `modelterm`, `covariate`, `at`, `level`, `adjust`, `df_method`,
  `vcov_adjust`, `vcov`, `model_class`, `version`.

## Details

The slopes are on the scale of the linear predictor. They come from the
reference grid of
[`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md):
the change in the K matrix between `at - h` and `at + h`, divided by
`2h` with `h` = 0.0001 times the standard deviation of the covariate.
This is exact when the covariate enters the model linearly (also through
interactions) and a central difference for
[`poly()`](https://rdrr.io/r/stats/poly.html), splines and smooths,
whose slope depends on `at`. Standard errors and degrees of freedom
follow
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md):
for models fitted by REML (`lmer`, `lme`, `gls`) Kenward-Roger (with the
adjusted covariance when `vcov_adjust = TRUE`); with one df per term
(containment, between-within, nlme's df) the smallest df of the terms
with the covariate; otherwise the residual df.

## See also

[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats <- as.data.frame(nlme::Oats)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  sl <- pm_slopes(fm, "nitro", "Variety", adj = "tukey")
  sl
#> Slopes of ‘nitro’ at 0.3 (95% limits; Kenward-Roger df) 
#> 
#>      Variety estimate std.error df statistic   p.value conf.low conf.high
#>  Golden Rain    75.33     11.86 51     6.353 5.726e-08    51.53     99.14
#>   Marvellous    64.58     11.86 51     5.446 1.490e-06    40.78     88.39
#>      Victory    81.08     11.86 51     6.838 9.828e-09    57.28    104.89
#> 
#> Pairwise differences of slopes (p.adj: tukey)
#> 
#>       level1     level2 estimate std.error df statistic p.value  p.adj conf.low
#>  Golden Rain Marvellous    10.75     16.77 51    0.6410  0.5244 0.7982   -22.92
#>  Golden Rain    Victory    -5.75     16.77 51   -0.3429  0.7331 0.9373   -39.42
#>   Marvellous    Victory   -16.50     16.77 51   -0.9839  0.3298 0.5903   -50.17
#>  conf.high adj.conf.low adj.conf.high
#>      44.42       -29.73         51.23
#>      27.92       -46.23         34.73
#>      17.17       -56.98         23.98
  as.data.frame(sl, what = "pairs")
#>        level1     level2 estimate std.error       df  statistic   p.value
#> 1 Golden Rain Marvellous    10.75  16.77051 50.99999  0.6410061 0.5243889
#> 2 Golden Rain    Victory    -5.75  16.77051 50.99999 -0.3428637 0.7331102
#> 3  Marvellous    Victory   -16.50  16.77051 50.99999 -0.9838699 0.3298242
#>       p.adj conf.low conf.high adj.conf.low adj.conf.high
#> 1 0.7982489 -22.9182   44.4182    -29.73367      51.23367
#> 2 0.9373365 -39.4182   27.9182    -46.23367      34.73367
#> 3 0.5903123 -50.1682   17.1682    -56.98367      23.98367
```
