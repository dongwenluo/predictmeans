# Slopes of a Covariate by Group

The slope (rate of change) of the fitted response in a numeric covariate
within each level of a model term, with standard errors, degrees of
freedom, confidence limits and tests, and pairwise comparisons of the
slopes, as
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html)
and SAS `ESTIMATE` statements on slope terms.

## Usage

``` r
pm_slopes(model, covariate, modelterm, at = NULL, level = 0.05, adj = "none",
          df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-SAS"),
          vcov_adjust = TRUE,
          weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL)
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

- level:

  Significance level for the confidence limits, the default is 0.05.

- adj:

  Adjustment of the p-values of the pairwise comparisons of slopes, as
  in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- df_method, vcov_adjust, weights, vcov:

  As in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

## Value

A list of two data frames:

- slopes:

  one row per level of `modelterm`: the factor columns, `slope`,
  `std.error`, `df`, `conf.low`, `conf.high`, `statistic` (t) and
  `p.value` (of slope = 0).

- pairs:

  one row per pair of levels: `level1`, `level2`, `estimate` (slope of
  `level1` minus slope of `level2`), `std.error`, `df`, `statistic`,
  `p.value` and `p.adj`.

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
for `lmer` models fitted by REML Kenward-Roger (with the adjusted
covariance when `vcov_adjust = TRUE`), for `lme` the smallest nlme df of
the terms with the covariate, otherwise the residual df.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats <- as.data.frame(nlme::Oats)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  pm_slopes(fm, "nitro", "Variety")
#> Slopes of ‘nitro’ at 0.3 
#> 
#>      Variety slope std.error df conf.low conf.high statistic   p.value
#>  Golden Rain 75.33     11.86 51    51.53     99.14     6.353 5.726e-08
#>   Marvellous 64.58     11.86 51    40.78     88.39     5.446 1.490e-06
#>      Victory 81.08     11.86 51    57.28    104.89     6.838 9.828e-09
#> 
#> Pairwise differences of slopes (p.adj: none)
#> 
#>       level1     level2 estimate std.error df statistic p.value  p.adj
#>  Golden Rain Marvellous    10.75     16.77 51    0.6410  0.5244 0.5244
#>  Golden Rain    Victory    -5.75     16.77 51   -0.3429  0.7331 0.7331
#>   Marvellous    Victory   -16.50     16.77 51   -0.9839  0.3298 0.3298
```
