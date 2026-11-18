# Curves and Their Differences as Data Frames

The results of
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
and
[`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
hold their tables in the layout of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md):
`means` (one row per point of a curve or surface), `slopes`
(`semipred(deriv = )`) and, with `contr`, `pairs` (the difference of two
curves at each covariate value), with the same column names (`estimate`,
`std.error`, `df`, `conf.low`, `conf.high`, ...) and the covariates
under their own names. `tidy()` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
them,
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns the plot and [`print()`](https://rdrr.io/r/base/print.html)
shows a short summary.

## Usage

``` r
# S3 method for class 'pm_curves'
tidy(x, what = c("means", "pairs", "slopes"), ...)

# S3 method for class 'pm_curves'
as.data.frame(x, row.names = NULL, optional = FALSE,
              what = c("means", "pairs", "slopes"), ...)

# S3 method for class 'pm_curves'
autoplot(object, ...)

# S3 method for class 'pm_curves'
print(x, n = 10, digits = 4, ...)
```

## Arguments

- x, object:

  A result of
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  or
  [`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md).

- what:

  Which table: `"means"`, `"pairs"` or `"slopes"`. The default is the
  one the result has (`"pairs"` for a result with `contr`, `"slopes"`
  for `semipred(deriv = 1)`).

- row.names:

  Optional row names for the result.

- optional:

  Not used.

- n:

  The number of rows printed.

- digits:

  Significant digits printed.

- ...:

  Not used.

## Value

`tidy()` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html): a data
frame with unrounded values. In every table `conf.low` and `conf.high`
are pointwise limits (t with the row's `df`); with
`band = "simultaneous"` the simultaneous band is in `adj.conf.low` and
`adj.conf.high`, with `p.adj` where there is a test (Scheffe over the
covariate values for
[`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
the simulated maximum for
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md);
the adjusted p-value is below `level` exactly where the band leaves out
zero, up to the resolution of the 10000 draws for
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)).
Columns beginning `response` are the same quantities on the response
scale (back-transformed limits for means, delta-method limits for
differences).

- `means`:

  The factor columns of `modelterm`, the covariate(s), `estimate`,
  `std.error`, `df`, `conf.low`, `conf.high` (the model scale: the link,
  or the transformed response) and, when the curves are
  back-transformed, `response`, `response.std.error` (delta method),
  `response.conf.low`, `response.conf.high`.

- `slopes`:

  (`semipred(deriv = 1, 2)`) The same columns for the slopes
  (curvatures) on the model scale, with `statistic` (t) and `p.value` of
  slope = 0.

- `pairs`:

  `level1`, `level2`, the covariate(s), `estimate` (`level1` minus
  `level2` on the model scale), `std.error`, `df`, `statistic` (t),
  `p.value`, `conf.low`, `conf.high`; `ratio` or `odds.ratio` with their
  (pointwise) limits for a log or logit scale; and `response`,
  `response.std.error`, `response.conf.low`, `response.conf.high`, the
  difference of the back-transformed curves.

`settings` (a list in the result) records `modelterm`, `covariate`,
`level`, `contrast`, `band`, `df_method`, `vcov_adjust`, `vcov`,
`back_transform`, `test`
([`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md):
the joint F test that the two curves are the same, with the columns of
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md):
`df1`, `df2`, `F.ratio`, `p.value`), for
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
also `deriv`, `unconditional` and `boundary`, then `component`
(zero-inflated `glmmTMB` fits), `model_class` and `version`.
`significant` (in the result) holds the ranges of the covariate where
the band asked for leaves out zero.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  data(Oats, package="nlme")
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  cm <- covariatemeans(fm, "Variety", covariate="nitro", at=c(0, 0.3, 0.6))
  cm
#> Curves of 'Variety' along nitro
#> 95% limits; Kenward-Roger df
#>      Variety nitro estimate std.error     df conf.low conf.high
#>  Golden Rain   0.0    81.90     8.571 12.849    63.36    100.44
#>  Golden Rain   0.3   104.50     7.797  8.869    86.82    122.18
#>  Golden Rain   0.6   127.10     8.571 12.849   108.56    145.64
#>   Marvellous   0.0    90.42     8.571 12.849    71.88    108.95
#>   Marvellous   0.3   109.79     7.797  8.869    92.11    127.47
#>   Marvellous   0.6   129.17     8.571 12.849   110.63    147.70
#>      Victory   0.0    73.30     8.571 12.849    54.76     91.84
#>      Victory   0.3    97.62     7.797  8.869    79.95    115.30
#>      Victory   0.6   121.95     8.571 12.849   103.41    140.49
#> autoplot(x) draws the plot.
  as.data.frame(cm)
#>       Variety nitro  estimate std.error        df  conf.low conf.high
#> 1 Golden Rain   0.0  81.90000  8.570716 12.849059  63.36196 100.43804
#> 2 Golden Rain   0.3 104.50000  7.797494  8.869179  86.82110 122.17890
#> 3 Golden Rain   0.6 127.10000  8.570716 12.849059 108.56196 145.63804
#> 4  Marvellous   0.0  90.41667  8.570716 12.849059  71.87862 108.95471
#> 5  Marvellous   0.3 109.79167  7.797494  8.869179  92.11277 127.47056
#> 6  Marvellous   0.6 129.16667  8.570716 12.849059 110.62862 147.70471
#> 7     Victory   0.0  73.30000  8.570716 12.849059  54.76196  91.83804
#> 8     Victory   0.3  97.62500  7.797494  8.869179  79.94610 115.30390
#> 9     Victory   0.6 121.95000  8.570716 12.849059 103.41196 140.48804
  cd <- covariatemeans(fm, "Variety", covariate="nitro", contr=c("Golden Rain", "Victory"),
                       band="simultaneous")
  head(as.data.frame(cd))
#>        level1  level2      nitro estimate std.error       df statistic
#> 1 Golden Rain Victory 0.00000000 8.600000  8.684676 21.57462 0.9902500
#> 2 Golden Rain Victory 0.01224490 8.529592  8.567347 20.58273 0.9955931
#> 3 Golden Rain Victory 0.02448980 8.459184  8.453379 19.63654 1.0006867
#> 4 Golden Rain Victory 0.03673469 8.388776  8.342909 18.73710 1.0054976
#> 5 Golden Rain Victory 0.04897959 8.318367  8.236079 17.88507 1.0099912
#> 6 Golden Rain Victory 0.06122449 8.247959  8.133031 17.08076 1.0141310
#>     p.value     p.adj  conf.low conf.high adj.conf.low adj.conf.high
#> 1 0.3330287 0.6191042 -9.431527  26.63153    -14.22282      31.42282
#> 2 0.3310091 0.6162890 -9.309200  26.36838    -14.06393      31.12311
#> 3 0.3291509 0.6136450 -9.195202  26.11357    -13.91598      30.83435
#> 4 0.3274667 0.6111909 -9.089732  25.86728    -13.77922      30.55677
#> 5 0.3259693 0.6089456 -8.992966  25.62970    -13.65382      30.29056
#> 6 0.3246712 0.6069286 -8.905059  25.40098    -13.53996      30.03588
  if (requireNamespace("ggplot2", quietly = TRUE)) ggplot2::autoplot(cd)
```
