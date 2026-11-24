# Predicted Means and Pairwise Comparisons as Data Frames

Return the `means` or `pairs` table of a
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
result (`glance()`: one row describing it). The tables have fixed column
names (see the Value section of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md))
and unrounded values, which makes them easy to use in further analysis,
reports, Shiny apps or AI tools. `tidy()` is the generic from package
generics (also used by broom);
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) needs no
other package. [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) return the means (on the
model scale, as in the `estimate` column) and their covariance matrix
(Kenward-Roger adjusted when the means are), named by the mean labels.
With them, `multcomp::glht()` or
[`marginaleffects::hypotheses()`](https://rdrr.io/pkg/marginaleffects/man/hypotheses.html)
can test any contrast of the means, e.g.
`hypotheses(pm, hypothesis = "n2 - n0 = 0")` (normal-based). Recent
versions of marginaleffects need the class listed first:
`options(marginaleffects_model_classes = "predictmeans")`.

## Usage

``` r
# S3 method for class 'predictmeans'
tidy(x, what = c("means", "pairs"), ...)

# S3 method for class 'predictmeans'
glance(x, ...)

# S3 method for class 'predictmeans'
as.data.frame(x, row.names = NULL, optional = FALSE,
              what = c("means", "pairs"), ...)

# S3 method for class 'predictmeans'
coef(object, ...)

# S3 method for class 'predictmeans'
vcov(object, ...)
```

## Arguments

- x:

  A result of
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- object:

  A result of
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- what:

  Which table: `"means"` (the default) or `"pairs"`.

- row.names:

  Optional row names for the result.

- optional:

  Not used.

- ...:

  Not used.

## Value

A data frame; for [`coef()`](https://rdrr.io/r/stats/coef.html) a named
vector and for [`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix.
`glance()` gives one row: `modelterm`, `n.means`, `n.pairs`, `nobs`
(observations in the cells), `df.method`, `adjust`, `level`, `vcov`,
`vcov.adjust`, `sed.min`, `sed.mean`, `sed.max` (standard errors of the
differences), `lsd.mean`, `nperm` and `model.class`.

## See also

[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`autoplot.predictmeans`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  as.data.frame(pm)
#>    nitro     Variety  estimate std.error       df  conf.low conf.high
#> 1      0 Golden Rain  80.00000   9.10701 16.08144  60.70194  99.29806
#> 2    0.2 Golden Rain  98.50000   9.10701 16.08144  79.20194 117.79806
#> 3    0.4 Golden Rain 114.66667   9.10701 16.08144  95.36861 133.96472
#> 4    0.6 Golden Rain 124.83333   9.10701 16.08144 105.53528 144.13139
#> 5      0  Marvellous  86.66667   9.10701 16.08144  67.36861 105.96472
#> 6    0.2  Marvellous 108.50000   9.10701 16.08144  89.20194 127.79806
#> 7    0.4  Marvellous 117.16667   9.10701 16.08144  97.86861 136.46472
#> 8    0.6  Marvellous 126.83333   9.10701 16.08144 107.53528 146.13139
#> 9      0     Victory  71.50000   9.10701 16.08144  52.20194  90.79806
#> 10   0.2     Victory  89.66667   9.10701 16.08144  70.36861 108.96472
#> 11   0.4     Victory 110.83333   9.10701 16.08144  91.53528 130.13139
#> 12   0.6     Victory 118.50000   9.10701 16.08144  99.20194 137.79806
#>    comparison.low comparison.high group n estimable
#> 1        69.75210        90.24790   A   6      TRUE
#> 2        88.25210       108.74790   AB  6      TRUE
#> 3       104.41877       124.91456    BC 6      TRUE
#> 4       114.58544       135.08123     C 6      TRUE
#> 5        76.41877        96.91456    A  6      TRUE
#> 6        98.25210       118.74790     B 6      TRUE
#> 7       106.91877       127.41456     B 6      TRUE
#> 8       116.58544       137.08123     B 6      TRUE
#> 9        61.25210        81.74790    A  6      TRUE
#> 10       79.41877        99.91456    A  6      TRUE
#> 11      100.58544       121.08123     B 6      TRUE
#> 12      108.25210       128.74790     B 6      TRUE
  head(as.data.frame(pm, what = "pairs"))
#>            level1          level2          by  estimate std.error       df
#> 1   0:Golden Rain 0.2:Golden Rain Golden Rain -18.50000  7.682948 45.00014
#> 2   0:Golden Rain 0.4:Golden Rain Golden Rain -34.66667  7.682948 45.00014
#> 3   0:Golden Rain 0.6:Golden Rain Golden Rain -44.83333  7.682948 45.00014
#> 4 0.2:Golden Rain 0.4:Golden Rain Golden Rain -16.16667  7.682948 45.00014
#> 5 0.2:Golden Rain 0.6:Golden Rain Golden Rain -26.33333  7.682948 45.00014
#> 6 0.4:Golden Rain 0.6:Golden Rain Golden Rain -10.16667  7.682948 45.00014
#>   statistic      p.value        p.adj  conf.low   conf.high      lsd
#> 1 -2.407930 2.020357e-02 9.004248e-02 -33.97425  -3.0257499 15.47425
#> 2 -4.512157 4.576126e-05 2.600347e-04 -50.14092 -19.1924165 15.47425
#> 3 -5.835434 5.476774e-07 3.213104e-06 -60.30758 -29.3590832 15.47425
#> 4 -2.104227 4.097729e-02 1.673400e-01 -31.64092  -0.6924165 15.47425
#> 5 -3.427504 1.312508e-03 6.939258e-03 -41.80758 -10.8590832 15.47425
#> 6 -1.323277 1.924277e-01 5.532559e-01 -25.64092   5.3075835 15.47425
#>   adj.conf.low adj.conf.high p.method
#> 1    -38.99579      1.995794    tukey
#> 2    -55.16246    -14.170873    tukey
#> 3    -65.32913    -24.337540    tukey
#> 4    -36.66246      4.329127    tukey
#> 5    -46.82913     -5.837540    tukey
#> 6    -30.66246     10.329127    tukey
  if (requireNamespace("generics", quietly = TRUE)) generics::glance(pm)
#>       modelterm n.means n.pairs nobs     df.method adjust level  vcov
#> 1 nitro:Variety      12      66   72 Kenward-Roger  tukey  0.05 model
#>   vcov.adjust  sed.min sed.mean sed.max lsd.mean nperm model.class
#> 1        TRUE 7.682948 9.160819 9.71502 18.64525    NA     lmerMod

  # Any contrast of the means from coef() and vcov(): Victory - Golden Rain at nitro 0.6
  b <- coef(pm)
  L <- (names(b) == "0.6:Victory") - (names(b) == "0.6:Golden Rain")
  c(estimate = sum(L*b), SE = sqrt(drop(L %*% vcov(pm) %*% L)))
#>  estimate        SE 
#> -6.333333  9.715020 
```
