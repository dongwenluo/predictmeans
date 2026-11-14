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

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
  as.data.frame(pm)
#>    nitro     Variety  estimate std.error       df  conf.low conf.high group n
#> 1      0 Golden Rain  80.00000   9.10701 16.08144  60.70194  99.29806   A   6
#> 2      0  Marvellous  86.66667   9.10701 16.08144  67.36861 105.96472    A  6
#> 3      0     Victory  71.50000   9.10701 16.08144  52.20194  90.79806    A  6
#> 4    0.2 Golden Rain  98.50000   9.10701 16.08144  79.20194 117.79806   AB  6
#> 5    0.2  Marvellous 108.50000   9.10701 16.08144  89.20194 127.79806     B 6
#> 6    0.2     Victory  89.66667   9.10701 16.08144  70.36861 108.96472    A  6
#> 7    0.4 Golden Rain 114.66667   9.10701 16.08144  95.36861 133.96472    BC 6
#> 8    0.4  Marvellous 117.16667   9.10701 16.08144  97.86861 136.46472     B 6
#> 9    0.4     Victory 110.83333   9.10701 16.08144  91.53528 130.13139     B 6
#> 10   0.6 Golden Rain 124.83333   9.10701 16.08144 105.53528 144.13139     C 6
#> 11   0.6  Marvellous 126.83333   9.10701 16.08144 107.53528 146.13139     B 6
#> 12   0.6     Victory 118.50000   9.10701 16.08144  99.20194 137.79806     B 6
#>    estimable
#> 1       TRUE
#> 2       TRUE
#> 3       TRUE
#> 4       TRUE
#> 5       TRUE
#> 6       TRUE
#> 7       TRUE
#> 8       TRUE
#> 9       TRUE
#> 10      TRUE
#> 11      TRUE
#> 12      TRUE
  head(as.data.frame(pm, what = "pairs"))
#>          level1          level2          by   estimate std.error       df
#> 1 0:Golden Rain    0:Marvellous        <NA>  -6.666667  9.715020 30.23080
#> 2 0:Golden Rain       0:Victory        <NA>   8.500000  9.715020 30.23080
#> 3  0:Marvellous       0:Victory        <NA>  15.166667  9.715020 30.23080
#> 4 0:Golden Rain 0.2:Golden Rain Golden Rain -18.500000  7.682948 45.00014
#> 5  0:Marvellous 0.2:Golden Rain        <NA> -11.833333  9.715020 30.23080
#> 6     0:Victory 0.2:Golden Rain        <NA> -27.000000  9.715020 30.23080
#>    statistic     p.value      p.adj   conf.low conf.high      lsd adj.conf.low
#> 1 -0.6862226 0.497802546         NA -26.501035 13.167702 19.83437           NA
#> 2  0.8749338 0.388508422         NA -11.334369 28.334369 19.83437           NA
#> 3  1.5611564 0.128895695         NA  -4.667702 35.001035 19.83437           NA
#> 4 -2.4079299 0.020203566 0.09004248 -33.974250 -3.025750 15.47425    -38.99579
#> 5 -1.2180451 0.232628340         NA -31.667702  8.001035 19.83437           NA
#> 6 -2.7792016 0.009281361         NA -46.834369 -7.165631 19.83437           NA
#>   adj.conf.high p.method
#> 1            NA    tukey
#> 2            NA    tukey
#> 3            NA    tukey
#> 4      1.995794    tukey
#> 5            NA    tukey
#> 6            NA    tukey
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
