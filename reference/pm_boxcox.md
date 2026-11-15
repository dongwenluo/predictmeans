# Box-Cox and Yeo-Johnson Transformations of a Response

Power transformations to use for the response in a model formula, e.g.
`lm(pm_boxcox(y, 0.25) ~ A)`.
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
and
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
recognise them (as they recognise `log(y)`, `log(y + c)`, `log10(y)`,
`log2(y)`, `log1p(y)`, `sqrt(y)`, `y^p`, `1/y`, `qlogis(y)` and
`asin(sqrt(y))`) and give the means back-transformed to the original
scale as well; `bias_adjust = TRUE` there gives the mean rather than the
median of the original scale.

## Usage

``` r
pm_boxcox(y, lambda)
pm_yeojohnson(y, lambda)
```

## Arguments

- y:

  A numeric vector; positive for `pm_boxcox()`, any values for
  `pm_yeojohnson()`.

- lambda:

  The power: one number. Box-Cox: 1 leaves the shape of `y` as it is,
  0.5 is like a square root, 0 is the logarithm, -1 like a reciprocal.
  [`pm_check`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md)
  reports the value the data suggest for an `lm` fit.

## Value

The transformed values.

## Details

Box-Cox: \\(y^\lambda - 1)/\lambda\\, and \\\log y\\ for \\\lambda =
0\\. Yeo-Johnson (2000): the same for \\y + 1\\ when \\y \ge 0\\, and
\\-((1 - y)^{2 - \lambda} - 1)/(2 - \lambda)\\ (\\-\log(1 - y)\\ for
\\\lambda = 2\\) when \\y \< 0\\, so zero and negative values are
allowed.

Choose a transformation from the residuals of the model, not from the
raw response, and prefer a generalized linear (mixed) model for counts
and proportions.

## References

Box, G. E. P. and Cox, D. R. (1964), *An analysis of transformations*,
Journal of the Royal Statistical Society B, 26(2), 211-252.

Yeo, I.-K. and Johnson, R. A. (2000), *A new family of power
transformations to improve normality or symmetry*, Biometrika, 87(4),
954-959.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lm(pm_boxcox(yield, 0.5) ~ Block + nitro, data = Oats)
  predictmeans(fm, "nitro", plot = FALSE)$means
#> The response is pm_boxcox(yield, 0.5): the means are also given back-transformed to the original scale (trans = I for the model scale only).
#>   nitro estimate std.error df conf.low conf.high comparison.low comparison.high
#> 1     0 15.69863 0.3751289 63 14.94899  16.44826       15.16855        16.22870
#> 2   0.2 17.77431 0.3751289 63 17.02468  18.52395       17.24424        18.30438
#> 3   0.4 19.28102 0.3751289 63 18.53139  20.03066       18.75095        19.81110
#> 4   0.6 20.12632 0.3751289 63 19.37668  20.87595       19.59624        20.65639
#>    response response.std.error response.conf.low response.conf.high  n
#> 1  78.31034           3.319633          71.81707           85.08459 18
#> 2  97.75584           3.708958          90.48457          105.30809 18
#> 3 113.22048           3.991564         105.38447          121.33748 18
#> 4 122.39347           4.150111         114.24062          130.82729 18
#>   estimable
#> 1      TRUE
#> 2      TRUE
#> 3      TRUE
#> 4      TRUE
```
