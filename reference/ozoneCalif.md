# Ozone Concentration and Meteorology in California

Daily ozone concentration with three meteorological variables (345
days), a version of the Los Angeles ozone data of Breiman and Friedman
(1985); the data of M. P. Wand's additive model examples (fitted by
[`lme()`](https://rdrr.io/pkg/nlme/man/lme.html) and by `gam()`).

## Usage

``` r
data(ozoneCalif)
```

## Format

A data frame with 345 rows and 4 columns:

- ozone:

  ozone concentration.

- inversion.base.height:

  inversion base height.

- daggett.pressure.gradient:

  pressure gradient (Los Angeles to Daggett).

- inversion.base.temp:

  inversion base temperature.

## Source

Breiman, L. and Friedman, J. H. (1985). Estimating optimal
transformations for multiple regression and correlation. *Journal of the
American Statistical Association*, 80, 580-598.

## Examples

``` r
data(ozoneCalif)
# an additive model: three O'Sullivan splines
fit <- semireg(ozone ~ s(inversion.base.height, bs = "os", k = 15) +
                 s(daggett.pressure.gradient, bs = "os", k = 15) +
                 s(inversion.base.temp, bs = "os", k = 15), data = ozoneCalif, prt = FALSE)
anova(fit)$smooth     # is each curve more than a straight line?
#>                          term  variance    logLik logLik.without      LRT
#> 1     s_inversion_base_height 12.286694 -1017.003      -1024.002 13.99700
#> 2 s_daggett_pressure_gradient  4.627735 -1017.003      -1034.664 35.32149
#> 3       s_inversion_base_temp  1.291126 -1017.003      -1025.642 17.27717
#>        p.value                      test
#> 1 9.155105e-05 0.5 chi2(0) + 0.5 chi2(1)
#> 2 1.397650e-09 0.5 chi2(0) + 0.5 chi2(1)
#> 3 1.615225e-05 0.5 chi2(0) + 0.5 chi2(1)
# the effect of one covariate, the others at their means
semipred(fit, covariate = "inversion.base.temp", prt = FALSE)$plt
```
