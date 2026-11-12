# Effect Sizes and Model Checks from Permutations

Uses the permuted statistics kept by
[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
to give, for each term and coefficient, a standardised effect size and a
check of the usual parametric test against the permutation distribution.

## Usage

``` r
permstats(x, level = 0.05)
```

## Arguments

- x:

  A result of
  [`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  with Freedman-Lane permutations (`method = "auto"`, the default, or
  `"freedman-lane"`).

- level:

  The significance level of the check of the parametric test (default
  0.05).

## Value

A data frame with one row per term and coefficient (rows
`"coef:<name>"`):

- statistic:

  the observed statistic: F for terms, \|t\| for coefficients (robust
  ones with `robust = TRUE`).

- df1, df2:

  numerator and denominator df of the usual parametric test: the
  residual df for `lm`, `aov`, `glm` and `gls`, Satterthwaite df for
  `lmer`, nlme's df for `lme`, and `Inf` (Wald chi-square) for `glmer`
  and `glmmTMB`.

- p.param:

  the parametric p-value of the observed statistic with these df.

- p.perm:

  the permutation p-value (as in the tables of
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)).

- ses:

  the standardised effect size, (observed - mean of the permuted
  statistics) / their standard deviation (Gotelli and McCabe 2002): how
  far the observed statistic lies from what the permutations give when
  the effect is absent, in standard deviations.

- param.size:

  the proportion of permuted data sets in which the parametric test
  rejects at `level`. The permuted data sets have no effect for the
  tested term, so this is the type I error rate of the parametric test
  for data like these: close to `level` when the parametric reference
  distribution suits the data, clearly above it when the parametric test
  is too liberal here (for example with unequal variances,
  overdispersion or few groups), below it when it is conservative. With
  few permutations it is rough (its standard error is about
  `sqrt(level (1 - level) / nperm)`).

## Details

Every permuted data set comes from the model without the tested term
(Freedman-Lane), so the permuted statistics show how the statistic
behaves when that term has no effect but everything else in the data is
kept. Comparing the parametric p-value with the permutation p-value, and
`param.size` with `level`, shows whether the usual test can be trusted
for these data. For mixed models the permuted statistics are Wald
statistics (unadjusted covariance), so `param.size` checks the Wald F
with the given df, not the Kenward-Roger test.

## References

Gotelli, N. J. and McCabe, D. J. (2002). Species co-occurrence: a
meta-analysis of J. M. Diamond's assembly rules model. *Ecology*, 83,
2091-2096.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
mtcars$cyl <- factor(mtcars$cyl)
fm <- lm(mpg ~ wt + cyl, data = mtcars)
pfl <- permmodels(fm, nperm = 999, type = "II", seed = 1, prt = FALSE)
permstats(pfl)
#>           statistic df1 df2      p.param p.perm       ses param.size
#> wt        18.080056   1  28 0.0002130435  0.001 11.695665 0.04904905
#> cyl        7.285567   2  28 0.0028353022  0.005  4.622323 0.05705706
#> coef:wt    4.252065   1  28 0.0002130435  0.001  5.599425 0.04904905
#> coef:cyl6  3.070244   1  28 0.0047178339  0.009  3.291703 0.05905906
#> coef:cyl8  3.674214   1  28 0.0009991893  0.003  4.186066 0.05405405
```
