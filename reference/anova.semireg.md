# Tests for a Semi-Parametric Model

Tests of the fixed effects and of the smooth terms of a
[`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
or
[`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
fit, with the package's own code (no lmerTest needed).

## Usage

``` r
# S3 method for class 'semireg'
anova(object, type = c("III", "II", "I"), smooth = TRUE, ...)
```

## Arguments

- object:

  A
  [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  or
  [`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  fit.

- type:

  `"III"` (default), `"II"` or `"I"` (also `3`, `2`, `1`): the type of
  the tests of the fixed effects. Type I is sequential (Kenward-Roger F
  for gaussian fits, Wald chi-square otherwise); see
  [`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md)
  for which type to use.

- smooth:

  Logical: test the smooth terms (refits the model once per smooth
  block)? The default is TRUE.

- ...:

  Not used.

## Value

A list of class `"pm_semi_anova"` with `fixed` (the table of the fixed
effects) and `smooth` (a data frame: the variance of each smooth block,
the log-likelihoods of the model with and without it, the likelihood
ratio statistic, its p-value and the reference distribution), printed
together. Advice on the type of tests is kept as for
[`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md)
(see
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).

## Details

**Fixed effects**: Type III (or II) tests from `pm_anova()`: for
gaussian fits Wald F tests with Kenward-Roger df and the adjusted
covariance (the df equal those of the same model fitted by
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html)); for generalized
fits Wald chi-square tests. The fixed parts of the smooths (e.g. `x` and
`f:x` for `s(x, by = f)`) are tested here: they are the straight-line
parts of the curves.

**Smooth terms**: for each block of smooth random effects (one
variance), the likelihood ratio test of the variance being zero, i.e. of
the curve being its straight-line part: the model is refitted without
the block (with the same fixed effects, so REML likelihoods can be
compared for gaussian fits). The tests are those of
`varcomp(fit, test = TRUE)`: as the variance is on the boundary of its
range under the null hypothesis, the reference is the mixture
\\\frac12\chi^2_0 + \frac12\chi^2_1\\ (Self and Liang, 1987), and for
the only random effect of a gaussian fit the exact null distribution of
the restricted likelihood ratio (Crainiceanu and Ruppert, 2004), where
the mixture is conservative;
[`lmerTest::ranova()`](https://rdrr.io/pkg/lmerTest/man/ranova.html)
used the full chi-square, more conservative still. Column `test` names
the distribution. For an `s(x, f)` smooth the block of the deviation
curves tests whether the curves of the levels of `f` have the same
shape. For small samples
[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
or
[`bootmodels`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
on `fit$semer` give resampling tests.

## References

Crainiceanu, C. M. and Ruppert, D. (2004). Likelihood ratio tests in
linear mixed models with one variance component. *Journal of the Royal
Statistical Society B*, 66, 165-185.

Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82,
605-610.

## See also

[`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md),
[`varcomp`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
[`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)

## Examples

``` r
library(predictmeans)
data(Dialyzer, package = "nlme")
fm <- semireg(rate ~ s(pressure, QB, k = 6) + (1 | Subject), data = Dialyzer)
anova(fm)
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: rate
#> 
#> Advice:
#> ! A numeric covariate is in an interaction: Type III tests of the other
#>   terms are at covariate = 0; centre the covariate if 0 is not a meaningful
#>   value.
#> 
#>                    F Df  Df.res  Pr(>F)    
#> (Intercept)   6.5438  1 116.395 0.01181 *  
#> pressure    183.8452  1 109.641 < 2e-16 ***
#> QB            6.8948  1  33.301 0.01296 *  
#> pressure:QB 219.5102  1 114.064 < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);
#> reference: exact restricted LRT for the only random effect of a gaussian fit, otherwise 0.5 chi2(0) + 0.5 chi2(1)
#>           term variance logLik logLik.without   LRT  p.value
#>     s_pressure      176 -361.9         -368.4 13.06 0.000151
#>  s_pressure_QB        0 -361.9         -361.9  0.00 1.000000
#>                       test
#>  0.5 chi2(0) + 0.5 chi2(1)
#>  0.5 chi2(0) + 0.5 chi2(1)
```
