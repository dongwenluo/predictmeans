# Variance Components with Standard Errors and Confidence Intervals

The variance components of a fitted model, with standard errors and
confidence intervals, for every model class that
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
handles: the random-effect variances and covariances and the residual
variance of `lmer`, `lme` and `glmmTMB` models (also with correlation
and variance structures and random slopes), the random-effect
(co)variances of `glmer` models, the correlation and variance-function
parameters of `gls` and `lme` models, the stratum variance components of
`aov` models with `Error()`, the smoothing variance components of
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) models and
[`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
fits, and the error term of `lm`, `aov` and `glm` models (residual
variance, dispersion, or the negative binomial \\\theta\\).

## Usage

``` r
varcomp(model, ci = TRUE, level = 0.95, ci_method = c("Satterthwaite", "Wald", "profile"),
        information = c("observed", "average", "expected"))
```

## Arguments

- model:

  A fitted model (see Description).

- ci:

  `TRUE` (the default): add confidence intervals.

- level:

  Confidence level of the intervals; the default is 0.95.

- ci_method:

  How the intervals of variances are made: `"Satterthwaite"` (the
  default) chi-square intervals with \\\nu = 2 (estimate/SE)^2\\ df, as
  SAS PROC MIXED (`COVTEST CL`); `"Wald"` Wald intervals on the log
  scale; `"profile"` (`lmer`, `glmer`) lme4's profile-likelihood
  intervals, which for a REML fit are those of the ML fit (lme4 profiles
  the ML deviance), so they do not go with the REML estimates. The other
  components always get the Wald intervals of the Details.

- information:

  For REML fits of `lmer`, `lme` and `gls` models: the information
  matrix of the variance parameters behind the SEs: `"observed"` (the
  default, as SAS PROC MIXED), `"average"` (the average information of
  Genstat's and ASReml's REML: Genstat's s.e.) or `"expected"`; the same
  choices as for the Kenward-Roger df (`df_method` of
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).
  Other fits use the observed information.

## Value

A data frame with one row per component (random-effect variances and
covariances by group, named `group.term` and `group.term1.term2`;
correlation and variance-function parameters named `cor.name` and
`var.name`; `residual`) and columns `vcov` (the estimate), `SE` and,
with `ci = TRUE`, the limits. Values are rounded to 4 decimals;
attribute `"ci_method"` says how the intervals were made.

## Details

Standard errors are from the observed information of the variance
parameters (the inverse Hessian of the REML, or ML, deviance at the
estimates, as SAS PROC MIXED and the Kenward-Roger df of predictmeans),
carried to the variance components by the delta method: `lmer` from the
deviance in (theta, sigma), `glmer` from the Hessian lme4 stores for
(theta, beta) (`nAGQ = 0` fits: the deviance in theta), `lme` and `gls`
from the REML (ML) log-likelihood in their unconstrained `modelStruct`
parameters (any `pdMat`, `corStruct` and `varFunc`), `glmmTMB` from
`vcov(full = TRUE)`.

Intervals of variances: Satterthwaite's chi-square interval \\\[\nu
\hat\sigma^2/\chi^2\_{\nu, 1-\alpha/2}, \nu \hat\sigma^2/\chi^2\_{\nu,
\alpha/2}\]\\ with \\\nu = 2 (\hat\sigma^2/SE)^2\\ (positive and skewed
as the sampling distribution of a variance; for the residual variance of
an `lm` it is the exact interval, \\\nu\\ the residual df; for the
stratum components of an `aov` Satterthwaite's interval from the mean
squares). Other components: Wald intervals on the scale where the
parameter is unbounded: atanh for correlations (a covariance: the
interval of its correlation times the two SDs), log for other positive
parameters (`varIdent` ratios, ranges, \\\theta\\), and the natural
scale otherwise (e.g. the power of `varPower`). A variance estimated at
its bound, zero, has no SE (the interval is (0, NA)), as SAS prints it;
an analysis of variance estimate can be negative.

`lm`, `aov`, gaussian `glm`: the residual variance \\s^2\\ with SE \\s^2
\sqrt{2/df}\\ and the exact chi-square interval. Other `glm` families
with an estimated dispersion (Gamma, inverse gaussian, quasi): the
Pearson dispersion, with the same (approximate) SE and interval. Poisson
and binomial: the dispersion is fixed at 1 (no SE). Negative binomial
([`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html)): \\\theta\\
with its SE.

`aov` with `Error()`: the analysis of variance (method of moments)
estimates of the variance of each error stratum, from the stratum mean
squares and their expectations (as Genstat's stratum variances; they
equal the REML estimates in balanced designs when positive, and can be
negative), with SEs from \\var(MS) = 2 MS^2/df\\ and Satterthwaite's
chi-square intervals; with missing values, those of the design completed
as in
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

`gam`:
[`mgcv::gam.vcomp()`](https://rdrr.io/pkg/mgcv/man/gam.vcomp.html)
(variances of the smooths, and the scale), its intervals squared, and
the SEs from their width on the log scale; it needs `method = "REML"` or
`"ML"`.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
varcomp(fm)
#>                               vcov       SE    2.5 %    97.5 %
#> Variety:Block.(Intercept) 106.0618  67.8754  40.9773  657.7703
#> Block.(Intercept)         214.4810 168.8381  70.8628 2575.3156
#> residual                  177.0831  37.3324 121.8272  280.9241
fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
varcomp(fm1)
#>                         vcov       SE    2.5 %    97.5 %
#> Block.(Intercept)   214.4745 168.8318  70.8610 2575.1691
#> Variety.(Intercept) 106.0618  67.8756  40.9772  657.7745
#> residual            177.0835  37.3325 121.8275  280.9250
# the same components by the analysis of variance
varcomp(aov(yield ~ nitro*Variety + Error(Block/Variety), data=Oats))
#>                   vcov       SE    2.5 %    97.5 %
#> Block         214.4771 168.8340  70.8619 2575.2100
#> Block:Variety 106.0618  67.8755  40.9773  657.7726
#> residual      177.0833  37.3324 121.8274  280.9246

data(Orthodont, package="nlme")
mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
varcomp(mod1)
#>                            vcov     SE   2.5 %   97.5 %
#> Subject.(Intercept)      5.4151 4.9583  1.5819 128.6996
#> Subject.age.(Intercept) -0.3211 0.4244 -0.4943   0.1570
#> Subject.age              0.0513 0.0415  0.0166   0.6844
#> residual                 1.7162 0.3303  1.2163   2.6042
# correlation and variance structures
mod2 <- gls(distance ~ age*Sex, correlation=corAR1(form=~1|Subject),
            weights=varIdent(form=~1|Sex), data=Orthodont)
varcomp(mod2)
#>              vcov     SE  2.5 %  97.5 %
#> cor.Phi    0.7473 0.0681 0.5815  0.8535
#> var.Female 0.5647 0.0908 0.4121  0.7738
#> residual   9.5035 2.8443 5.7032 18.9234
# the error term of a linear model
varcomp(lm(yield ~ Block + nitro*Variety, data=Oats))
#>              vcov      SE    2.5 %   97.5 %
#> residual 254.2192 48.4777 180.6923 384.1423
```
