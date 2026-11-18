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
varcomp(model, ci = TRUE, conf.level = 0.95, ci_method = c("Satterthwaite", "Wald", "profile"),
        information = c("observed", "average", "expected"), test = FALSE,
        test_method = c("auto", "mixture", "bootstrap"), nsim = 10000, nboot = 999, ncore = 1, ...)
```

## Arguments

- model:

  A fitted model (see Description).

- ci:

  `TRUE` (the default): add confidence intervals.

- conf.level:

  Confidence level of the intervals; the default is 0.95 (`level` before
  2.0).

- ...:

  Not used: `level` (before 2.0) stops with its 2.0 spelling.

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

- test:

  `TRUE`: add a likelihood ratio test of each component being zero
  (`lmer`, `glmer`, `glmmTMB` and `semireg` fits; see ‘Tests’). The
  default is `FALSE`: each test refits the model once.

- test_method:

  `"auto"` (the default) or `"mixture"`: the reference distribution of a
  single variance (see ‘Tests’); `"bootstrap"`: a parametric bootstrap
  of each likelihood ratio, for small samples, few levels, or where the
  asymptotic references are conservative.

- nsim:

  Number of draws of the exact null distribution of the restricted
  likelihood ratio (default 10000, with a fixed seed, so the p-values
  repeat).

- nboot:

  Maximum number of bootstrap samples per test (default 999); a test
  stops early once 20 simulated statistics reach the observed one.

- ncore:

  Number of cores for the bootstrap (default 1; not on Windows).

## Value

A data frame with one row per component (random-effect variances and
covariances by group, named `group.term` and `group.term1.term2`;
correlation and variance-function parameters named `cor.name` and
`var.name`; `residual`) and columns `vcov` (the estimate), `SE` and,
with `ci = TRUE`, the limits. Values are rounded to 4 decimals;
attribute `"ci_method"` says how the intervals were made. With
`test = TRUE` also `LRT` and `p.value`, and class `"pm_varcomp"`
(printed with the hypotheses).

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

## Tests

With `test = TRUE` each row gets the likelihood ratio test of the model
one would fit without it (columns `LRT` and `p.value`; the hypotheses
and reference distributions are printed below the table and kept in
`attr(, "tests")`):

|  |  |  |
|----|----|----|
| **Row** | **Model refitted** | **Reference distribution** |
| variance of a one-variance term `(1 | g)`, smooth block | the term dropped | \\\frac12\chi^2_0 + \frac12\chi^2_1\\; exact restricted LRT when it is the only random effect of a gaussian `lmer`/`semireg` REML fit |
| slope variance in `(1 + x | g)` (k components) | `(1 | g)` | \\\frac12\chi^2\_{k-1} + \frac12\chi^2_k\\ |
| intercept variance of `(1 + x | g)`, variance of `ar1()` etc. | the whole term dropped | \\\chi^2\\ with the term's number of parameters (conservative) |
| covariance in `(1 + x | g)` | `(1 + x || g)` (glmmTMB: [`diag()`](https://rdrr.io/r/base/diag.html)) | \\\chi^2\\ with the number of covariances |
| correlation of `ar1()`, `ou()`, `cs()` | `homdiag()` / [`diag()`](https://rdrr.io/r/base/diag.html) | \\\chi^2_1\\ |
| residual variance; a structure that is the residual covariance (`dispformula = ~0`) | no test |  |

A variance is on the boundary of its range under the hypothesis, so its
likelihood ratio has the mixture distributions of Self and Liang (1987)
and Stram and Lee (1994), not \\\chi^2_1\\ (which is conservative, as
[`lmerTest::ranova()`](https://rdrr.io/pkg/lmerTest/man/ranova.html)).
When the tested variance is the only random effect of a gaussian model
fitted by REML, the exact finite-sample null distribution of Crainiceanu
and Ruppert (2004) is used (simulated); with other random effects
present that distribution is no longer exact (in simulations it was too
liberal), and the mixture is used. The mixture p-values equal SAS PROC
GLIMMIX `COVTEST`, the likelihood ratios SAS PROC MIXED and Genstat REML
(`tests/testthat/test-vc-test.R`). Where a p-value may be inaccurate in
either direction (a term with fewer than 6 levels; binary data with
fewer than 5 observations per level, where the Laplace approximation is
poor), a warning names the rows and suggests the bootstrap.

`test_method = "bootstrap"`: for each test, responses are simulated from
the model without the component (new random effects and errors), both
models are refitted to each, and the p-value is \\(\\\\LRT^\* \ge
LRT\\ + 1)/(b + 1)\\. It needs no distributional approximation, so it is
the check for the rows the warning names. A test stops once 20 simulated
statistics reach the observed one (Besag and Clifford 1991): a clearly
non-significant row needs some 40-400 refits, a significant one all
`nboot` pairs. Each draw has its own seed, so the p-values repeat, with
or without `ncore`. Run time: a few seconds to a minute for usual `lmer`
and `glmer` models; `glmmTMB` fits are slower (each refit rebuilds the
model). For `lme`, `gls`, `gam` and models without random effects
`test = TRUE` is not available.

Which method: in a simulation of seven designs (`dev/sim-vc-defaults.R`
in the package sources) the default rejected at most 4.5% of the data
sets without the component at the 5% level: it was never liberal. For
GLMMs it was conservative (binary data in clusters of four 0.5%, Poisson
4%, negative binomial 0%), and there the bootstrap rejected 1.5-6% and
found 5-7% more of the real components. Use the default for gaussian
models and as a quick first look; use the bootstrap for GLMMs, small
data sets, the rows a warning names, and a p-value near the cut-off.

## References

Besag, J. and Clifford, P. (1991). Sequential Monte Carlo p-values.
*Biometrika*, 78, 301-304.

Crainiceanu, C. M. and Ruppert, D. (2004). Likelihood ratio tests in
linear mixed models with one variance component. *Journal of the Royal
Statistical Society B*, 66, 165-185.

Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82,
605-610.

Stram, D. O. and Lee, J. W. (1994). Variance components testing in the
longitudinal mixed effects model. *Biometrics*, 50, 1171-1177.

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
# is each variance needed? (one refit per row)
varcomp(fm, test = TRUE)
#>                               vcov       SE    2.5 %    97.5 %    LRT  p.value
#> Variety:Block.(Intercept) 106.0618  67.8754  40.9773  657.7703 7.6615 0.002821
#> Block.(Intercept)         214.4810 168.8381  70.8628 2575.3156 4.9782 0.012830
#> residual                  177.0831  37.3324 121.8272  280.9241     NA       NA
#> 
#> Tests (likelihood ratio, each row's model refitted):
#>   Variety:Block.(Intercept): drop (1 | Variety:Block); 0.5 chi2(0) + 0.5 chi2(1)
#>   Block.(Intercept): drop (1 | Block); 0.5 chi2(0) + 0.5 chi2(1)
# \donttest{
# the same by parametric bootstrap (a few seconds; more draws and ncore for real use)
varcomp(fm, test = TRUE, test_method = "bootstrap", nboot = 99)
#>                               vcov       SE    2.5 %    97.5 %    LRT p.value
#> Variety:Block.(Intercept) 106.0618  67.8754  40.9773  657.7703 7.6615    0.01
#> Block.(Intercept)         214.4810 168.8381  70.8628 2575.3156 4.9782    0.02
#> residual                  177.0831  37.3324 121.8272  280.9241     NA      NA
#> 
#> Tests (likelihood ratio, each row's model refitted):
#>   Variety:Block.(Intercept): drop (1 | Variety:Block); parametric bootstrap (99 draws)
#>   Block.(Intercept): drop (1 | Block); parametric bootstrap (99 draws)
# }
# a GLMM: the boundary mixture is conservative here, the bootstrap nearer 5%
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = lme4::cbpp)
varcomp(gm, test = TRUE)
#>                    vcov     SE  2.5 % 97.5 %     LRT   p.value
#> herd.(Intercept) 0.4125 0.2294 0.1757 1.8476 14.0058 9.112e-05
#> 
#> Tests (likelihood ratio, each row's model refitted):
#>   herd.(Intercept): drop (1 | herd); 0.5 chi2(0) + 0.5 chi2(1)
fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
varcomp(fm1)
#>                         vcov       SE    2.5 %    97.5 %
#> Block.(Intercept)   214.4745 168.8318  70.8610 2575.1692
#> Variety.(Intercept) 106.0618  67.8756  40.9772  657.7746
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
