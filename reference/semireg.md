# Semiparametric Mixed Models Fitted with lme4

Fit a semiparametric (additive mixed) model with lme4 (`lmer` for
gaussian, `glmer` for other families): smooth terms written in the
formula, as in [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), are
fitted as random effects. For negative binomial, zero-inflated or
correlated data use
[`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md).

## Usage

``` r
semireg(formula, data, family = NULL, REML = TRUE,
        control = lmerControl(optimizer = "nloptwrap"),
        start = NULL, verbose = FALSE, drop_unused_levels = TRUE, subset, weights,
        offset, contrasts = NULL, ...)
```

## Arguments

- formula:

  A two-sided linear formula object describing both the fixed-effects
  and random-effects part of the model, with the response on the left of
  a ~ operator and the terms, separated by + operators, on the right.
  Random-effects terms are distinguished by vertical bars ("\|")
  separating expressions for design matrices from grouping factors.

- data:

  A data frame or list containing the model response variable and
  covariates required by the formula. By default the variables are taken
  from environment(formula), typically the environment from which
  semireg is called.

- family:

  A GLM family (see [`family`](https://rdrr.io/r/stats/family.html));
  `NULL` (the default) is gaussian.

- REML:

  Logical scalar - Should the estimates be chosen to optimize the REML
  criterion (as opposed to the log-likelihood)?

- control:

  A list (of correct class, resulting from lmerControl() or
  glmerControl() respectively) containing control parameters, including
  the nonlinear optimizer to be used and parameters to be passed through
  to the nonlinear optimizer, see the \*lmerControl documentation for
  details. The default is `lmerControl(optimizer = "nloptwrap")`
  (`"bobyqa"` in 1.x), also used for generalized fits.

- start:

  Starting values, as for `lmer` or `glmer`.

- verbose:

  Passed on to the lme4 fitting functions.

- drop_unused_levels:

  By default unused levels are dropped from factors before fitting. For
  some smooths involving factor variables you might want to turn this
  off. Only do so if you know what you are doing.

- subset:

  An optional expression indicating the subset of the rows of data that
  should be used in the fit. This can be a logical vector, or a numeric
  vector indicating which observation numbers are to be included, or a
  character vector of the row names to be included. All observations are
  included by default.

- weights:

  An optional vector of \<e2\>\<80\>\<98\>prior
  weights\<e2\>\<80\>\<99\> to be used in the fitting process. Should be
  NULL or a numeric vector.

- offset:

  This can be used to specify an a priori known component to be included
  in the linear predictor during fitting. This should be NULL or a
  numeric vector of length equal to the number of cases. One or more
  offset terms can be included in the formula instead or as well, and if
  more than one is specified their sum is used. See model.offset.

- contrasts:

  An optional list, as `contrasts.arg` of
  [`model.matrix.default`](https://rdrr.io/r/stats/model.matrix.html).

- ...:

  Further arguments for passing on to model setup routines. Arguments
  removed in 2.0 (`smoothZ`, `ncenter`, `nscale`, `resp_scale`, `prt`,
  `predict_info`, `ngbinomial`) stop with their 2.0 spelling; nothing is
  printed during the fit (`print(fit$semer)`). For a negative binomial
  model use
  [`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)`(family = glmmTMB::nbinom2)`:
  the 1.x `ngbinomial = TRUE` estimated theta from a Poisson fit, which
  is not the maximum-likelihood fit.

## Value

A list of class `semireg`.
[`print()`](https://rdrr.io/r/base/print.html) shows a short account of
the fit: the formula as written, the smooth terms (one row per
random-effect block: the smooth, its random-effect term, its covariate,
basis, basis dimension, smoothing variance and whether it is estimated
as a straight line), the fixed effects, the other variance components
and any advice. The `smooth`, `term` and `covariate` columns are the
values
[`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
takes for `sm_term` and `covariate`.
[`summary()`](https://rdrr.io/r/base/summary.html) adds the summary of
the fitted mixed model (`summary(fit)$model`); `tidy(fit)` returns the
smooth table and `tidy(fit, "fixed")` the fixed effects with standard
errors. [`coef()`](https://rdrr.io/r/stats/coef.html),
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`family()`](https://rdrr.io/r/stats/family.html),
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) and
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) work as
for the fitted mixed model, and
[`formula()`](https://rdrr.io/r/stats/formula.html) gives the formula as
written. [`anova()`](https://rdrr.io/r/stats/anova.html) tests the terms
and smooths; notes on the fit are in `attr(, "advice")`
([`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).
Its elements:

- semer:

  The fitted `lmer` or `glmer` model (class `semerMod` for gaussian
  fits).

- call:

  The call of `semireg()`.

- data:

  A data.frame with generated variables in the fitting.

- fomul_vars:

  Name of variables in the formula of semireg model.

- sm_vars:

  Names of the variables of the smooths.

- cov_lst:

  Covariance matrix list for each smooth term.

- u_lst:

  Random effects list for each smooth term.

- CovMat:

  Covariance matrix for all smooth terms.

- Cov_ind:

  Covariance matrix index for each smooth term.

- Cov_indN:

  Covariance matrix index for each smooth term of several blocks (e.g.
  `s(x, by = f)`).

- df:

  Degree of freedom of all random terms.

- lmerc:

  Call used in the mer model in the fitting.

- smooths:

  For smooths written in the formula: the mgcv smooth objects and the
  maps that give the random-effect design of each block for new data
  (used by
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  and
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).

## Details

A semiparametric model is fitted as a linear (or generalized linear)
mixed model whose random effects are the penalised parts of smooth
functions of covariates (Wand and Ormerod, 2008); each smooth's variance
is its smoothing parameter, estimated by (RE)ML.

## Smooth terms in the formula

Smooths can be written in the model formula, as for
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), next to fixed and
lme4 random terms: `y ~ A + s(x, by = A) + (1 | Block)`. Each smooth is
set up by mgcv (bases, knots, identifiability constraints) and split, as
in gamm4, into an unpenalised part, which is added to the fixed effects
as ordinary terms, and penalised parts, each a block of random effects
with one variance. The fixed terms are added automatically (terms
already written are not doubled), so a smooth by a factor cannot lose
its group slopes:

|  |  |  |
|----|----|----|
| **Written** | **Meaning** | **Fixed terms added** |
| `s(x)`, `s(x, bs = "cr", k = 12)` | one curve | `x` |
| `s(x, by = f)` | a curve per level of `f`, a variance each (`id = 1`: one shared) | `f + f:x` |
| `s(x, by = z)`, `z` numeric | varying coefficient \\f(x) z\\ | `z + z:x` |
| `s(x, f, bs = "fs")` | random curves per level, one variance | none |
| `(s(x, k = 10) | g)` | random curve per level of `g`: `(x | g)` (correlated intercepts and slopes) plus a spline deviation per level, one variance | none |
| `s(x1, x2)` | isotropic surface (thin plate) | `x1 + x2` |
| `s(x, bs = "os", k = 10)` | O'Sullivan spline with `k` interior knots (1.x `smZ(x, k = 10)`); `xt = list(knots = , range = )` sets the interior knots and the range | `x` |
| `s(x1, x2, bs = "ztps", k = 6)` | low-rank thin plate with `k` knots (1.x `smZ(cbind(x1, x2), type = "Ztps", k = 6)`); `xt = list(knots = )` sets the knots (a two-column matrix) | `x1 + x2` |
| `s(x, bs = "user", xt = list(basis = f))` | a basis of your own: `f(x)` (`f(x1, x2)` for two covariates) returns a matrix with a row per value, whose columns get iid random coefficients with one variance (as a matrix in `smoothZ`); checked before the fit (see below); `xt = list(basis = f, null = "constant")` or `"none"` for a smaller fixed part | `x` (`x1 + x2`) |
| `t2(x1, x2)`, `t2(x1, x2, bs = c("cr", "ps"), k = c(6, 8))` | tensor-product surface (covariates on different scales), a variance per penalty; `te()` is fitted as `t2()` | `x1 + x2 + x1:x2` |
| `s(x, bs = "cs")`, `"ts"` | shrinkage smooths (no unpenalised part) | none |
| `s(x, f)` | common curve + a deviation curve per level of `f`, one variance for the deviations | `x + f + f:x` |
| `s(x, f1:f2)` | the same with the cells of `f1:f2` | `x + f1:f2 + f1:f2:x` |
| `s(x, f1*f2)` | common curve + deviation curves for `f1`, `f2` and `f1:f2`, a variance each | `x + f1*f2 + f1:x + f2:x + f1:f2:x` |
| `s(x1, x2, f)` | common surface + a deviation surface per level | `x1 + x2 + f + f:x1 + f:x2` |

A factor term among the arguments of `s()` gives this ANOVA-like split
(not an mgcv form): the deviation curves of a factor term are random
effects with one variance, so
[`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
shows how much the curves differ by `f1`, by `f2` and by their
interaction, and a variance of zero means the curves of that term have
the same shape (their levels still differ by the fixed terms). The
pieces are set up without constraints between them: the mixed model
separates them by their variances, as nested random effects. Each
level's curve is the fixed part plus the common curve plus its
deviations.

A random curve `(s(x) | g)` (O'Sullivan basis unless `bs` is given) is
the group-specific curve model of Durban et al. (2005) and Harezlak,
Ruppert and Wand (2018, section 9.2): each level of `g` (a subject) has
its own straight line, `(x | g)`, and its own spline deviation, all
levels sharing one variance; e.g. heights of boys over age,
`height ~ s(age, bs = "os", k = 20) + (s(age, k = 10) | id)`. Unlike
`s(x, g, bs = "fs")`, the intercepts and slopes are correlated, as
`pdSymm(~ x)` in [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html).
The random curves are random effects:
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
and
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
give the population curves without them, and
[`anova()`](https://rdrr.io/r/stats/anova.html) tests whether they are
needed (their variance is zero: each subject's curve is a straight
line). The fits equal `gam(..., method = "REML")` for gaussian models
(with `s(g, bs = "re")` for a random intercept), and
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
and
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
give the curves, means and their standard errors of `predict.gam()`. The
variables are used as they are (not centred or scaled). Tensor products:
`t2(x1, x2)` (Wood, Scheipl and Faraway 2013), a surface for covariates
on different scales, whose penalties (one per product of the marginal
penalised and unpenalised parts; `full = TRUE` for more) are independent
random effects, each with its own variance; the fixed terms
`x1 + x2 + x1:x2` are added, and `by = f` gives a surface per level. It
is set up without mgcv's sum-to-zero constraint on the whole product
(which leaves a little of the penalised space in the unpenalised one),
so its fixed part is exactly those terms; the fit equals `gam()` with
that basis and its penalties, and `gam()`'s own `t2()` when there are no
other random effects. `te()` is fitted as `t2()` (with a message): its
overlapping penalties cannot be written as random effects; `ti()` is not
available. Objects named in `s()` (e.g. `xt = list(knots = kn)`) are
taken from the environment of the formula.

A basis of your own (`bs = "user"`) must give the same row for the same
covariate value whatever data it is given, as prediction evaluates it on
new values: knots and range belong in the function, e.g.
`f <- function(x) ZOSull(x, intKnots = c(2, 4, 6, 8), range.x = c(0, 10))`.
It is checked before the fit, and the fit stops with what is wrong and
how to put it right when it fails on the data, does not return a numeric
matrix with a row per value (also for a single value), gives missing or
infinite values, gives different values on two calls, gives other rows
for the same values on part of the data (knots or range taken from the
data), or has linearly dependent columns or columns in the span of its
fixed part; a warning says when it fails just beyond the ends of the
data (the slopes of `semipred(deriv = )` there need it). Its values are
checked again at prediction.

The penalised-spline examples fitted by
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) in Wand's course
scripts and in HRW (additive models, additive mixed models, thin-plate
surfaces with given knots, group-specific curves with a contrast curve)
give the same REML log-likelihood and fitted values with `semireg()`
(`tests/testthat/test-semireg-golden.R`); their data are in the package:
[`growthIndiana`](https://dongwenluo.github.io/predictmeans/reference/growthIndiana.md),
[`SBMDfem`](https://dongwenluo.github.io/predictmeans/reference/SBMDfem.md),
[`ozoneCalif`](https://dongwenluo.github.io/predictmeans/reference/ozoneCalif.md),
[`scallop`](https://dongwenluo.github.io/predictmeans/reference/scallop.md),
[`wageTUnion`](https://dongwenluo.github.io/predictmeans/reference/wageTUnion.md),
[`copper`](https://dongwenluo.github.io/predictmeans/reference/copper.md),
whose help pages show the analyses. The optimizer is `nloptwrap` (lme4's
default for `lmer`) since 2.0: on these examples it reached the same
optimum as `bobyqa` (1.x's default;
`control = lmerControl(optimizer = "bobyqa")`) faster, also for binomial
and Poisson fits. Each smooth block of a formula smooth is scaled as a
random intercept (sum of squares equal to the number of observations),
so that the optimizer starts near the optimum; the variance of a smooth
block reported by
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) or
[`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
is on that scale.

Moving from 1.x:
[`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
and `smoothZ` are superseded and were removed in 2.0 (they stop with a
pointer to the formula syntax);
[`vignette("migrating")`](https://dongwenluo.github.io/predictmeans/articles/migrating.md)
shows each `smoothZ = list(sm = smZ(...))` form beside its formula form,
which gives the same fit (checked in `dev/check-smz-migration.R`). The
variables are no longer centred (1.x's `ncenter = TRUE`), which put
knots given in the data's units (`smZ(x, intKnots = )`) in the wrong
place.

## References

Durban, M., Harezlak, J., Wand, M.P. and Carroll, R.J. (2005). Simple
fitting of subject-specific curves for longitudinal data. *Statistics in
Medicine*, **24**, 1153-1167.

Harezlak, J., Ruppert, D. and Wand, M.P. (2018). *Semiparametric
Regression with R*. Springer, New York.

Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with
O'Sullivan penalized splines. *Australian and New Zealand Journal of
Statistics.* **50**, 179-198.

Wood, S.N., Scheipl, F. and Faraway, J.J. (2013). Straightforward
intermediate rank tensor product smoothing in mixed models. *Statistics
and Computing*, **23**, 341-360.

## See also

[`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md),
[`anova.semireg`](https://dongwenluo.github.io/predictmeans/reference/anova.semireg.md),
[`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md),
[`varcomp`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
## 1. Curves in a mixed model (nlme::Dialyzer): a curve of the ultrafiltration rate over
##    pressure for each blood flow rate, a random intercept and slope per dialyzer; the
##    fixed terms QB + QB:pressure are added automatically
data(Dialyzer, package = "nlme")
fm <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
              data = Dialyzer)
anova(fm)       # fixed effects (Kenward-Roger df) and a likelihood ratio test of each smooth
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: rate
#> 
#> Advice:
#> ! A numeric covariate is in an interaction: Type III tests of the other
#>   terms are at covariate = 0; centre the covariate if 0 is not a meaningful
#>   value.
#> 
#>                    F Df Df.res    Pr(>F)    
#> (Intercept) 212.0332  1 19.128 8.354e-12 ***
#> QB            7.4995  1 19.128     0.013 *  
#> QB:pressure 850.3842  2 19.129 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);
#> reference: exact restricted LRT for the only random effect of a gaussian fit, otherwise 0.5 chi2(0) + 0.5 chi2(1)
#>              term variance logLik logLik.without   LRT p.value
#>  s_pressure_QB200    28.92 -355.6         -476.4 241.5  <2e-16
#>  s_pressure_QB300    24.53 -355.6         -468.3 225.3  <2e-16
#>                       test
#>  0.5 chi2(0) + 0.5 chi2(1)
#>  0.5 chi2(0) + 0.5 chi2(1)
varcomp(fm)     # variance components, the smooths' among them
#>                                 vcov      SE   2.5 %   97.5 %
#> Subject.(Intercept)           5.5166  2.7181  2.5401  19.7492
#> Subject.pressure.(Intercept) -1.4806  1.3641 -2.8448   0.6615
#> Subject.pressure              2.4138  1.0703  1.1844   7.3436
#> s_pressure_QB200.(Intercept) 28.9181 21.3160 10.0569 273.5581
#> s_pressure_QB300.(Intercept) 24.5258 18.6606  8.3226 259.2940
#> residual                      4.8277  0.7119  3.6879   6.5945
predictmeans(fm, "QB", covariate = c(pressure = 1))$means
#>    QB estimate std.error       df conf.low conf.high group
#> 1 200 36.30263 0.9244804 36.54632 34.42867  38.17659     A
#> 2 300 38.81293 0.9262667 36.79558 36.93578  40.69007     A
semipred(fm, "QB", "pressure", contr = c("300", "200"), band = "simultaneous")$plt


## 2. A designed experiment: treatments on the whole plots of six blocks, curves over time;
##    treatment means at a time get Satterthwaite df from the whole-plot stratum
set.seed(11)
sp <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
sp$y <- 2 + sin(sp$t/2) + 0.08*(sp$A == "a2")*sp$t + rnorm(6, 0, 0.5)[sp$Block] +
  rnorm(18, 0, 0.3)[interaction(sp$Block, sp$A)] + rnorm(nrow(sp), 0, 0.3)
fsp <- semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A),
               data = sp)
predictmeans(fsp, "A", covariate = c(t = 5), pairwise = TRUE)$pairs
#>   level1 level2   estimate std.error       df  statistic    p.value      p.adj
#> 1     a1     a2 -0.2449559 0.1611056 25.50375 -1.5204680 0.14069552 0.14069552
#> 2     a1     a3  0.1291968 0.1602868 25.01663  0.8060353 0.42781848 0.42781848
#> 3     a2     a3  0.3741527 0.1615742 25.78521  2.3156705 0.02879782 0.02879782
#>      conf.low  conf.high       lsd p.method
#> 1 -0.57642702 0.08651528 0.3314711        t
#> 2 -0.20090884 0.45930242 0.3301056        t
#> 3  0.04189739 0.70640793 0.3322553        t

## 3. Two crossed factors: a common curve and deviation curves for f1, f2 and f1:f2, each
##    with its own variance; which deviations are needed?
set.seed(1)
dx <- expand.grid(x = seq(0, 1, length = 15), f1 = factor(c("p", "q")), f2 = factor(c("u", "v")),
                  rep = 1:3)
dx$y <- sin(2*pi*dx$x) + (dx$f1 == "q")*sin(pi*dx$x)^2 + rnorm(nrow(dx), 0, 0.2)
fx <- semireg(y ~ s(x, f1*f2, bs = "os", k = 6), data = dx)
varcomp(fx, test = TRUE)
#>                         vcov     SE  2.5 %      97.5 %     LRT   p.value
#> s_x_f1_f2.(Intercept) 0.0000     NA 0.0000          NA  0.0000 1.0000000
#> s_x_f2.(Intercept)    0.0005 0.0011 0.0001 191055.5255  0.3757 0.2700000
#> s_x_f1.(Intercept)    0.0532 0.0426 0.0173      0.6851 11.6162 0.0003269
#> s_x.(Intercept)       0.3776 0.2198 0.1559      1.8653  7.1297 0.0037910
#> residual              0.0345 0.0039 0.0280      0.0435      NA        NA
#> 
#> Tests (likelihood ratio, each row's model refitted):
#>   s_x_f1_f2.(Intercept): smooth s_x_f1_f2 dropped; 0.5 chi2(0) + 0.5 chi2(1)
#>   s_x_f2.(Intercept): smooth s_x_f2 dropped; 0.5 chi2(0) + 0.5 chi2(1)
#>   s_x_f1.(Intercept): smooth s_x_f1 dropped; 0.5 chi2(0) + 0.5 chi2(1)
#>   s_x.(Intercept): smooth s_x dropped; 0.5 chi2(0) + 0.5 chi2(1)
#> 
#> Advice:
#> ! A variance is estimated as zero (or a correlation as +-1): its interval
#>   and SE are not reliable, and 'no variation' is not shown. Test it (test =
#>   TRUE) or simplify the random effects.
#> i To test a variance against zero use varcomp(test = TRUE): it uses the
#>   boundary (mixture or bootstrap) reference, not the naive chi-square,
#>   which is conservative.

## 4. An additive model with a varying coefficient: a smooth in x1, and the effect of z
##    changing smoothly with x2 (numeric 'by')
set.seed(2)
da <- data.frame(x1 = runif(200), x2 = runif(200), z = rnorm(200))
da$y <- sin(3*da$x1) + da$z*sin(2*pi*da$x2) + rnorm(200, 0, 0.3)
fa <- semireg(y ~ s(x1, bs = "os") + s(x2, by = z, bs = "os"), data = da)
anova(fa)$smooth
#>     term   variance   logLik logLik.without      LRT      p.value
#> 1   s_x1 0.03875716 -65.4386      -126.8535 122.8299 7.596575e-29
#> 2 s_x2_z 0.33715120 -65.4386      -159.4806 188.0839 4.164359e-43
#>                        test
#> 1 0.5 chi2(0) + 0.5 chi2(1)
#> 2 0.5 chi2(0) + 0.5 chi2(1)

## 5. A basis of your own (checked before the fit): a seasonal cycle on a straight trend
set.seed(3)
ds <- data.frame(day = 1:730)
ds$y <- 10 + 0.002*ds$day + 2*sin(2*pi*ds$day/365) + rnorm(730)
season <- function(d) cbind(sin(2*pi*d/365), cos(2*pi*d/365), sin(4*pi*d/365), cos(4*pi*d/365))
fs <- semireg(y ~ day + s(day, bs = "user", xt = list(basis = season, null = "none")),
              data = ds)
semipred(fs, covariate = "day")$plt

## 6. A tensor-product surface for two covariates on different scales (data(ozoneCalif)),
##    and a curve in a third
data(ozoneCalif)
fo <- semireg(ozone ~ t2(inversion.base.height, inversion.base.temp) +
                s(daggett.pressure.gradient, bs = "os", k = 15), data = ozoneCalif)
anova(fo)$smooth          # a variance per penalty of the surface
#>                                             term  variance    logLik
#> 1                    s_daggett_pressure_gradient  4.591269 -1022.055
#> 2 t2_inversion_base_height_inversion_base_temp_1  0.000000 -1022.055
#> 3 t2_inversion_base_height_inversion_base_temp_2  1.124284 -1022.055
#> 4 t2_inversion_base_height_inversion_base_temp_3 12.919219 -1022.055
#>   logLik.without          LRT      p.value                      test
#> 1      -1040.213 3.631529e+01 8.392096e-10 0.5 chi2(0) + 0.5 chi2(1)
#> 2      -1022.055 2.169145e-10 4.999941e-01 0.5 chi2(0) + 0.5 chi2(1)
#> 3      -1022.438 7.662851e-01 1.906844e-01 0.5 chi2(0) + 0.5 chi2(1)
#> 4      -1026.259 8.407200e+00 1.868689e-03 0.5 chi2(0) + 0.5 chi2(1)
semipred(fo, covariate = c("inversion.base.height", "inversion.base.temp"))$plt

# \donttest{
## 7. Binary responses with a random intercept per child (HRW::indonRespir): the probability
##    of a respiratory infection over age
if (requireNamespace("HRW", quietly = TRUE)) {
  data(indonRespir, package = "HRW")
  fb <- semireg(respirInfec ~ s(age, bs = "os", k = 10) + vitAdefic + female + height + stunted +
                  visit2 + visit3 + visit4 + visit5 + visit6 + (1 | idnum),
                family = binomial, data = indonRespir)
  semipred(fb, "age", "age")$plt
}


## 8. Subject-specific curves: each subject's own straight line, (Days | Subject), and its
##    own spline deviation, one variance for all; are the subjects' curves straight lines?
fr <- semireg(Reaction ~ s(Days, bs = "os", k = 4) + (s(Days, k = 4) | Subject),
              data = lme4::sleepstudy)
anova(fr)$smooth
#>             term     variance    logLik logLik.without      LRT      p.value
#> 1         s_Days 3.055837e-07 -862.9684      -862.9684  0.00000 1.000000e+00
#> 2 s_Days_Subject 2.251268e+02 -862.9684      -871.7951 17.65343 1.325185e-05
#>                        test
#> 1 0.5 chi2(0) + 0.5 chi2(1)
#> 2 0.5 chi2(0) + 0.5 chi2(1)

## 9. Curves by group, the price of apartments in Warsaw by district (HRW::WarsawApts)
if (requireNamespace("HRW", quietly = TRUE)) {
  data(WarsawApts, package = "HRW")
  WarsawApts$district <- factor(WarsawApts$district)
  fw <- semireg(areaPerMzloty ~ s(construction.date, by = district, bs = "os", k = 15),
                data = WarsawApts)
  semipred(fw, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"))$plt
  predictmeans(fw, "district", covariate = c(construction.date = 1980))$means
}
#>      district estimate std.error  df conf.low conf.high group
#> 1     Mokotow 130.9225  3.510809 401 124.0207  137.8244    A 
#> 2 Srodmiescie 138.2202  8.320780 401 121.8624  154.5780    A 
#> 3        Wola 118.9150  3.389961 401 112.2507  125.5794     B
#> 4    Zoliborz 131.6012  5.959483 401 119.8854  143.3169    AB
# }
```
