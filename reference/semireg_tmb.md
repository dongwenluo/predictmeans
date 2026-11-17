# Fitting Semi Parametric Models Using glmmTMB

Fit a semi parametric model based on glmmTMB.

## Usage

``` r
semireg_tmb(formula, data, family = gaussian(), smoothZ = list(), ziformula = ~0, 
        dispformula = ~1, weights = NULL, offset = NULL, contrasts = NULL, na.action, 
    se = TRUE, verbose = FALSE, doFit = TRUE, control = glmmTMBControl(), 
    REML = TRUE, start = NULL, map = NULL, sparseX = NULL, prt=TRUE, 
    predict_info=TRUE)
```

## Arguments

- formula:

  A two-sided linear formula object describing both the fixed-effects
  and random-effects part of the model, with the response on the left of
  a ~ operator and the terms, separated by + operators, on the right.
  Random-effects terms are distinguished by vertical bars ("\|")
  separating expressions for design matrices from grouping factors.

- data:

  data frame (tibbles are OK) containing model variables. Not required,
  but strongly recommended; if `data` is not specified, downstream
  methods such as prediction with new data
  (`predict(fitted_model, newdata = ...)`) will fail. If it is necessary
  to call `glmmTMB` with model variables taken from the environment
  rather than from a data frame, specifying `data=NULL` will suppress
  the warning message.

- family:

  a family function, a character string naming a family function, or the
  result of a call to a family function (variance/link function)
  information. See [`family`](https://rdrr.io/r/stats/family.html) for a
  generic discussion of families or `family_glmmTMB` for details of
  `glmmTMB`-specific families.

- smoothZ:

  Superseded since 2.0 (kept so that 1.x code runs; no new features):
  write the smooths in the formula instead (see ‘Smooth terms in the
  formula’). The older form of the smooths: a list includes a set of
  smooth Z matrixs (called 'smooth term') used in the mixed effects
  model, the name of 'smooth term' should be different any variables in
  the model, each 'smooth term' is the result of function `smZ`. e.g.
  smoothZ=list(sm1=smZ(x1), sm2=smZ(x2, by=f1), sm3=smZ(x3, by=f2,
  group=TRUE), ...) where 'sm1' to 'sm3' should be new variable names in
  the `data`, and x1 to x3 are covariates, and f1, f2 are factors.

- ziformula:

  a *one-sided* (i.e., no response variable) formula for zero-inflation
  combining fixed and random effects: the default `~0` specifies no
  zero-inflation. Specifying `~.` sets the zero-inflation formula
  identical to the right-hand side of `formula` (i.e., the conditional
  effects formula); terms can also be added or subtracted. **When using
  `~.` as the zero-inflation formula in models where the conditional
  effects formula contains an offset term, the offset term will
  automatically be dropped**. The zero-inflation model uses a logit
  link.

- dispformula:

  a *one-sided* formula for dispersion containing only fixed effects:
  the default `~1` specifies the standard dispersion given any family.
  The argument is ignored for families that do not have a dispersion
  parameter. For an explanation of the dispersion parameter for each
  family, see [`sigma`](https://rdrr.io/pkg/lme4/man/sigma.html). The
  dispersion model uses a log link. In Gaussian mixed models,
  `dispformula=~0` fixes the residual variance to be 0 (actually a small
  non-zero value), forcing variance into the random effects. The precise
  value can be controlled via
  `control=glmmTMBControl(zero_dispval=...)`; the default value is
  `sqrt(.Machine$double.eps)`.

- weights:

  weights, as in `glm`. Not automatically scaled to have sum 1.

- offset:

  offset for conditional model (only).

- contrasts:

  an optional list, e.g., `list(fac1="contr.sum")`. See the
  `contrasts.arg` of
  [`model.matrix.default`](https://rdrr.io/r/stats/model.matrix.html).

- na.action:

  a function that specifies how to handle observations containing `NA`s.
  The default action (`na.omit`, inherited from the 'factory fresh'
  value of `getOption("na.action")`) strips any observations with any
  missing values in any variables. Using `na.action = na.exclude` will
  similarly drop observations with missing values while fitting the
  model, but will fill in `NA` values for the predicted and residual
  values for cases that were excluded during the fitting process because
  of missingness.

- se:

  whether to return standard errors.

- verbose:

  whether progress indication should be printed to the console.

- doFit:

  whether to fit the full model, or (if FALSE) return the preprocessed
  data and parameter objects, without fitting the model.

- control:

  control parameters, see `glmmTMBControl`.

- REML:

  whether to use REML estimation rather than maximum likelihood.
  Gaussian fits use REML by default; for other families the default is
  maximum likelihood (Laplace), as `glmer`, unless `REML = TRUE` is
  given. (Before 2.0 the default never reached glmmTMB, so gaussian fits
  were maximum likelihood.)

- start:

  starting values, expressed as a list with possible components `beta`,
  `betazi`, `betad` (fixed-effect parameters for conditional,
  zero-inflation, dispersion models); `b`, `bzi` (conditional modes for
  conditional and zero-inflation models); `theta`, `thetazi`
  (random-effect parameters, on the standard deviation/Cholesky scale,
  for conditional and z-i models); `psi` (extra family parameters, e.g.,
  shape for Tweedie models).

- map:

  a list specifying which parameter values should be fixed to a constant
  value rather than estimated. `map` should be a named list containing
  factors corresponding to a subset of the internal parameter names (see
  `start` parameter). Distinct factor values are fitted as separate
  parameter values, `NA` values are held fixed: e.g.,
  `map=list(beta=factor(c(1,2,3,NA)))` would fit the first three
  fixed-effect parameters of the conditional model and fix the fourth
  parameter to its starting value. In general, users will probably want
  to use `start` to specify non-default starting values for fixed
  parameters. See
  [`MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html) for more
  details.

- sparseX:

  a named logical vector containing (possibly) elements named "cond",
  "zi", "disp" to indicate whether fixed-effect model matrices for
  particular model components should be generated as sparse matrices,
  e.g. `c(cond=TRUE)`. Default is all `FALSE`

- prt:

  Logical scalar - Should the info to be print on screen in the middle
  of the process or not?

- predict_info:

  Logical scalar - Should provide the info for function semipred or not?
  In case of there is a correlation theta parameter appearing, you may
  set predict=FALSE.

## Value

- semer:

  A glmmTMB model used in the fitting.

- data:

  A data.frame with generated variables in the fitting.

- fomul_vars:

  Name of variables in the formula of semireg_tmb model.

- sm_vars:

  Name of variables in the smoothZ list.

- smoothZ_call:

  A call used to produce smooth terms in the fitting.

- knots_lst:

  Knots used in each smooth term in the fitting.

- range_lst:

  Range of covariate used in each smooth term in the fitting.

- cov_lst:

  Covariance matrix list for each smooth term.

- u_lst:

  Random effects list for each smooth term.

- type_lst:

  Smooth type list of smooth terms.

- CovMat:

  Covariance matrix for all smooth terms.

- Cov_ind:

  Covariance matrix index for each smooth term.

- Cov_indN:

  Covariance matrix index for each smooth term when `group=TRUE` in
  `smoothZ` argument.

- df:

  Degree of freedom of all random terms.

- tmbf:

  The glmmTMB model result using doFit=FALSE.

## Details

A semi parametric model can be parameterized as a linear (or generalized
linear) mixed model in which its random effects are smooth functions of
some covariates (named \<e2\>\<80\>\<98\>smooth term\<e2\>\<80\>\<99\>).
`semireg_tmb` follows the approach suggested by Wand and Ormerod (2008)
and represents the 'smooth term' using O'Sullivan-type of Z.

## Correlation structures

glmmTMB's covariance structures can be written next to the smooths and
are kept as written: `ar1()`, `ou()` (unequally spaced times, with
[`glmmTMB::numFactor()`](https://rdrr.io/pkg/glmmTMB/man/numFactor.html)),
`cs()`, `toep()`, [`diag()`](https://rdrr.io/r/base/diag.html), and the
spatial [`exp()`](https://rdrr.io/r/base/Log.html), `gau()`, `mat()`.
Residuals correlated in time within subjects are written as
`ar1(time_factor + 0 | subject)` with `dispformula = ~0` (the structure
is then the residual covariance; with a residual variance as well it is
AR1 plus a nugget, and `cs()` plus a residual variance is not
identified), e.g.
`semireg_tmb(y ~ s(time, by = trt) + ar1(tf + 0 | id), dispformula = ~0, data = d)`.
[`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
(lme4) has no such structures and refuses them.

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
| `s(x1, x2)` | isotropic surface (thin plate) | `x1 + x2` |
| `s(x, bs = "os", k = 10)` | O'Sullivan spline with `k` interior knots, as `smZ(x, k = 10)` | `x` |
| `s(x1, x2, bs = "ztps", k = 6)` | thin plate of `Ztps()`, as `smZ(cbind(x1, x2), type = "Ztps", k = 6)` with `k = 6` | `x1 + x2` |
| `s(x, bs = "user", xt = list(basis = f))` | a basis of your own (see [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)) | `x` |
| `t2(x1, x2)` | tensor-product surface, a variance per penalty | `x1 + x2 + x1:x2` |
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
deviations. The fits equal `gam(..., method = "REML")` for gaussian
models (with `s(g, bs = "re")` for a random intercept), and
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
and
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
give the curves, means and their standard errors of `predict.gam()`. The
variables of formula smooths are not centred or scaled (`ncenter` and
`nscale` apply to the other numeric covariates). Tensor products:
`t2(x1, x2)`, as in
[`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
(`te()` is fitted as `t2()`; `ti()` is not available). Smooths are given
either in the formula or in `smoothZ`, not both.

## Zero inflation (glmmTMB)

For a `glmmTMB` fit with zero inflation (`ziformula`) the means are
those of the conditional (count) part, \\\mu\\, as `emmeans` and
`predict(type = "conditional")` give, not the overall mean \\(1 -
\pi)\mu\\ (\\\pi\\ the probability of an extra zero); ratios are of
these means, and equal the ratios of overall means only when \\\pi\\ is
the same. A message says so and `settings$component` is `"conditional"`.
The two can order treatments differently when \\\pi\\ differs between
them. `mean_type = "overall"`
([`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md))
gives the overall means and their ratios.

## References

Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with
O'Sullivan penalized splines. *Australian and New Zealand Journal of
Statistics.* **50**, 179-198.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
## Counts with extra zeros: curves of two treatments, a random site effect, a negative
## binomial with zero inflation; the rate ratio of the curves along x and where it differs
## from 1
set.seed(5)
cn <- expand.grid(x = seq(0, 3, length = 30), trt = factor(c("ctl", "trt")), site = factor(1:8))
cn$y <- rnbinom(nrow(cn), size = 3,
                mu = exp(1 + sin(cn$x) + (cn$trt == "trt")*0.5*cn$x + rnorm(8, 0, 0.3)[cn$site]))
cn$y[runif(nrow(cn)) < 0.1] <- 0
fc <- semireg_tmb(y ~ s(x, by = trt, bs = "os", k = 6) + (1 | site), ziformula = ~ 1,
                  family = glmmTMB::nbinom2, data = cn, prt = FALSE)
rr <- semipred(fc, "trt", "x", contr = c("trt", "ctl"), band = "simultaneous", prt = FALSE)
#> The model has zero inflation (ziformula ~1): the means are those of the conditional (count) part, mu, as emmeans and predict(type = "conditional") give; the overall mean is (1 - pi) mu, with pi the probability of an extra zero. Ratios are of these means (the same for overall means only when pi is the same). mean_type = "overall" gives the overall means.
attr(rr$pred_df, "significant")
#>        from to direction
#> 1 0.7241379  3         >
predictmeans(fc, "trt", covariate = c(x = 1.5))$means
#> The model has zero inflation (ziformula ~1): the means are those of the conditional (count) part, mu, as emmeans and predict(type = "conditional") give; the overall mean is (1 - pi) mu, with pi the probability of an extra zero. Ratios are of these means (the same for overall means only when pi is the same).
#>   trt estimate std.error  df conf.low conf.high  response response.std.error
#> 1 ctl 1.842115 0.1565336 469 1.534521  2.149709  6.309871           0.987707
#> 2 trt 2.852923 0.1517429 469 2.554743  3.151103 17.338385           2.630977
#>   response.conf.low response.conf.high group
#> 1          4.639104           8.582363    A 
#> 2         12.867989          23.361818     B
# \donttest{
## Repeated measures: curves over time, AR(1) correlation within subjects (the residual
## covariance, dispformula = ~0); means at a time with Satterthwaite df
set.seed(9)
d <- expand.grid(t = 1:8, id = factor(1:12))
d$A <- factor(rep(c("a1", "a2", "a3"), each = 32))
d$y <- sin(d$t/3) + rnorm(12, 0, 0.4)[d$id] +
  unlist(lapply(1:12, function(i) as.numeric(stats::arima.sim(list(ar = 0.6), 8, sd = 0.4))))
d$tf <- factor(d$t)
fa <- semireg_tmb(y ~ s(t, by = A, bs = "os", k = 4) + (1 | id) + ar1(tf + 0 | id),
                  dispformula = ~0, data = d, prt = FALSE)
#> Warning: Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
predictmeans(fa, "A", covariate = c(t = 4), pairwise = TRUE)$pairs
#>   level1 level2   estimate std.error       df  statistic   p.value     p.adj
#> 1     a1     a2 -0.2146406 0.3486284 20.51249 -0.6156717 0.5448827 0.5448827
#> 2     a1     a3 -0.3028874 0.3191386 17.43094 -0.9490782 0.3555549 0.3555549
#> 3     a2     a3 -0.0882468 0.3358240 19.17360 -0.2627769 0.7955254 0.7955254
#>   p.method
#> 1        t
#> 2        t
#> 3        t
varcomp(fa)
#>                       vcov     SE  2.5 %       97.5 %
#> id.(Intercept)      0.0000     NA 0.0000           NA
#> id.1.var            0.3314 0.0850 0.2122 5.893000e-01
#> cor.id.1            0.7043 0.0773 0.5191 8.264000e-01
#> s_t_Aa1.(Intercept) 0.0056 0.0097 0.0009 1.807443e+02
#> s_t_Aa2.(Intercept) 0.0274 0.0307 0.0067 2.432300e+00
#> s_t_Aa3.(Intercept) 0.0016 0.0057 0.0002 2.658432e+16
#> residual            0.0000     NA 0.0000           NA

## Binary responses with a random intercept per child (HRW::indonRespir)
if (requireNamespace("HRW", quietly = TRUE)) {
  data(indonRespir, package = "HRW")
  fb <- semireg_tmb(respirInfec ~ s(age, bs = "os", k = 10) + vitAdefic + female + height +
                      stunted + visit2 + visit3 + visit4 + visit5 + visit6 + (1 | idnum),
                    family = binomial, data = indonRespir, prt = FALSE)
  semipred(fb, "age", "age", prt = FALSE)$plt
}

# }
```
