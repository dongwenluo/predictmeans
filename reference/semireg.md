# Fitting Semi Parametric Models Using lme4 Ecosystem

Fit a semi parametric model based on lme4 ecosystem including `lmer`,
`glmer` and `glmer.nb`.

## Usage

``` r
semireg(formula, data, family = NULL, ngbinomial=FALSE, REML = TRUE, 
        smoothZ = list(), ncenter=TRUE, nscale=FALSE, resp_scale=FALSE, 
        control = lmerControl(optimizer="nloptwrap"), start = NULL, 
        verbose = FALSE, drop.unused.levels=TRUE, subset, weights, 
        offset, contrasts = NULL,  prt=TRUE, predict_info=TRUE, ...)
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
  from environment(formula and smoothZ), typically the environment from
  which semireg is called.

- family:

  A GLM family, see glm and family.

- ngbinomial:

  Logical scalar - Should a negative binomial GLMMs be used? .

- REML:

  Logical scalar - Should the estimates be chosen to optimize the REML
  criterion (as opposed to the log-likelihood)?

- smoothZ:

  The older form of the smooths (see ‘Smooth terms in the formula’ for
  the newer one): a list includes a set of smooth Z matrixs (called
  'smooth term') used in the mixed effects model, the name of 'smooth
  term' should be different any variables in the model, each 'smooth
  term' is the result of function `smZ`. e.g. smoothZ=list(sm1=smZ(x1),
  sm2=smZ(x2, by=f1), sm3=smZ(x3, by=f2, group=TRUE), ...) where 'sm1'
  to 'sm3' should be new variable names in the `data`, and x1 to x3 are
  covariates, and f1, f2 are factors.

- ncenter:

  Logical scalar - Should the numeric predictors to be centered or not?

- nscale:

  Logical scalar - Should the numeric predictors to be scaled or not?

- resp_scale:

  Logical scalar - Should the response be involved in the scaling action
  or not?

- control:

  A list (of correct class, resulting from lmerControl() or
  glmerControl() respectively) containing control parameters, including
  the nonlinear optimizer to be used and parameters to be passed through
  to the nonlinear optimizer, see the \*lmerControl documentation for
  details. The default is `lmerControl(optimizer = "nloptwrap")`
  (`"bobyqa"` in 1.x), also used for generalized fits.

- start:

  Starting value list as used by lmer or glmer.

- verbose:

  Passed on to fitting lme4 fitting routines.

- drop.unused.levels:

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

  An optional list. See the contrasts.arg of model.matrix.default.

- prt:

  Logical scalar - Should the info to be print on screen in the middle
  of the process or not?

- predict_info:

  Logical scalar - Should provide the info for function semipred or not?

- ...:

  Further arguments for passing on to model setup routines.

## Value

- semer:

  A mer model used in the fitting.

- data:

  A data.frame with generated variables in the fitting.

- fomul_vars:

  Name of variables in the formula of semireg model.

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

A semi parametric model can be parameterized as a linear (or generalized
linear) mixed model in which its random effects are smooth functions of
some covariates (named \<e2\>\<80\>\<98\>smooth term\<e2\>\<80\>\<99\>).
`semireg` follows the approach suggested by Wand and Ormerod (2008) and
represents the 'smooth term' using O'Sullivan-type of Z.

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
| `s(x, bs = "os", k = 10)` | O'Sullivan spline with `k` interior knots, as `smZ(x, k = 10)`; `xt = list(knots = , range = )` sets the interior knots and the range | `x` |
| `s(x1, x2, bs = "ztps", k = 6)` | thin plate of `Ztps()`, as `smZ(cbind(x1, x2), type = "Ztps")`; `xt = list(knots = )` sets the knots (a two-column matrix) | `x1 + x2` |
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
variables of formula smooths are not centred or scaled (`ncenter` and
`nscale` apply to the other numeric covariates). Tensor products
(`te()`, `t2()`, `ti()`) are not available yet. Smooths are given either
in the formula or in `smoothZ`, not both. Objects named in `s()` (e.g.
`xt = list(knots = kn)`) are taken from the environment of the formula.

The penalised-spline examples fitted by
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) in Wand's course
scripts and in HRW (additive models, additive mixed models, thin-plate
surfaces with given knots, group-specific curves with a contrast curve)
give the same REML log-likelihood and fitted values with `semireg()`
(data in `system.file("extdata", "semireg", package = "predictmeans")`;
`tests/testthat/test-semireg-golden.R`). The optimizer is `nloptwrap`
(lme4's default for `lmer`) since 2.0: on these examples it reached the
same optimum as `bobyqa` (1.x's default;
`control = lmerControl(optimizer = "bobyqa")`) faster, also for binomial
and Poisson fits. Each smooth block of a formula smooth is scaled as a
random intercept (sum of squares equal to the number of observations),
so that the optimizer starts near the optimum; the variance of a smooth
block reported by
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) or
[`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
is on that scale.

Moving from the `smoothZ` form:
`y ~ x, smoothZ = list(sm = smZ(x, k = 10))` is
`y ~ s(x, bs = "os", k = 10)`; `smZ(x, by = f, group = TRUE)` is
`s(x, by = f, bs = "os")`; `smZ(x, by = f)` (one variance) is
`s(x, by = f, bs = "os", id = 1)`; `smZ(cbind(x1, x2), type = "Ztps")`
is `s(x1, x2, bs = "ztps")`; `smZ(x, intKnots = kn, range.x = r)` is
`s(x, bs = "os", xt = list(knots = kn, range = r))`. The fits are the
same
([`vignette("smooth")`](https://dongwenluo.github.io/predictmeans/articles/smooth.md)).

## References

Durban, M., Harezlak, J., Wand, M.P. and Carroll, R.J. (2005). Simple
fitting of subject-specific curves for longitudinal data. *Statistics in
Medicine*, **24**, 1153-1167.

Harezlak, J., Ruppert, D. and Wand, M.P. (2018). *Semiparametric
Regression with R*. Springer, New York.

Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with
O'Sullivan penalized splines. *Australian and New Zealand Journal of
Statistics.* **50**, 179-198.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
if (requireNamespace("HRW", quietly = TRUE)) {
  data(WarsawApts, package = "HRW")
  WarsawApts$district <- factor(WarsawApts$district)
  # smooths in the formula: a curve per district; district and its slopes are added as fixed terms
  fit0 <- semireg(areaPerMzloty ~ s(construction.date, by = district), data = WarsawApts, prt = FALSE)
  semipred(fit0, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"), prt = FALSE)$plt
  predictmeans(fit0, "district", covariate = c(construction.date = 1980))$means
}
#>      district estimate std.error        df conf.low conf.high group
#> 1     Mokotow 131.5985  3.327568 15.395491 124.5218  138.6752    A 
#> 2 Srodmiescie 120.0818  6.036074  7.781159 106.0941  134.0694    AB
#> 3        Wola 118.8170  3.473279 10.389220 111.1172  126.5169     B
#> 4    Zoliborz 131.4446  6.091278 34.210731 119.0684  143.8207    AB
if (requireNamespace("HRW", quietly = TRUE)) {
  data(WarsawApts, package = "HRW")
  fit1 <- semireg(areaPerMzloty ~ construction.date,
                  smoothZ=list(
                    grp=smZ(construction.date, k=25)
                  ),
                  data = WarsawApts)
  sp_out1 <- semipred(fit1, "construction.date", "construction.date")

  WarsawApts$district <- factor(WarsawApts$district)
  fit2 <- semireg(areaPerMzloty ~ construction.date*district, resp_scale = TRUE,
                  smoothZ=list(group=smZ(construction.date, k=15,
                                         by = district, group=TRUE)),
                  data=WarsawApts)
  sp_out2_1 <- semipred(fit2, "district", "construction.date")
  sp_out2_2 <- semipred(fit2, "district", "construction.date", contr=c(2,1))
}
#> Linear mixed model fit by REML ['semerMod']
#> Formula: areaPerMzloty ~ construction.date + (1 | grp)
#>    Data: WarsawApts
#> REML criterion at convergence: 3544.905
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  grp      (Intercept)  1.191  
#>  Residual             17.812  
#> Number of obs: 409, groups:  grp, 27
#> Fixed Effects:
#>       (Intercept)  construction.date  
#>         109.38694            0.06692  

#> Linear mixed model fit by REML ['semerMod']
#> Formula: areaPerMzloty ~ construction.date * district + (1 | group.Mokotow) +  
#>     (1 | group.Srodmiescie) + (1 | group.Wola) + (1 | group.Zoliborz)
#>    Data: WarsawApts
#> REML criterion at convergence: 3501.211
#> Random effects:
#>  Groups            Name        Std.Dev.
#>  group.Mokotow     (Intercept)  1.0336 
#>  group.Srodmiescie (Intercept)  2.9070 
#>  group.Wola        (Intercept)  0.1686 
#>  group.Zoliborz    (Intercept)  0.5914 
#>  Residual                      16.4102 
#> Number of obs: 409, groups:  
#> group.Mokotow, 17; group.Srodmiescie, 17; group.Wola, 17; group.Zoliborz, 17
#> Fixed Effects:
#>                           (Intercept)                      construction.date  
#>                               1.02674                                0.09594  
#>                   districtSrodmiescie                           districtWola  
#>                              -9.67368                                2.89944  
#>                      districtZoliborz  construction.date:districtSrodmiescie  
#>                              -5.35968                               -0.21367  
#>        construction.date:districtWola     construction.date:districtZoliborz  
#>                              -0.35512                               -0.17735  

#> Srodmiescie differs from Mokotow (pointwise band, level 0.05) for construction.date: 1944 to 1951 (<); 1960 to 1975 (<); 1990 to 1998 (<); 2003 to 2007 (<).

# \donttest{
if (requireNamespace("HRW", quietly = TRUE)) {
  data(indonRespir, package = "HRW")
  fit3 <- semireg(respirInfec ~ age+vitAdefic + female + height
                  + stunted + visit2 + visit3 + visit4  + visit5 + visit6+(1|idnum),
                  smoothZ=list(
                    grp=smZ(age)
                  ),
                  family = binomial,
                  data = indonRespir)
  sp_out3 <- semipred(fit3, "age", "age")
  library(ggplot2)
  print(sp_out3$plt +
    geom_rug(data = subset(indonRespir, respirInfec==0), sides = "b", col="deeppink") +
    geom_rug(data = subset(indonRespir, respirInfec==1), sides = "t", col="deeppink") +
    coord_cartesian(ylim = c(0, 0.2)))
}
#> Generalized linear mixed model fit by maximum likelihood (Laplace
#>   Approximation) [glmerMod]
#>  Family: binomial  ( logit )
#> Formula: respirInfec ~ age + vitAdefic + female + height + stunted + visit2 +  
#>     visit3 + visit4 + visit5 + visit6 + (1 | idnum)
#>    Data: indonRespir
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>  673.9676  740.1386 -323.9838  647.9676      1187 
#> Random effects:
#>  Groups Name        Std.Dev.
#>  idnum  (Intercept) 0.7937  
#>  grp    (Intercept) 0.4288  
#> Number of obs: 1200, groups:  idnum, 275; grp, 8
#> Fixed Effects:
#> (Intercept)          age    vitAdefic       female       height      stunted  
#>   -3.157550    -0.500128     0.666216    -0.505965    -0.033001     0.454787  
#>      visit2       visit3       visit4       visit5       visit6  
#>   -1.107518    -0.573018    -1.276192     0.472608    -0.005654  


# }
```
