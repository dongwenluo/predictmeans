# Predicted Curves and Surfaces of a Semiparametric Model

Predicted curves (or surfaces) of a semiparametric model along a
covariate, with pointwise or simultaneous bands; their slopes (`deriv`),
or the difference between two curves (`contr`). The result holds the
tables in the layout of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
(see
[`tidy.pm_curves`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md))
and the plot, which
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns.

## Usage

``` r
semipred(model, modelterm=NULL, covariate, sm_term=NULL, contr=NULL, at=NULL,
         boundary=NULL, conf_level=0.95, trans=NULL, contr_scale=c("model", "response"),
         band=c("pointwise", "simultaneous"),
         unconditional=getOption("predictmeans.semi_unconditional"), deriv=0,
         plot_args=list(), mean_type=c("conditional", "overall"),
         bias_adjust=FALSE, sigma=NULL, plot=getOption("predictmeans.plot", FALSE), ...)
```

## Arguments

- model:

  A fit of
  [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  or
  [`semireg_tmb`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  (`semireg` before 2.0).

- modelterm:

  Name (in "quotes") of the factor term whose curves are predicted, as
  in the model formula, e.g. "A" or "A:B" (any number of factors: every
  combination is computed; `autoplot(x, group = , panel = )` chooses the
  layout). With `NULL` or the name of the covariate, one curve without
  grouping. Factors that are not in `modelterm`, including `by` factors
  of the smooth terms, are averaged over with equal weights, in both the
  parametric and the smooth parts of the model.

- covariate:

  Name (in "quotes") of the covariate, or two names for a surface, e.g.
  `c("x1", "x2")`.

- sm_term:

  The smooths used in the prediction; the default `NULL` uses all of
  them. `print(model)` (or `tidy(model)`) lists the smooths of the fit
  with the values to use here. Each entry may be a term as written in
  the formula (`"s(pressure, by = QB)"`, all its pieces), a smooth as
  labelled (`"s(pressure):QB200"`, column `smooth`), a covariate
  (`"pressure"`: every smooth of it, column `covariate`) or a
  random-effect term (`"s_pressure_QB200"`, column `term`). Smooths of
  other covariates enter the curves as a constant; the advice
  `smooth_other_covariate` says when the prediction includes such
  smooths. `settings$sm_term` records the blocks used.

- contr:

  Two levels of `modelterm` (numbers or names, e.g. `c(4, 1)` or
  `c("d", "a")`): the difference of their curves (level 4 minus level
  1), with limits, along `covariate`.

- at:

  The values of `covariate` to predict at (`covariateV` before 2.0): a
  numeric vector (one covariate) or a matrix with one named column per
  covariate (two covariates, e.g. for a `Ztps` smooth); each row of the
  matrix is one prediction point.

- boundary:

  For two covariates (a surface): the region where predictions are
  shown. A matrix or data frame of two columns (the covariates in the
  order of `covariate`), the vertices of a polygon, or of several rings
  separated by a row of `NA` (a point is inside when it is inside an odd
  number of them: separate pieces and holes), such as
  [`semiboundary()`](https://dongwenluo.github.io/predictmeans/reference/semiboundary.md)
  returns; or `TRUE` for
  [`semiboundary()`](https://dongwenluo.github.io/predictmeans/reference/semiboundary.md)
  of the model's data (the region the data cover, following their
  shape). `NULL` (default): the whole grid.

- conf_level:

  Confidence level of the two-sided limits, default 0.95 (`level = 0.05`
  before 2.0; before 2.0 the limits used `qnorm(1 - level)`, i.e. 90%).
  For gaussian fits the limits are t limits with Satterthwaite df at
  each point (column `df`; see Details); for generalized fits t limits
  with the df of the design's strata (inner-outer, as for GLMMs), normal
  limits when the model has no random term besides the smooths.

- trans:

  A function that back-transforms the curves, e.g. `trans = exp`.
  Without it the curves are back-transformed by the inverse link of a
  generalized fit, or by the inverse of a response transformed in the
  formula (`log(y)`, `sqrt(y)`, `log(y + 1)`, ...), as in
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  (one rule for both); `trans = I` keeps the model scale.

- bias_adjust, sigma:

  For back-transformed curves, as in
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md):
  by default a curve is the back-transformed fitted value (the median
  for a transformed response, a typical cluster for a GLMM);
  `bias_adjust = TRUE` gives the mean of the response, over the residual
  variance (transformed response) and the variances of the random
  intercepts other than the smooths (the smooths' random effects are the
  curves themselves), or `sigma^2`.

- mean_type:

  Only for `glmmTMB` fits with zero inflation (`ziformula`; ignored,
  with a message, for other models): `"conditional"` (default), curves
  of the means of the count part, \\\mu\\; `"overall"`, curves of the
  overall means \\(1 - \pi)\mu\\ (\\\pi\\ the probability of an extra
  zero), their bands \\\exp(\log m \mp crit\\SE)\\ with the SE of \\\log
  m\\ from the joint covariance of the count and zero-inflation
  coefficients (delta method), and with `contr` the ratio and the
  difference of the overall means. Log-link count families only.

- plot_args:

  Superseded by the arguments of
  [`autoplot.pm_curves`](https://dongwenluo.github.io/predictmeans/reference/autoplot.pm_curves.md),
  kept for 1.x scripts (before 2.0 arguments of their own; given so they
  stop with this spelling). A named list: `trellis` (`TRUE`, default:
  one panel per level), `scales` (`"fixed"`, default, `"free"`,
  `"free_x"`, `"free_y"`: the axes of the panels), `plotord` (the order
  of the two factors of a two-way term, e.g. `c(2, 1)`; giving it draws
  the plot), `ci` (draw the bands, `TRUE`), `point` (draw the data,
  `TRUE`), `jitterv` (jitter of the points, 0), `threeD` (a 3-D surface
  with plotly, `FALSE`); and colours and sizes, any of: `line` and
  `band` (one curve; dark and light blue), `diff_line` and `diff_band`
  (a difference of two curves, `contr`; dark and light green), `ref`
  (the dashed reference line at 0, or 1 for ratios, and the marks under
  a difference where its band leaves out the reference; red), `alpha`
  (of the bands, 0.35), `linewidth` (1), `points`, `point_alpha`,
  `point_size` (the data: grey, 0.6, 1.3), `palette` (curves by group:
  eight colours that colour-blind readers can tell apart; ggplot2's own
  beyond eight groups), `surface_points` (the data on surfaces: red).
  Bands are drawn first and their lines on top. For example
  `plot_args = list(line = "black", band = "grey70")`; the plot is a
  ggplot, so anything else can be added to it.

- plot:

  Draw the plot now? The default is `FALSE` (since 2.0; 1.x drew it):
  the plot is in the result
  ([`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  returns it). `plot = TRUE`, or `options(predictmeans.plot = TRUE)` for
  a whole session, draws it; so does `plot_args = list(plotord = )`.

- ...:

  Not used: arguments removed in 2.0 (`semireg`, `covariateV`, `level`,
  `prt`, `simultaneous`, the plot settings) stop with their 2.0 spelling
  (see
  [`vignette("migrating")`](https://dongwenluo.github.io/predictmeans/articles/migrating.md)).

- band:

  `"pointwise"` (the default) or `"simultaneous"`: simultaneous
  confidence bands over the whole curve(s) instead of pointwise limits;
  the critical value is the `1 - level` quantile of the largest
  standardised deviation over the predicted points, from 10000 draws of
  the estimates' distribution (Ruppert, Wand and Carroll, 2003); for
  gaussian fits a multivariate t with the df of each point (Details).
  The same argument as in
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
  where lines of a parametric model use Scheffe's critical value
  instead; for smooths that would be far too wide.

- unconditional:

  If `TRUE`, the covariance of the coefficients includes the uncertainty
  of the smoothing parameters (the variances of the smooth terms), as
  mgcv's `unconditional = TRUE`: the first-order correction \\J
  V\_\theta J'\\ is added, \\J\\ the derivative of the estimated
  coefficients in the variance parameters and \\V\_\theta\\ the
  covariance of their estimates (from the Hessian of the deviance; for
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  fits TMB's joint precision of all parameters). The default (`NULL`) is
  `TRUE` for generalized models and `FALSE` for gaussian ones: in a
  simulation of the difference between two curves
  (`dev/sim-curve-difference.R`, 150 data sets, nominal 0.95) the bands
  without it covered 0.93 pointwise and 0.89-0.91 over the whole curve
  for Poisson and binary data, and 0.955-0.957 and 0.94-0.95 with it;
  for gaussian data 0.97 without it (0.98 with it).
  `options(predictmeans.semi_unconditional = TRUE)` (or `FALSE`) sets it
  for `semipred()` and
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  on `semireg` fits.

- deriv:

  `0` (the default): the curves; `1`: their slopes, `2`: their
  curvatures (first and second derivatives along `covariate`, one
  covariate), with bands as for the curves, and with `contr` the
  difference of the slopes (curvatures) of two curves. See ‘Slopes and
  curvature’.

- contr_scale:

  With `contr` and a generalized model: `"model"` (the default) plots
  the difference on the link scale, as the ratio (log link) or the odds
  ratio (logit link) of the two curves; `"response"` the difference of
  the back-transformed curves with a delta-method band. See ‘Comparing
  two lines’ in
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md):
  the test is on the model scale.

## Value

An object of class `c("semipred", "pm_curves")`, a list. `tidy()` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
its table,
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
the plot, and [`print()`](https://rdrr.io/r/base/print.html) shows a
summary (see
[`tidy.pm_curves`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md));
notes in `attr(, "advice")`
([`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).

- means:

  The curves: the factor columns of `modelterm`, the covariate(s),
  `estimate`, `std.error`, `df`, `conf.low`, `conf.high` on the model
  scale (pointwise limits), with a simultaneous band `adj.conf.low`,
  `adj.conf.high`; for generalized fits and with `trans` also
  `response`, `response.std.error`, `response.conf.low`,
  `response.conf.high` (and `response.adj.conf.*`). `NULL` with `contr`
  or `deriv`.

- slopes:

  With `deriv = 1` (`2`): the slopes (curvatures) in the same columns,
  with `statistic`, `p.value` and, with a simultaneous band, `p.adj`.

- pairs:

  With `contr`: the difference of the two curves, as
  `covariatemeans()$pairs`; with a simultaneous band `p.adj` is the
  share of the simulated maxima at least as large (below 0.05 exactly
  where the band leaves out zero; the smallest value is 1/10000).

- settings:

  The options: `df_method`, `band`, `deriv`, `unconditional`,
  `back_transform`, ...; with `contr`, `test`: whether the curves differ
  anywhere along the covariate, the max-\|t\| test (`statistic` the
  largest \|t\| over the `points`, `p.value` the share of simulated
  maxima at least as large, from the covariance of the fixed and smooth
  coefficients: the smallest level at which the simultaneous band leaves
  out zero somewhere). A Wald F test of the smooth coefficients is not
  used: they are shrunken, and the statistic would not have an F
  distribution. In a simulation of equal curves (200 data sets) it
  rejected 3.0% at the 5% level.

- significant:

  With `contr` or `deriv`: the ranges of the covariate where the band
  asked for leaves out zero.

- plt:

  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  of the result with `plot_args`: a ggplot (one plotly widget for
  `threeD`); drawn only when `plot = TRUE`.

- pred_df:

  The 1.x table: a data.frame with the predicted values, their limits
  `LL` and `UL` (the band asked for), `SE` and `df`; for generalized
  fits the means and limits are on the response scale (limits
  `h(eta -+ crit SE)`), `SE` on the link scale and `SE_mu` on the
  response scale (delta method), `df` when the limits are t; with
  `contr` the columns described in Details.

## Details

Degrees of freedom (gaussian fits, both engines): a point of a curve, or
of a difference of curves, is a combination of fixed and smooth (random)
coefficients. Its df are Satterthwaite's for that combination (as
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
on the fit and SAS for estimates that involve random effects), counting
the variance components of the design (blocks, plots, ...) and the
residual; the smoothing parameters are held at their estimates (their
uncertainty is what `unconditional` adds), and components estimated at
zero are left out. So the df are the residual df when there is no other
random term, and few where a design stratum dominates (curves of
treatments applied to whole plots). Pointwise limits use
`qt(1 - level/2, df)` at each point. Simultaneous bands use a
multivariate t: the simulated deviations are divided by variance
estimates that move together over the points (each point with its own
df), the points are put on one scale by their t probabilities, and each
point gets the t quantile of its df at the simulated `1 - level`
quantile (with one df for all points: the exact multivariate t).

Generalized fits (both engines) take the df that
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
gives GLMMs: the inner-outer df of the stratum of the design's random
effects in which the fixed coefficients of a point vary (a curve of a
treatment applied to whole plots: the whole-plot df), the smooth blocks
left out (they are not strata of the design); with no random term
besides the smooths, normal limits (as for a `glm` with a known
dispersion). In a simulation (`dev/sim-semireg-glmm-df.R`; Poisson
counts in the split plot above, 1000 data sets per setting) normal
limits rejected 5.8-9.2% of true equalities of two treatments at a time
and covered a treatment's curve 0.91-0.93, the t limits (10 df) 3.1-5.9%
and 0.955-0.970; bands of differences covered 0.92-0.96 (normal) and
0.95-0.996 (t).

In a simulation (`dev/sim-semipred-df.R`; split plot, 6 blocks, 3
treatments on whole plots, curves over 11 times; 500 data sets per
setting) normal limits of a treatment's curve covered 0.91 pointwise and
0.92 over the whole curve, these t limits (about 8-11 df) 0.94-0.95 and
0.96-0.97; for differences of curves (15-31 df) normal 0.94-0.96, t
0.96-0.97. With curves by group and no design strata both covered
0.95-0.98.

With `contr` the difference between two curves is computed as in
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
from the contrast of the two curves and the covariance of the fixed and
smooth coefficients. For a linear model `pred_df` has the difference (in
the column of the response), `LL`, `UL`, `SE_contr` and the indicators
`includes0` (the limits include zero) and `significant`. For a
generalized model it has the columns of `covariatemeans(contr = )`:
`diff`, `ses`, `LL`, `UL` and `includes0` on the link scale; the `ratio`
or `odds.ratio` with its limits and `..includes1`; `response.diff`,
`response.ses`, `response.LL`, `response.UL` and `response.includes0` on
the response scale; and `significant` (the link scale). Before 2.0 the
band of a generalized model's difference used the standard error of the
link-scale difference multiplied by a weight evaluated at the difference
of the means, which is not a standard error of any difference. The same
difference is in the table `pairs` in the layout of
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).
The ranges of the covariate where the curves differ are printed with the
result and kept in `significant` and `attr(pred_df, "significant")`;
with `band = "simultaneous"` they hold for the whole curve, and the
curves differ somewhere at level `1 - conf_level` when there is such a
range. (A joint F test, as
[`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
gives, is not given: Wald tests of penalised smooth coefficients are not
reliable.)

## Slopes and curvature

With `deriv = 1` (`2`) each row of the prediction design is replaced by
its first (second) derivative along the covariate, by central
differences of the whole design (the fixed part and every smooth basis,
on the same grid), so the slopes are those of the fitted curves for any
basis; for the O'Sullivan basis (`bs = "os"`, 1.x
[`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md))
they equal the exact derivatives of the basis (`ZOSull(drv = )`), to
1e-8. Standard errors, df and pointwise or simultaneous bands are
computed as for the curves; generalized models give the slopes of the
linear predictor (link scale). `pred_df` holds the slopes in the column
of the response, with `includes0` and `significant`; the ranges of the
covariate where the band lies above or below zero (the curve rises or
falls) are printed with the result and kept in
`attr(pred_df, "significant")`. A covariate scaled by
`semireg(nscale = TRUE)` is differentiated in its own units. The table
`slopes` of the result has the slopes with `statistic` and `p.value` of
slope = 0 (and `p.adj` with a simultaneous band), and `significant` the
ranges. In a simulation (`dev/sim-semipred-deriv.R`; \\y = \sin(x/1.5)\\
plus noise, 150 points, 500 data sets) the bands covered the true slope
0.97 pointwise and over the whole range with the default 6 knots,
0.99-1.00 with 15; the curvature 0.98-1.00.

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
`semipred()`) gives the overall means and their ratios.

## References

Ruppert, D., Wand, M. P. and Carroll, R. J. (2003), *Semiparametric
Regression*, Cambridge University Press.

## See also

[`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md),
[`tidy.pm_curves`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md),
[`semiboundary`](https://dongwenluo.github.io/predictmeans/reference/semiboundary.md),
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
data(Dialyzer, package = "nlme")
library(ggplot2)
ggplot(Dialyzer, aes(x = pressure, y = rate, col = QB)) +
  geom_line() +
  facet_wrap(vars(Subject))


# a curve per blood flow rate, a random intercept and slope per dialyzer
fm <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
              data = Dialyzer)
residplot(fm$semer, group = "QB")

anova(fm)   # fixed effects (Kenward-Roger) and a test of each smooth term
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

# the curves with pointwise bands (t limits with Satterthwaite df), then simultaneous ones
out <- semipred(fm, "QB", "pressure")
out
#> Curves of 'QB' along pressure
#> 95% limits; Satterthwaite (fixed + smooth coefficients) df
#>   QB pressure estimate std.error    df conf.low conf.high
#>  200    0.235    2.720    1.0383 40.53   0.6225     4.818
#>  200    0.240    3.030    1.0224 38.63   0.9612     5.099
#>  200    0.245    3.340    1.0080 36.95   1.2971     5.382
#>  200    0.250    3.650    0.9950 35.47   1.6305     5.669
#>  200    0.255    3.960    0.9833 34.17   1.9616     5.957
#>  200    0.260    4.270    0.9727 33.02   2.2906     6.248
#>  200    0.280    5.508    0.9396 29.67   3.5886     7.428
#>  200    0.285    5.818    0.9332 29.07   3.9093     7.726
#>  200    0.290    6.127    0.9275 28.55   4.2286     8.025
#>  200    0.292    6.253    0.9252 28.35   4.3586     8.147
#> ... 220 more rows: as.data.frame(x) gives the table.
#> autoplot(x) draws the plot (layouts: ?autoplot.pm_curves).
autoplot(out)

# \donttest{
autoplot(semipred(fm, "QB", "pressure", band = "simultaneous"), panel = character(0))

# the smooths of one term only, as written in the formula
semipred(fm, "QB", "pressure", sm_term = "s(pressure, by = QB)")$settings$sm_term
#> [1] "s_pressure_QB200" "s_pressure_QB300"
# at given pressures
semipred(fm, "QB", "pressure", at = c(0.5, 1, 1.5, 2))
#> Curves of 'QB' along pressure
#> 95% limits; Satterthwaite (fixed + smooth coefficients) df
#>   QB pressure estimate std.error    df conf.low conf.high
#>  200      0.5    18.28    0.9329 34.54    16.38     20.17
#>  200      1.0    36.30    0.9245 36.55    34.43     38.18
#>  200      1.5    44.34    0.9639 29.54    42.37     46.31
#>  200      2.0    46.06    1.1326 28.62    43.74     48.38
#>  300      0.5    18.14    0.9227 33.10    16.27     20.02
#>  300      1.0    38.81    0.9263 36.80    36.94     40.69
#>  300      1.5    50.58    0.9626 29.39    48.61     52.55
#>  300      2.0    56.64    1.1303 28.41    54.32     58.95
#> autoplot(x) draws the plot (layouts: ?autoplot.pm_curves).

# the difference of the two curves, and where they differ
dif <- semipred(fm, "QB", "pressure", contr = c("300", "200"), band = "simultaneous")
dif$significant
#>       from   to direction
#> 1 1.204694 3.03         >
head(as.data.frame(dif))
#>   level1 level2 pressure   estimate std.error       df  statistic   p.value
#> 1    300    200    0.235 -1.0449653  1.527457 45.91274 -0.6841208 0.4973376
#> 2    300    200    0.240 -1.0298933  1.501755 43.65666 -0.6857933 0.4964694
#> 3    300    200    0.245 -1.0147676  1.478069 41.61607 -0.6865494 0.4961745
#> 4    300    200    0.250 -0.9995865  1.456254 39.77410 -0.6864096 0.4964372
#> 5    300    200    0.255 -0.9843485  1.436166 38.11365 -0.6854001 0.4972393
#> 6    300    200    0.260 -0.9690517  1.417671 36.61803 -0.6835519 0.4985610
#>    p.adj  conf.low conf.high adj.conf.low adj.conf.high
#> 1 0.9806 -4.119735  2.029805    -5.310165      3.220235
#> 2 0.9802 -4.057153  1.997367    -5.233305      3.173518
#> 3 0.9801 -3.998448  1.968912    -5.161730      3.132195
#> 4 0.9801 -3.943306  1.944133    -5.094995      3.095822
#> 5 0.9806 -3.891430  1.922733    -5.032675      3.063979
#> 6 0.9810 -3.842538  1.904434    -4.974370      3.036266

# slopes of the curves, and where they differ from zero over the whole range
slp <- semipred(fm, "QB", "pressure", deriv = 1, band = "simultaneous")
slp$significant
#>   level  from       to direction
#> 1   200 0.235 1.661020       > 0
#> 2   300 0.235 2.060306       > 0

# a surface: copper grade over the plane and a smooth in depth (data(copper))
data(copper)
fcu <- semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6),
               data = copper)
autoplot(semipred(fcu, covariate = "zcoord", trans = exp))

autoplot(semipred(fcu, covariate = c("xcoord", "ycoord"), trans = exp))


# a binary response: the probability scale, and the odds ratio of two curves
set.seed(4)
db <- expand.grid(x = seq(0, 1, length = 25), g = factor(c("a", "b")), cl = factor(1:12))
db$y <- rbinom(nrow(db), 1, plogis(-1 + 2*sin(pi*db$x) + 0.8*(db$g == "b")*db$x +
                                     rnorm(12, 0, 0.5)[db$cl]))
fbin <- semireg(y ~ s(x, by = g, bs = "os", k = 6) + (1 | cl), family = binomial, data = db)
autoplot(semipred(fbin, "g", "x"))

head(as.data.frame(semipred(fbin, "g", "x", contr = c("b", "a"))))
#>   level1 level2          x    estimate std.error  df   statistic   p.value
#> 1      b      a 0.00000000 -0.45544553 0.6498042 585 -0.70089661 0.4836462
#> 2      b      a 0.02040816 -0.33512989 0.5511548 585 -0.60805035 0.5433901
#> 3      b      a 0.04081633 -0.22158170 0.4722011 585 -0.46925283 0.6390638
#> 4      b      a 0.04166667 -0.21701785 0.4693481 585 -0.46238142 0.6439797
#> 5      b      a 0.06122449 -0.11618598 0.4138089 585 -0.28077204 0.7789845
#> 6      b      a 0.08163265 -0.02032772 0.3767891 585 -0.05394987 0.9569935
#>     conf.low conf.high odds.ratio odds.ratio.conf.low odds.ratio.conf.high
#> 1 -1.7316787 0.8207876  0.6341654           0.1769871             2.272289
#> 2 -1.4176131 0.7473533  0.7152452           0.2422917             2.111404
#> 3 -1.1489976 0.7058342  0.8012505           0.3169543             2.025536
#> 4 -1.1388304 0.7047947  0.8049156           0.3201933             2.023431
#> 5 -0.9289181 0.6965461  0.8903096           0.3949808             2.006809
#> 6 -0.7603519 0.7196965  0.9798775           0.4675019             2.053810
#>       response response.std.error response.conf.low response.conf.high
#> 1 -0.085658253         0.11998578        -0.3213136          0.1499971
#> 2 -0.068798133         0.11203980        -0.2888474          0.1512511
#> 3 -0.048869676         0.10372380        -0.2525861          0.1548467
#> 4 -0.047989164         0.10338662        -0.2510433          0.1550650
#> 5 -0.027076246         0.09632543        -0.2162620          0.1621095
#> 6 -0.004923197         0.09124626        -0.1841334          0.1742870
# }
```
