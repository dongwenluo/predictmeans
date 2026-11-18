# Predicted Means of a Linear Model with Covariate Variable(s)

Predicted means along a covariate for each level of a factor term
(curves) with their confidence limits, or the difference between two of
the curves with its band. The result holds the tables in the layout of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
and the plot, which `autoplot()` returns.

## Usage

``` r
covariatemeans(model, modelterm=NULL, covariate, as.is=FALSE, at=NULL, data=NULL,
  conf.level=0.95, trans=NULL, cov_reduce=mean, exclude=NULL,
  weights=c("equal", "proportional", "outer", "cells", "flat"),
  df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
  vcov_adjust=TRUE, contr=NULL, contr_scale=c("model", "response"),
  band=c("pointwise", "simultaneous"), bias_adjust=FALSE, sigma=NULL,
  mean_type=c("conditional", "overall"), plot=getOption("predictmeans.plot", FALSE),
  plot_args=list(), ...)
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

- modelterm:

  Name (in "quotes") for indicating which factor term's predicted mean
  to be calculated. The `modelterm` must be given exactly as it appears
  in the printed model, e.g. "A" or "A:B".

- covariate:

  Name (in "quotes") of one the covariate variables in the `model`.

- as.is:

  `TRUE`: each line only over the range of the covariate observed in its
  group (rows outside it are left out of the result); the default
  `FALSE` uses the values of `at` for every group.

- at:

  The values of `covariate` at which the lines are predicted
  (`covariateV` before 2.0); by default 50 values over the observed
  range. See
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).

- data:

  In some cases, you need to provide the data set used in model fitting,
  especially when you have applied some variable trnasformation in the
  model.

- conf.level:

  The confidence level of the limits, the default is 0.95
  (`level = 0.05` before 2.0). The limits are for each line on its own:
  whether the bands of two groups overlap is not a test of their
  difference (the lines are estimated from the same model and are
  correlated; see `letterCI` in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).
  Use `contr` for the band of the difference between two groups, or
  compare groups at chosen covariate values with
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  (`covariate = `) and their slopes with
  [`pm_slopes`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md).

- trans:

  A function object for calculating the back transformed means, e.g.
  `trans=exp`. When it is not given, the inverse link of the model is
  used, or the inverse of a transformation of the response in the model
  formula (`log(y)`, `sqrt(y)`, ...; see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)),
  with a message; `trans = I` keeps the model scale. The
  back-transformed line is the back-transformed fitted value. (Before
  2.0 half the standard error of the fitted value was added inside
  [`exp()`](https://rdrr.io/r/base/Log.html) for a log.)

- plot_args:

  A named list of plot settings (names may be abbreviated):

  `trellis`

  :   `TRUE` (default): one panel per level.

  `plotord`

  :   The order of the factors for two or three way terms (e.g.
      `c(2, 1)` or `c("B", "A")`): the first gives the colours, the
      second (and third) the panels. Giving it draws the plot.

  `mtitle`

  :   The main title.

  `ci`

  :   Draw the confidence band (default `TRUE`).

  `point`

  :   Draw the data points (default `TRUE`).

  `jitterv`

  :   Jitter of the points in x and y (default 0).

  Before 2.0 these were arguments of their own; given so they stop with
  this spelling.

- ...:

  Not used: arguments removed in 2.0 (`level`, `Df`, `covariateV`,
  `transOff`, `responsen`, `newwd`, the plot settings) stop with their
  2.0 spelling (see
  [`vignette("migrating")`](https://dongwenluo.github.io/predictmeans/articles/migrating.md)).

- cov_reduce:

  How covariates that are not varied are handled when building the
  predicted means: a function such as `mean` (the default) or `median`,
  or `"average"` to average the predictions over the observed covariate
  values (e.g. means adjusted for a spatial smooth in a `gam`). See
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).

- exclude:

  For `gam` models: smooth terms to leave out of the predictions. The
  default excludes random-effect smooths. See
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).

- weights:

  Passed to
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md):
  how the cells of the factors that are not in `modelterm` are weighted
  when averaging over them: `"equal"` (the default; equal weights),
  `"proportional"` (weights proportional to the frequencies of those
  cells in the data, the same for every mean, as SAS `LSMEANS / OM`),
  `"outer"` (products of each factor's own frequencies), `"cells"` (each
  mean weights the cells by their own frequencies) or `"flat"` (equal
  weights over the non-empty cells). The definitions are those of
  emmeans.

- df_method, vcov_adjust:

  As in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  (a number in `df_method` is a fixed df, `Df` before 2.0): by default
  Kenward-Roger df (each point of the curves has its own) and adjusted
  standard errors for models fitted by REML (`lmer`, `lme`, `gls`).
  predictmeans 1.x used the unadjusted covariance and, for `lmer`, one
  df per term; `vcov_adjust = FALSE` and, for `lme` and `gls`,
  `df_method = "containment"` give its results.

- contr:

  Two levels of `modelterm`, by their names as in the data (e.g.
  `c("Golden Rain", "Victory")`) or by their numbers (`c(1, 3)`; for a
  term with several factors the combined labels such as `"a:x"`): the
  difference between the first and the second along the covariate, with
  its confidence band, is returned and plotted instead of the lines. See
  ‘Comparing two lines’.

- contr_scale:

  With `contr`: what is plotted. `"model"` (the default) is the
  difference on the scale of the model, where the test is made; for a
  log link or a log response it is shown as the ratio of the two lines
  and for a logit link as the odds ratio (exact: exp of the difference
  and of its limits). `"response"` is the difference of the
  back-transformed lines, in the units of the response, with a
  delta-method band; it describes the size of the difference, the test
  stays on the model scale.

- band:

  With `contr`: `"pointwise"` (the default) gives limits at the
  confidence level `conf.level` for each covariate value on its own;
  `"simultaneous"` gives Scheffe limits that hold for all covariate
  values together, for saying over which range of the covariate the
  lines differ.

- mean_type:

  Only for `glmmTMB` fits with zero inflation (`ziformula`; ignored,
  with a message, for other models): `"conditional"` (default), curves
  of the means of the count part, \\\mu\\; `"overall"`, curves of the
  overall means \\(1 - \pi)\mu\\ (\\\pi\\ the probability of an extra
  zero), their bands \\\exp(\log m \mp crit\\SE)\\ with the SE of \\\log
  m\\ from the joint covariance of the count and zero-inflation
  coefficients (delta method), and with `contr` the ratio and the
  difference of the overall means. Log-link count families only.

- plot:

  Draw the plot now? The default is `FALSE` (since 2.0; 1.x drew it):
  the plot is in the result (`autoplot()` returns it). `plot = TRUE`, or
  `options(predictmeans.plot = TRUE)` for a whole session, draws it as
  1.x did; `plot_args = list(plotord = )` also draws it.

- bias_adjust, sigma:

  For back-transformed lines (a link, a transformed response or
  `trans`), as in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).
  By default a line is the back-transformed fitted value: the median of
  the response for a transformed response (the geometric mean for a
  log), and the response of a typical cluster (random effects at zero)
  for a generalized linear mixed model. `bias_adjust = TRUE` gives the
  mean of the response instead (for a mixed model: averaged over the
  clusters), computed over the variation on the model scale: the
  residual variance for a transformed response plus the variances of the
  random intercepts, or `sigma^2` when `sigma` is given (needed with
  random slopes). It is exact (for a log the factor \\exp(\sigma^2/2)\\,
  for a square root \\+ \sigma^2\\); `"taylor"` is the second-order
  formula of emmeans. Use it for totals and expected amounts; ratios and
  tests do not need it. The limits treat the variance as known. With
  `contr` the response-scale difference is that of the adjusted lines.

## Value

An object of class `c("covariatemeans", "pm_curves")`, a list. Its
tables have the layout of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
(columns `estimate`, `std.error`, `df`, `conf.low`, `conf.high`, ...,
the covariate under its own name; see
[`tidy.pm_curves`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md));
`tidy()` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
them, `autoplot()` the plot, and
[`print()`](https://rdrr.io/r/base/print.html) shows a summary.

- means:

  The curves: one row per level of `modelterm` and covariate value
  (`NULL` with `contr`).

- pairs:

  With `contr`: the difference of the two curves at each covariate
  value.

- settings:

  The options that produced them (df method, level, band,
  back-transformation, ...).

- plt:

  The plot (a `ggplot`); drawn only with `plot = TRUE`.

- pltdf:

  The data of the plot as in 1.x: the fitted values with their limits,
  or with `contr` the differences on the model scale and, when there is
  a back-transformation, on the response scale.

- test, significant:

  With `contr`: the joint F test that the two lines are the same
  (numerator and denominator df, F, p-value), and the ranges of the
  covariate where the band leaves out zero.

## Comparing two lines

Whether the bands of two lines overlap does not say whether the lines
differ: the lines are estimated from the same model and are correlated.
The band of their difference does, and it is computed from the contrast
of the two lines, not from their limits.

*The scale of the test.* The difference is tested on the scale of the
model (the link, or the transformed response): there the model is
additive, the estimates are closest to normal and the standard error
does not depend on the mean. \\d = (K_a - K_b)\hat\beta\\, \\SE^2 =
(K_a - K_b) V (K_a - K_b)'\\, limits \\d \pm c\\SE\\ (columns `diff`,
`ses`, `Df`, `LL`, `UL`). With a log the difference back-transforms
exactly to a ratio and with a logit to an odds ratio (columns `ratio` or
`odds.ratio` with their limits): the limits leave out 1 exactly when
those of `diff` leave out 0, so nothing is approximated. For other
transformations (square root, probit, a power) the back-transformed
difference has no meaning and the model scale is shown.

*Why back-transform.* Only to read the size of the difference in the
units of the response; it is not needed to decide whether the lines
differ. The response-scale difference \\h(\eta_a) - h(\eta_b)\\ is
another quantity: its standard error is an approximation (the delta
method, gradient \\h'(\eta_a) K_a - h'(\eta_b) K_b\\; columns
`response.diff`, `response.ses`, `response.LL`, `response.UL`), it
depends on the values of the other variables of the model, and lines
that are parallel on the model scale are not parallel on the response
scale. Its band can disagree with the test near the boundary, so with
`contr_scale = "response"` the covariate values where the lines differ
are still those of the model scale (marked along the axis). For a
generalized linear mixed model it is the difference for a typical
cluster (random effects at zero), as the lines are.

*Significance at a level.* Three questions have three answers:

1.  At one covariate value chosen beforehand: the pointwise band
    (`band = "pointwise"`); the lines differ at the level
    `1 - conf.level` when it leaves out zero (one for a ratio).

2.  Over which range of the covariate: looking at many values needs
    wider limits. `band = "simultaneous"` uses Scheffe's critical value
    \\c = \sqrt{r F\_{1-\alpha;\\r,\\df}}\\, \\r\\ the rank of the
    contrasts along the covariate (2 for straight lines: the
    Working-Hotelling band). It holds for every covariate value at once
    and is slightly conservative over a limited range.

3.  Whether the lines differ anywhere: the joint F test that all the
    differences along the covariate are zero (`test` in the result).

Each set of limits in `pltdf` has a column that says whether it includes
the value of no difference: `includes0` for `diff`, `ratio.includes1` or
`odds.ratio.includes1` (always equal to `includes0`) and
`response.includes0` for the response-scale difference (which can differ
from `includes0` near the boundary). Column `significant` (`!includes0`)
and `significant` in the result (the ranges of the covariate) follow the
band asked for, on the model scale.

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
`covariatemeans()`,
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md))
gives the overall means and their ratios.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  data(Oats, package="nlme")
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  cm <- covariatemeans(fm, "Variety", covariate="nitro")
  ggplot2::autoplot(cm)

  covariatemeans(fm, "Variety", covariate="nitro", at=seq(0, 0.6, 0.3))
#> Curves of 'Variety' along nitro
#> 95% limits; Kenward-Roger df
#>      Variety nitro estimate std.error     df conf.low conf.high
#>  Golden Rain   0.0    81.90     8.571 12.849    63.36    100.44
#>  Golden Rain   0.3   104.50     7.797  8.869    86.82    122.18
#>  Golden Rain   0.6   127.10     8.571 12.849   108.56    145.64
#>   Marvellous   0.0    90.42     8.571 12.849    71.88    108.95
#>   Marvellous   0.3   109.79     7.797  8.869    92.11    127.47
#>   Marvellous   0.6   129.17     8.571 12.849   110.63    147.70
#>      Victory   0.0    73.30     8.571 12.849    54.76     91.84
#>      Victory   0.3    97.62     7.797  8.869    79.95    115.30
#>      Victory   0.6   121.95     8.571 12.849   103.41    140.49
#> autoplot(x) draws the plot.
# the difference between two varieties along nitro, with a simultaneous band
# (the two levels by name, or by number: contr=c(1, 3))
  cd <- covariatemeans(fm, "Variety", covariate="nitro", contr=c("Golden Rain", "Victory"),
                       band="simultaneous")
  cd$settings$test
#>   df1      df2   F.ratio   p.value
#> 1   2 21.59705 0.5146026 0.6048849
  head(as.data.frame(cd))
#>        level1  level2      nitro estimate std.error       df statistic
#> 1 Golden Rain Victory 0.00000000 8.600000  8.684676 21.57462 0.9902500
#> 2 Golden Rain Victory 0.01224490 8.529592  8.567347 20.58273 0.9955931
#> 3 Golden Rain Victory 0.02448980 8.459184  8.453379 19.63654 1.0006867
#> 4 Golden Rain Victory 0.03673469 8.388776  8.342909 18.73710 1.0054976
#> 5 Golden Rain Victory 0.04897959 8.318367  8.236079 17.88507 1.0099912
#> 6 Golden Rain Victory 0.06122449 8.247959  8.133031 17.08076 1.0141310
#>     p.value     p.adj  conf.low conf.high adj.conf.low adj.conf.high
#> 1 0.3330287 0.6191042 -9.431527  26.63153    -14.22282      31.42282
#> 2 0.3310091 0.6162890 -9.309200  26.36838    -14.06393      31.12311
#> 3 0.3291509 0.6136450 -9.195202  26.11357    -13.91598      30.83435
#> 4 0.3274667 0.6111909 -9.089732  25.86728    -13.77922      30.55677
#> 5 0.3259693 0.6089456 -8.992966  25.62970    -13.65382      30.29056
#> 6 0.3246712 0.6069286 -8.905059  25.40098    -13.53996      30.03588
```
