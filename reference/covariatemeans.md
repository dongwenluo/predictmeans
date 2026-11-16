# Predicted Means of a Linear Model with Covariate Variable(s)

This function obtains predicted means with graph for a new set of
covariate values.

## Usage

``` r
covariatemeans(model, modelterm=NULL, covariate, as.is=FALSE, covariateV=NULL, 
  data=NULL, level=0.05, Df=NULL, trans=NULL, transOff=0, responsen=NULL, trellis=TRUE, 
  plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=FALSE,
  cov_reduce=mean, exclude=NULL,
  weights=c("equal", "proportional", "outer", "cells", "flat"),
  df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
  vcov_adjust=TRUE, contr=NULL, contr_scale=c("model", "response"),
  band=c("pointwise", "simultaneous"), bias_adjust=FALSE, sigma=NULL)
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

  A logic value to specify wheather or not using original covariate
  values' rage for graph, the default is FALSE.

- covariateV:

  A numeric vector when as.is is FALSE, then covariatemeans will produce
  the result for `covariate` at value of `covariateV`.

- data:

  In some cases, you need to provide the data set used in model fitting,
  especially when you have applied some variable trnasformation in the
  model.

- level:

  A significant level for calculating confident interval. The default
  value is 0.05. The limits are for each line on its own: whether the
  bands of two groups overlap is not a test of their difference (the
  lines are estimated from the same model and are correlated; see
  `letterCI` in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).
  Use `contr` for the band of the difference between two groups, or
  compare groups at chosen covariate values with
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  (`covariate = `) and their slopes with
  [`pm_slopes`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md).

- Df:

  A degree of freedom for calculating CI of predicted means (you can
  manually specified ddf here). For the above models, ddf is obtained
  from the function automatically.

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

- transOff:

  When you use `trans=exp(x+1)`, then `transOff=1`, the default is 0.

- responsen:

  Name (in "quotes") of the back transformed response variable in the
  `model`.

- trellis:

  A logical variable. If set to TRUE (default), a trellis plots of
  predicted means with CI will be drawn.

- plotord:

  A numeric (or character) vector specifying the order of plotting for
  two or three way interaction (e.g. `plotord = c(2, 1, 3)` or
  `plotord = c("B", "A", "C")` will put the second variable 'B' in
  `modelterm` on the `X` axis, the first variable 'A' as the grouping
  variable, and the third one 'C' as the panel variable). The defaults
  are `c(1, 2)` and `c(1, 2, 3)` for two and three way interactions.

- mtitle:

  The main title in the graph.

- ci:

  A logical variable to indicate whether to print confidence interval.
  The default value is TRUE.

- point:

  A logical variable to indicate whether to print raw data points. The
  default value is TRUE.

- jitterv:

  A degree of jitter in x and y direction in the graph. The default is
  zero.

- newwd:

  A logical variable to indicate whether to print graph in a new window.
  The default value is FALSE.

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
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md):
  by default Kenward-Roger df (each point of the curves has its own) and
  adjusted standard errors for models fitted by REML (`lmer`, `lme`,
  `gls`). predictmeans 1.x used the unadjusted covariance and, for
  `lmer`, one df per term; `vcov_adjust = FALSE` and, for `lme` and
  `gls`, `df_method = "containment"` give its results.

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

  With `contr`: `"pointwise"` (the default) gives limits at the level
  `level` for each covariate value on its own; `"simultaneous"` gives
  Scheffe limits that hold for all covariate values together, for saying
  over which range of the covariate the lines differ.

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

- plt:

  The plot (a `ggplot`).

- pltdf:

  The data of the plot: the fitted values with their limits, or with
  `contr` the differences on the model scale and, when there is a
  back-transformation, on the response scale.

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
    (`band = "pointwise"`); the lines differ at level `level` when it
    leaves out zero (one for a ratio).

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

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  data(Oats, package="nlme")
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  covariatemeans(fm, "Variety", covariate="nitro")

  covariatemeans(fm, "Variety", covariate="nitro", covariateV=seq(0, 0.6, 0.1))$pltdf

#>        Variety xvar      ses      Mean        LL        UL     factors
#> 1  Golden Rain  0.0 8.570716  81.90000  63.36196 100.43804 Golden Rain
#> 2  Golden Rain  0.1 8.150210  89.43333  71.40550 107.46117 Golden Rain
#> 3  Golden Rain  0.2 7.887152  96.96667  79.20735 114.72598 Golden Rain
#> 4  Golden Rain  0.3 7.797494 104.50000  86.82110 122.17890 Golden Rain
#> 5  Golden Rain  0.4 7.887152 112.03333  94.27402 129.79265 Golden Rain
#> 6  Golden Rain  0.5 8.150210 119.56667 101.53883 137.59450 Golden Rain
#> 7  Golden Rain  0.6 8.570716 127.10000 108.56196 145.63804 Golden Rain
#> 8   Marvellous  0.0 8.570716  90.41667  71.87862 108.95471  Marvellous
#> 9   Marvellous  0.1 8.150210  96.87500  78.84716 114.90284  Marvellous
#> 10  Marvellous  0.2 7.887152 103.33333  85.57402 121.09265  Marvellous
#> 11  Marvellous  0.3 7.797494 109.79167  92.11277 127.47056  Marvellous
#> 12  Marvellous  0.4 7.887152 116.25000  98.49069 134.00931  Marvellous
#> 13  Marvellous  0.5 8.150210 122.70833 104.68050 140.73617  Marvellous
#> 14  Marvellous  0.6 8.570716 129.16667 110.62862 147.70471  Marvellous
#> 15     Victory  0.0 8.570716  73.30000  54.76196  91.83804     Victory
#> 16     Victory  0.1 8.150210  81.40833  63.38050  99.43617     Victory
#> 17     Victory  0.2 7.887152  89.51667  71.75735 107.27598     Victory
#> 18     Victory  0.3 7.797494  97.62500  79.94610 115.30390     Victory
#> 19     Victory  0.4 7.887152 105.73333  87.97402 123.49265     Victory
#> 20     Victory  0.5 8.150210 113.84167  95.81383 131.86950     Victory
#> 21     Victory  0.6 8.570716 121.95000 103.41196 140.48804     Victory
# the difference between two varieties along nitro, with a simultaneous band
# (the two levels by name, or by number: contr=c(1, 3))
  covariatemeans(fm, "Variety", covariate="nitro", contr=c("Golden Rain", "Victory"),
                 band="simultaneous")$test

#> Lines equal along nitro: F(2, 21.6) = 0.515, p = 0.605.
#> The simultaneous band (level 0.05) includes no difference at every value of nitro shown.
#>      numDF      denDF          F    p.value 
#>  2.0000000 21.5970538  0.5146026  0.6048849 
```
