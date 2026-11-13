# Predicted Means of a Linear Model

This function obtains predicted means, SE of means, SED of means, LSDs
and plots of means with SE bar or LSD bar for parametric models such as
`aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`. The function also perfomrs
pairwise comparisons and permutation tests.

## Usage

``` r
predictmeans(model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL,
  lsd_bar=TRUE, level=NULL, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = I,
  transOff = 0, responsen=NULL, count=FALSE, plotord=NULL, lineplot=TRUE, plottitle=NULL, 
  plotxlab=NULL, plotylab=NULL, which.plots = c("mean", "pval", "back"),
  mplot=NULL, barplot=NULL, pplot=NULL, bkplot=NULL, plot=getOption("predictmeans.plot", TRUE),
  jitterv=0.2, basesz=12, prtnum=TRUE, prtplt=TRUE, newwd=FALSE, permlist=NULL, 
  ncore=1, ndecimal=4, cov_reduce=mean, exclude=NULL, vcov_adjust=TRUE,
  weights=c("equal", "proportional", "outer", "cells", "flat"),
  df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"), vcov=NULL,
  equivalence=NULL, bias_adjust=FALSE, sigma=NULL, average=c("link", "response"))
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`,
  `glmer`, `glmmTMB` and `gam` (package mgcv), or a
  [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fit. An `aov` model with `Error()` strata is analysed by strata, as
  Genstat's ANOVA: each comparison uses the mean squares of the strata
  it draws on, with Satterthwaite's df across strata
  (`settings$df_method` `"multistratum"`); standard errors of single
  means include the block stratum. With missing values in an orthogonal
  design (rows whose response is `NA`, or a few combinations of the
  design factors absent) they are estimated as Genstat's ANOVA does (the
  residual sum of squares of the lowest stratum made smallest; the
  residual df reduced) and the completed design is analysed by strata,
  with exact standard errors (Genstat's are those of the complete
  design); a message says that REML
  ([`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) with the `Error()`
  terms as random effects) uses the information in all strata and is the
  better analysis when plots are missing. For a `semireg` fit the means
  include the smooth terms (as
  [`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)),
  at the covariates' means or at the values of `covariate` (a named
  vector, e.g. `c(construction.date = 1980)`); the result holds only the
  `means`, `pairs` and `settings` tables, with normal-based tests and
  any `adj`.

- modelterm:

  Name (in "quotes") for indicating which factor term's predicted mean
  to be calculated. The `modelterm` must be factors and given exactly as
  it appears in the printed model, e.g. "A" or "A:B".

- data:

  In some cases, you need to provide the data set used in model fitting,
  especially when you have applied some variable trnasformation in the
  model.

- pairwise:

  An option for showing pair-wise LSDs and p-values, or not. The default
  is FALSE.

- atvar:

  When `pairwise = TRUE`, a quoted name indicating within levels of
  which variable in `modelterm` the multiple comparison will be
  performed.

- adj:

  Name (in "quote") for indicating a method for adjusting p-values of
  pairwise comparisons. The choices are "none", "tukey", "holm",
  "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "sidak",
  "scheffe" and "mvt". "sidak", "scheffe" (rank of the family: number of
  means - 1) and "mvt" (single-step multivariate t, the exact
  counterpart of SAS `ADJUST=SIMULATE`; needs package mvtnorm) are as in
  emmeans. Adjusted p-values use each pair's own df (as SAS
  `ADJDFE=ROW`; SAS's default `ADJDFE=SOURCE` uses the df of the term's
  F test instead, which differs when the pairs' df differ). With
  `permlist`, "tukey" and "mvt" give the single-step max-\|t\|
  adjustment of the permutation p-values (the permutation version of
  Tukey's test: each pair's \|t\| is compared with the largest \|t\| of
  the family in each permuted data set), and "westfall" the step-down
  max-\|t\| of Westfall and Young (1993), which is uniformly more
  powerful; both use the correlation between the comparisons. "westfall"
  needs `permlist`; "scheffe" is not available with it. The default
  method is "none"; letter case does not matter ("Tukey" works).
  Confidence limits that agree with the adjusted p-values are in the
  columns `adj.conf.low` and `adj.conf.high` of `pairs`. The printed LSD
  is adjusted only by "bonferroni".

- Df:

  A degree of freedom for calculating LSD. For the above models, Df is
  obtained from the function automatically.

- lsd_bar:

  A logical variable to indicate to print an average LSD or SED bar on
  the means plot. The default is TRUE.

- level:

  A significant level for calculating LSD, CI etc. The default value is
  0.05.

- covariate:

  A numerical vector to specify values of covariates for calculating
  predicted means. The default values are the means of the associated
  covariates.

- meandecr:

  A logical variable to indicate whether to print letters for multiple
  comparisons by decreasing order of means in the mean_table. The
  default is NULL which indicates the mean order follows the associated
  factor levels.

- letterCI:

  A logical variable to indicate printed letters for multiple
  comparisons by whether or not CI overlap in the mean_table. The
  default is FALSE. With `level = NULL` the CIs are 83.4% CIs: for two
  uncorrelated means with equal SEs they fail to overlap when the means
  differ at about the 5% level. Otherwise overlap only approximates the
  pairwise test: it finds fewer differences when the means are
  positively correlated, as in most mixed models, and more when their
  SEs differ (Maghsoodloo and Huang, 2010). With `letterCI = TRUE`,
  `level` sets the CIs (e.g. `level = 0.05` gives 95% CIs), `adj` is
  "none", and with `atvar` the letters compare the means within each
  `atvar` group. The note under the mean table gives the CI level and
  the approximate test level.

- trans:

  A function object for calculating the back transformed means, e.g.
  `trans=exp`.

- transOff:

  When you use `trans=exp(x+1)`, then `transOff=1`, the default is 0.

- responsen:

  Name (in "quotes") of the back transformed response variable in the
  `model`.

- count:

  An option for indicating the back transformed mean values are counts
  or not. The default is FALSE.

- plotord:

  A numeric (or character) vector specifying the order of plotting for
  two or three way interaction (e.g. `plotord = c(2, 1, 3)` or
  `plotord = c("B", "A", "C")` will put the second variable 'B' in
  `modelterm` on the `X` axis, the first variable 'A' as the grouping
  variable, and the third one 'C' as the panel variable). The defaults
  are `c(1, 2)` and `c(1, 2, 3)` for two and three way interactions.

- lineplot:

  An option for drawing a line chart, or dot chart. The default is TRUE.

- plottitle:

  A character vector specifying the main title for plot(s). The default
  is NULL.

- plotxlab:

  A character vector specifying the x label for plot(s). The default is
  NULL.

- plotylab:

  A character vector specifying the y label for plot(s). The default is
  NULL.

- which.plots:

  Which plots to draw: any of `"mean"` (means plot), `"bar"` (bar
  chart), `"pval"` (p-value plot, when there are more than six p-values)
  and `"back"` (back-transformed plot). The default is
  `c("mean", "pval", "back")`.

- mplot:

  Deprecated, use `which.plots`. An option for drawing a means plot, or
  not.

- barplot:

  Deprecated, use `which.plots`. An option for drawing a bar chart, or
  not.

- pplot:

  Deprecated, use `which.plots`. An option for drawing a p-values plot,
  or not.

- bkplot:

  Deprecated, use `which.plots`. An option for drawing back transformed
  plot, or not.

- plot:

  An option for drawing plots, or not. The default is TRUE, unless
  `options(predictmeans.plot = FALSE)` is set (e.g. for scripts, reports
  and pipelines that use the `means` and `pairs` tables only).

- jitterv:

  A degree of jitter in x and y direction in the back transformed means
  graph. The default is zero.

- basesz:

  The base font size. The default is 12.

- prtnum:

  An option for printing covariate information on the screen, or not.
  The default is TRUE.

- prtplt:

  An option for printing plots on the screen, or not. The default is
  TRUE.

- newwd:

  A logical variable to indicate whether to print graph in a new window.
  The default is FALSE.

- permlist:

  A model parameter list produced by the function `permmodels`. When
  `permlist != NULL`, the option `Df` will be non-functional. This is a
  key option for pairwise comparisons via permutation tests. With the
  result of `permmodels(method = "freedman-lane")` (`lm` and `aov`
  models), each pair is tested by Freedman-Lane permutations of the
  model in which that pair's means are equal. With a `permlist`,
  `vcov_adjust` is not used: the observed t statistics are compared with
  those of the permuted refits, which use the unadjusted covariance.

- ncore:

  Number of cores for parallel computing when `permlist != NULL`, the
  default value is 1 (before predictmeans 2.0: 3).

- ndecimal:

  An option for specifying number of decimal point to be print at
  predicted means table. The default is 4.

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

- vcov_adjust:

  For models fitted by REML (`lmer`, `lme`, `gls`) with Kenward-Roger
  df: if `TRUE` (the default), standard errors use the Kenward-Roger
  adjusted covariance matrix of the fixed effects, as SAS `DDFM=KR` and
  emmeans do, to go with the Kenward-Roger degrees of freedom. `FALSE`
  uses the unadjusted covariance and reproduces results of predictmeans
  1.x. It has no effect for other models.

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
  emmeans. Empty cells (combinations of interacting factors without
  data, not those absent by design, as for nested factors): with
  `"equal"`, a mean that averages over them is the average of its
  observed cells (as predictmeans 1.x, and ASReml and Genstat with
  `PRESENT`), with a warning naming those means. Two such means average
  over different cells, so their difference would also carry the effects
  of the other factors; a pair with such a mean is therefore compared on
  the cells both means have (column `over` of `pairs`; NA,
  `over = "none"`, when they have none). Only factors whose effects the
  two means share are matched: with `Sample_date/Farm` in the model, two
  dates are not matched on farms (each date has its own farm effects),
  but two treatments at one date are. In simulations comparing the
  averages rejected 8-78% of true null hypotheses at the 5% level when
  the other factor had an effect, the common cells 4-6%
  (`dev/review-defaults.md` in the package sources).
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  and
  [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
  do the same.

- df_method:

  Denominator degrees of freedom for `lmer` models. `"auto"` (the
  default) uses Kenward-Roger for REML fits and Satterthwaite for ML
  fits. Kenward-Roger uses the observed information of the variance
  parameters, as SAS PROC MIXED `DDFM=KR`; it stays valid with missing
  values, where the expected information (`"KR-expected"`, as pbkrtest
  and emmeans) is not: in simulations (`dev/sim-df/` in the package
  sources) F tests of repeated measures with missing values rejected 6%
  of true hypotheses at the 5% level with the observed information and
  8.5% with the expected. The two agree in balanced designs.
  `"KR-average"` uses the average information of the REML iterations of
  Genstat and ASReml (Gilmour, Thompson and Cullis 1995) and gives
  Genstat's REML df; Genstat's SEDs are the unadjusted ones
  (`vcov_adjust = FALSE`). `"containment"` gives the SAS `DDFM=CONTAIN`
  df: the smallest rank contribution to \\\[X Z\]\\ of the random
  effects that contain `modelterm`, or the residual df when none does;
  all means and comparisons of the term then share that df, and the
  unadjusted covariance is used unless `vcov_adjust = TRUE` is given.
  `"between-within"` gives the SAS `DDFM=BETWITHIN` df, for `lmer`,
  `lme` and `gls` models: the subject is the coarsest random-effect
  grouping (`lmer`, `lme`) or the grouping of the correlation structure
  (`gls`); an effect constant within subjects gets the between-subject
  df, and the others the df of their subject-by-effect interaction (for
  example Block x Variety in a split plot); with an unstructured
  covariance (`gls` with `corSymm`) and no random effects every effect
  gets the between-subject df, as in SAS. Checked against SAS PROC
  MIXED. As for containment, the unadjusted covariance is used unless
  `vcov_adjust = TRUE` is given. For `lme` and `gls` models fitted by
  REML, `"auto"` is Kenward-Roger too (with adjusted standard errors),
  for any correlation and variance structure: in simulations the
  residual df (`gls`) rejected 16% of true hypotheses in F tests of
  unstructured repeated measures (24% with missing values),
  Kenward-Roger 4-6%. `"Satterthwaite"` is also available (the
  covariance of the variance parameters from the Hessian of the REML
  log-likelihood, as SAS; its F tests are liberal when the covariance
  has several parameters), and `"containment"` gives nlme's df (`lme`)
  or the residual df (`gls`), as predictmeans 1.x. All are checked
  against SAS PROC MIXED. ML fits of `lme` and `gls` use Satterthwaite
  from the ML log-likelihood, as ML `lmer` fits (lmerTest); before 2.0
  they used nlme's or the residual df (`df_method = "containment"`). In
  simulations of ML `gls` fits with an unstructured covariance the
  residual df rejected 6.9% (t) and 9.3-11.7% (F) of true null
  hypotheses, Satterthwaite 5.5% and 6.3-8.8%; fit by REML for
  Kenward-Roger, which held 5%. For GLMMs (`glmer`, `glmmTMB`) `"auto"`
  is `"inner-outer"`: the df of the level of the random effects at which
  each effect varies, nlme's rule (Pinheiro and Bates 2000), computed
  from the design and extended to crossed random effects; a contrast
  gets the smallest df of the coefficients it uses, and the means of a
  term the df of the term. In simulations it rejected 6-8% of true
  hypotheses for effects that vary between clusters, where the residual
  df (`"residual"`, predictmeans 1.x) rejected 8-15%;
  [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  is the check with few clusters. `"containment"` and `"between-within"`
  are also available for `glmer`, and `"inner-outer"` for `lmer` and
  `lme`. `glm` models with a known dispersion (Poisson, binomial,
  negative binomial) get z tests (`"asymptotic (z)"`), as emmeans and
  SAS GLIMMIX; `"residual"` gives the t tests with the residual df of
  predictmeans 1.x, which were conservative in simulations (0.3-4.5% of
  true null hypotheses rejected at the 5% level, z tests 3.6-5.7%).
  Other choices are ignored for other model classes.

- vcov:

  Optional covariance of the coefficients to use instead of the model's.
  For `lm`, `aov` and `glm` models: `"HC0"`, `"HC1"`, `"HC2"` or `"HC3"`
  (heteroscedasticity-consistent, as
  [`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)),
  or a one-sided formula naming a cluster variable of the data, e.g.
  `~ Block` (cluster-robust CR1, as
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)).
  For any model: a covariance matrix of the coefficients, or a function
  that returns one from the model (for mixed models e.g. from the
  clubSandwich package). The degrees of freedom stay those of the model
  (the residual df; for `lmer` Kenward-Roger or Satterthwaite df
  computed from the model), and `vcov_adjust` is not used. With a
  `permlist`, `vcov` is not used (for robust permutation tests see
  `permmodels(robust = TRUE)`). The default `NULL` uses the model's
  covariance.

- equivalence:

  Optional equivalence margins for equivalence tests (two one-sided t
  tests, TOST) of the differences: a positive number \\d\\ (equivalence
  region \\-d\\ to \\d\\) or `c(lower, upper)`, on the scale of the
  comparisons (the link scale for GLMs: e.g. log ratios). The TOST
  p-value is the larger of the two one-sided p-values; below `level`,
  the difference lies within the margins (as emmeans'
  `test(..., side = "equivalence")`). Not adjusted for multiplicity;
  with a `permlist` it is t-based (normal-based in `contrastmeans`).

- bias_adjust:

  Logical: adjust the back-transformed means (non-identity link or
  `trans`) for the random effects (and, for a transformed response of a
  linear model, the residual error): \\h(\eta) + \sigma^2 h''(\eta)/2\\
  instead of \\h(\eta)\\, the second-order approximation of the mean
  over the population (as emmeans' `bias.adjust`); the response-scale
  SEs follow by the delta method. Without it (the default) the
  back-transformed mean is that of a unit with all random effects zero.

- sigma:

  With `bias_adjust = TRUE`: the SD to adjust for. The default is the
  square root of the summed random-intercept variances (`glmer`,
  `glmmTMB`, `lme`), plus the residual variance for a transformed
  response of a linear (mixed) model; models with random slopes need it
  given.

- average:

  How means average over the other factors when they are
  back-transformed: `"link"` (the default; average on the link scale,
  then back-transform) or `"response"` (back-transform each cell of the
  reference grid, then average with the same weights, as emmeans'
  `regrid()`: the mean of the response, e.g. a mean count). Both change
  only the `response` columns of `means`; the printed tables and the
  `pairs` stay on the link scale.

## Value

- Predicted Means:

  A table of predicted means.

- Standard Error of Means:

  A table of standard errors of predicted means.

- Standard Error of Differences:

  Standard errors of differences between predicted means.

- LSD:

  Least significant differences between predicted means.

- Pairwise p-value:

  A matrix with t-values above the diagonal and p-values below the
  diagonal, or matrix of pairwise comparison p-values for each level of
  `atvar`.

- mean_table:

  A summary of predicted means result including 'Predicted means',
  'Standard error', 'Df' and 'CIs'. When `trans!=NULL` or `trans!=I`, a
  table of back transformed means with CIs are also shown.

- predictmeansPlot:

  ggplot of predicted means.

- predictmeansBKPlot:

  ggplot of back transformed means.

- predictmeansBarPlot:

  gg bar plot of predicted means.

- p_valueMatrix:

  p_value matrix for pairwise comparison.

- means:

  A data frame with one row per predicted mean and fixed column names:
  the factor (and covariate) columns, `estimate`, `std.error`, `df`,
  `conf.low`, `conf.high`; when the means are back-transformed also
  `response`, `response.std.error` (delta method, as emmeans with
  `type = "response"`), `response.conf.low` and `response.conf.high`;
  `group` (letters, when pairwise comparisons are made), `n` (number of
  observations) and `estimable`. Values are not rounded. Also returned
  by `tidy()` and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- pairs:

  A data frame with one row per pairwise comparison: `level1`, `level2`,
  `by` (the level of `atvar`, when given; `NA` for pairs across levels),
  `estimate` (mean of `level1` minus mean of `level2`), `std.error`,
  `df`, `statistic` (t), `p.value`, `p.adj` (adjusted by `adj`, within
  each level of `atvar` when given), `p.equiv` (with `equivalence`: the
  TOST p-value), `conf.low` and `conf.high` (limits of each difference
  at `level`, not adjusted), `lsd` (the half-width of those limits),
  `adj.conf.low` and `adj.conf.high` (simultaneous limits within each
  family, which exclude 0 exactly when `p.adj` is below `level`, as SAS
  `AdjLower` and `AdjUpper` and emmeans; for `adj` "tukey" (each pair's
  own df, as SAS `ADJDFE=ROW`), "bonferroni", "sidak", "scheffe" and
  "mvt"; the step-down and false discovery rate methods have none) and
  `p.method`; with `permlist` the limits are permutation-t limits
  (estimate +- c SE, c the critical value of the permutation test), and
  the adjusted limits use the permutation distribution of the family's
  largest \|t\| ("tukey", "mvt", "westfall") or each pair's at level/m
  ("bonferroni", "sidak"); for a log link or `trans = exp` also `ratio`
  (of the back-transformed means, `exp(estimate)`), `ratio.conf.low` and
  `ratio.conf.high`, and for a logit link `odds.ratio` with its limits.

- settings:

  A list with the options used: model term, significance level,
  adjustment method, `atvar`, degrees-of-freedom method, `vcov_adjust`,
  `vcov` (`"model"`, `"HC3"`, ...), `singular` (a singular lmer/glmer
  fit; a message says so too), number of permutations, covariate values,
  model class and package version.

The result has class `c("predictmeans", "pdmlist")`. Printing shows the
same output as earlier versions; the `means`, `pairs` and `settings`
slots are not printed.

## Note

The `predictmeans` function becomes confused if a factor or covariate is
changed to the other in a model formula. Consequently, formulae that
include calls `as.factor`, `factor`, or `numeric` (e.g.
`as.factor(income)`) will cause errors. Instead, create the modified
variables outside of the model formula (e.g.,
`fincome <- as.factor(income)`) and then use them in the model formula.

Factors cannot have colons in level names (e.g., `"level:A"`) for the
legacy tables; the `predictmeans` function will confuse the colons with
interactions; rename levels to avoid colons. The `means` and `pairs`
tables are not affected.

For `predictmeans` function, it is assumed that methods `coef`, `vcov`,
`model.matrix`, `model.frame` and `terms` are available for `model`.

## References

Torsten Hothorn, Frank Bretz and Peter Westfall (2008), *Simultaneous
Inference in General Parametric Models. Biometrical*, Journal 50(3),
346–363.

Gilmour, A. R., Thompson, R. and Cullis, B. R. (1995), *Average
information REML: an efficient algorithm for variance parameter
estimation in linear mixed models*, Biometrics, 51(4), 1440-1450.

Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
*Prediction in linear mixed models*, Australian and New Zealand Journal
of Statistics, 46(3), 325-347.

Vanessa, C. (05 October 2022), *Confidence tricks: the 83.4% confidence
interval for comparing means*,
https://vsni.co.uk/blogs/confidence_trick.

Maghsoodloo, S. and Huang, C.-Y. (2010), *Comparing the overlapping of
two independent confidence intervals with a single confidence interval
for two normal population parameters*, Journal of Statistical Planning
and Inference, 140(11), 3295-3305.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data=Oats)
  ## The results as tables (columns described under Value)
  pm <- predictmeans(fm, "nitro", plot=FALSE)
  pm$means
#>   nitro  estimate std.error       df  conf.low conf.high  n estimable
#> 1     0  79.38889  7.174754 6.791824  62.31730  96.46047 18      TRUE
#> 2   0.2  98.88889  7.174754 6.791824  81.81730 115.96047 18      TRUE
#> 3   0.4 114.22222  7.174754 6.791824  97.15064 131.29381 18      TRUE
#> 4   0.6 123.38889  7.174754 6.791824 106.31730 140.46047 18      TRUE
  pm$pairs
#>   level1 level2   estimate std.error       df statistic      p.value
#> 1      0    0.2 -19.500000  4.435752 45.00014 -4.396098 6.656699e-05
#> 2      0    0.4 -34.833333  4.435752 45.00014 -7.852859 5.648849e-10
#> 3    0.2    0.4 -15.333333  4.435752 45.00014 -3.456761 1.205045e-03
#> 4      0    0.6 -44.000000  4.435752 45.00014 -9.919400 6.696622e-13
#> 5    0.2    0.6 -24.500000  4.435752 45.00014 -5.523302 1.583303e-06
#> 6    0.4    0.6  -9.166667  4.435752 45.00014 -2.066542 4.456070e-02
#>          p.adj  conf.low   conf.high      lsd p.method
#> 1 6.656699e-05 -28.43406 -10.5659375 8.934062        t
#> 2 5.648849e-10 -43.76740 -25.8992708 8.934062        t
#> 3 1.205045e-03 -24.26740  -6.3992708 8.934062        t
#> 4 6.696622e-13 -52.93406 -35.0659375 8.934062        t
#> 5 1.583303e-06 -33.43406 -15.5659375 8.934062        t
#> 6 4.456070e-02 -18.10073  -0.2326042 8.934062        t
  pm$settings
#> $modelterm
#> [1] "nitro"
#> 
#> $level
#> [1] 0.05
#> 
#> $adjust
#> [1] "none"
#> 
#> $atvar
#> NULL
#> 
#> $df_method
#> [1] "Kenward-Roger"
#> 
#> $vcov_adjust
#> [1] TRUE
#> 
#> $vcov
#> [1] "model"
#> 
#> $equivalence
#> NULL
#> 
#> $bias_adjust
#> [1] NA
#> 
#> $singular
#> [1] FALSE
#> 
#> $average
#> [1] "link"
#> 
#> $nperm
#> [1] NA
#> 
#> $covariates
#> NULL
#> 
#> $model_class
#> [1] "lmerModLmerTest"
#> 
#> $version
#> [1] "1.1.1.9000"
#> 
  summary(pm)
#> Predicted means for 'nitro' (lmerModLmerTest model; df: Kenward-Roger, Kenward-Roger adjusted SEs)
#> Confidence limits: 95%
#> 
#>  nitro estimate std.error    df conf.low conf.high
#>      0    79.39     7.175 6.792    62.32     96.46
#>    0.2    98.89     7.175 6.792    81.82    115.96
#>    0.4   114.22     7.175 6.792    97.15    131.29
#>    0.6   123.39     7.175 6.792   106.32    140.46
#> 
#> 6 pairwise comparisons.
#> SED: min 4.436, mean 4.436, max 4.436
#> LSD (5%): min 8.934, mean 8.934, max 8.934
#> p-values: t tests, adjusted by 'none'.
  ## Means within each variety, Tukey-adjusted pairs, plots
  predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="tukey",
    which.plots=c("mean", "bar", "pval", "back"), line=FALSE)




  names(predictout)
#>  [1] "Predicted Means"                                         
#>  [2] "Standard Error of Means"                                 
#>  [3] "Standard Error of Differences"                           
#>  [4] "LSD"                                                     
#>  [5] "Pairwise comparison p-value (adjusted by 'tukey' method)"
#>  [6] "Golden Rain"                                             
#>  [7] "Marvellous"                                              
#>  [8] "Victory"                                                 
#>  [9] "mean_table"                                              
#> [10] "predictmeansPlot"                                        
#> [11] "predictmeansBarPlot"                                     
#> [12] "p_valueMatrix"                                           
#> [13] "means"                                                   
#> [14] "pairs"                                                   
#> [15] "settings"                                                
  print(predictout$predictmeansPlot)
#> $meanPlot

#> 
#> $ciPlot

#> 
  print(predictout$predictmeansBarPlot)

  ## lme fits work the same way
  fl <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  predictmeans(fl, "nitro", adj="BH", plot=FALSE)$pairs
#>   level1 level2   estimate std.error       df statistic      p.value
#> 1      0    0.2 -19.500000  4.435758 44.99995 -4.396092 6.656859e-05
#> 2      0    0.4 -34.833333  4.435758 44.99995 -7.852848 5.649197e-10
#> 3    0.2    0.4 -15.333333  4.435758 44.99995 -3.456756 1.205064e-03
#> 4      0    0.6 -44.000000  4.435758 44.99995 -9.919388 6.697195e-13
#> 5    0.2    0.6 -24.500000  4.435758 44.99995 -5.523295 1.583360e-06
#> 6    0.4    0.6  -9.166667  4.435758 44.99995 -2.066539 4.456099e-02
#>          p.adj  conf.low   conf.high      lsd p.method
#> 1 9.985289e-05 -28.43407 -10.5659250 8.934075        t
#> 2 1.694759e-09 -43.76741 -25.8992584 8.934075        t
#> 3 1.446077e-03 -24.26741  -6.3992584 8.934075        t
#> 4 4.018317e-12 -52.93407 -35.0659250 8.934075        t
#> 5 3.166719e-06 -33.43407 -15.5659250 8.934075        t
#> 6 4.456099e-02 -18.10074  -0.2325917 8.934075        t
```
