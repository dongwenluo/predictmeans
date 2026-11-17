# predictmeans 1.1.1.9000 (development version, to be released as 2.0.0)

## In short

predictmeans 2.0 keeps the 1.x functions and printed output, and changes
defaults where simulations or the literature showed a more reliable method.
Plots are drawn only when asked for (`plot = TRUE` or `autoplot()`), and
lmerTest is no longer needed (`anova()` of `lmer` fits has F tests with
Kenward-Roger df from the package's own code).
Results checked against SAS PROC MIXED, Genstat and emmeans are in the tests.

* **Degrees of freedom** are chosen by the model in one place (`df_method =
  "auto"`): Kenward-Roger with adjusted SEs for REML `lmer`, `lme` and `gls`
  (observed information, as SAS `DDFM=KR`); Satterthwaite for ML fits; the df
  of the grouping level for GLMMs; z tests for Poisson, binomial and negative
  binomial `glm`. Options: `"Satterthwaite"`, `"containment"`,
  `"between-within"`, `"KR-expected"` (pbkrtest, emmeans), `"KR-average"`
  (Genstat), `"residual"` (1.x). `vcov_adjust = FALSE` gives 1.x SEs.
* **Stable result tables**: `means`, `pairs` and `settings` with fixed column
  names (`tidy()`, `as.data.frame()`, `autoplot()`), including adjusted
  confidence limits that agree with the adjusted p-values.
* **Empty cells**: means keep the 1.x average of the observed cells, with a
  warning naming them; their comparisons use the cells both means have.
* **`aov()` with `Error()`** is analysed by strata, as Genstat (missing
  values estimated as Genstat's ANOVA does), not refitted as `lmer`.
* **Type III tests** do not depend on the model's contrasts (as SAS Type 3).
* **Permutation tests**: Freedman-Lane by default, also for mixed models,
  GLMs and GLMMs; permutation confidence limits, Tukey and Westfall-Young
  adjustments, robust (HC3) statistics.
* **Model checking**: `pm_check()` and a ggplot2 `residplot()` for every
  model class (`vignette("checking")`).
* More adjustments (`"sidak"`, `"scheffe"`, `"mvt"`), contrast families,
  joint tests, slopes, equivalence tests, robust covariances (`vcov =`),
  response-scale means with bias adjustment, and an adapter for new model
  classes (`?pm_coef`).

## Changes that affect results or scripts

* **`semipred()` bands of gaussian `semireg()` / `semireg_tmb()` fits are t
  limits with Satterthwaite df** at each point (new column `df`), as
  `predictmeans()` on the fit; simultaneous bands use a multivariate t with
  each point's df; differences of curves (`contr`) their own df.
  Generalized fits keep normal limits. In a split-plot design with curves
  over time (`dev/sim-semipred-df.R`, 500 data sets per setting) the normal
  band of a treatment's curve covered 0.91 (pointwise) and 0.92 (whole
  curve) at the 95% level; with t (8-11 df) 0.94-0.95 and 0.96-0.97.
  Differences of curves: normal 0.94-0.96, t 0.96-0.97. Curves without
  design strata: 0.95-0.98 either way (the df are then the residual df).
* **Generalized `semireg()` / `semireg_tmb()` fits get the df of the
  design's strata** in `predictmeans()` and `semipred()`, the inner-outer
  df that GLMMs get (the smooth blocks are not strata), instead of z tests;
  without random terms besides the smooths they stay z (as `glm` with a
  known dispersion). Poisson counts in a split-plot design with curves over
  time (`dev/sim-semireg-glmm-df.R`, 1000 data sets per setting): z tests of
  two treatments at a time rejected 5.8-9.2% of true equalities and normal
  bands covered a treatment's curve 0.91-0.93; with t on the 10 whole-plot
  df 3.1-5.9% and 0.955-0.970. `settings$df_method` says
  `"inner-outer (design strata; smooths within)"`. The rule for a row (the
  smallest df of the coefficients it uses; means: the df of their term) is
  now one function for every model (`pm_io_row_df()`).
* **The Satterthwaite df of `semireg()` fits hold the smoothing parameters
  at their estimates** (`predictmeans()` and `semipred()`): they count the
  design's variance components and the residual. With the smoothing
  variances free, a point where one of them dominates got 2 df (a
  WarsawApts fit: 1.85-369 df along the curves, now 401, the residual df);
  their uncertainty is what `unconditional = TRUE` adds. In the split-plot
  designs the df hardly change (the level of `predictmeans()`'s tests stays
  4.6%). The glmmTMB engine holds the same quantity (the ratio of each
  smoothing variance to the residual variance), so both engines agree.
* **Kenward-Roger df with a variance component estimated as zero** (a
  singular fit): the component is left out of the adjustment, as in SAS
  PROC MIXED. Its near-singular information gave nonsense df before (in a
  `semireg()` fit whose deviation curves had variance 0: 0.54 and 0.04 df
  for the factor terms; now 33 and 114, those of the model without the
  component). Fits without a zero variance are unchanged.
* **`semireg()` uses the optimizer `nloptwrap`** (lme4's default for `lmer`)
  instead of `bobyqa`: the same optimum (log-likelihoods equal to 4e-6 on
  the golden examples and on binomial and Poisson fits) 1.3-4 times faster
  (growth curves of 2257 rows: 53 s -> 18 s; `lme` takes 31 s). Estimates
  move by optimizer tolerance only. `control = lmerControl(optimizer =
  "bobyqa")` gives the 1.x optimizer.
* **`semipred()` and `predictmeans()` on generalized `semireg()` /
  `semireg_tmb()` fits include the uncertainty of the smoothing parameters**
  (`unconditional`, as mgcv's): the first-order correction J V_theta J' is
  added to the covariance of the coefficients (lme4: J from the deviance
  function, V_theta from its Hessian; glmmTMB: TMB's joint precision of all
  parameters). In a simulation of the difference between two curves
  (`dev/sim-curve-difference.R`) the bands of Poisson and binary models
  covered 0.93 pointwise and 0.89-0.91 over the whole curve without it,
  0.955-0.957 and 0.94-0.95 with it. Gaussian bands already covered 0.97, so
  they stay conditional by default. `unconditional = TRUE`/`FALSE` (or
  `options(predictmeans.semi_unconditional =)`) chooses. The corrected SEs are
  within 1% of mgcv's `unconditional = TRUE` (which adds a second-order term).
  The prediction design of the `smoothZ` form now takes coefficients, columns
  and covariance rows in one order (that of the fit's blocks).
* **`semireg()` gaussian fits: `anova()` df and `varcomp()` SEs were
  wrong.** lmerTest's Satterthwaite df and `varcomp()` used the deviance
  function of the placeholder model (plain random intercepts where the spline
  design goes), not of the fitted one. Now they equal the same model fitted
  by `nlme::lme()` (Warsaw apartments: 150 df instead of 1871 for the linear
  term; SE of the smoothing variance 1.45 instead of 7.7, with a sensible
  upper limit instead of 5e45). Fitted values, coefficients and `semipred()`
  results are unchanged.
* **`semireg_tmb()` with families whose `sigma()` is not a residual SD
  (negative binomial: theta; beta: phi): `semipred()` standard errors were
  scaled by it.** The covariance now uses the information weights of the
  working model (`pm_working_tmb()`), so negative binomial SEs from
  `semireg_tmb()` agree with those of `semireg(ngbinomial = TRUE)`; with
  nearly Poisson data (large theta) they were up to 10^6 times too large.
  Poisson and binomial fits are unchanged.
* `semireg()` and `semireg_tmb()` share one fitting engine
  (`R/semireg-engine.R`):
  - The model structure is built without fitting (`lme4::lFormula()`,
    `glFormula()`, `glmmTMB(doFit = FALSE)`) and fitted once; before, a full
    preliminary fit was made first. 13 example fits: 19 s -> 14 s.
  - The data are no longer written into the caller's workspace, and `data`
    can be an expression (`data = subset(d, x > 0)`).
  - `cbind(successes, failures)` responses work (the second column was
    centred as a covariate: "Response is constant").
  - A smooth whose variance is estimated as zero (a straight line) no longer
    stops the fit ("Matrix seems negative semi-definite"): the covariance of
    the coefficients is computed without inverting Lambda (the same matrix
    otherwise), and those coefficients get zero covariance.
  - Weighted gaussian fits use their prior weights in that covariance.
* `covariatemeans(contr = c(1, 3))` gives the difference between two levels
  of the term along the covariate with its confidence band (as `semipred()`
  does): the two lines differ where the band leaves out zero, which the
  overlap of the lines' own bands does not show. The band comes from the
  contrast of the two lines.
  - The test is on the model scale, and that is what is plotted
    (`contr_scale = "model"`, the default): the difference, shown as the
    ratio for a log link or log response and as the odds ratio for a logit
    link (exact; a log axis with a line at 1).
  - `contr_scale = "response"` plots the difference of the back-transformed
    lines with a delta-method band (gradient h'(eta_a) K_a - h'(eta_b) K_b):
    the size of the difference in the units of the response, not a second
    test; the values where the lines differ on the model scale are marked.
    Both scales equal emmeans for GLMs and GLMMs (`pairs()` and `regrid()`).
    In a simulation with a Poisson GLMM the delta band covered the true
    difference in 94% of data sets; the difference of the two lines' limits
    in 100% (far too wide).
  - `band = "simultaneous"` gives Scheffe limits that hold for all covariate
    values together (rank 2 for straight lines: the Working-Hotelling band),
    for saying over which range the lines differ; `"pointwise"` (default) is
    for one value chosen beforehand.
  - The result has `test`, the joint F test that the two lines are the same,
    and `significant`, the ranges of the covariate where the band leaves out
    zero; both are printed. The table has an indicator column for each set
    of limits: `includes0`, `ratio.includes1` or `odds.ratio.includes1`,
    and `response.includes0`. `?covariatemeans`, section "Comparing two
    lines", and `vignette("response-scale")` give the reasoning.
* `covariatemeans(bias_adjust = TRUE)` (also `"taylor"`, and `sigma`), as in
  `predictmeans()` and with the same code: the back-transformed lines are
  the mean of the response, not its median (a transformed response) or the
  typical cluster (a GLMM). Exact: the factor exp(sigma^2/2) for a log, +
  sigma^2 for a square root. The default is unchanged. With `contr` the
  response-scale difference is that of the adjusted lines; ratios and tests
  do not change.
* `semipred(contr = )` uses the same calculation as
  `covariatemeans(contr = )` (one function, `pm_diff_table()`).
  - Linear models: the same numbers as before; `pred_df` also has
    `includes0` and `significant`.
  - Generalized models: **the band was wrong** (the standard error of the
    link-scale difference times a weight evaluated at the difference of the
    means). Now the link-scale difference with its limits, the ratio (log
    link) or odds ratio (logit link), and the response-scale difference by
    the delta method, each with its indicator column. The plot is the ratio
    or odds ratio; `contr_scale = "response"` plots the response-scale
    difference.
  - The ranges of the covariate where the curves differ are printed and
    kept in `attr(pred_df, "significant")`; with `simultaneous = TRUE` they
    hold for the whole curve.
* `bootmodels()` takes `lm`, `aov`, `glm` (gaussian, binomial, Poisson,
  Gamma, `MASS::glm.nb`), `gls`, `lme` and `glmmTMB` models as well as
  `lmer` and `glmer`. Responses are simulated from the model without its
  fixed effects (intercept, offsets, the same random effects and
  covariance structure); `lme` and `gls` responses from the null model's
  covariance of the data (tested: the simulated covariance is the fitted
  one). For `lm` the bootstrap p-values approach the exact F tests (tested
  with 999 samples). `glmmTMB` has no sequential tests: Type II is used, with
  a message. Quasi families, `gam` and `aov` with `Error()` are refused with
  a message that says what to use instead.
* `CookD()` takes `lm`, `aov`, `glm`, `glmer`, `glmmTMB` and `gam` models as
  well as `lme`, `lmer` and `gls`. `lm` and `aov` use the exact closed form
  (`stats::cooks.distance()`); other models are refitted for each deletion
  (exact also for `glm`, where `cooks.distance()` is a one-step
  approximation; for a `gam` p is the effective df, as mgcv). A deletion that
  leaves a coefficient without data gives NA. `aov` with `Error()` is refused
  with a message.
* New test `test-model-classes.R`: every main function on 13 model classes
  (lm, aov, aov with `Error()`, Poisson and binomial glm, gls, lme, lmer,
  Poisson and binomial glmer, Poisson and negative binomial glmmTMB, gam);
  each works or refuses the class with a message. It found that
  `pm_slopes()` and `CookD()` could not find the data of a `gam` fitted
  inside a function (mgcv records the global environment): fixed for every
  function that refits or rereads the data.
* `semipred()` means of generalized models: the limits are the
  back-transformed limits of the link scale, h(eta -+ c SE), as in
  `predictmeans()` and `covariatemeans()` (always inside the range of the
  response, not symmetric), and `SE_mu` is the delta-method SE |h'(eta)| SE.
  Before, the limits were mean -+ c SE_mu, clipped to [0, 1] for binomial,
  with a weight that is |h'(eta)| only for canonical links. `semipred()` on
  a `semireg()` fit now gives the same means, SEs and limits as
  `predictmeans()` on it.
* `semipred(band = "simultaneous")`, the argument of `covariatemeans()`;
  `simultaneous = TRUE` still works. The difference of two surfaces (two
  covariates, `contr`) uses the same calculation as one covariate.
* `R2_model(method = "nakagawa")` for `lme` models (equal to `lmer`; with a
  variance function the observation-level variance is the mean residual
  variance).
* **`semireg()` with a 0/1 response: standard errors were 15-25% too
  small.** The covariance of the fixed and smooth coefficients used the
  constant 3/pi^2 (logit; the latent variance that `R2_glmm()` uses for
  variance shares) instead of the working weights of the fit. It now uses
  the working weights for every family, and its fixed-effect block equals
  `vcov(model, use.hessian = FALSE)` for gaussian, Poisson and binary fits
  (tested). `semipred()` bands and `predictmeans()` on such fits are wider.
  A simulation (`dev/sim-curve-difference.R`) gives the coverage of the
  band of a difference of two curves: 0.96 pointwise and 0.96 simultaneous
  for a linear model; 0.92-0.94 pointwise and 0.87-0.88 simultaneous for
  Poisson and binary models (the covariance is conditional on the
  estimated smoothing parameters).
* `covariatemeans()` takes `aov` models with `Error()` that have a covariate
  (e.g. a baseline covariate per subject in repeated measures): analysed by
  strata, equal to the mixed model in balanced data (it failed before).
* **Fewer arguments in `predictmeans()`: `plot_args`.** The eleven plot
  settings (`plottitle`, `plotxlab`, `plotylab`, `plotord`, `lineplot`,
  `lsd_bar`, `responsen`, `jitterv`, `basesz`, `prtplt`, `newwd`) are one
  argument, `plot_args = list(plottitle = "Yield", basesz = 14)`. Scripts
  that give them by name as before (`plottitle = "Yield"`) keep working, also
  abbreviated. `mplot`, `barplot`, `pplot` and `bkplot` (deprecated for
  `which.plots`) are no longer arguments but still work, with the same
  warning. An argument name that does not exist is now an error. Only calls
  that passed arguments after `Df` by position need a change: `lsd_bar` was
  the eighth argument, where `level` is now.
* Inside `predictmeans()` the stages are ordinary functions with named
  results (`pm_setup()`, `pm_estimate()`, `pm_compare()`, `pm_tables()`)
  instead of code evaluated in the call's environment; the snapshot of 52
  calls is identical.

* **Back-transformation.**
  * *A transformed response in the formula is back-transformed by default*,
    as a link is: `log(y)`, `log(y + c)`, `log10(y)`, `log2(y)`, `log1p(y)`,
    `sqrt(y)`, `y^p`, `1/y`, `qlogis(y)`, `asin(sqrt(y))` and the new
    `pm_boxcox(y, lambda)` and `pm_yeojohnson(y, lambda)` (Box-Cox and
    Yeo-Johnson, to use in a formula). `predictmeans()` and
    `covariatemeans()` say so in a message; `trans = I` keeps the model scale
    only, and a `trans` that is given is used as before. For `log(y)` the
    pairs also have the ratios.
  * *What the back-transformed mean is* is now stated
    (`settings$back_transform`, `summary()`): for a GLM the mean; for a GLMM
    the response of a typical cluster (random effects at zero); for a
    transformed response the median of the original scale (the geometric
    mean for a log).
  * *`bias_adjust = TRUE` gives the mean exactly*: E h(eta + u) by
    Gauss-Hermite quadrature, which is exp(eta + sigma^2/2) for a log and
    eta^2 + sigma^2 for a square root. Before, it was the second-order
    formula h + sigma^2 h''/2 (as emmeans), which is too low for a log when
    sigma is large and was off by up to 0.018 for a logit with sigma = 1.2
    (0.389 against the true 0.371). `bias_adjust = "taylor"` keeps that
    formula.
  * `covariatemeans()` with a log link or `trans = exp`: the line is the
    back-transformed fitted value. 1.x drew exp(fitted + SE/2), adding half
    the standard error of the fitted value, which is not a bias adjustment
    and made the line depend on its own precision.
  * `pm_check()` reports the Box-Cox power the residuals of an `lm` fit
    favour, with its 95% interval, and names the usual transformation within
    it; it points to a generalized linear model for counts and proportions.

* **Comparison intervals: overlap that answers the tests.** Confidence
  limits of means are not for comparing means: with random blocks (oats) two
  means with p < 0.001 have overlapping 95% limits, and the limits that just
  touch at p = 0.05 are 49% limits, not the 83.4% of the rule for
  independent means. Every result now has comparison intervals, mean +- c SE
  with c = crit x SED/(SE1 + SE2) averaged over the pairs compared (all
  pairs, or those within a level of `atvar`; crit is t, or the critical
  value of `adj` when it has one): two means differ at `level` when their
  intervals do not overlap. This is exact when the SEDs are equal (balanced
  designs); otherwise the number of pairs judged differently by their own
  tests is reported. For independent means with equal SEs c is t/sqrt(2)
  (83.4% for a 5% test), with SE ratio K it is t sqrt(1 + K^2)/(1 + K), as
  in Maghsoodloo and Huang (2010). They are in `means` (`comparison.low`,
  `comparison.high`), `settings$comparison`, `summary()` and
  `autoplot(interval = "comparison")`.
* **`letterCI = TRUE` reads its letters from the comparison intervals**,
  which the mean table shows as LL and UL with their level, and `level` is
  the significance level of the comparisons (default 0.05), as everywhere
  else. Before, the intervals were 83.4% confidence intervals whatever the
  design (`level` set their level, 0.166 by default), which finds too few
  differences when the means are correlated. A script with `letterCI = TRUE,
  level = 0.166` now compares at the 16.6% level: drop `level`.
* `autoplot(interval = "lsd")` and the LSD line of `summary()` use the LSDs
  that go with the adjusted p-values.
* **Long plot titles are broken into lines** in every plot of the package
  (the mean, bar and back-transformed plots of `predictmeans()`, `PMplot()`,
  `autoplot()`, `covariatemeans()`, `semipred()`, `residplot()`), also a
  title given by the user (`plottitle`, `mtitle`); line breaks in it are
  kept. The line length follows the font size;
  `options(predictmeans.title_width = )` sets a number of characters.

* **The printed LSDs go with the adjusted p-values.** `predictmeans()` printed
  the unadjusted LSD (t x SED) with every adjustment except Bonferroni, so
  with `adj = "tukey"` a difference could be larger than the LSD and have an
  adjusted p-value above 0.05 (oats, nitrogen: LSD 10.65, but 14.08 is the
  value that matches the Tukey p-values). Now the LSD is the SED times the
  critical value of the method: t (`"none"`), the studentized range over
  sqrt(2) (`"tukey"`, the HSD), t at the adjusted level (`"bonferroni"`,
  `"sidak"`), and the Scheffe and multivariate t values (`"scheffe"`,
  `"mvt"`), so a difference is beyond its LSD exactly when its adjusted
  p-value is below the level. The family is that of the p-values: with
  `atvar`, the pairs within a level (Bonferroni used all pairs: oats,
  `atvar = "Variety"`: 32.8, now 25.2). The step-wise methods (`"holm"`,
  `"hochberg"`, `"hommel"`, `"BH"`, `"BY"`) have no critical value for all
  pairs: their LSDs stay unadjusted, and a note printed with the LSDs says
  so. The LSD bar of the mean plot follows. `adj = "none"` is unchanged.

* **`Kmatrix()` reviewed** (the 52-call snapshot of `predictmeans()` is
  unchanged; 549 of 592 `Kmatrix()` cases are identical, the others are the
  changes below).
  * *Covariate values by name.* `covariate = list(x = c(-1, 0, 1), z = 0.3)`:
    several values give rows, one value fixes the covariate, the others stay
    at their means; a data frame gives joint points. Before, the names of
    `covariate = c(z = 0.2, x = 0.5)` and of a `covariateV` list were ignored
    (values went by position, without a warning) and one value was recycled
    to every covariate; unnamed values must now be one per covariate.
  * *Any combination of variables.* `modelterm` need not be a term of the
    model: `"A:B"` for `y ~ A + B`, and a covariate may be in it (`"A:x"`).
    `predictmeans()` accepts combinations of factors that are not a term
    (with a note: such means are sums of the factors' effects), shows a
    combination without data when its mean can be estimated, and uses the
    smallest df of the terms of its factors where df go by term. The
    factors of a term in another order (`"Variety:nitro"`) now give the
    means in that order instead of an error.
  * *Faster for factors that do not interact.* Each group of variables that
    share a model term is averaged on its own grid: `y ~ Line + Block + Site
    + Year` with 300, 90, 4 and 3 levels took 3.6 s (324,000 grid rows) and
    takes 0.14 s; `cov_reduce = "average"` with 200,000 observations 20 s ->
    1.1 s. The model matrix of the data is no longer built to read the
    coefficient names. `refgrid = TRUE` returns the grid of all combinations
    (`$refgrid` is otherwise `NULL` unless it was needed).
  * *Empty cells are those of each interaction term.* Before, the cells were
    the combinations of all factors that interact with anything, so in
    `y ~ date/farm + date*trt` a date-farm-treatment combination without data
    made the mean of that date and treatment an average over its observed
    farms (with a warning), though farm and treatment do not interact and the
    mean over all farms can be estimated. Now it is that mean. The same for a
    factor that interacts only with a covariate. The exception "every cell
    has 0 or 1 observation" is gone: unreplicated interacting factors with an
    empty cell get the average of the observed cells and the warning, as
    replicated ones.
  * *Nested factors* declared in the formula (`G/T`) are also left out of
    the averages with `weights = "proportional"` and `"outer"` (gave NA).
    Nesting that is only in the data (`y ~ G + T`): see "Confounded
    factors" below.
  * *Confounded factors.* Two factors that do not interact in the model and
    do not occur in all combinations (a factor nested in another only in the
    data, `y ~ G + T`; a spray volume used with one droplet size; each
    pasture type on one farm at a date) have no mean over all their
    combinations. The means now average over the combinations that exist,
    with a warning that names the two factors and the ways to say it in the
    model: write the nesting (`G/T`) or leave out the factor that adds
    nothing (the same means, equal weights), or `weights = "cells"` (the
    combinations weighted by their numbers of observations). The pair is
    found by trying, pairs among the factors averaged over first. Without
    this the per-term cells above turned such means into NA (the
    maintainer's pasture trial: 9 of 12 date x forest means and 60 of 66
    pairs); 1.x gave the same averages without saying why. Means and pairs
    that still cannot be estimated (no common cell) are NA.
    `vignette("missing-cells")` has the examples and the table of the three
    ways.
  * `Kmatrix()` accepts `aov` models with `Error()`; an `aov` fit with
    aliased coefficients gets the estimability check; clearer errors for
    unknown variables and covariates. New adapter generic `pm_contrasts()`
    (a new model class needs `coef()`, `vcov()` and `terms()` only).

* **Several hundred means** (e.g. the lines of a breeding trial). With more
  than 100 means and no comparison asked for, `predictmeans()` gives the
  means, their SEs and the SED and LSD summaries (maximum, minimum, average;
  within the levels of each factor for a term with several factors) from the
  covariance of the means, without a row for every pair, and `pairs` is
  `NULL` (a message says so; `summary()` prints the SED and LSD ranges). 200
  lines in an `lmer` model took 47 s and now take 0.7 s with the same
  results; 300 lines (44,850 pairs) take about a second. `pairwise = TRUE`
  (or `adj`) still compares every pair (200 lines: 5 s). Letter groups are
  not made for more than 100 means: they cannot be read and the letter
  algorithm takes minutes. `options(predictmeans.max_means = )` changes the
  limit. Up to 100 means nothing changes.

* **`PMplot()` is drawn with ggplot2** and returns the plot: a square level
  plot of the lower triangle, like a correlation plot, with the legend in
  the empty corner. New: `values = TRUE` prints the p-values on the tiles,
  `basesz` sets the font size, a list of matrices (`atvar`) is drawn as
  blocks down the diagonal with lines between them, and a matrix with as few
  as two levels is plotted. `margin` is no longer used. Inside a function or
  a loop wrap the call in `print()`. `autoplot(type = "pvalue")` draws the
  same plot (before: tiles of unequal sides, one panel per `atvar` level) and
  gains `values`.

* Faster with many pairs: the SEDs no longer build the covariance matrix of
  all pairs (96 means, 4,560 pairs: 1.4 s -> 0.2 s), and Tukey-adjusted limits
  compute one quantile for each distinct df (60 means: 1.8 s -> 0.2 s).

* **`varcomp()` for every model class**, with one rule for SEs and
  intervals. New: `lm`, `aov` and `glm` (the error term: residual variance,
  dispersion, negative binomial theta), `aov` with `Error()` (analysis of
  variance estimates of the stratum variances, also with missing values),
  `gls` (correlation and variance-function parameters), `gam` (any mgcv
  version: `gam.vcomp()` returns a list since mgcv 1.9-3) and `semireg`
  fits; `glmer` with several grouping factors (failed) and with `nAGQ = 0`;
  `lme` with correlation and variance structures (failed, e.g. AR(1)); lmeInfo
  is no longer used. SEs are from the observed information (as SAS PROC
  MIXED) by the delta method; `information = "average"` gives Genstat's s.e.
  (equal to the printed digits in 17 of 17 values) and `"expected"` the
  expected information. Checked against real runs (2026-11-14): SAS's
  standard errors and limits for 10 designs (missing plots, random slope,
  AR(1), unstructured repeated measures with the covariances' SEs, a
  component at zero, crossed factors) and Genstat's s.e. for 5 (a negative
  component equals the `aov` estimate), all to the printed digits; SAS's
  limits for covariances and correlations are Wald on the natural scale,
  predictmeans' on the Fisher z scale. Intervals of variances are Satterthwaite's
  chi-square intervals with nu = 2 (estimate/SE)^2 (SAS's `COVTEST CL`; the
  exact interval for an `lm` residual variance); before, `lmer` had lme4's
  profile intervals, which for a REML fit are those of the ML fit, and `lme`
  nlme's. A variance at zero has no SE. One naming for every class
  (`group.term`, `group.term1.term2`, `cor.`, `var.`, `residual`); glmmTMB
  `ar1`, `cs` and `toep` terms give their variance and correlation instead of
  every element. Much faster for `lmer`: 1,200 rows 7.6 s -> 0.06 s; 6,000
  rows over 20 minutes -> 0.2 s.

* **`predictmeans()` no longer draws plots unless asked** (`plot = FALSE` is
  the default). Drawing took most of the time of a call (0.22 s against
  0.03 s for the oats example) and put graphics pages into scripts, reports,
  Shiny apps and loops. The result holds the same tables, and `autoplot()`
  draws from it. `plot = TRUE`, or `options(predictmeans.plot = TRUE)` for a
  session, draws the 1.x plots; so does asking for plots by name
  (`which.plots`, `plotord`, or the old `mplot`, `barplot`, `pplot`,
  `bkplot`).
* **lmerTest is no longer needed**: it moved from Depends to Suggests (used
  only in the tests, to compare against). predictmeans' own df and tests
  never used it. `anova()` of an `lmer` fit now gives F tests with df by the
  package's own code (see "New features"); `summary()` with t tests and
  `ranova()` need `library(lmerTest)`. Gaussian `semireg()` fits are class
  `semerMod` (an `lmerMod`), not lmerTest's: `anova(fit)` tests their fixed
  effects and smooth terms (replacing `lmerTest::ranova(fit$semer)`).
* **Adjusted confidence limits of the pairs.** With `adj` "tukey",
  "bonferroni", "sidak", "scheffe" or "mvt", `pairs` has the columns
  `adj.conf.low` and `adj.conf.high`: simultaneous limits within each family
  (each level of `atvar`) that exclude 0 exactly when `p.adj` is below
  `level`, as SAS `AdjLower`/`AdjUpper` and emmeans' `confint()` (Tukey with
  each pair's own df, as SAS `ADJDFE=ROW`; equal to SAS's Tukey-Kramer limits
  in `test-sas-adjust.R`). The step-down and false discovery rate methods
  (holm, hochberg, hommel, BH, BY) have no simultaneous limits and no such
  columns. With a `permlist`, "tukey", "mvt" and "westfall" use the
  permutation distribution of the family's largest |t|. `conf.low`,
  `conf.high` and `lsd` are now always the limits of each comparison on its
  own: with `adj = "bonferroni"` they were Bonferroni limits over all pairs
  (not per family), and with a `permlist` and a max-|t| adjustment they were
  the simultaneous limits now in `adj.conf.*`. The printed LSD is unchanged.
  `contrastmeans()` gives the same limits as `adj.lower.CL` and
  `adj.upper.CL` ("bonferroni", "sidak", "scheffe" with the rank of
  `ctrmatrix`, "mvt"; with a `permlist` also "westfall"), equal to emmeans'
  `confint()`; its `lower.CL` and `upper.CL` are per contrast, also with a
  `permlist` and a max-|t| adjustment.
* `permmodels()` for lmer models prints its tables (coefficients with
  Satterthwaite df and p-values; the type I ANOVA) from predictmeans' own
  Satterthwaite code instead of converting the fit to an lmerTest object.
  The numbers are lmerTest's: the type I hypotheses now use lmerTest's
  unit-diagonal rows, on which its Satterthwaite df of several rows depend
  (unbalanced data: 0.03-0.04 df closer to lmerTest's in `permstats()`).
* Review of the week's additions, with edge cases (model classes, odd
  designs):
  * `aov()` with `Error()` and missing values: rows left out by `subset =`
    were taken as missing values to estimate (13 instead of 1 in a test);
    they no longer are. The estimation no longer forms n x n matrices: 2,000
    plots took 11 s, now 0.2 s; results identical. A check that observed
    plots are in the fixed part of the model never ran; it does now.
  * The empty-cell warning named the means with covariate values when the
    grid had a covariate (`pm_slopes()`, `covariatemeans()`: "b2:0.1494...",
    many times over); it names the factor levels, once. `pm_slopes()` gives
    it once instead of twice.
  * `residplot()`: with fewer than 7 distinct fitted values (a model with
    factors only) the trend line joins the mean residual at each fitted
    value instead of a loess smooth, which warned. `pm_check()`: no warning
    from the equal-variance test when the deviations are equal within every
    level.
* **One computation of the pairs.** The `pairs` table, the printed p-value
  matrices and the letter groups now read one table of pair statistics, and
  `predictmeans()`, `contrastmeans()`, `pm_slopes()` and the means of
  `semireg` fits use one set of adjustment functions (`R/p-adjust.R`). The
  printed output is unchanged except where the two computations disagreed:
  * With a `permlist` and `level` other than 0.05, the letters among all
    means used 0.05 (the letters within `atvar` levels and their note used
    `level`). Both now use `level`, and `settings$level` is `level`.
  * With means that cannot be estimated (confounded factors), Tukey's k
    within an `atvar` level counted them in `pairs$p.adj` but not in the
    printed matrices; both now count the means that take part in the
    estimable comparisons. Scheffe's rank is the rank of those comparisons
    (as emmeans), and the adjusted limits use the same sizes as `p.adj`.
* Fixed: `predictmeans()` with `level` given (and `letterCI = FALSE`)
  stopped with "object 'slevel' not found" when making pairwise comparisons,
  as in 1.1.1.
* `permmodels(type = "III")` computes the map from the model's coefficients
  to the sum-to-zero hypotheses once instead of in every refit: `lm` with 999
  raw permutations 1.9 s (3.3 s before, 2.3 s before the contrast-free Type
  III tests); results identical.

* **`aov()` with `Error()` and missing values: estimated as Genstat's ANOVA.**
  In an orthogonal design with some plots missing (rows with an `NA`
  response, or a few combinations of the design factors absent), the missing
  values are estimated so that the residual sum of squares of the lowest
  stratum is smallest (the next stratum up for what that leaves open, e.g. a
  whole plot), the residual df are reduced, and the completed design is
  analysed by strata. Means, stratum mean squares and the missing value equal
  Genstat's (oats with one sub-plot missing: 113.0; Victory x 0.6 120.7; the
  strata combined in 2.0-dev gave 121.8), Yates' formula for randomized
  blocks is reproduced, and the SEDs are exact (Genstat's are those of the
  complete design, "not adjusted for missing values"; they agree for the
  comparisons the missing plot does not touch). A message recommends REML,
  which uses the information in all strata. Designs that cannot be completed
  to an orthogonal one keep the combination of the strata, with a message.

* **`residplot()` is rewritten with ggplot2 and works for every model class
  predictmeans handles**: `lm`, `aov` (also with `Error()`), `glm`, `gls`,
  `lme`, `lmer`, `glmer`, `glmmTMB`, `gam`, `semireg` fits and other classes with
  `fitted()` and `residuals()` methods. Four plots: residuals vs fitted (with
  a loess smooth), normal Q-Q of the residuals, then by model the random
  effects' Q-Q plot, the autocorrelation of the normalized residuals (`gls`)
  or scale-location, and Cook's distances (`lm`, `aov`, `glm`) or fitted vs
  observed (a calibration plot for a 0/1 response). Poisson, binomial and
  negative binomial models use randomized quantile residuals (Dunn and Smyth
  1996), normal when the model is right. The axis labels of the 1.x plots of
  mixed models were wrong (swapped in the Q-Q plot; raw residuals called
  standardized). `interactive = TRUE` shows the plots with plotly (hover for
  each point's row); `plot = FALSE` returns the ggplot objects to change or
  combine. `id` now labels the most extreme points (3 by default) instead of
  asking for clicks.

* **Type III tests no longer depend on the model's contrasts** (`permmodels(type
  = "III")` and its tables). They are the tests of sum-to-zero contrasts, as SAS
  Type 3, computed from the model's own coefficients without refitting. Before,
  as `car::Anova()`, a model fitted with R's default treatment contrasts had
  each main effect tested at the first level of the other factors. Oats with a
  whole plot missing: Variety F 1.28 on 8.26 df and Nitro 37.92, as SAS (were
  1.17 on 22.1 and 12.41). Interactions, and models fitted with sum-to-zero
  contrasts, are unchanged.

* Comparisons of means over empty cells match only the factors whose effects
  the two means share: with `Sample_date/Farm` in the model, two dates are no
  longer matched on farms (each date has its own farm effects), which gave NA
  for pairs that are estimable; two treatments at one date still are. Found
  on the maintainer's pasture data (Sample_date:ExForest: 0 NA pairs, was 11).
* Fixed: `predictmeans()` stopped ("missing value where TRUE/FALSE needed")
  when some pairs cannot be estimated (confounded factors, e.g. a spray volume
  used only with fine droplets) and the df vary by pair (Kenward-Roger,
  Satterthwaite, inner-outer). Those pairs are NA, with NA df; the others are
  given. The same error stopped predictmeans 1.1.1.

* **`glm` models with a known dispersion (Poisson, binomial, negative
  binomial) use z tests** (normal; chi-square for joint tests), as emmeans and
  SAS GLIMMIX; they used t with the residual df, as for `lm`, which is
  conservative: in simulations (`dev/sim-df-glm/`) it rejected 0.3–4.5% of
  true null hypotheses at the 5% level (F tests 0.3–3.0%), z tests 3.6–5.7%.
  `df_method = "residual"` gives the 1.x t tests; `settings$df_method` is
  `"asymptotic (z)"`. Quasi families and Gaussian/Gamma `glm` keep the
  residual df (their dispersion is estimated).

* The Kenward-Roger F tests in the `permmodels()` tables of `lmer` models
  (`test.statistic = "F"`) use the observed information, as `pm_joint()` and
  the rest of the package; they used the expected one, as `car::Anova()`.
  The two agree for balanced data.

* **ML fits of `lme` and `gls` get Satterthwaite df** (from the ML
  log-likelihood), as ML `lmer` fits; they used nlme's df or the residual df,
  which are still `df_method = "containment"`. In simulations of ML `gls`
  fits (`dev/sim-df-ml/`), unstructured covariance: t tests 5.5% (residual df
  6.9%), F tests 6.3–8.8% (9.3–11.7%); AR(1) the same either way. Kenward-Roger
  asked for an ML fit gives Satterthwaite with a warning, as for `lmer`. An ML
  `lme` fit now gives the same df as the ML `lmer` fit (lmerTest).

* **Means over empty cells are compared on the cells both means have.** With
  `weights = "equal"`, a mean that averages over cells without data keeps the
  1.x rule (the average of its observed cells), but two such means average
  over different cells, so their difference also carried the effects of the
  other factors: in simulations (`dev/review-defaults.md`,
  `dev/sim-empty-cells.R`) it rejected 8–78% of true null hypotheses when the
  other factor had an effect. Now a pair with such a mean is compared on the
  cells both means have (4–6% in the same simulations), and the new column
  `over` of `pairs` says which (`"none"`: no cell in common, the pair is NA).
  Example: oats with Victory lost at N = 0.6: 0.6 − 0 is 42.5 (Golden Rain and
  Marvellous) instead of 46.4, and 0.6 vs 0.4 has p = 0.08 instead of 0.03.
  `contrastmeans()` and `pm_joint()` (column `over`) do the same; means,
  and pairs of means without empty cells, are unchanged. The warning now
  names the empty cells and the means concerned, and is no longer given
  where no mean averages over an empty cell (e.g. for the cell means
  themselves).

* New `df_method = "KR-average"` (in `predictmeans()`, `contrastmeans()`,
  `pm_joint()`, `pm_slopes()`, `covariatemeans()` and `df_term()`, for
  `lmer`, `lme` and `gls`): Kenward-Roger with the average information of
  the variance parameters, y'P G_i P G_j P y / 2 (Gilmour, Thompson and
  Cullis 1995), the matrix of Genstat's REML iterations. It gives Genstat's
  REML df (oats with a whole plot missing: Variety 8.75, Genstat 8.7; AR(1)
  with missing values: Group x Time 25.76, Genstat 25.8). Genstat's SEDs
  are the unadjusted ones (`vcov_adjust = FALSE`). The default stays
  Kenward-Roger with the observed information.

* `pm_joint()` without `ctrmatrix` writes each set of contrasts as each
  level against the last (SAS's Type 3 basis). Wald and Kenward-Roger F
  tests are unchanged; the Satterthwaite df of a test with several df depend
  on the basis (Fai and Cornelius 1996) and now equal SAS's (AR(1) repeated
  measures: 26.1 df for Time, was 25.5). New tests against a SAS run with
  AR(1) and missing values: predictmeans' Kenward-Roger equals SAS
  `DDFM=KR(FIRSTORDER)` everywhere and has the df of `KR` and `KR2`;
  Satterthwaite and between-within agree everywhere.

* **`aov()` models with `Error()` strata are analysed by strata**, as
  Genstat's ANOVA does, in `predictmeans()`, `contrastmeans()`, `pm_joint()`
  and `pm_slopes()` (new `R/aovlist.R`). They were refitted as `lmer`
  models, which changed the analysis: when a stratum variance is negative by
  ANOVA, REML fixes it at zero and the treatments were compared with the
  wrong error (tree-improvement Exam 3.1, diameter: SED 1.5 on 75 df, where
  the plot stratum gives 1.074 on 4 df, as Genstat). Now each comparison
  uses the mean squares of the strata it draws on, with Satterthwaite's df
  across strata (Genstat's rule: split plot, 0.583 on 14.69 df); the F tests
  are the stratum F tests; `settings$df_method` is `"multistratum"`.
  Standard errors of single means include the block stratum (Genstat's
  e.s.e., which leaves it out, is the SED divided by sqrt(2)). `residplot()`
  shows the residuals of the lowest stratum (`proj()`). To get the 1.x
  numbers fit the mixed model yourself, e.g.
  `lmer(y ~ trt + (1 | block/plot))`.

* New tests from the examples of *Experimental Design and Analysis for Tree
  Improvement*, against Genstat and SAS runs: nested treatments (Genstat's
  means, which average over the treatments that exist; SAS and emmeans
  weight by the nesting), a split plot (Genstat's SEDs, and its 14.69 df for
  comparisons across strata, equal to predictmeans' Kenward-Roger df), an
  incomplete site-by-seedlot table (SAS LS-means), families as random
  effects (SAS PROC MIXED), trees within plots (Genstat, by strata),
  inoculation on whole replicates crossed with provenances within countries
  (Genstat, by strata: 2, 116 and 31.61-117.16 df) and a clone trial with
  standards and random rows, columns and clones (SAS `DDFM=KR`: means, SEs,
  df and F tests).
* The warning "Missing treatments' combination appeared" is no longer given
  for combinations absent by design (a factor nested in another, such as
  provenances within countries); the means are the same. It is still given
  for combinations missing from the data.

* **GLMMs (`glmer`, `glmmTMB`) get the df of the grouping level at which
  each effect varies** (`df_method = "inner-outer"`, the new default), in
  place of the residual df (n minus parameters). This is nlme's rule
  (Pinheiro and Bates 2000), computed from the design (`R/df-innerouter.R`;
  equal to nlme's `fixDF` for nested random effects) and extended to
  crossed random effects: a treatment applied to whole clusters gets the
  number of clusters minus the treatments, an effect within clusters the
  residual stratum's df. A contrast gets the smallest df of the coefficients
  it uses; the means of a term share the term's df. In simulations
  (`dev/sim-df/README.md`) the residual df rejected 8-15% of true hypotheses
  for between-cluster effects, inner-outer 6-8% (nested and crossed); the
  rest is the Wald standard error with few clusters, for which
  `bootmodels()` is the check. `df_method = "residual"` gives the previous
  df. `"containment"` and `"between-within"` now work for `glmer` too, and
  `"inner-outer"` for `lmer` and `lme`.

* **Kenward-Roger df and adjusted SEs are the default for `lme` and `gls`
  models fitted by REML** (as for `lmer`), in `predictmeans()`,
  `contrastmeans()`, `pm_joint()`, `pm_slopes()` and `covariatemeans()`.
  They were nlme's df (`lme`) and the residual df (`gls`), which in
  simulations rejected 16% of true hypotheses in F tests of unstructured
  repeated measures and 24-27% with 10% of values missing; Kenward-Roger
  4-6% (`dev/sim-df/README.md`). `df_method = "containment"` gives the
  previous df (and, by default, the unadjusted covariance), so 1.x results.
  ML fits keep nlme's or the residual df (with a message); a covariance
  structure that cannot be rebuilt falls back to them with a message.
* `covariatemeans()` gains `df_method` and `vcov_adjust` and follows the same
  defaults: Kenward-Roger df for each point of the curves and adjusted SEs
  for REML fits. It used the unadjusted covariance and, for `lmer`, one df
  per term; `vcov_adjust = FALSE` (and `df_method = "containment"` for
  `lme`, `gls`) give its previous results.

* **Kenward-Roger uses the observed information of the variance parameters**
  (as SAS PROC MIXED `DDFM=KR`), for `lmer` (the default for REML fits),
  `lme` and `gls`. It stays valid with missing values, where the expected
  information (pbkrtest, emmeans) is not (Kenward and Molenberghs 1998). In
  simulations (`dev/sim-df/README.md`) it was never worse, and with
  unstructured repeated measures and 10% of values missing its F tests
  rejected 6% of true hypotheses at the 5% level where the expected
  information rejected 8.5%. Balanced designs give the same results as
  before; with missing values the df and SEs change a little (split plot with
  a whole plot missing: 8.35 df, was 9.25). `df_method = "KR-expected"` gives
  the previous results (and pbkrtest's, emmeans'). The option `"KR-SAS"` of
  the development version is gone: it is now `"Kenward-Roger"`. The tables of
  `permmodels()` keep pbkrtest's version, as `car::Anova()`.

* **One set of inference rules for `predictmeans()`, `contrastmeans()`,
  `pm_joint()`, `pm_slopes()` and `covariatemeans()`** (new internal
  `pm_inference()`, `R/inference.R`): which covariance of the coefficients is
  used, which denominator df, and how they are reported, decided in one
  place. The 39-call `predictmeans()` snapshot and the 33 permutation calls
  are unchanged apart from `settings$vcov_adjust`. Where the functions used
  to disagree they now follow the rule `predictmeans()` used:
  * `lme` models with nlme's df: `contrastmeans()` and `pm_joint()` used the
    smallest df of the term and its marginal terms (and printed a note); they
    now use the term's df, as `predictmeans()` and nlme's `anova()` do. One df
    per term is liberal for some comparisons (e.g. varieties within a nitrogen
    level of a split plot) and conservative for others; Satterthwaite or
    Kenward-Roger df (`df_method =`) give each contrast its own df.
  * `contrastmeans()` on an `lmer` model with a single fixed term gave every
    contrast the term's average df; each contrast now has its own
    Kenward-Roger or Satterthwaite df.
  * `pm_slopes()` with `df_method = "containment"`: the df of the terms with
    the covariate (the smallest), as with `"between-within"`; it used the
    factor's df.
  * The legacy p-value matrix of `predictmeans()` (`pairwise = TRUE`) used
    the mean df of the pairs for `lme`/`gls` with Satterthwaite or
    Kenward-Roger df; it now uses each pair's df, as the `pairs` table (and
    the letters) do.
  * `settings$vcov_adjust` is `TRUE` only when the Kenward-Roger adjusted
    covariance was used (it was `TRUE` for `lm`, `glm` and ML fits too).
  * `pm_joint()` reports nlme's df as `"containment (nlme)"` (was `"nlme"`),
    as `settings$df_method` does; its result has plain row names.

* **`ncore = 1` by default** in `permmodels()`, `permlmer()` (and
  `permglmer()`, `permglmmTMB()`), `permanova.lmer()` and `predictmeans()`
  (it was 3; `bootmodels()` already used 1). Results are the same; scripts
  that relied on three cores for speed should set `ncore`. CRAN and many
  shared machines allow at most 2 cores.

* **`permmodels()` uses Freedman-Lane permutations by default**
  (`method = "auto"`, which chooses the more robust scheme for the model
  class; at present Freedman-Lane for every class). Each term, coefficient
  or pair is then tested with the model without it, so the other terms keep
  their effects; permuting the response (1.x) tests the global null. In a
  simulation study (study 1b, `dev/sim-study/`), with a covariate in the
  model, permuting a GLM's response rejected only 0.6-1.9% of true null
  hypotheses at the 5% level (little power), Freedman-Lane 4.2-6.1%; for
  linear mixed models Freedman-Lane with refits kept 3.6-5.4% in every
  case, and pairs of means 5.2%. The result has no `permlist` of refits but
  can still be given to `predictmeans()` and `contrastmeans()`
  (`permlist = `). Use `method = "raw"` for the 1.x scheme and results.

* **Permutation tests of mixed models decorrelate residuals with the
  symmetric square root of their covariance** instead of the Cholesky
  factor: `permmodels()` for `lme`, `gls` and `lmer` models and
  `permlmer()` (so `permanova.lmer()`). The Cholesky factor depends on the
  order of the data (it treats the first unit of each block differently), so
  results changed when the data were sorted differently, and tests of terms
  that vary within blocks rejected too often: in simulations 6-8% of true
  null hypotheses at the 5% level, increasing with the number of blocks;
  5.2-5.3% with the symmetric root (`dev/permutation-plan.md`, §3.3).
  Permutation p-values of these functions therefore differ from earlier
  versions. The root is computed from q x q matrices (q random effects),
  without forming any n x n matrix, so it is also much faster with crossed
  random effects: 0.11 s instead of 3.9 s for n = 3000, and 21 s for
  n = 20000, where the Cholesky factor took over 30 minutes and 2.9 GB.
  Confirmed with the package functions (1000 data sets each): `permlmer()`
  split plot 5.9%, crossed design 5.4% (6.4% before); `permmodels()` with no
  fixed effects 4.8% (`lmer`, 6.2-7.1% before) and 5.9% (`lme`).

* With a `permlist`, `predictmeans()` and `contrastmeans()` no longer use
  the Kenward-Roger adjusted covariance (`vcov_adjust`): the permutation
  p-value compares the observed t statistic with those of the permuted
  refits, which use the unadjusted covariance, so the observed one must use
  it too. Before, unbalanced REML `lmer` fits compared an adjusted observed
  t (slightly smaller) with unadjusted permuted ones, so p-values were
  slightly too large; balanced designs are unchanged. `vcov_adjust = TRUE`
  given explicitly with a `permlist` gives a warning.

* **Kenward-Roger adjusted standard errors are now the default** for `lmer`
  models fitted by REML, in `predictmeans()` and `contrastmeans()`
  (`vcov_adjust = TRUE`). This matches SAS `DDFM=KR` and emmeans. Standard
  errors, confidence limits, SEDs, LSDs and p-values can differ slightly
  from predictmeans 1.x (e.g. 7.896 instead of 7.887 for an unbalanced Oats
  example); balanced designs barely change. Use `vcov_adjust = FALSE` to
  reproduce 1.x results.
* **glmmTMB is no longer attached** by `library(predictmeans)`; it moved
  from Depends to Imports. Scripts that fit models with `glmmTMB()` need
  `library(glmmTMB)`. lme4, nlme and lmerTest are still attached.
* `permmodels()` gives the same tables (including Satterthwaite df and
  p-values) whether an lmer model was fitted with `lme4::lmer()` or
  `lmerTest::lmer()`.
* `glmer` Poisson models now give back-transformed means as counts (whole
  numbers, never below zero), as `glm` and `glmmTMB` Poisson models already
  did.
* `semipred()` limits are now two-sided at the stated level: the default
  `level = 0.05` gives 95% limits. Earlier versions used `qnorm(1 - level)`,
  i.e. 90% limits with the default; `level = 0.1` gives those limits.
* `predictmeans(letterCI = TRUE)` now works with `atvar`: the CI-overlap
  letters compare the means within each `atvar` group (as in predictmeans
  1.1.2). Before, `atvar` was dropped and the letters compared all means.
  The note under the mean table now says the letters come from the overlap
  of, e.g., 83.4% CIs, and gives the approximate test level (about 0.05),
  instead of "significant level '0.166'". `level` keeps its meaning (it sets
  the CIs). With two `atvar` factors the note is one line
  ("at each level of A:B") instead of two.


## New features

* **Tensor products in `semireg()` and `semireg_tmb()`: `t2(x1, x2)`**
  (Wood, Scheipl and Faraway 2013), a surface for covariates on different
  scales, each penalty a random effect with its own variance, the fixed
  terms `x1 + x2 + x1:x2` added, `by = f` for a surface per level; with
  `semipred()`, `predictmeans()`, `anova()` and `varcomp()` as for other
  smooths. It is set up without mgcv's sum-to-zero constraint on the whole
  product, which leaves part of the penalised space in the unpenalised one
  (so the fixed part would not be those terms, and with a random intercept
  the population surface shifted by 0.04 on a test surface); the fit equals
  `gam()` with that basis and its penalties, and `gam()`'s own `t2()`
  without other random effects (tested). `te()` is fitted as `t2()` with a
  message (its overlapping penalties cannot be random effects); `ti()` is
  refused with a message.

* **`smZ()` and `smoothZ` are superseded** (maintainer's decision of
  2026-11-17): they keep working, so that 1.x code runs, but get no new
  features; `?smZ` has the table of their formula equivalents. Every help
  page, vignette and article now writes smooths in the formula, with more
  examples of the semi-parametric system: curves in a mixed model, a
  split-plot design with curves over time, crossed factors (`s(x, f1*f2)`),
  an additive model with a varying coefficient, a basis of your own, binary
  responses, subject-specific curves, counts with zero inflation, AR(1)
  residuals, surfaces (`?semireg`, `?semireg_tmb`, `?semipred`, section
  "More situations" of `vignette("smooth")`). The migration table said
  `smZ(cbind(x1, x2), type = "Ztps")` is `s(x1, x2, bs = "ztps")`: the
  defaults differ (6 and 10 knots per axis), it is `bs = "ztps", k = 6`.
* **New data sets: the data of M. P. Wand's semiparametric regression
  examples** (the maintainer's golden examples): `growthIndiana`
  (subject-specific growth curves, black-white contrast curve), `SBMDfem`
  (bone density by ethnic group, additive mixed model), `ozoneCalif`
  (additive model), `scallop` with `scallopKnots` and `scallopBdry`
  (thin-plate surface inside a boundary), `wageTUnion` (penalised spline
  logistic regression) and `copperKnots` (for `copper`), each with its
  analysis in the formula syntax; `semireg()` reproduces the original
  `lme()` fits (`tests/testthat/test-semireg-golden.R`).

* **`s(x, bs = "user", xt = list(basis = f))`: a basis of your own in the
  formula of `semireg()` / `semireg_tmb()`**, the formula form of a matrix
  in `smoothZ`: `f(x)` (or `f(x1, x2)`) returns a matrix with a row per
  value, its columns get iid random coefficients with one variance, and the
  constant and the covariate are the fixed part, added to the formula
  (`xt$null = "constant"` or `"none"` for less). Prediction evaluates `f` on
  the grid, so `semipred()` (also `deriv`), `predictmeans()`, `by`,
  `s(x, f)` and random curves work. The function is checked before the fit,
  which stops with what is wrong and how to put it right when it fails on
  the data, does not return a numeric matrix with a row per value (also for
  a single value, as at prediction), gives NA or infinite values, gives
  other values on a second call, gives other rows for the same values on
  part of the data (knots or a range taken from the data: prediction would
  use other ones), or has dependent columns or columns in the fixed part; a
  warning when it fails just beyond the ends of the data. With the
  O'Sullivan or Ztps basis it reproduces the `bs = "os"` / `"ztps"` fits
  exactly (tested, both engines).

* **`semipred(deriv = 1)` and `deriv = 2`: slopes and curvatures of the
  fitted curves** along the covariate, with pointwise or simultaneous bands
  and the same df as the curves; with `contr` the difference of two slopes.
  The ranges where a slope is above or below zero (the curve rises or
  falls) are printed and kept in `attr(pred_df, "significant")`. Central
  differences of the whole prediction design, so any basis works; for the
  O'Sullivan basis (`bs = "os"`, `smZ()`) they equal the exact derivatives
  of `ZOSull(drv = )` to 1e-8 (tested). Generalized models: link scale.
  The bands covered the true slope and curvature 0.97-1.00 at the 95% level
  (`dev/sim-semipred-deriv.R`). This replaces `smZ(drv = )` for formula
  smooths.

* **`anova()` of `lmer` fits with df and p-values, without lmerTest**: Type
  III F tests with Kenward-Roger df by default (as SAS `DDFM=KR`),
  `ddf = "Satterthwaite"` (lmerTest's numbers), `type = "II"` or `"I"`
  (sequential, lmerTest's), `ddf = "lme4"` for lme4's own table.
  `anova(m1, m2)` is lme4's likelihood ratio comparison, unchanged. Also for
  gaussian `semireg()` fits (`anova(fit$semer)`). With lmerTest attached as
  well there is no conflict: `lmerTest::lmer()` fits keep lmerTest's
  `anova()` and `summary()`, `lme4::lmer()` fits get this method, and
  `lmerTest::ranova(fit$semer)` still refits a `semireg()` fit with its
  smooth bases (tested).
* **`anova()` of `semireg()` / `semireg_tmb()` fits**: the fixed effects
  (Kenward-Roger F for gaussian, Wald chi-square for generalized fits) and,
  for each smooth term, the likelihood ratio test of its variance being zero
  (refit without it; p-value halved for the boundary, Self and Liang 1987;
  `lmerTest::ranova()` used the full chi-square). For `s(x, f)` the
  deviation block tests whether the levels' curves have the same shape.
* **Formula smooths: each smooth block is scaled as a random intercept**
  before the fit (the same model; only the reported variance of a smooth
  block is on another scale). lme4 then starts near the optimum: an
  `s(x, by = f, id = 1)` fit had stopped at a worse local optimum, and fits
  are faster (growth curves of 2257 rows: 8 s; `lme` takes 31 s).
* **`vignette("smooth")`** is written with the smooths in the formula, with
  a table of the forms, a section on what predictmeans adds (means with
  Satterthwaite df, tests of smooth variances, rate ratios between count
  curves, subject-specific curves) and a table for moving from the
  `smoothZ = list(smZ())` form, which still runs.
* **`varcomp(model, test = TRUE)`: a test of each variance component** for
  `lmer`, `glmer`, `glmmTMB` and `semireg()`/`semireg_tmb()` fits: columns
  `LRT` and `p.value`, each row's hypothesis and reference distribution
  printed below the table. Each row is tested by refitting the model one
  would fit without it: a term dropped, `(1 + x | g)` reduced to `(1 | g)`
  for the slope, covariances set to 0 (`||`, glmmTMB `diag()`), the
  correlation of `ar1()`/`ou()`/`cs()` set to 0 (`homdiag()`). The
  likelihood ratios equal SAS PROC MIXED and Genstat REML runs; the
  boundary mixtures (0.5 chi2(0) + 0.5 chi2(1) for a variance, 0.5
  chi2(k-1) + 0.5 chi2(k) for a slope) equal SAS PROC GLIMMIX `COVTEST`.
  When the tested variance is the only random effect of a gaussian REML fit
  the exact null distribution of the restricted likelihood ratio
  (Crainiceanu and Ruppert 2004) is used: there the mixture was
  conservative (2-3% at the 5% level, the exact test 4.3%); with other
  random effects the "exact" distribution was too liberal (7%) and the
  mixture is used (5.0%) (`dev/sim-vc-test.R`). Rows where the reference
  may be inaccurate (fewer than 6 levels; binary data with fewer than 5
  observations per level) are named in a warning with how to check them.
  `anova()` of a `semireg()` fit uses the same tests for its smooth terms
  (before: always the halved chi-square). `test_method = "bootstrap"`
  gives a parametric bootstrap of each likelihood ratio (responses
  simulated from the model without the component, both models refitted;
  `nboot`, `ncore`), with sequential stopping (a clearly non-significant
  test stops after 20 simulated statistics reach the observed one) and one
  seed per draw (the same p-values with any number of cores); the warning
  points to it. In a simulation of seven designs (`dev/sim-vc-defaults.R`)
  the default was never liberal (0-4.5% at the 5% level); for GLMMs it was
  conservative (binary, clusters of 4: 0.5%; Poisson 4%; negative binomial
  0%), where the bootstrap was nearer 5% (1.5-6%) and added 5-7 points of
  power. For `s(x, f1*f2)` deviation curves every reference was conservative
  (0%), with full power for a real difference in shape. `?varcomp` and
  `vignette("resampling")` say which to use: the default for gaussian models
  and a first look; the bootstrap for GLMMs, small data, flagged rows and
  p-values near the cut-off.
* **`semireg_tmb()` gaussian fits are REML, as documented.** The default
  `REML = TRUE` never reached glmmTMB (the call passed only the arguments
  given), so since 1.x gaussian `semireg_tmb()` fits were maximum
  likelihood, with variance components biased downwards. They are now REML
  and equal `semireg()` (lme4) on the same model. Other families keep
  maximum likelihood (Laplace), as `glmer`, unless `REML = TRUE` is given.
* **`predictmeans()` on gaussian `semireg_tmb()` fits: Satterthwaite df**,
  as for `semireg()`, also with glmmTMB covariance structures such as AR(1)
  residuals (`ar1(tf + 0 | id)`, `dispformula = ~0`). The df equal those of
  the lme4 engine on the same model (to 1e-4). In a repeated-measures
  simulation with AR(1) residuals and no treatment effect, t tests rejected
  5.0% at the 5% level (df 18-49), z tests 5.5% (`dev/sim-semireg-tmb-df.R`). In both
  engines, variance components estimated at zero are left out of the
  Satterthwaite df, as in Kenward-Roger (a log-SD far out on a flat
  likelihood had given df below 1).
* **`predictmeans()` on gaussian `semireg()` fits: t tests with
  Satterthwaite df** instead of normal-based tests. The means at a
  covariate value include the smooths' random coefficients, so the df are
  those of the combination of fixed and smooth coefficients: the prediction
  error variance as a function of the variance parameters (as SAS
  `DDFM=SATTERTHWAITE` for estimates with random effects;
  `pm_sr_satt_ddf()`). In a split-plot design with curves over time the
  treatment means at a time get about 9 df and their differences about 17,
  not infinity. With no treatment effect these t tests rejected 4.6% at the
  5% level, z tests 5.2-6.8% (`dev/sim-semireg-df.R`). For fixed effects
  alone the df equal lmer's Satterthwaite df. Generalized fits keep
  normal-based tests. `anova()`'s tests of the smooth variances are
  conservative: 1-3% rejections at the 5% level with no difference in
  shape (`dev/sim-semireg-lrt.R`).
* **Random curves `(s(x) | g)` in `semireg()`**: the group-specific curve
  model of Durban et al. (2005) and Harezlak, Ruppert and Wand (2018, section
  9.2): `(x | g)` (correlated intercept and slope per subject) plus a spline
  deviation per subject, one variance for all subjects, e.g.
  `height ~ s(age, bs = "os", k = 20) + (s(age, k = 10) | id)`. The curves
  are random effects, so `semipred()` and `predictmeans()` give population
  curves without them; `anova()` tests whether they are needed.
  `s(x, bs = "os", xt = list(knots = , range = ))` and
  `s(x1, x2, bs = "ztps", xt = list(knots = ))` set knots and range.
* **Golden examples** (the maintainer's lme scripts of penalised splines;
  data in `system.file("extdata", "semireg", package = "predictmeans")`):
  additive models, an additive mixed model, thin-plate surfaces with given
  knots, group-specific curves with a black-white contrast curve (as HRW's
  `maleGrowthIndiananlme`). `semireg()` with formula smooths gives the same
  REML log-likelihood and fitted values (to 2e-5 SD) as `nlme::lme`, and
  equals `gam(method = "REML")` where the script used `gamm()` (whose `lme`
  fit stopped short of the REML optimum) or `gam()`
  (`tests/testthat/test-semireg-golden.R`, `dev/semireg-golden/`).
* **`semireg_tmb()` with glmmTMB covariance structures** (`ar1()`, `ou()`,
  `cs()`, `toep()`, `diag()`, `exp()`, ...) next to formula smooths: the
  formula was rebuilt from `findbars()`, which keeps only the bar inside
  such a term, so `ar1(tf + 0 | id)` was fitted as an unstructured
  covariance without a message (78 parameters instead of 2 for 12 times).
  The terms are now kept as written, also in `anova()`'s refits and in
  `bootmodels()`' null model for glmmTMB fits. The covariance of the
  coefficients is computed with the matrix scaled to unit diagonal, so AR1
  residuals with `dispformula = ~0` (or a nugget variance near 0) no longer
  fail as "computationally singular". `semireg()` (lme4, which has no such
  structures) refuses them with a message.
* `semireg_tmb(family = nbinom1)` (a glmmTMB family named without glmmTMB
  attached) failed with "object 'nbinom1' not found" after the engine
  change; the family is found in glmmTMB again, as in 1.x.
* Formula smooths whose covariate has fewer distinct values than basis
  columns (e.g. a depth measured every 5 m), and factor smooths (`bs = "fs"`)
  with subjects of few observations, failed with "singular matrix"; the map
  from the basis to the random effects is now `smooth2random()`'s own
  transformation, exact in these cases.
* **Smooths in the formula of `semireg()` and `semireg_tmb()`**, as for
  `mgcv::gam`, next to fixed and lme4 random terms:
  `semireg(y ~ A + s(x, by = A) + (1 | Block), data = d)`. Each smooth is set
  up by mgcv and split, as in gamm4, into fixed terms and blocks of random
  effects. The fixed terms are added automatically (`x` for `s(x)`,
  `f + f:x` for `s(x, by = f)`, `z + z:x` for a numeric `by`), so they cannot
  be forgotten; terms already written are not doubled.
  - Types: any one-penalty mgcv basis (`tp`, `cr`, `ps`, shrinkage `cs`/`ts`),
    `by` factors (`id =` for one shared variance) and numeric `by`, factor
    smooths `bs = "fs"`, isotropic `s(x1, x2)`.
  - Factor terms inside `s()` give an ANOVA-like split of the curves:
    `s(x, f)` is a common curve plus a deviation curve per level of `f` (one
    variance for the deviations); `s(x, f1*f2)` adds deviation curves for
    `f1`, `f2` and `f1:f2`, a variance each, so `varcomp()` shows which factor
    terms change the shape of the curves; `s(x, f1:f2)` uses the cells;
    `s(x1, x2, f)` does the same for surfaces. Fixed terms `x + f1*f2 + f1:x
    + f2:x + f1:f2:x` are added. `s(x, f)` equals the same model written for
    `lme` with `pdIdent` blocks (Harezlak, Ruppert and Wand).
  - The package's own bases as `bs = "os"` (O'Sullivan, `k` interior knots)
    and `bs = "ztps"`: the same fits as `smZ()` in `smoothZ`.
  - Fits equal `gam(method = "REML")` (fitted values within 3e-5 in 11
    cases), and `semipred()` / `predictmeans()` give `predict.gam()`'s curves,
    means, SEs and curve differences (to 1e-5).
  - Tensor products (`te()`, `t2()`) are not available yet: mgcv's `t2()`
    design could not be split reliably (`dev/review-semireg.md`).
  - mgcv moves from Suggests to Imports.
* New vignette "Checking the model" (`vignette("checking")`): `pm_check()`
  and `residplot()` on a split plot and on overdispersed counts (Poisson, then
  negative binomial), with what the checks do and do not decide. The main
  vignette has a table of the default df by model class (Kenward-Roger for
  REML fits, Satterthwaite for ML fits, inner-outer for GLMMs, z tests for
  Poisson/binomial/negative binomial `glm`, residual df otherwise, strata for
  `aov` with `Error()`), with `"KR-expected"` and `"KR-average"` among the
  options. `dev/benchmarks.R` also times permutation tests and means over
  empty cells.

* New tests against Genstat 23.1 REML (`test-genstat-reference.R`; runs by
  the maintainer, output files in `dev/genstat/`): six analyses of the oats
  split plot with missing plots, an empty cell and scattered missing sub-plots,
  and AR(1) repeated measures with missing values. Variance components,
  predicted means, SEDs (`vcov_adjust = FALSE`), Wald statistics and 17 of 18
  denominator df (`df_method = "KR-average"`) agree; Genstat's F is
  Kenward-Roger's scaling factor times Wald / df. Genstat's means over the
  combinations present equal predictmeans' means over empty cells.

* `pm_check(model)`: a table of checks of a fitted model with a plain note on
  each, for any model predictmeans handles: the fit (singular, convergence),
  outliers (beyond 3, flagged when improbable for normal residuals),
  influence (Cook's distance), normality of the residuals and of each
  random-effect term, equal variances across the levels of each factor
  (Brown-Forsythe), overdispersion of Poisson and binomial models,
  interactions (Type III tests), whether a factor with numeric levels follows
  a straight line (curvature contrasts with the actual spacing, by
  `pm_joint()`), and aliased coefficients. The checks inform; they do not
  simplify the model. `plot = TRUE` also draws `residplot()`.

## Stable result tables

* `predictmeans()` results gain three slots with fixed names and unrounded
  values. `tidy()` (generics/broom) and `as.data.frame()` return them:
  * `means`: one row per mean, with `estimate`, `std.error`, `df`,
    `conf.low`, `conf.high`, back-transformed `response` columns, `group`
    letters, `n` and `estimable`;
  * `pairs`: one row per comparison, with `level1`, `level2`, `by` (for
    `atvar`), `estimate`, `std.error`, `df`, `statistic`, `p.value`,
    `p.adj`, `lsd` and `p.method`;
  * `settings`.
  They are checked against emmeans in the tests. The legacy slots and the
  printed output are unchanged, and the result has class
  `c("predictmeans", "pdmlist")`.
* `coef()` and `vcov()` on a `predictmeans()` result give the means (model
  scale) and their covariance matrix, Kenward-Roger adjusted when the means
  are (equal to emmeans' `vcov()`). Any contrast of the means can then be
  tested with `multcomp::glht()` or `marginaleffects::hypotheses()`, e.g.
  `hypotheses(pm, hypothesis = "pairwise")` (marginaleffects is in
  Suggests; predictmeans registers its `get_coef`/`set_coef`/`get_vcov`
  methods; recent marginaleffects versions also need
  `options(marginaleffects_model_classes = "predictmeans")`).

## Plots and summaries from the result tables

* `autoplot()` for `predictmeans()` results gives ggplot objects that work
  in Shiny, reports and `plotly::ggplotly()`:
  * `type = "means"`, with confidence limits, half-LSD bars or no bars;
  * `type = "pvalue"`, a heatmap of the adjusted pairwise p-values in the
    colour bands of `PMplot()`, one panel per `atvar` level;
  * `type = "response"`, the back-transformed means.
* `summary()` prints the means with a footer stating the model and df
  method, the SED and LSD ranges, the p-value method and adjustment, what the
  letters mean, and the covariate values used.
* Both read only the `means` and `pairs` tables: they compute nothing new,
  open no graphics device and print nothing unless printed.

## Supporting new model classes

* `Kmatrix()`, `predictmeans()`, `contrastmeans()` and `covariatemeans()`
  now get what they need from a model through six exported generics:
  * `pm_coef()`, `pm_vcov()`: estimates and their covariance;
  * `pm_basis()`: design matrix for new data;
  * `pm_data()`: the fitting data;
  * `pm_df()`: df for each contrast;
  * `pm_linkinv()`: inverse link.

  The defaults work for any class with the usual `coef`, `vcov`, `terms`,
  `df.residual` and `family` methods. Supporting another class means writing
  methods for these generics only (see `?pm_coef`), and model classes
  predictmeans doesn't know are no longer rejected. Built-in classes behave
  exactly as before.
* Ordered factors keep their contrasts in the reference grid even when a
  model class stores no contrasts.

## Named contrasts

* New `specify_contrast()` writes contrasts with level names instead of
  counted columns, e.g.
  `specify_contrast(Kmatrix(fm, "A:B"), a = (`x:b1` + `y:b1`)/2 - `z:b1`)`.
  The result goes straight into `contrastmeans(ctrmatrix = )`. It checks
  that the weights sum to zero (unless `.zero_sum = FALSE`) and that no
  missing or non-estimable mean is used.
* `contrastmeans()`:
  * accepts a contrast matrix whose columns are named by level, in any order
    and for a subset of levels;
  * keeps the contrast names in its result;
  * gains `weights`.

## New options

* `permmodels(method = "freedman-lane")` for `lm` and `aov` models: each
  term (and coefficient) is tested by permuting the residuals of the model
  without it and adding them back to that model's fitted values (Freedman
  and Lane 1983). Use it when a covariate with extreme values is related to
  the tested term: in a simulation study (`dev/sim-freedman-lane.md`)
  permuting the response then rejected 6-7% of true null hypotheses at the
  5% level, Freedman-Lane 4-5%; in other designs the two agreed. The
  statistics come from QR projections of blocks of permuted responses, with
  no refits (4999 permutations of a 1000-row model with 16 tests take about
  11 seconds on one core). Weights and offsets are allowed. `permlist` is
  `NULL`; instead the result has a component `perm` (what is needed to test
  any comparison with the same permutations), and the whole result can be
  given to `predictmeans()` and `contrastmeans()` (argument `permlist`):
  each pair or contrast is then tested by Freedman-Lane permutations of the
  model in which that comparison is zero, so its p-value does not depend on
  the other means. The default, `method = "raw"`, is unchanged. For mixed models with such a covariate
  the help page recommends the Kenward-Roger or Satterthwaite F test: the
  raw scheme was slightly liberal there (6.3%), the F test kept 5%.

* **`permmodels(method = "freedman-lane")` for mixed models, GLMs and
  GLMMs** (`lmer`, `lme`, `gls`, `glm`, `glmer`, `glmmTMB`). Linear mixed
  models: the model without the tested term is fitted, its decorrelated
  residuals permuted and added back, and the model refitted to each
  permuted data set (one set of refits per term; pairs and contrasts get
  one set for their family); the statistics are the model's Wald t and F
  (unadjusted covariance). `refit = FALSE` fixes the covariance of the data
  at the fit instead (decorrelated least squares, no refits, much faster,
  but liberal when the variance parameters rest on few units: for the
  whole-plot factor of a split plot 8.6%, 7.2%, 6.0% and 5.6% with 6, 12,
  24 and 48 blocks, against 5.0%, 6.2%, 5.4% and 4.6% with refits; 6.2% for
  whole-plot pairs with 6 blocks). The help of `permmodels()` says when
  `refit = FALSE` is enough (error df of the tested terms above about 100). GLMs: the V-fixed version on
  the working (link) scale, with the dispersion estimated from the working
  residuals; GLMMs: the working linear mixed model, refitted to each
  permuted working response (or V fixed with `refit = FALSE`). Terms (types I, II, III), coefficients, and pairs and
  contrasts in `predictmeans()`/`contrastmeans()` all use the same
  permutations. `glmer` and `glmmTMB` models can now be given to
  `permmodels()` (Freedman-Lane only).

* **Confidence limits of pairs and contrasts**, also by permutation. The
  `pairs` table of `predictmeans()` gains `conf.low` and `conf.high` (limits
  of each difference at `level`), and `contrastmeans()` tables gain
  `lower.CL` and `upper.CL` (new argument `level`). With a `permlist` they
  are permutation-t limits: estimate +- c SE with c the critical value of
  the permutation test (the k-th largest permuted |t|, k = floor(level
  (nperm + 1))), from the same permutations; with a max-|t| adjustment
  (`adj = "tukey"`, `"mvt"`, `"westfall"`) c comes from the largest |t| of
  the family, so the limits are simultaneous (a permutation version of
  Tukey's intervals). The `lsd` column with a `permlist` is now this
  permutation LSD (it was 2 SED), and so are the ratio limits for log and
  logit models.

* `contrastmeans(joint = TRUE)` also tests all the contrasts together: the
  F test of `pm_joint()` (Kenward-Roger or Satterthwaite df for `lmer`,
  nlme's df for `lme`, the residual df otherwise; contrasts that are linear
  combinations of others are dropped), and with a `permlist` a permutation
  p-value of the joint Wald F from the same permutations. The result gains
  a third component, `Joint F test of the contrasts`.

* `permmodels(robust = TRUE)` for `lm`, `aov` and `glm` models whose error
  variances are unequal while the model assumes them equal:
  heteroscedasticity-robust statistics (HC3 sandwich t and Wald F;
  DiCiccio and Romano 2017), with the signs of the reduced model's residuals
  flipped at random instead of permuted (a wild bootstrap), so each residual
  keeps its own variance. In a simulation with group SDs 3, 1, 1 the usual
  and permutation F tests rejected about 20% of true null hypotheses at the
  5% level; robust permutations 6-9%, sign flips less (see
  `dev/sim-perm-robust.R`). Type II F values equal
  `car::Anova(white.adjust = "hc3")`. Pairs and contrasts from the result
  use the robust t, and their confidence limits use robust standard errors.
  The `perm` component of `permmodels()` results now prints a one-line
  summary instead of its contents.

* New `permstats()`: from the permutations kept by `permmodels()`
  (Freedman-Lane), for each term and coefficient a standardised effect size
  (observed statistic against the mean and SD of the permuted ones) and a
  check of the parametric test: its p-value and `param.size`, the proportion
  of permuted data sets (no effect of the term) in which it rejects at
  `level` (residual, Satterthwaite or nlme df; Wald chi-square for GLMMs).

* New `permcombine()`: nonparametric combination of permutation tests made
  with the same permutations, e.g. one model fitted to several responses on
  the same units (same `seed`, `nperm`): for each term, a combined p-value
  (Fisher or Tippett) and min-p adjusted p-values over the responses, with
  their correlation kept.

* New argument `vcov` in `predictmeans()`, `contrastmeans()`,
  `pm_slopes()` and `pm_joint()`: another covariance of the coefficients
  for the standard errors, tests and limits. For `lm`, `aov` and `glm`,
  `vcov = "HC0"` to `"HC3"` (heteroscedasticity-consistent; equal to
  `car::hccm()` and to `sandwich::vcovHC()`) or `vcov = ~ cluster`
  (cluster-robust CR1, as `sandwich::vcovCL()`), with no new dependency;
  for any model a covariance matrix or a function of the model (for mixed
  models e.g. from clubSandwich). The df stay the model's. Means, pairs and
  slopes equal emmeans with the same covariance (`vcov. =`). The `settings`
  table records it (`vcov`). With a `permlist` it is not used (see
  `permmodels(robust = TRUE)`).

* Singular `lmer`/`glmer` fits (a variance component estimated as zero, or
  a correlation as +-1) are noted by `predictmeans()` and `contrastmeans()`
  (SAS's "G matrix not positive definite") and recorded in
  `settings$singular`. Adjusted p-values use each pair's own df, as SAS
  `ADJDFE=ROW` (now tested and documented).

* `glance()` (generics/broom) of a `predictmeans()` result: one row with
  the term, numbers of means and pairs, df method, adjustment, covariance,
  SED range and mean LSD. `options(predictmeans.plot = FALSE)` turns the
  plots off by default (scripts, reports, pipelines). The examples of
  `?predictmeans` start from the `means`, `pairs` and `settings` tables.

* SAS validation: the tests compare every variance component, LS-mean,
  difference and Type 3 test of 16 SAS PROC MIXED runs (split plot,
  sub-plot or whole plot missing, unstructured repeated measures; containment,
  between-within, Satterthwaite and Kenward-Roger df) with predictmeans
  (`tests/testthat/test-sas-reference.R`).

* `df_method = "between-within"` (SAS `DDFM=BETWITHIN`) in
  `predictmeans()`, `contrastmeans()`, `pm_joint()`, `pm_slopes()` and
  `df_term()`, for `lmer`, `lme` and `gls` models: effects constant within
  subjects get the between-subject df, the others the df of their
  subject-by-effect interaction; with an unstructured covariance
  (`gls` with `corSymm`) every effect gets the between-subject df. Equal to
  SAS PROC MIXED on the split plot (Variety 10, Nitro 15, Variety x Nitro 30,
  29 with a plot missing) and the unstructured repeated measures (10). As
  for containment, the unadjusted covariance is used unless
  `vcov_adjust = TRUE`. `pm_slopes()` with containment df now also uses the
  unadjusted covariance by default, as `predictmeans()` does.

* Kenward-Roger with the observed information of the variance parameters,
  as SAS PROC MIXED `DDFM=KR` (computed in closed form); since 2026-11-13
  this is `df_method = "Kenward-Roger"` and the expected information
  (pbkrtest, emmeans) is `"KR-expected"` (see the first entry). Same results
  in balanced designs; with one plot missing in the split plot of the SAS
  reference run the observed information reproduces SAS's df and SEs (e.g.
  9.79 and 31.3, where the expected gives 9.87 and 31.46). With `"between-within"`, `"containment"` and
  `"Satterthwaite"`, every df of that SAS run now agrees with SAS.

* `df_method = "Satterthwaite"` for `lme` and `gls` models, with any
  correlation and variance structure, in `predictmeans()`,
  `contrastmeans()`, `pm_joint()` and `pm_slopes()`: the covariance of the
  data is rebuilt from the fit (checked against `vcov()`), and the
  Satterthwaite df use the Hessian of the REML log-likelihood, as SAS does.
  For the unstructured repeated measures of the SAS reference run, the df
  are 10 to 10.2 as SAS gives, where the residual df (the default for
  `gls`) are 30. `"Kenward-Roger"` and `"KR-expected"` give Kenward-Roger df,
  F tests and adjusted standard errors for these models too (numerical
  derivatives of the data covariance; the second-derivative term is left
  out, so the result does not depend on nlme's parameterisation and equals
  SAS DDFM=KR for variance components and unstructured covariances): the
  SAS run's KR F tests (Time 1616.03 with 9 df) are reproduced, and `lme`
  equals `lmer`'s Kenward-Roger. The default for `lme` and `gls` is
  unchanged.

* Equivalence tests: `predictmeans(equivalence = d)` (or
  `c(lower, upper)`) adds `p.equiv` to the `pairs` table, the TOST p-value
  (two one-sided t tests; below `level`, the difference lies within the
  margins), equal to emmeans' `test(..., side = "equivalence")`;
  `contrastmeans(equivalence =)` adds it to the contrast table. The margins
  are on the scale of the comparisons and are recorded in `settings`.

* Back-transformed means for GLMMs and transformed responses:
  `predictmeans(bias_adjust = TRUE)` gives h(eta) + sigma^2 h''(eta)/2,
  the mean over the population of random effects (default sigma: the SD of
  the random intercepts, plus the residual SD for a transformed response),
  and `average = "response"` back-transforms the cells of the reference
  grid before averaging them (the mean of the response, e.g. a mean count),
  both with delta-method SEs; equal to emmeans' `bias.adjust` and
  `regrid()`. They change only the `response` columns of `means`; the
  settings record them (`bias_adjust`, `average`).

* **Permutation Tukey and Westfall-Young adjustments.** With a `permlist`,
  `predictmeans()` and `contrastmeans()` take `adj = "tukey"` (or `"mvt"`):
  each comparison's |t| is compared with the largest |t| of its family in
  each permuted data set (single-step max-|t|, the permutation version of
  Tukey's test), and the new `adj = "westfall"`: the same step down
  (Westfall and Young 1993), uniformly more powerful. Families are all the
  pairs, or the pairs within each `atvar` level, or the rows of a contrast
  matrix. They use the correlation between the comparisons, so they are
  more powerful than Holm or Bonferroni. For the nitrogen pairs of the Oats
  example the permutation Tukey p-values agree with Tukey's to within 0.005.
  Before, these adjustments were not available with a `permlist`, and
  `contrastmeans()` ignored `adj` when given a `permlist`; now it applies it.

* **`permlmer()` tests GLMMs** (`glmer` or `glmmTMB` models; also as
  `permglmer()` and `permglmmTMB()`), random effects and fixed effects. By
  default counts and binomial data with several trials are permuted on the
  working (link) scale: the H0 model is linearised (working response and
  variances, the covariance of the data also having the random effects,
  as in `R2_glmm()`), the decorrelated working residuals are permuted, and
  the statistic is the likelihood-ratio statistic of weighted linear mixed
  models of the working responses. 0/1 responses use randomised quantile
  residuals mapped back to 0/1 and GLMM refits (`method = "quantile"`). In
  the simulation study (`dev/sim-study/`) the working scale kept the 5%
  level for counts and binomial data with trials (4.6-5.9%) and was much
  more robust than the usual likelihood-ratio test with overdispersed counts
  (9.8% against 63%); for 0/1 responses with rare events it failed (36%),
  and the quantile method kept the level (5.6-6.2%). glmmTMB families:
  poisson, binomial, nbinom1, nbinom2, Gamma, gaussian, beta (no zero
  inflation or dispersion model).

* `permanova.lmer()` also tests the terms of `glmer` models (each term by
  `permlmer()` on the working scale), and builds the reduced models where
  the model was fitted (it used `update()`, which looks for the data in the
  wrong place when the model was fitted inside a function).

* `permmodels()` and `permlmer()` (so `permanova.lmer()`) take a matrix of
  permutations, `perms`, e.g. from `permindex()`: permutations within
  blocks or of whole plots, for tests that permute only as the design
  allows (randomisation tests). With `method = "freedman-lane"` the same
  permutations are used for pairs and contrasts in `predictmeans()`.

* New `bootmodels()`: parametric bootstrap p-values for the terms and
  coefficients of `lmer` and `glmer` models, simulating from the model
  without its fixed effects (global null) and refitting. Its result can be
  passed to `predictmeans()` and `contrastmeans()` as `permlist`, giving
  bootstrap p-values for pairwise comparisons and contrasts (`p.method`
  `"bootstrap"` in the `pairs` table).

* `semipred(simultaneous = TRUE)` gives simultaneous confidence bands over
  the whole curve(s), and over difference curves with `contr` (Ruppert,
  Wand and Carroll): the critical value is the quantile of the largest
  standardised deviation from 10,000 draws, with a fixed seed so results
  repeat (the user's random numbers are left alone).

* `predictmeans()` accepts `semireg()` fits: the means of a factor term
  include the smooth terms (fixed part plus the smooth's random effects,
  as `semipred()`), at the covariates' means or at named `covariate`
  values, with pairwise comparisons (normal-based, any `adj`), letters,
  and `coef()`/`vcov()`. The result holds the `means`, `pairs` and
  `settings` tables. `semipred()`'s design step is now an internal
  function shared by both (its output is unchanged).

* New `R2_model()`: R-squared of a mixed model as numbers (class
  `"pm_r2"`, printed as percentages). `method = "piepho"` (default) gives
  the values of `R2_glmm()` as proportions; `method = "nakagawa"` gives
  Nakagawa and Schielzeth's marginal and conditional R-squared, with
  random slopes as Johnson (2014), for gaussian and binomial models
  (sleepstudy: 0.279 and 0.799, as `performance::r2_nakagawa()`).
  `R2_glmm()` is unchanged.

* New `pm_slopes()`: the slope of a numeric covariate within each level of
  a factor term, with SEs, df, confidence limits and tests, and pairwise
  comparisons of the slopes (any `adj`). Equal to emmeans' `emtrends()`
  (Kenward-Roger for `lmer`); for `poly()` or splines the slope is taken
  at `at` (default: the covariate's mean) by a central difference.

* Response scale: when the means are back-transformed (non-identity link
  or `trans`), the `means` table has `response.std.error` (delta method,
  equal to emmeans with `type = "response"`). For a log link or
  `trans = exp`, the `pairs` table has the ratio of the means with its
  confidence limits (`ratio`, `ratio.conf.low`, `ratio.conf.high`); for a
  logit link, the odds ratio (`odds.ratio`, ...). The limits use the same
  df as the differences.

* New `pm_joint()`: F tests of sets of contrasts among the means of a term.
  With `by`, the means are compared within each level of the `by`
  factor(s), as SAS `SLICE` (e.g. nitrogen within each variety); with
  `ctrmatrix`, any set of contrasts is tested together. For `lmer`: the
  Kenward-Roger F test equal to `pbkrtest::KRmodcomp()`, or Satterthwaite
  equal to `lmerTest::contest(joint = TRUE)`, or containment df; nlme df
  for `lme`; residual df otherwise.

* New `contrast_family()` builds the standard families of contrasts among
  the means of a term for `contrastmeans()`: `"pairwise"`,
  `"revpairwise"`, `"trt.vs.ctrl"` (Dunnett; exact with `adj = "mvt"`),
  `"consec"`, `"poly"` (integer orthogonal polynomials), `"eff"` and
  `"del.eff"`, optionally within each level of `by` factors. Coefficients
  and names are emmeans'; results through `contrastmeans()` equal
  emmeans' `contrast()`.

* New `adj` methods in `predictmeans()` and `contrastmeans()`: `"sidak"`,
  `"scheffe"` and `"mvt"`, as in emmeans (equal to emmeans' p-values for
  sidak and scheffe; within 1e-3 for mvt). `"mvt"` is the single-step
  multivariate t adjustment, the exact counterpart of SAS
  `ADJUST=SIMULATE`; it needs mvtnorm (Suggests), uses a fixed seed so
  results repeat, and leaves the user's random numbers as they were. The
  Scheffe rank is the number of means - 1 for pairs (within each `atvar`
  group) and the rank of `ctrmatrix` in `contrastmeans()`. With `permlist`
  only `"sidak"` of the three is allowed.

* `Kmatrix()`, `predictmeans()` and `covariatemeans()` gain
  `weights = c("equal", "proportional", "outer", "cells", "flat")`, which
  sets how the cells of the factors not in `modelterm` are weighted.
  `"proportional"` uses the observed frequencies, the same for every mean,
  like SAS `LSMEANS / OM`. The definitions are emmeans', and the results
  match emmeans. The default `"equal"` is unchanged.

* `predictmeans(vcov_adjust = TRUE)` and `contrastmeans(vcov_adjust = TRUE)`:
  for `lmer` models fitted by REML, standard errors use the Kenward-Roger
  adjusted covariance (as SAS `DDFM=KR` and emmeans do) to go with the
  Kenward-Roger df. The default keeps the earlier, unadjusted standard errors.
* `predictmeans()` and `covariatemeans()` pass `cov_reduce` and `exclude` on
  to `Kmatrix()`. For example, `cov_reduce = "average"` gives `gam` means
  adjusted for a spatial smooth.

## Bug fixes

* **Gaussian `semireg()` fits moved a little each time their df or
  variance components were computed.** lme4 maps the fitted model's theta
  into the deviance function that the fit keeps; the Satterthwaite df
  (`predictmeans()` on the fit, `anova()`), `varcomp()` and `unconditional
  = TRUE` evaluate that function at nearby points and then reset it to the
  model's theta, which by then was the last of those points. Each call
  moved the variance parameters by about 1-3% (Block's theta 1.66 -> 1.60
  -> 1.56 in a split-plot example) and with them `vcov()` of the fit, so
  the second `predictmeans()` call gave other df than the first. The fit
  now keeps its own copy of theta and the deviance function is reset from
  copies made before it is evaluated; repeated calls give the same results
  (tested). The first call on a fit was right before too.
* Refits inside `permmodels()`, `permlmer()` and the other functions that
  refit a model failed for models fitted inside a function with the formula
  or family given as a variable (e.g. `fit <- function(f, d) lmer(f, data =
  d)`): the variable was looked up where the model was fitted, and not
  found (or, for `family`, the function `family()` was found). The refits
  now use the model's formula and family objects. Also, raw permutations of
  `lmer` models took the variable names from the call (`lmer(f, ...)` gave
  a column `f`).

* `permmodels()` for `lmer` models, type I: the observed table was
  lmerTest's type III table (its `anova()` default) while the permuted
  statistics were sequential F statistics from lme4, so in unbalanced data
  the p-values compared different statistics (e.g. F 3.64 against 3.99 for
  the same data). Now the table is lmerTest's type I table and the permuted
  refits compute the same type I F (Wald tests of the same hypotheses).
  Balanced data are unchanged.

* `R2_glmm()` (and `R2_model()`) for `glmmTMB` models with families
  `nbinom1`, `nbinom2` and `beta` used the wrong observation-level variance
  (it multiplied by `sigma()^2`, which for these families is not a
  variance), and ignored offsets in all families. For example, a Poisson
  model with an offset gave 67% with glmmTMB and 82% with glmer; now both
  give 82%. A beta model gave 0.08%, now 37%. The variances are those of
  the new working-scale code, checked against `simulate()`.

* `permindex()` returns indices for the rows of the data as given. Before,
  with `block` or `group` it sorted the data first and returned indices
  for the sorted copy, which were wrong unless the data were already sorted
  by block and group (as in its examples). For sorted data and the same
  `seed` the permutations are unchanged. Also: duplicated permutations are
  found with a separator between indices (`1, 12` and `11, 2` were taken
  as the same), whole groups of unequal size can no longer change places
  (that mixed rows of different groups; now an error), and `seed` leaves
  the user's random number stream as it was.

* `permmodels()` for `glm`, `lm` and `aov` models permutes all the variables
  of the response together: `cbind(R, F)` (successes, failures) permuted only
  `F`, so the totals changed; a transformed response such as `log(y)` failed
  even with `data`. `cbind(R, N - R)` was already right and gives the same
  results. When the response's variables are not in the data (e.g. a model
  of `log(y)` without `data`), the error now says so, and a model with
  several responses gets a clear error.

* `permlmer()` (and so `permanova.lmer()`): when the tested variance is
  estimated at 0, the likelihood-ratio statistic is 0 up to optimizer noise
  (e.g. 3e-14). Permuted statistics of 0 then counted as smaller than the
  observed one, giving a small p-value instead of 1. Statistics are now
  rounded to 4 decimals before they are compared. In a simulation study
  (`dev/sim-study/`) this affected 0.1-0.8% of data sets under the null.

* `predictmeans()` and `covariatemeans()` ignored `trans` for `glm` and
  `glmmTMB` models with an identity link: a gaussian `glmmTMB(log(y) ~ ...)`
  with `trans = exp` showed no back-transformed means. `trans` is now used
  whenever the link is the identity; models with another link still use its
  inverse.
* In the back-transformed plots of `predictmeans()` and `covariatemeans()`,
  the data points were not shifted by `transOff` when the response was
  transformed before modelling (e.g. a column holding `log(y + 1)`, with
  `transOff = 1`), so they sat `transOff` above the means.
* `R2_glmm()`:
  * for `glmmTMB` models with correlated random effects, e.g.
    `(Days | Subject)`, the random-effect variance was wrong (84.32% instead
    of the 78.71% that `lmer` gives for `sleepstudy`) and the slope was
    missing from the breakdown. glmmTMB puts correlations on its own scale,
    so the covariance is now built from `VarCorr()`;
  * the breakdown by term split a correlated term's columns into two
    halves instead of by component (`lmer` too). Each component (intercept,
    slope) now gets its own share; the covariance between them is only in
    the random-effect total. Models without correlated terms are unchanged;
  * failed for `glmmTMB` models with random effects in the zero-inflation
    part, or with coefficients dropped as aliased;
  * `over_disp = TRUE` finds the observation-level random effect wherever
    it is in the formula (glmmTMB keeps the formula's order).
* `varcomp()` for `glmer` models:
  * the SEs of the variances were wrong: converting from SDs used the
    derivative 0.5/sqrt(sd) instead of 0.5/sd (0.2845 instead of 0.2280 for
    the `cbpp` example; glmmTMB gives 0.2294), and covariances ignored their
    dependence on the two SDs. `vcov_glmerMod(ranpar = "var")` now applies
    the delta method with the full Jacobian;
  * failed for models with more than one variance parameter, e.g.
    `(Days | Subject)` ("invalid 'row.names' length").
* The `letterCI` help says when CI overlap matches the pairwise tests and
  when it finds fewer or more differences, with a reference (Maghsoodloo and
  Huang, 2010; from the 1.1.2 help).
* Clearer messages for two common mistakes (ideas from James Curran's
  `isValidTerm()` and `adj` check on the `master` branch):
  * a `modelterm` with its factors in another order or letter case:
    `'Variety:nitro' is not a term of the model. Did you mean
    'nitro:Variety'?`; otherwise the message lists the model's terms;
  * a misspelled `adj` lists the valid choices (before, `p.adjust()` gave an
    error that did not mention `"tukey"`). `adj` may be given in any letter
    case, e.g. `"Tukey"` or `"bh"`, in `predictmeans()` and
    `contrastmeans()`. The `contrastmeans()` help no longer lists `"tukey"`,
    which it never supported.
* `residplot()`:
  * failed for `glmmTMB` models with coefficients dropped as aliased;
  * in the "random effects vs fitted" panel of `glmer` and `glmmTMB`
    models with a correlated term, e.g. `(Days | Subject)`, the random
    effects were matched to the wrong columns of Z (`ranef()` lists them by
    component, Z by level).
* `ci_mcp()` stopped when a lower and upper limit were equal (a mean with
  SE 0).
* The internal `model.matrix()` method for `lme` models accepts `data =`.
* `varcomp()` for `glmmTMB` models gave wrong confidence limits (the upper
  limit could be below the estimate), read correlation parameters as log
  standard deviations, and had no residual variance. It now gives the same
  table as for `lmer`: the variances and covariances of each term
  (zero-inflation terms prefixed `zi.`), the residual variance of gaussian
  models, delta-method SEs (equal to the `lmer` ones for `sleepstudy`), and
  Wald CIs on the log scale for variances and the Fisher z scale for
  correlations. `ci = FALSE` now works for these models.
* `semireg_tmb()` with `predict_info = TRUE`: the random-effect covariance
  was wrong for correlated terms (same cause as in `R2_glmm()`). A singular
  covariance is detected without computing all eigenvalues, and a system
  that is not positive definite falls back to `nearPD()` instead of failing.
* `contrastmeans()` failed with "non-conformable arguments" for every model
  with a missing treatment combination, even for contrasts that do not
  involve it. It now works, and stops with a clear message only when a
  contrast uses a mean that cannot be estimated.
* Tukey-adjusted p-values within `atvar` levels used the number of pairs
  instead of the number of means. With 4 means per level, Tukey was computed
  for 6 means, so p-values were too large (0.350 instead of 0.197 in the Oats
  example). For `lmer` models they also used one averaged df instead of each
  pair's df. Both now match emmeans.
* `predictmeans()` failed for every `lmer` model with a missing treatment
  combination ("replacement has 12 rows, data has 11").
* `predictmeans()` failed for `lmer` models fitted by ML (`REML = FALSE`).
  Satterthwaite df are now used for such fits, with a message; `df_term()`
  falls back to Satterthwaite with a warning instead of stopping.
* `permmodels()` and `CookD()` refitted `lme`/`gls` models by re-evaluating
  the model call from inside the package. Data fitted inside a function were
  not found, and a data set with the same name elsewhere on the search path
  could be used silently. For example, a local `Oats` with `nitro` as a factor
  was replaced by `nlme::Oats`, so the permutation tests were for a different
  model. Refits now use the model's own data and check that they reproduce
  its coefficients.
* `Kmatrix()`, `predictmeans()`, `permmodels()` and `CookD()` now find the
  data of `gls` models fitted inside a function; `nlme::getData()` does not.
* `CookD()` failed for every `lmer` model: it refitted with ML and then asked
  for Kenward-Roger df, which needs REML. It now works and agrees with `lme`.
* `permmodels()` failed for `gls` models without a correlation structure.
* `R2_glmm()`:
  * for `glmmTMB` models it now finds data that are not in the global
    environment;
  * `over_disp = TRUE` gives a clear message when the model has no
    observation-level random effect.
* `predictmeans()`:
  * passing two or more of the deprecated `mplot`, `barplot`, `pplot`,
    `bkplot` arguments gave an error instead of a deprecation warning;
  * `trans = NULL` gave an error.
* `contrastmeans()` with `permlist` counted permuted statistics strictly
  greater than the observed one. It now counts "greater or equal" (after
  rounding), as `predictmeans()` does, so both give the same permutation
  p-value for the same contrast.
* `predictmeans()`, `contrastmeans()` and `permmodels()` no longer leave
  `options(scipen = 6)` set after they return; printed results look the same.
* `contrastmeans()` without `ctrmatrix` only opens the contrast editor in an
  interactive session.
* `semipred()` failed ("replacement has 0 rows") on `semireg()` fits whose
  data had a numeric variable not in the prediction grid, such as a numeric
  grouping id in `(1 | idnum)`, because it tried to undo the centring for
  it in the grid too.

## Speed

* Kenward-Roger and Satterthwaite df of `gls` and `lme` fits (the default
  for these models) and `varcomp()` of them: the observed information is
  computed from its formula instead of numerical second derivatives of the
  log-likelihood (numDeriv's Hessian, 1 + 4k(k + 1) likelihood evaluations
  for k variance parameters), and the derivative of the coefficients'
  covariance exactly. Only the pieces of the covariance (`corMatrix()`,
  `varWeights()`, `pdMatrix()`) are differenced. Same df (to 1e-7) on 11
  gls/lme models with every kind of structure, REML and ML; an unstructured
  `gls` with 6 times, 2000 observations and 21 variance parameters: 13.5 s
  -> 2.7 s; small models 3-10 times faster.

* The Satterthwaite df of the curves of `semireg()` fits (`semipred()`,
  `predictmeans()`): central differences instead of Richardson
  extrapolation, and only the block of the inverse mixed-model matrix that
  the curves use (the fixed and smooth coefficients, not the subjects'
  random effects), and the Hessian of the variance parameters computed
  once per fit (cached). On the growth curves of 2257 rows with 130
  subject curves: the difference of two curves with a simultaneous band
  46 s -> 6 s, slopes 23 s -> 3 s. `semipred(contr = )` no longer computes
  the df of the curves it does not report.

* Satterthwaite and Kenward-Roger df for `lme` and `gls` models build the
  covariance of the data as a sparse matrix, block diagonal by subject, so the
  time grows about linearly with the data: an unstructured `gls` with 2000
  observations takes 3 s (was about 70 s), an `lme` random-intercept model
  with 8000 observations under 1 s. Results are unchanged (to 1e-9, the
  accuracy of the numerical derivatives).

* Satterthwaite df (used for `lmer` models fitted by ML, and by
  `df_term(type = "Satterthwaite")`) are computed for all contrasts at once
  from ingredients computed once per model. They equal lmerTest's values,
  about 25x faster for many contrasts (780 pairs: 0.8 s -> 0.03 s), and no
  longer convert the model to lmerTest's class.

* Kenward-Roger computations no longer form or invert the n x n covariance
  of the data. They use the model's low-rank structure (Woodbury identity),
  with the same results as pbkrtest. n = 3,000 takes under 0.05 s instead of
  21.5 s; n = 20,000 takes 0.3 s (previously out of reach).

* Kenward-Roger df are computed for all contrasts at once (identical to
  `pbkrtest::Lb_ddf()` row by row), and the adjusted covariance is computed
  once per `predictmeans()` call instead of once per model term, for the
  means, and again for the pairs. `predictmeans()` on `lmer` models is about
  10-50x faster: Oats split-plot 0.96 s -> 0.09 s; 40 means with 780 pairwise
  comparisons 5.3 s -> 0.1 s.
* `permmodels()` no longer computes Kenward-Roger df for every permutation
  refit (the df were discarded). About 10x faster for `lmer` models
  (99 permutations of the Oats split-plot: 21.8 s -> 2.2 s), with identical
  permutation p-values.
* `permmodels()` and `permlmer()` use a sparse Cholesky factor and
  triangular solves instead of a dense n x n inverse.
* Letter groups (`multcompLetters()`, used for all letters) are much faster
  for many means, with identical letters (checked on 375 random cases):
  100 means 4.1 s -> 0.3 s; 200 means 103 s -> 2.1 s. The "absorb" step now
  only checks the letter columns that changed.
* `permmodels()` on `lmer` models takes the t values of each refit from the
  covariance matrix it already has, instead of calling `summary()`, which
  computed it again: about 15-20% faster (300 permutations of the Oats
  split-plot, type III: 4.3 s -> 3.5 s), with identical results. The nine
  copies of the per-refit code are now one internal function.
* `R2_glmm()` computes its traces without n x n matrices. For n = 4,000 it
  goes from 4.4 s and 529 MB to 0.01 s with identical results; n = 100,000
  now takes 0.02 s.

## Non-estimable means

* Means that cannot be estimated (they depend on aliased coefficients) are
  now NA in every table, with a message naming them; their pairwise
  differences are still shown when estimable. Before, a number that
  depended on the arbitrary choice of aliased coefficient was printed.
  Means over missing cells keep the 1.x rule (average of the observed
  cells, with the "Missing treatments' combination" warning).

## Containment df

* New argument `df_method` in `predictmeans()` and `contrastmeans()` (and
  `type = "containment"` in `df_term()`) for lmer models: `"auto"` (the
  default, unchanged: Kenward-Roger for REML fits, Satterthwaite for ML fits),
  `"Kenward-Roger"`, `"Satterthwaite"` or `"containment"`. Containment gives
  the SAS `DDFM=CONTAIN` df (the smallest rank contribution to [X Z] of the
  random effects containing the term, else the residual df), with the
  unadjusted covariance unless `vcov_adjust = TRUE` is given. For the oats
  split plot this is 10 for Variety and 45 for Nitro and Variety:Nitro; with
  a whole plot missing, 9 and 42. The ranks use M'M, so a split plot with
  n = 6000 and 750 whole plots takes about 3 s.

## Bug fixes found by the refactor

* `predictmeans()` with both `permlist` and `atvar` failed with "undefined
  columns selected". It now gives permutation p-values within each `atvar`
  level.
* Factor levels containing `":"` (e.g. `"Golden:Rain"`) broke the legacy
  pairwise tables: `atvar` failed and the within-factor SEDs and LSDs were
  computed from misread labels (B4). Level labels are now taken from the
  label table, never split from row names.
* A user-supplied `Df` for an lmer model gave NA p-values for all pairs but
  the first ("missing value where TRUE/FALSE needed" when letters were
  computed). The given `Df` is now used for every pair.

## Documentation

* Letter groups, overall and within `atvar` groups, now come from one
  function (`pm_letter_groups()`) that handles unavailable p-values
  (missing levels or non-estimable differences) the same way everywhere.
  Fixes repeated errors in the letters with missing p-values, and an
  error in the mean/SE tables when non-estimable means meet empty cells.

* Pairwise differences are now checked for estimability, not only the
  means: with confounded factors (e.g. one volume used only with fine
  droplets) a difference that depends on a dropped coefficient is NA
  instead of being computed with that coefficient taken as 0. Letters are
  given to the means whose comparisons can be estimated, and the message
  for non-estimable means suggests `weights = "cells"`. Fixes an error
  ("missing value where TRUE/FALSE needed") in the letters.

* `lmerTest::ranova()` on a gaussian `semireg()` fit failed ("object
  'grp_geo' not found"), and its refits would have treated the remaining
  smooth terms as ordinary random intercepts. Gaussian `semireg()` fits now
  have class `semerLmerTest`, whose `update()` refits through `semireg()`
  with the same smooth bases. `smoothZ` may also be a list given by name.
* New data set `copper` (used in the overview article).

* The deprecated plot arguments (`mplot`, `barplot`/`bar`, `pplot`,
  `bkplot`) turned off every plot not given, so `bar = TRUE` gave a NULL
  mean plot. Plots not given now keep their `which.plots` setting, and the
  deprecation warning shows the equivalent `which.plots` call.

* Fixed (regression in this development version): with nested treatments
  and one observation per block cell (e.g. `treat/contcomp` with
  `repl`), means were averaged over structurally empty cells and were
  wrong. Empty cells are now judged on the factors that interact only.
* `atvar` with several factors failed ("Names required") when a group of
  `atvar` levels had only one or no observed level of the other factors,
  or ("replacement has 8 rows, data has 12") when a group had one level.
* Within-`atvar` comparisons: adjusted p-values could be assigned to the
  wrong group or pair when the groups were not in row order (for example
  two `atvar` factors, or empty groups). This came from 1.x. P-values and
  letters are now matched to groups and levels by name.

* The examples of `permmodels()`, `permlmer()`, `permanova.lmer()`,
  `semireg()`, `semireg_tmb()` and `semipred()`, and the `lmer` example of
  `residplot()`, run again (they were commented out). Permutation examples
  use 99 permutations on one core; slow parts are in `\donttest{}`, and
  examples with HRW data run only when HRW is installed. The `semireg_tmb()`
  example no longer passes `resp_scale`, which only `semireg()` has; the
  `sleepstudy` plot uses ggplot2 instead of latticeExtra.
* For AI assistants and scripted use, `inst/llms-context.md` describes the
  result tables (every column of `means`, `pairs` and `settings`), the
  defaults and common recipes (a test checks that it lists every column).
  The README, which pkgdown turns into the site's home page and `llms.txt`,
  now describes the package, its main functions and its vignettes.
* New vignette "Missing cells and non-estimable means": an empty cell in
  the oats split plot (the 1.x averaging rule and its warning, comparisons
  within the observed cells) and confounded factors in a spray trial (NA
  means with estimable differences, `weights = "cells"`, and writing the
  nesting into the model).
* New vignette "Smooth curves in mixed models": `semireg()` and
  `semipred()` on `nlme::Dialyzer`, `ranova()` of the smooth terms,
  pointwise and simultaneous bands, a difference curve, and `predictmeans()`
  on the fit. `semipred()`'s plot title says when the band is simultaneous.
* New vignette "Permutation and bootstrap tests": `permmodels()` and its
  use as `permlist` in `predictmeans()` and `contrastmeans()`, `permlmer()`
  for a random slope, `permanova.lmer()`, and `bootmodels()` for a binomial
  GLMM, with advice on which to use.
* New vignette "Link and response scale": binomial and Poisson GLMs and
  GLMMs and a log-transformed response (`trans = exp`), with response-scale
  means and SEs, odds ratios and ratios, what a back-transformed GLMM mean
  means, and Piepho's vs Nakagawa's R² (`R2_model()`).
* New vignette "Contrasts, joint tests and slopes": the `adj` methods
  (including `"mvt"`), `contrast_family()` (Dunnett, trends), `pm_joint()`
  (slices), `pm_slopes()` and `coef()`/`vcov()` with marginaleffects, on the
  oats split plot. The getting-started and migrating vignettes link to it.
* Three vignettes: "Getting started with predictmeans", "Adjusted R² for
  generalized linear mixed models" and "Migrating from predictmeans 1.x".
* A pkgdown website (built by GitHub Actions and published to GitHub
  Pages), with two longer worked-example articles: an overview of the
  package and experimental design and analysis for tree improvement.
* `semireg_tmb()` works with glmmTMB >= 1.1.12, which no longer stores
  `Lambdat`, `Lind` and `lower` (they are rebuilt when missing).
* `R2_glmm()` failed for gaussian glmmTMB fits (including `semireg_tmb()`)
  with glmmTMB >= 1.1.12, whose dispersion parameter is named `betadisp`;
  the fixed effects now come from `fixef(model)$cond`.
* New data sets `Body_length` and `Poisson_ind` (used in the R² vignette).

## Internal

* `smZ()`: the code that puts back rows with a missing covariate (7
  copies) and splits a basis by the levels of `by` (4 copies) is one helper
  each (`pm_smz_na()`, `pm_smz_by()`); the results of 114 argument
  combinations (every type, `by`, `group`, missing values, knots) are
  identical to before.

* The permutation and bootstrap tests share one engine (`R/perm-engine.R`):
  permutation indices, refits (in parallel if asked), failure counts and
  p-values from statistics rounded before they are compared. `permmodels()`,
  `permlmer()` (so `permanova.lmer()`), `bootmodels()` and the permutation
  p-values of pairs and contrasts in `predictmeans()` and `contrastmeans()`
  use it. Results are unchanged (checked with `dev/snapshot-perm.R` on 33
  calls and `dev/snapshot-predictmeans.R`), with three small differences:
  with `seed`, the user's random number stream is left as it was (as
  `bootmodels()` already did); refits that fail with an error are left out
  and counted instead of stopping the function, the count shown in the note
  under the tables and in attribute `"failures"`; `bootmodels()` also uses
  several cores on Windows.

* The confidence-interval plot uses `geom_errorbar(orientation = "y")`
  instead of `geom_errorbarh()`, which ggplot2 4.0 deprecates; the plot is
  unchanged.

* `1:length(x)`, `1:nrow(x)` and `1:ncol(x)` are now `seq_along(x)` /
  `seq_len()` throughout (idea from James Curran's code in 1.1.2); results
  are unchanged, and empty inputs no longer loop over `c(1, 0)`.

* Help pages are now written with roxygen2: each page's source is the
  block above its function (data sets in `R/data.R`, the package page in
  `R/predictmeans-package.R`), and `man/` is generated with
  `roxygen2::roxygenise(roclets = "rd")`. NAMESPACE is still maintained by
  hand. The rendered pages are unchanged.

* `predictmeans()` is split into stages (setup, estimate, compare, mean table,
  back-transform, plots, assemble). All statistics are computed before any
  plotting, so `plot = FALSE` skips all graphics work. The work of each
  stage is done by functions that take their inputs as arguments and return
  their results, so they can be tested and changed one at a time. Results
  are unchanged.

## Packaging

* James Curran is credited as a contributor (DESCRIPTION now uses `Authors@R`).
* `which.plots` is documented. DESCRIPTION no longer carries CRAN-generated
  fields, and `NeedsCompilation` is "no".
* Each pull request is timed against its base branch on the same runner
  (`.github/workflows/benchmarks.yaml`, `dev/benchmarks-ci.R`): it fails when
  a key timing (means and comparisons of a 1,200-row lmer, empty cells,
  Satterthwaite df, contrasts, aov with missing values, permutation tests,
  residual plots) is more than 1.5 times slower and at least 0.1 s slower.
* GitHub Actions R CMD check on Linux, macOS and Windows, with `--as-cran`,
  and one job that runs the tests as CRAN does (without `NOT_CRAN`: 140 s).
  Release preparation (no release yet): `cran-comments.md` (draft), and a
  reverse dependency check run by hand (`.github/workflows/revdep-check.yaml`,
  `dev/revdep-check.R`: each CRAN package that uses predictmeans checked with
  the CRAN version and with this one).

## Kmatrix() and semipred()

* `Kmatrix()` now builds predicted means from a reference grid evaluated with
  the model's own basis information. It handles covariates that:
  * enter through transformations: `poly()`, `ns()`, `bs()`, `log()`,
    `scale()`, `I()`;
  * appear only inside interactions (`A:x` without `x`, `x:z`);
  * are mgcv smooth terms (`s()`, `te()`, `s(x, by = A)`).
  These cases previously stopped with "subscript out of bounds" or
  "object not found".
* `Kmatrix()` and `predictmeans()` now support `mgcv::gam()`/`bam()` models,
  using the lpmatrix of the smooths. Random-effect smooths (`bs = "re"`,
  `bs = "fs"`) are excluded by default (`exclude` argument).
* New `Kmatrix()` argument `cov_reduce`: fix covariates with any summary
  function (default `mean`, e.g. `median`), or use `"average"` to average the
  predictions over the observed covariate values.
* `Kmatrix()` now also handles:
  * character and logical predictors;
  * `factor(x)` terms in the formula (`modelterm` may be `"x"` or
    `"factor(x)"`);
  * models fitted without a `data` argument;
  * a transformed response together with `data =`.
* `Kmatrix()` returns the reference grid (`refgrid`) and flags non-estimable
  rows in `attr(K, "estimable")`.
* A matrix `covariateV` for several covariates is now used as a set of joint
  prediction points, as documented; previously it was ignored.
* `semipred()` uses `Kmatrix()` on the `semireg` object directly, without
  refitting an auxiliary `lm()`. The smooth basis is averaged over the same
  reference grid as the fixed part. As a result, factors not in `modelterm`,
  including `by` factors of the smooths, are averaged correctly instead of
  causing an error.
* The covariate summary printed by `prtnum = TRUE` now goes through
  `message()`, so it can be suppressed.
* Kmatrix() is faster on large data (about 4-15x on a 20,000-row, 360-cell
  model).
* New tests check that results for every case supported by 1.1.1 are
  unchanged, and compare the new cases with emmeans and mgcv.

## Fewer dependencies

* `pbkrtest` moved from Imports to Suggests. Kenward-Roger adjusted
  covariances, df and F tests are computed by a port of pbkrtest's
  `vcovAdj()` and `KRmodcomp()` that gives the same results in the tests.
  pbkrtest alone brought 11 tidyverse packages (dplyr, broom, tidyr, ...).
  Installing predictmeans now pulls in 56 packages (112 for 1.1.1).

* `car` moved from Imports to Suggests. `permmodels()` computes its Type II
  and III tests itself: lm F; glm LR, Wald and F; lme and gls Wald
  chi-square; lmer Wald chi-square and Kenward-Roger F. They are ported from
  car and equal `car::Anova()` in the tests. Installing predictmeans now
  pulls in 68 packages (112 for 1.1.1).
* `permmodels()` with `type = "II"` now also works for `lm` models with
  aliased coefficients. 1.1.1 stopped with an error from `car::Anova()`.

* `plyr` is no longer used.
* `HRW`, `lmeInfo`, `lmeSplines`, `plotly` and `splines2` moved from Imports
  to Suggests. They are only needed for, respectively:
  * HRW: example data;
  * lmeInfo: `varcomp()` on `lme` models;
  * lmeSplines: `smZ(type = "smspline")`;
  * plotly: `semipred(threeD = TRUE)`;
  * splines2: the splines2 spline types in `smZ()`.
  A clear message asks for the package when such a feature is used.
  Installing predictmeans now pulls in 76 packages instead of 112. This also
  avoids current HRW's requirement of R >= 4.5.
* `HRW::pointsInPoly()` is replaced by an internal copy. It also closes open
  polygons correctly; HRW's version adds a spurious vertex when the last
  vertex does not repeat the first.
