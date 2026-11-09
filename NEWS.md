# predictmeans 1.1.1.9000 (development version)

## Changes that affect results or scripts

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
* `predictmeans(letterCI = TRUE)` now works with `atvar`: the CI-overlap
  letters compare the means within each `atvar` group (as in predictmeans
  1.1.2). Before, `atvar` was dropped and the letters compared all means.
  The note under the mean table now says the letters come from the overlap
  of, e.g., 83.4% CIs, and gives the approximate test level (about 0.05),
  instead of "significant level '0.166'". `level` keeps its meaning (it sets
  the CIs). With two `atvar` factors the note is one line
  ("at each level of A:B") instead of two.

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
  methods).

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
* GitHub Actions R CMD check on Linux, macOS and Windows.

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
