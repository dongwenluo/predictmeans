# predictmeans 2.0.0 development log

## In short

predictmeans 2.0 keeps the 1.x functions and printed output, and changes
defaults where simulations or the literature showed a more reliable method.
Plots are drawn from results with `autoplot()` or `plot()`, and
lmerTest is no longer needed (`anova()` of `lmer` fits has F tests with
Kenward-Roger df from the package's own code).
Results checked against SAS PROC MIXED, Genstat and emmeans are in the tests.

* **`covariatemeans()` and `semipred()` plots come from `autoplot()`**
  (`?autoplot.pm_curves`, new `R/autoplot-curves.R`), built from the result
  tables and the data rows kept in `attr(x, "plot_data")`: `type = "lines"`,
  `"surface"`, `"3d"` (one plotly widget, a scene per level, instead of a list),
  `"difference"`; `group`/`panel` for any number of factors (two panel factors:
  a grid); `scale`, `data`, `jitter`, `interval`, `title`, `style`. The 1.x look is
  kept. `plot = TRUE` prints it, `$plt` holds it, and `plot_args` (superseded)
  maps to these arguments. Numbers unchanged (tables compared with tolerance 0).
  The curve and surface plots map `y` (and `z`) in each layer, not for the
  whole plot, so layers with other data can be added (e.g.
  `+ geom_rug(data = mydata)`, which only has the covariate).
* `covariatemeans()$pltdf` is now `$pred_df`, as in `semipred()`; `x$pltdf` and
  `x[["pltdf"]]` still give it (a message once per session). With `as_is = TRUE`
  `pred_df` no longer has the helper columns `min_value`, `max_value`.
* `semipred(sm_term = )` accepts a term as written (`"s(pressure, by = QB)"`),
  a smooth as labelled (`"s(pressure):QB200"`) or a covariate (`"pressure"`)
  besides block names; an unknown name lists the choices. New advice
  `smooth_other_covariate` when a smooth in the prediction does not involve
  `covariate` (it adds a constant). `settings$sm_term` records the blocks.
* Fixed: `covariatemeans()` with four or more factors (no plot was possible,
  the call stopped); `plotord` check (`unique()` in `if`); `as_is` rows chosen by
  position, not by row name; two covariates need a matrix `at` (clear error).
  `semipred()` with three or more factors drew no plot; the 2-D `contr` flag of a
  generalized fit now also checks `deriv == 0` (as the 1-D one). The joint F test of `covariatemeans(contr = )` is not added to
  `semipred()`: Wald tests of penalised smooth coefficients are not reliable
  (as `?semipred` says).
* Tests: `R2_model()` reproduces Table 2 of Piepho (2023) (normal, Poisson with
  `over_disp = TRUE`, binary) to the published two decimals (`test-r2-model.R`).
* `semipred(contr = )` reports whether the two curves differ anywhere along the
  covariate: a max-|t| test (`settings$test`, printed and in the plot subtitle),
  simulated from the covariance of the fixed and smooth coefficients; its p-value
  is the smallest level at which the simultaneous band leaves out zero. A Wald F
  test of the shrunken smooth coefficients is not used (it would not have an F
  distribution). In 200 simulated data sets with equal curves it rejected 3.0%
  at the 5% level.
* `semipred()` back-transforms its curves with the same rule as
  `covariatemeans()`: `trans`, else the inverse link, else the inverse of a
  response transformed in the formula (`log(y)` now gives curves on the y scale;
  before it stayed on the log scale unless `trans` was given). New arguments
  `bias_adjust` and `sigma` give the mean of the response; the smooths' own
  random effects are not counted as variation to average over. Results of
  existing calls are unchanged (checked at tolerance 0).
* `CookD()` draws its plot with ggplot2, as `residplot()`; the distances it
  returns are unchanged.
* `autoplot()` of `covariatemeans()` / `semipred()` results uses the plot
  settings of the call (`plot_args`, kept in `attr(x, "plot_args")`) as its
  defaults: `trellis = FALSE` and `ci = FALSE` were ignored. One legend for the
  curves, bands and data, complete and in the data's level order (with many
  groups the colours repeated and the keys were split over two legends);
  distinct hues past 8 groups, a compact legend past 12, no legend when each
  panel is one group; axis labels that would overlap are dropped.
* `semireg()` / `semireg_tmb()` fits print a short account: the formula as
  written, a table of the smooth terms (smooth, random-effect term, covariate,
  basis, basis dimension, smoothing variance, straight line or not) with how to
  pass them to `semipred(sm_term = )`, the fixed effects, the other variance
  components and advice. `summary()` adds the summary of the mixed model
  (`$model`); `tidy(fit)` gives the smooth table, `tidy(fit, "fixed")` the fixed
  effects. `coef()`, `fixef()`, `vcov()`, `fitted()`, `residuals()`, `logLik()`,
  `AIC()`, `nobs()`, `family()`, `VarCorr()`, `ranef()` and `formula()` work on the
  fit. The fit keeps its `call`.
* Difference plots of `covariatemeans()` and `semipred()` (`contr`, `deriv`)
  again mark where the band leaves out 0 (or 1): on the log axis of ratio and
  odds-ratio plots the marks were dropped (y = -Inf). They are now a rug along
  the x axis.
* `residplot()`: `gam` fits show fitted vs observed as the fourth plot instead
  of Cook's distances (only approximate for a penalized fit).
* `residplot()`: the fitted-vs-observed plot takes the response from `pm_data()`
  (rows used in the fit) and labels its axis with the response as modelled
  (`log(yield + 1)`, `Observed proportion of incidence`) instead of "Observed
  values"; a two-level factor response of glmmTMB now gets the calibration plot
  (was "response not available"); points are labelled with their row of the data for every model class (before, lme, gam and glmmTMB count models renumbered the rows when rows with missing values were dropped).
* `permmodels()` says when it does not use the `test_statistic` asked for
  (`method = "raw"`: gls/lme use Chisq, Type I uses F, glm LR/Wald/F;
  Freedman-Lane permutes F or t). Before, it was replaced silently.
* **Documentation reviewed** (2026-11-23): help pages, vignettes, README and
  DESCRIPTION checked against the 2.0 API (snake_case arguments, `autoplot()`,
  `pairwise`); each page states its defaults, df method and scale, links to
  `advice()`; `?predictmeans-package`, README and `vignette("predictmeans")`
  say what the package aims for and which tests check it against SAS, Genstat,
  pbkrtest, nlme and emmeans; `vignette("sas")` lists every reference run.
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

* **`se_ranef()` rewritten** (2026-11-23): a data frame (class `pm_ranef`; `group`, `level`, `term`,
  `estimate` (BLUP), `std.error`, `conf.low`, `conf.high`) with `print()`, `tidy()`, `as.data.frame()`; `as.list()`
  gives the 1.x list of SE matrices. New default `type = "prediction"`: the prediction error SE
  sqrt(Var(u_hat - u)) from the inverse of the mixed model equations (SAS "Std Err Pred"), larger than 1.x's
  (Oats `(1|Block/Variety)`: Block 8.34 vs 6.37, Block:Variety 7.78 vs 7.17); `type = "conditional"` gives the 1.x
  numbers (lme4 condVar). Now also `lme` (any correlation/variance structure) and `semireg` fits; gls/lm/aov stop
  clearly. `predict = TRUE` adds fixed part + BLUP with its SE (and the response scale for GLMMs). One joint
  covariance for all classes: `pm_ranef_cov()` (`R/ranef-cov.R`). Advice `variance_at_zero`, `few_re_levels`,
  new `ranef_conditional_fixed`. Tests `test-se-ranef.R`; SAS run to do: `dev/sas/se-ranef.sas`.

* **`lsd_varies` advice only when it matters** (2026-11-23): it now appears only
  when the average LSD (the printed one and the plot's bar) would judge at least
  one pair differently from that pair's own test, and says how many ("2 of 66
  pairs ..."); unequal LSDs alone no longer trigger it.

* **Results keep the observations, not the model** (2026-11-23): for
  `autoplot(data = TRUE)` a result keeps `attr(x, "plot_data")`, the factors of
  `modelterm` and the response (`.y`, and `.y_orig` when it is transformed), not
  the fitted model as before. A saved result is about the size of the means
  plus the data (Oats lmer: 5 kB, was 97 kB with the model; the model is
  92 kB). Without it, `autoplot(x, data = df)` takes the data.

* **Argument names are snake_case** (2026-11-23, step 6 of the predictmeans()
  review): `conf.level` is `conf_level` (predictmeans, contrastmeans,
  covariatemeans, pm_slopes, semipred, permstats, PMplot, varcomp),
  `letterCI` is `letter_ci`, `as.is` is `as_is` (covariatemeans),
  `test.statistic` is `test_statistic` (permmodels), `max.degree` is
  `max_degree` (contrast_family, which now checks its arguments) and
  `drop.unused.levels` is `drop_unused_levels` (semireg). The dotted names stop
  with the new spelling and the value given. Result columns keep the broom
  names (`conf.low`, `std.error`, `p.value`); arguments passed on to other
  packages keep theirs (`na.action`, `REML`).
* Fixed in this review: the back-transformed columns of `mean_table` (joined
  by level since step 1) failed when the back-transformed table lacks a factor
  column (a GLMM with a treatment constant within clusters); and a `gam` fitted
  inside a function found a data set of the same name in the global
  environment (mgcv records the global environment), so `pm_data()` and
  `CookD()` read the wrong rows.

* **`pairwise` only sets what is printed; one LSD only when all are equal**
  (2026-11-23, step 5 of the predictmeans() review). Comparisons are asked
  for by `adj`, `atvar`, `permlist`, `letterCI` or `pairwise = TRUE`; then the
  `pairs` table is made and adjusted whatever the number of means.
  `pairwise` now defaults to `NULL`: the p-value matrices and letters are
  printed when comparisons are asked for (as before); `pairwise = FALSE` hides
  them without changing any result (before, `adj` or `atvar` turned it back
  on). The 1.x "Pairwise LSDs" matrix is no longer printed (each pair's LSD is
  `pairs$lsd`; the LSD summary goes with the p-values). The LSD summary has
  attribute `"equal"`: when all LSDs compared are equal, `summary()` prints one
  value "(equal for all pairs)"; when they differ, the note of `LSD` says no
  single LSD tests every pair, `summary()` prints the range, and new advice
  `lsd_varies` (only when comparisons were asked for) says to use each pair's
  p-value or LSD. Tests: `test-pairwise.R`; snapshot `tidy` accepted.
* **`autoplot(scales = "free_y")`** (2026-11-23): each panel its own y range,
  with its LSD bar from its own lowest point; the caption says to read each bar
  within its panel. `"fixed"` (default) keeps one bar position for all panels.

* **`predictmeans()` draws no plots: `plot` and `plot_args` removed**
  (2026-11-23, step 4 of the predictmeans() review). Plots come from the
  result: `autoplot(pm)` or `plot(pm)` (`type = "means"`, `"bar"`, `"ci"`,
  `"pvalue"`). `plot = TRUE`, `plot = c(...)`, `plot_args` and the 1.x names
  (`mplot`, `plotord`, `plottitle`, `basesz`, ...) stop with the `autoplot()`
  call that replaces them; `plot = FALSE` is still accepted (it asked for
  nothing), and `options(predictmeans.plot)` no longer affects
  `predictmeans()` (it still does `covariatemeans()` and `semipred()`). The
  result has no `predictmeansPlot` and `predictmeansBarPlot` slots;
  `p_valueMatrix` stays (data, read by `PMplot()`). `pm_make_plots()` (about
  370 lines) is removed; `pm_plot_args()` moved to `R/plot-style.R` for
  `covariatemeans()` and `semipred()`. Vignettes (`predictmeans`,
  `migrating`, the overview and tree-improvement articles) use `autoplot()`.
  Tests: `test-plot-args.R`, `test-stages.R` and others updated.

* **`autoplot()` draws every plot of predictmeans(); `plot()` method**
  (2026-11-22, step 3 of the predictmeans() review). Types: `"means"` (points
  and lines, `lines = FALSE` for dots), `"bar"`, `"ci"` (each mean with its
  confidence limits, horizontally: the 1.x "back" plot, on the response scale
  when the means are back-transformed) and `"pvalue"`. `interval = "lsd"` (the
  new default) draws one average LSD bar beside the means of each panel, as
  1.x `lsd_bar`; it is the LSD that goes with the p-values (adjusted with
  `adj`), and with `atvar` each panel is a level of `atvar` with the LSD of
  the pairs compared within it (two `atvar` factors: a grid). The caption gives
  the range when the LSDs differ. `interval = "sed"` draws the average SED;
  `"ci"`, `"comparison"` and `"none"` as before. The half-LSD bars on each mean
  are gone: with equal SEDs they said what one bar says, and they looked like
  confidence limits. `data = TRUE` adds the observations (from `pm_data()`) on
  the scale shown, e.g. `log(y + 1)` on the model scale and `y` on the response
  scale (not on a link scale); `x`, `group`, `panel` place the factors
  (replacing 1.x `plotord`); `scale = "model"/"response"`. `type = "response"`
  is now `type = "ci"` (or `scale = "response"`), with an error saying so.
  `plot(pm)` prints `autoplot(pm)`. Results keep the model as attribute
  `"model"` for `data = TRUE`. The look is that of predictmeans 1.x: the
  LSD/SED bar is blue, labelled, at the lower left of each panel, half a
  category step left of the first level, at the same place in every panel:
  from the lowest point of the plot (the lowest mean, limit or observation;
  0 for bars); one group's means are red points
  joined by black lines, groups have coloured lines with line types; bar
  charts are light blue (groups: shades of blue) with axis lines and no grid,
  legend on top; the "ci" plot has red means and limits and the observations
  as open blue points. Tests: `test-autoplot.R`, helper `pm_layer()`.

* **`pm_data()`: the data of a fit with the response on both scales; plot
  labels and titles say what the means are** (2026-11-22, step 2 of the
  predictmeans() review). `pm_data(model, data, trans)` returns the rows used
  in the fit (subset, missing values) for every class predictmeans handles
  (lm, aov, aovlist, glm, gls, lme, lmer, glmer, glmmTMB, gam, semireg), with
  `.y` (the response as modelled, e.g. `log(y + 1)`; a proportion for
  `cbind()`), `.y_orig` (the original scale: the data's `y`, or the inverse of
  a transformation given as `trans`) and `.weights`, and attribute
  `"response"` (expression, variable, scale, link, what back-transforms the
  means). It stays the adapter generic for new classes (`?pm_coef`); the
  default method used to return NULL. Results have `settings$response`, and
  `autoplot()` labels axes and titles from it (`pm_scale_label()`):
  `log(y + 1)` on the model scale, "Median of y" back-transformed (the mean
  with `bias_adjust`), "Probability of y" for a logit link, "Mean of y,
  typical cluster" for a GLMM, `trans(ly)` for a user's `trans`;
  "back-transformed" appears only when the means are. Before, every plot
  said "Predicted mean" / "Back-transformed mean". Tests: `test-pm-data.R`.

* **`atvar`: checked, two factors, LSDs of the pairs compared, atvar order**
  (2026-11-22, step 1 of the predictmeans() review). `atvar` is checked at
  once: a name not in `modelterm`, or one naming every factor, is an error;
  `"A:B"` works as `c("A", "B")`; `""`/`"NULL"` mean none (the old check
  failed for two names). With `atvar` the printed LSDs (max, min, average)
  are now those of the pairs within its levels, each with its own critical
  value (its df; with `adj` its family's), not of all pairs with the largest
  critical value; attribute `"Within levels of atvar"` of `LSD` gives them per
  level (pairs, SED range, critical value, LSD range). `means`, `coef()`,
  `vcov()` and `pairs` are ordered atvar first (pairs across levels last), as
  `mean_table` was; `meandecr` now sorts within each level of `atvar`. Fixed:
  `mean_table` sorted only by the first `atvar` factor with two of them
  (`eval(parse())` of a vector), and its back-transformed columns were joined
  by row position, so a reordered table had them on the wrong rows. Tests:
  `test-atvar.R`.

* **Advice in outputs: one engine** (2026-11-22, step 1 of the plan "rule
  no. 1, Explained"): `R/advice.R` stores short advice triggered by the fit
  and data as `attr(x, "advice")` (columns id, level, text) and prints it as
  one "Advice:" block; `options(predictmeans.advice = FALSE)` turns printing
  off. No function uses it yet (`test-advice.R`).
* **anova() and pm_check() advice on the engine** (step 2): the advice of
  `anova()` (lmer, semireg) is stored as `attr(, "advice")` with stable ids
  (`anova_unbalanced_typeI`/`II`/`III`, `anova_empty_cells`,
  `anova_covariate_zero`) and printed as an "Advice:" block (`!` caution,
  `i` note) instead of "Note:" lines; `options(predictmeans.advice = FALSE)`
  hides it, `advice = FALSE` drops it. `pm_check()` stores the notes of the
  rows worth a look as advice (`check_<name>`), kept by `summary()`; its
  printout is unchanged.
* **Advice from predictmeans() and contrastmeans()** (step 3): stored as
  `attr(, "advice")` and printed at the end of `print()` and `summary()`
  (`R/advice.R`, thresholds in `pm_advice_limits`). Rules, each triggered
  by the fit and data only: `nonestimable_means`, `empty_cells`,
  `aliased_covariate`, `large_se` (SE over 3 times the median),
  `bias_adjust` (back-transformed means that are a typical cluster's or a
  median, not the mean; mentions 'sigma' when there are random slopes),
  `singular_fit`, `small_df` (df < 5), `many_tests_unadjusted` (6 or more
  requested pairs/contrasts with `adj = "none"`); for contrasts also
  `contrast_common_cells`. `Kmatrix()` returns `level_covariates`.
* **Advice from permmodels()/bootmodels(), R2_model() and varcomp()**
  (step 4): `few_permutations` (fewer than 999 resamples and a p-value at
  the floor or within 2 Monte Carlo SEs of 0.05); `r2_latent_scale`
  (non-gaussian R2 is on the link scale) and `singular_fit` for
  `R2_model()`; `variance_at_zero`, `few_re_levels` (fewer than 5 levels)
  and `variance_test_boundary` (points to `varcomp(test = TRUE)`, lmer,
  glmer and glmmTMB) for `varcomp()`. A `varcomp()` table with advice has
  class `pm_varcomp` so that it prints.
* **Advice from semireg() and semipred()** (step 5): `smooth_linear` (a
  smooth estimated as a straight line, variance 0) and `singular_fit` for
  `semireg()` and `semireg_tmb()` fits, printed by `print()`;
  `extrapolation` (predictions outside the range of the data) and
  `pointwise_regions` (regions of difference read from pointwise limits;
  `band = "simultaneous"` controls the error over the curve) for
  `semipred()`, printed by `print()`.
* **New `advice()` and help page `?advice`** (step 6): `advice(x)` returns
  the advice of a result as a data frame (id, level, text); `?advice` lists
  every id, the functions that give it and when. A test checks that every
  id in the code is on the page.
* **Advice from covariatemeans() and pm_slopes()**: `extrapolation` (lines
  or `at` beyond the data), `pointwise_regions` (`contr` without a
  simultaneous band) and `bias_adjust` for `covariatemeans()`;
  `slopes_equal_by_model` (no covariate x term interaction),
  `slopes_at_value` (curved term or link: slopes depend on `at`),
  `extrapolation`, and the shared `small_df` / `many_tests_unadjusted` /
  `singular_fit` for `pm_slopes()`, printed by `print()`
  (`test-advice-lines.R`).

* **No partial-match warnings from mgcv** (2026-11-22): with
  `options(warnPartialMatchDollar = TRUE)`, `semireg()` with `t2()` and
  `R2_model()` on `gam` fits warned from inside mgcv (`er$value`, `ed$val`).
  `t2()` bases are built with the option off; `R2_model()` takes a gam's
  scale from `model$sig2` (equal to `summary()$dispersion`). Results unchanged.

* **Covariates constant within a factor's levels** (2026-11-22): in a fit
  such as `lm(ht ~ site + sitemean * seedlot)`, where `sitemean` is one value
  per site and so aliased with `site`, every mean from `predictmeans()` was
  `NA` (the covariate was set at its overall mean, which no site has). Each
  grid row now takes its own level's value, with a message; the means equal
  the averages of the fitted values over the cells averaged over. (Version
  1.1.1 set the aliased coefficient to 0, so its means for such fits depended
  on the reference level of the factor and were not estimates.) Checked
  against Williams et al., *Experimental Design and Analysis for Tree
  Improvement*, Table 5.7 (`test-level-covariate.R`). New data set
  `DataExam5.2` (the book's Table 5.2) and a section "Covariates aliased with
  a factor" in `?predictmeans`.

* **Printed `permmodels()` / `bootmodels()` tables** (2026-11-21): model
  p-values below 1e-4 print as `<1e-04` instead of `0` (the stored tables
  still hold the rounded numbers), and the notes give the smallest attainable
  permutation p-value, 1/(nperm + 1) (1/nperm with `exact = TRUE`), so a
  permutation p-value at that floor beside a model p-value near 0 is not read
  as disagreement.
* **`semipred()` / `predictmeans()` on `semireg` fits with `unconditional =
  TRUE`** (2026-11-21): the smoothing-parameter correction is computed once
  per fit and kept in `fit$cache` (`pm_sr_uncond_cached()`), not on every
  call. Fits made before this change have no cache and compute it each time,
  as before. Results are unchanged.
* **`permmodels()` on `lme` fits** no longer warns "partial match of
  'reStruc' to 'reStruct'" (`R/build-variance-matrices.R`).
* **`varcomp(test_method = "bootstrap", ncore = )`** (2026-11-21): on Windows
  the parametric bootstrap ran on one core whatever `ncore` was; it now uses a
  socket cluster, started once for all batches (`pm_perm_mapper()`,
  `R/perm-engine.R`), about twice as fast with 4 cores in a 199-draw example.
  The early-stopping checks are made at the same draws for any `ncore` (a
  fixed batch of 40; before, 10 x `ncore`), so with the per-draw seeds the
  p-values are now identical for every `ncore` and platform (before they
  could differ when the bootstrap stopped early).
* **`permlmer()` with different fixed effects** (2026-11-21): the permuted
  statistic was `logLik(REML = FALSE)` of the REML fits (the ML likelihood at
  the REML estimates), while the table printed the ML likelihood ratio of
  `anova()` (refitted by ML): two statistics (66.43 against 66.33 in the Oats
  example). Both models are now refitted by ML once before permuting, so the
  null model that generates the permuted data, every refit and the observed
  statistic are ML, and the permuted statistic is the printed `Chisq`. Tests
  of random effects (same fixed effects) stay REML.
* **Kenward-Roger with an information matrix that is not positive definite**
  (2026-11-21). The observed (default) and average information of the
  variance parameters need not be positive definite away from a clear
  maximum. For `lmer` fits such a matrix was inverted anyway (only near
  singularity was checked, not the sign); for `lme`/`gls` the observed one
  silently dropped those directions and the average one was inverted as it
  was. Now `pm_kr_info_check()` (`R/kr.R`) checks the eigenvalues for both
  model classes and, when one is not positive, the expected information is
  used with a message; `attr(, "information")` of the adjusted covariance
  says which was used. Fits with a positive definite information are
  unchanged.
* **`bias_adjust = TRUE` for `aov()` with `Error()`** (2026-11-21):
  `predictmeans()` used the residual variance of the fixed-effects `lm` as
  sigma^2, while `covariatemeans()` used the sum of the stratum variance
  components (rounded to 4 decimals). Both now use the unrounded sum of the
  stratum components (equal on balanced data; 0.0442 against 0.0427 in an
  unbalanced Oats example). The check of `bias_adjust` and the choice of
  sigma^2 are made in one place, `pm_bias_setup()` (`R/response-scale.R`).
* **`semireg_tmb()` refits with prior weights or an offset** (2026-11-21):
  the refits behind `anova()` (likelihood ratio tests of the smooths) and
  `varcomp(test = TRUE)` left out `weights` and `offset`, so a weighted fit
  was compared with an unweighted null model. They are now kept (as for
  `semireg()`); a refit with nothing dropped reproduces the fit's
  log-likelihood. Also: Type III Wald tests of a model with an `offset()`
  term failed for glmmTMB fits ("object not found"); the offset is now left
  out of the rebuilt design (`pm_type3_map()`), as it is part of no hypothesis.
* **`predictmeans(trans = )` for models with a link** (2026-11-21): a given
  `trans` was silently replaced by the inverse link. It is now used, with a
  message, and `trans = I` keeps the link scale, as in `covariatemeans()`;
  `settings$back_transform` says the means are the user's `trans` of the
  link-scale mean. Without `trans` nothing changes.
* **O'Sullivan smooths (`s(x, bs = "os")`)** (2026-11-21): `ZOSull()` returned
  a vector for one value of x, so a one-point prediction grid failed; and
  predictions outside the basis range (the data range widened by 1%) silently
  fell back to the linear part of the curve, because the B-splines are 0
  there. One value now gives a one-row matrix; values outside the range stop
  with the range and how to widen it (`xt = list(range = )`). The instability
  check of `ZOSull()` is a `warning()`, not `print()`.
* **Kenward-Roger for `lmer` fits with prior weights or an offset** (fixed
  2026-11-21). The adjusted covariance and df ignored the weights (residual
  covariance taken as sigma^2 I, not sigma^2 W^-1: adjusted SEs smaller than
  unadjusted ones and negative df were possible), and the observed and average
  information left the offset out of the residuals. Now equal to `lme()` with
  `varFixed(~ 1/w)` for weights, and to the fit of `y - offset` for offsets
  (`tests/testthat/test-kr.R`). pbkrtest and lmerTest also ignore lmer
  weights, so they are not the reference there.
* **`semipred()` standard errors of Gamma and inverse-Gaussian `semireg()`
  fits** (lme4 engine; fixed 2026-11-21) were too small by a factor sigma: the
  working weights already hold 1/sigma^2 and `semipred()` multiplied by
  sigma^2 again. `CovMat` is now relative to sigma^2 as for the other families
  (fixed-effect block times sigma^2 equals lme4's `(RX'RX)^-1`); binomial and
  Poisson fits are unchanged.
* **A failed Kenward-Roger adjustment is no longer silent** (2026-11-21).
  Before, `df_method = "auto"` (or `"Kenward-Roger"`) on an `lmer` fit where
  the adjustment failed kept the Kenward-Roger label with unadjusted SEs. Now a
  message gives the reason, and Satterthwaite df with unadjusted SEs are used.
  A Kenward-Roger df that is not positive and finite stops with a clear
  message instead of giving a test with a negative or infinite df.
* **`permlmer()` and `permmodels(method = "raw")` with weighted `lmer` fits**
  (2026-11-21). `permlmer()` decorrelated the residuals as if V = Z Lambda
  Lambda' Z' + I, ignoring the prior weights (the residual part is
  diag(1/w)), so the permuted residuals were not exchangeable; it also left
  the offset out of the fitted values. `permmodels(method = "raw")` refitted
  the model to its formula's variables and stopped ("object 'w' not found").
  Both now whiten with V = Z Lambda Lambda' Z' + diag(1/w) (`pm_mer_resvar()`,
  `R/whiten.R`), keep the offset in the fitted values, and permute through
  `refit()`, which keeps weights, offset and the rows used.
* **Confidence limits of the means with a `permlist`** (2026-11-21) were
  mean +/- 2 SE, labelled 95% whatever `conf.level` was. The permutation tests
  are of comparisons; the limits of each mean are now the model's t limits at
  `conf.level` with the model's df (the same limits as without a `permlist`),
  in `means`, the printed table and the back-transformed limits. `df` stays
  NA in the tables. `summary()` says which limits these are.

  even with random permutations, so p = 0 was possible. It now applies only
  when `perms` holds distinct permutations (the user's complete enumeration);
  otherwise it is ignored with a warning and p = (b + 1)/(B + 1).

* **Arguments removed or renamed in 2.0 (steps 2a, 2b: `predictmeans()`,
  `contrastmeans()`, `pm_slopes()`, `Kmatrix()`, `df_term()`,
  `covariatemeans()`).** A 1.1.1
  spelling stops with the 2.0 one, e.g. `'level' was removed in predictmeans
  2.0: use conf.level = 0.95.` (`R/defunct.R`); abbreviated 1.x names too.
  - `level` (a significance level) -> `conf.level = 1 - level` (default 0.95).
  - `Df` -> `df_method = <number>`.
  - `which.plots`, `mplot`, `barplot`, `pplot`, `bkplot` -> `plot = c("mean",
    "bar", "pval", "back")`; plot settings as arguments of their own
    (`plottitle`, `plotord`, `lsd_bar`, ...) -> `plot_args = list(...)`;
    `prtplt`, `newwd` and `plot_args$responsen` removed.
  - `prtnum`, `ndecimal` removed: nothing is printed during the call; the
    `means`/`pairs` tables are unrounded.
  - `transOff` -> the offset inside `trans` (`trans = function(x) exp(x) - 1`);
    the pairs then have no `ratio` columns (with an offset they were ratios of
    y + 1, labelled as ratios of y). `count` removed (Poisson and negative
    binomial means are detected).
  - `Kmatrix(covariateV =)` -> `at =`; `df_term(type =)` -> `df_method =`.
  - `covariatemeans()`: `level`, `Df`, `covariateV` (-> `at`), `transOff`
    as above; `trellis`, `plotord`, `mtitle`, `ci`, `point`, `jitterv` ->
    `plot_args = list(...)`; `responsen`, `newwd` removed (the response
    variable is read from the formula). `as.is` stays (each line over its
    group's observed range). 15 calls of `dev/snapshot-covariatemeans.R`
    written both ways give identical results, except `transOff` (as above).
  - `semipred()` (step 3a): the first argument is `model` (was `semireg`);
    `covariateV` -> `at`; `level` -> `conf.level`; `simultaneous = TRUE` ->
    `band = "simultaneous"`; `prt` -> `plot`; `trellis`, `scales`, `plotord`,
    `ci`, `point`, `jitterv`, `threeD` -> `plot_args = list(...)` (with the
    colours). 13 calls of `dev/snapshot-semipred.R` written both ways give
    identical results.
  - `semireg()`, `semireg_tmb()` (step 3b; both kept: their options differ by
    engine): `prt` and `predict_info` removed (nothing is printed during the
    fit; the prediction information is always kept). `ngbinomial` removed: it
    took theta from a Poisson fit, which is not the maximum-likelihood fit
    (on a test data set log-likelihood -1133.60 against -1132.98, theta 3.31
    against 2.80, true 3; means 1%, SEs 4% apart); use
    `semireg_tmb(family = glmmTMB::nbinom2)`. Six fits and their curves
    (`dev/snapshot-semipred.R`) are identical (tolerance 0).
  - `smZ()` and `smoothZ` removed (step 3c): smooths are written in the formula
    (`smZ(x)` is `s(x, bs = "os")`, `smZ(cbind(x, z), type = "Ztps")` is
    `s(x, z, bs = "ztps", k = 6)`, ...). `vignette("migrating")` shows each
    1.1.1 form beside its formula form, and the 1.1.1 models of the overview
    article beside their 2.0 code; `dev/check-smz-migration.R` fitted them both
    ways on the last build with `smZ()`: the same REML log-likelihoods and
    curves (to 1e-5). Where they differ 2.0 is right: knots given in the data's
    units were placed on the centred covariate in 1.x (`ncenter = TRUE`); the
    1.1.1 article's formulas for a curve per group left out the group slopes,
    which the formula form adds. `ncenter`, `nscale`, `resp_scale` removed: the
    variables are used as they are (curves unchanged; the coefficients of other
    numeric covariates are no longer centred). `fit$knots_lst`, `range_lst`,
    `type_lst`, `smoothZ_call` removed. splines2 and lmeSplines are no longer
    suggested.
  - Permutation functions (step 4): `permlmer(lmer0, lmer1)` ->
    `permlmer(model0, model1)`; `permglmer()`, `permglmmTMB()` (the same
    function) and `permanova.lmer()` (use `permmodels(model, type = )`,
    Freedman-Lane, other permutation p-values) removed; `permmodels()` and
    `bootmodels()` lose `prt` (their results are class `"pm_permmodels"` and
    print their tables; one print method for the five places that printed them)
    and `permmodels()` loses `fo` (a formula held in a variable is detected);
    `permstats(level =)` -> `conf.level`. 12 calls of `dev/snapshot-perm.R`
    written both ways give identical results.
  - Plots, diagnostics and small functions (step 5): `residplot(level =)` ->
    `re_level` (the random-effect term), `newwd`, `ask` removed;
    `CookD(idn =)` -> `id`, `newwd` removed; `PMplot(level =)` ->
    `conf.level`, `margin`, `legendx`, `newwd` removed; `varcomp(level =)` ->
    `conf.level` (the same value: it was a confidence level). `R2_glmm()`
    removed: `R2_model()` gives the same values as numbers. `ci_mcp()` removed:
    letters from overlapping confidence intervals do not answer the pairwise
    tests; `predictmeans(letterCI = TRUE)` gives letters from comparison
    intervals (its engine is now internal). The overview article's 83.4%
    letters along the dose became letters at a dose from `letterCI`. 12 calls
    of `dev/snapshot-diag.R` written both ways give identical results.
  - The migration guide in full (step 6): `vignette("migrating")` shows the
    1.1.1 and 2.0 code side by side for every changed call (means, covariates
    and curves, permutation tests, diagnostics, smooth models), one table of
    every changed argument and function, the defaults and the numbers that
    change. `dev/check-migrating.R` runs each 2.0 line and checks that each
    1.1.1 line stops with the 2.0 spelling.
  - Release checks (step 7): the migration guide also lists the packages
    `library(predictmeans)` no longer attaches (glmmTMB, lmerTest);
    `cran-comments.md` describes the removed functions and changed
    dependencies.
  - `print()` of a `semireg()` or `semireg_tmb()` fit shows the fitted model
    (as 1.x printed it during the fit with `prt = TRUE`) and the random-effect
    term of each smooth; it printed the whole list, data and matrices
    included. `summary()` gives the summary of the fitted model.
* James Curran (University of Auckland) is an author of the package (he was
  listed as a contributor).
* `permmodels()` keeps to one core for any model whose formula is not
  written in the call: R-devel's `lm()` records more than the name of a
  formula held in a variable, which the check missed. Help pages: package
  anchors for the glmmTMB links, the workflow figure sized in pixels.
* A glmmTMB fit whose design is rank deficient (aliased columns, e.g. from
  combinations of factor levels without data) stopped every means function
  with "logical subscript too long": glmmTMB drops the aliased columns from
  its model matrix but lists them in `fixef()` as NA. The coefficients are
  now matched to the model matrix (`mymodelparm()`); the means equal those
  of the lmer ML fit of the same model.
* The columns of a K matrix are matched to the model's coefficients by name
  in one place (`pm_K_cols()`, `R/mymodelparm.R`, used by `predictmeans()`,
  `contrastmeans()`, `covariatemeans()`, `pm_joint()` and `pm_slopes()`),
  not by position: a class that lists its coefficients in another order, or
  lists aliased ones its design leaves out (glmmTMB), gives the right means;
  names that do not match stop with a message naming the coefficients
  (`Kmatrix()` says the same, with both causes: empty factor levels, or a
  new class whose `pm_coef()` and `pm_basis()` names differ). Replaces the
  glmmTMB-only fix above. `covariatemeans(mean_type = "overall")` used
  `mp$coef[mp$estimable]`, wrong when coefficients are aliased. Test: a
  made-up class with reversed, partly aliased coefficients equals `lm()`;
  mismatched names stop (`test-model-classes.R`). All 107 results of
  `dev/snapshot-*.R` identical to main at tolerance 0.
* Checks on an older stack (R 4.3.3, lme4 1.1-35, ggplot2 3.4.4, glmmTMB
  1.1.8) after merging the local Copilot work (2026-11-24) found two faults
  and two version-dependent tolerances:
  - `autoplot()` with `data = TRUE` and a grouping factor used
    `position_jitterdodge()` with only `group`, which ggplot2 < 3.5 refuses;
    the points now also map `fill` (open circles draw no fill;
    `show.legend = FALSE`, so no extra legend).
  - lme4 1.1-35 `refit()` drops a fit's NA rows again from the new response
    unless it carries an "na.action" attribute, so every permutation,
    bootstrap and variance-component refit of an `lmer`/`glmer` fit with
    missing values failed ("permutation can't produce useful data!"). All
    such refits now go through `pm_mer_refit()` (`R/utils-refit.R`), which
    sets the attribute. Test: `test-permlmer.R`.
  - Tolerances: Piepho's Table 2 binomial GLMM to 0.01 (as the Poisson
    one); `se_ranef(type = "conditional")` against lme4's glmer `condVar`
    to 1e-4 (1e-6 with lme4 2.1-0, 3e-5 with 1.1-35).
* CRAN pre-tests of the 2.0.0 submission (2026-11-24): "Rd files without
  \usage: anova.lmerMod.Rd, anova.semireg.Rd, se_ranef.Rd" (r-devel Debian).
  Their roxygen blocks were detached (`@name x` + `NULL`); now on the
  functions, so roxygen writes `\usage` (`\method{anova}{lmerMod}` etc.).
* CRAN pre-test "Overall checktime 21 min > 10 min" (r-devel-windows; tests
  13 min). Tests timed as CRAN runs them (`NOT_CRAN` unset): 430 s here, 4
  tests took 166 s. The 114 tests over 0.5 s now start with
  `skip_on_cran()` (they still run on GitHub, locally and with `NOT_CRAN`):
  CRAN-mode tests 58 s, 0 failures. Every exported function is still
  called by a test CRAN runs (checked by tracing all exports); two quick
  tests added for the four that were not (`semireg_tmb()`, the removed
  `permanova.lmer()`/`permglmer()`/`permglmmTMB()` stubs).
* Vignettes `resampling` (31 s) and `smooth` (28 s), 74% of the vignette
  build time, are precomputed (maintainer's choice, 2026-11-24): sources in
  `vignettes/*.Rmd.orig`, knitted by `dev/precompute-vignettes.R` into the
  `.Rmd` with results and 16 figures (696 KB; the HTML size is unchanged,
  960 KB). Their build now takes 2.6 s. Rule in CLAUDE.md and AGENTS.md.
  - Checked: 52 calls of `dev/snapshot-predictmeans.R` written both ways give
    identical results (tolerance 0), except `trans = exp, transOff = 1`
    (above; `response.std.error` within 4e-10, from a numerical derivative).

* **`anova()` with empty cells.** Type III tests of `lmer`, `lm` and other
  fits (and `permmodels(type = "III")`) stopped with "there are aliased
  coefficients in the model" when a cell had no data. Each term is now
  tested on the cells with data (`pm_type3_map_cells()`): an interaction by
  its contrasts that use no empty cell (one df fewer per empty cell; for
  `lm` its sum of squares equals the comparison of the models with and
  without it), a main effect whose means average over empty cells on the
  cells its means share, the test of `pm_joint()`; the heading says so and
  the intercept is not tested. With no empty cell the hypotheses are the
  sum-to-zero ones (equal F and df). `glm` Type III likelihood ratio tests
  need a complete design and say so (`test.statistic = "Wald"` works).
  Checked against a SAS run (`dev/sas/pm-reference-emptycells.sas`,
  2026-11-18: PROC GLM contrasts, equal to SAS Type IV here, and PROC MIXED
  KR and Satterthwaite agree to the printed digits) and Genstat's REML Wald
  test of the interaction (run of 2026-11-14).
* **Means over empty cells: nothing read off them contradicts their tests.**
  Their pairs are compared on the cells both means have, so the difference
  of the shown means is not the difference tested. Before, `letterCI = TRUE`
  read letters from comparison intervals around these means and could call
  two means different when their test did not (oats with Victory's N = 0.6
  plots lost: Marvellous "A", Victory "B", p = 0.07), while
  `settings$comparison$disagree` reported 0 because it judged overlap from
  the pair estimates instead of the means shown. Now: column `over` of
  `means` names the cells such a mean averages (the factors with levels
  missing, e.g. `"N 0, 0.2, 0.4"`); when one of its pairs is not the
  difference of the shown means (compared on common cells, or NA) it has no
  comparison interval (`comparison.low`/`high` NA,
  `settings$comparison$left_out`), `letterCI` takes the letters from the
  p-values (with a message), the LSD note says to compare its pairs, not the
  differences of the means, with the LSD, and `autoplot()` draws no
  comparison or half-LSD bar for it; every mean over empty cells is an
  open point. One rule decides this (`pm_misread_means()`). `disagree`
  counts from the means shown. Means and pairs are unchanged; in the
  snapshot of 52 calls only the four with empty cells or confounded
  factors change (confounded factors: only `over` is added, their pairs
  being the differences of the means).
* **One layout for the comparison tables.** `contrastmeans()` returns class
  `"contrastmeans"`: the 1.x elements keep their places (`[[1]]` the
  rounded table with the adjusted p-value in `Pr(>|t|)`, `$K`, the joint
  test), followed by `contrasts` (unrounded: `contrast`, `estimate`,
  `std.error`, `df`, `statistic`, `p.value` unadjusted, `p.adj`,
  `p.equiv`, `conf.low`, `conf.high`, `adj.conf.low`, `adj.conf.high`
  instead of `adj.lower.CL`/`adj.upper.CL`, `p.method`, and `ratio` or
  `odds.ratio` on a log or logit scale when each contrast's weights sum to
  zero), `joint` (unrounded) and `settings`; it prints that table, and
  `tidy()`/`as.data.frame()` return it. `pm_slopes()`: the column `slope`
  is now `estimate` (columns in the order of `predictmeans()$pairs`); its
  `pairs` gain `conf.low`, `conf.high` and, with tukey, bonferroni, sidak,
  scheffe or mvt, `adj.conf.low`/`adj.conf.high`; `settings`, `tidy()` and
  `as.data.frame()` are added. The pairs of `predictmeans()` on `semireg`
  fits gain `conf.low`, `conf.high`, `lsd`, `adj.conf.*` and ratios or odds
  ratios, as for other models. F tests use the names of `pm_joint()`
  (`df1`, `df2`, `F.ratio`, `p.value`), also in
  `covariatemeans()$settings$test`.
* **`semipred()` returns a classed result and draws its plot only when
  asked** (`plot = TRUE`, `options(predictmeans.plot = TRUE)`, `plotord`,
  or `prt = TRUE`, the 1.x name; `prt = FALSE` still never draws). The
  result (class `c("semipred", "pm_curves")`) has the layout of
  `covariatemeans()` and `predictmeans()`: `means` (curves or surfaces),
  `slopes` (`deriv = 1, 2`, with the t test of slope = 0) or `pairs`
  (`contr`), with `estimate`, `std.error`, `df`, pointwise `conf.low`/
  `conf.high` and, with `band = "simultaneous"`, the simulated band in
  `adj.conf.low`/`adj.conf.high` with adjusted p-values `p.adj` (the
  share of the simulated maxima at least as large, so `p.adj < 0.05`
  exactly where the band leaves out zero); `response` columns for
  generalized fits and `trans`; `settings` (df method, band, `deriv`,
  `unconditional`, `boundary`, back-transformation); `significant`.
  `tidy()`, `as.data.frame()`, `autoplot()` and `print()` work on it. The
  1.x `pred_df` and `plt` are unchanged. The ranges where curves or slopes
  differ are printed with the result instead of as messages.
* **`covariatemeans()` returns a classed result and draws its plot only
  when asked** (`plot = TRUE`, `options(predictmeans.plot = TRUE)`,
  `plotord` or `newwd = TRUE`), as `predictmeans()`. The result (class
  `c("covariatemeans", "pm_curves")`) prints a short summary and holds
  the tables in the layout of `predictmeans()`: `means` (the curves:
  factor columns, the covariate under its own name, `estimate`,
  `std.error`, `df`, `conf.low`, `conf.high`, and `response` columns when
  back-transformed) and, with `contr`, `pairs` (`level1`, `level2`, the
  covariate, `estimate`, `std.error`, `df`, `statistic`, `p.value`,
  pointwise `conf.low`/`conf.high`; with `band = "simultaneous"` the
  Scheffé `p.adj` and `adj.conf.low`/`adj.conf.high`, which agree; ratio
  or odds-ratio and response-scale columns), plus `settings` (with the
  joint test in `settings$test`). `tidy()`, `as.data.frame()` and
  `autoplot()` return them. The 1.x fields `plt`, `pltdf`, `test` and
  `significant` are unchanged; the summary of the difference (the joint
  test and where the curves differ) is printed with the result instead of
  as messages. Scripts that relied on the plot appearing need
  `plot = TRUE` or `autoplot()`. (Before, `plot = FALSE` was taken as
  `plotord = FALSE` by partial matching.)
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
  session, draws the 1.x plots; `plot = c("mean", "bar")` chooses them.
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
  `library(glmmTMB)`. lme4 and nlme are still attached (lmerTest: see
  above).
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

* **Faster `R2_model(method = "nakagawa")`** for `lmer`/`glmer`/`glmmTMB`:
  the random-effect variance is computed with sparse matrices (the model's
  Z stays sparse; G = I (x) vcov as a sparse Kronecker product). 10,000
  observations with 2,000 subjects and random slopes: 136 s -> 0.02 s; same
  values (all.equal at 1e-12). Internal products `A diag(v) A'` are now
  column scalings (`varcomp()`, simultaneous intervals, nlme observed
  information); results unchanged to 3e-15. Standard errors of means and
  curves are `rowSums((K V) * K)` (no m x m product); the permutation and
  robust tests take (X'X)^-1 from the QR they already have (`pm_xtx_inv()`).
  The Wald statistic (L b)' (L V L')^-1 (L b) is computed in one place,
  `pm_wald_quad()`, by all tests (Type II/III tables, F tests,
  Kenward-Roger, permutations), without forming an inverse. For `lme` fits
  whose residual groups are not nested in the random-effect groups, the
  permutation tests build the covariance sparse and work on its independent
  blocks instead of one dense n x n matrix; unused block-matrix helpers
  removed. The covariance square roots of the permutation tests are set up
  in one place (`pm_mer_whiten()`, `pm_mer_perm_parts()`, `pm_nlme_whiten()`)
  for `permmodels()`, `permlmer()` and the Freedman-Lane tests.

* **`PMplot()` and `autoplot(type = "pvalue")`**: the legend now always lists
  every p-value band in order, each with its colour (a band with no tiles had
  been shown without a key, so keys and colours no longer lined up). With
  `atvar` each level is its own square panel titled by the level (before: one
  large grid, mostly empty, with level labels repeated on the axes and no
  name for the groups); `p_valueMatrix` of `predictmeans()` is named by the
  levels of `atvar`. The legend is anchored at its top-right corner inside a
  single matrix of 6 or more rows, else placed to the right, so it no longer
  runs off the plot or covers tiles.

* **Type of tests, one rule everywhere**: `anova()` of `lmer` and `semireg`
  fits, `permmodels()` and `bootmodels()` take `type` as `"I"`, `"II"`,
  `"III"` or `1`, `2`, `3` (any case); other values stop with a clear
  message. `anova()` of a `semireg` fit now has Type I tests (Satterthwaite F
  for gaussian fits, sequential Wald chi-square otherwise).
* **Advice on the type of tests**: with unbalanced data or empty cells,
  `anova()` of `lmer` and `semireg` fits adds a short note on what Types I, II
  and III test and when to use each (and to centre a covariate that is in an
  interaction for Type III); `advice = FALSE` turns it off. `?anova.lmerMod`
  has a new section "Which type of test".
* **`pm_check()`**: printed by section (Fit, Residuals, Variance,
  Distribution, Design, Terms) with a status per row (`!`, `i`, `OK`) and a
  summary line; `summary()` keeps the rows worth a look. New checks:
  dispersion in both directions for Poisson, binomial with trials and negative
  binomial `glm`/`glmer`/`glmmTMB` fits (the row is now "Dispersion (Pearson
  chi-square / df)"; for 0/1 responses and families with an estimated
  dispersion an information row), excess zeros (simulated, fixed seed),
  separation in binomial models, balance of cells, collinearity
  (GVIF^(1/(2 df))), random effects with fewer than 5 levels, and spread of the
  residuals against the fitted values. The result has new columns `section`
  and `status`.

* **`Caschool`**: the California test score data (420 school districts,
  1998-1999; Stock and Watson; as in package Ecdat), the data of the overview
  article's example with several smooth terms, now in the package
  (`data(Caschool)`). Its 1.1.1 and 2.0 models give the same fit
  (`dev/check-smz-migration.R`).

* `semipred()` on generalized fits: the table of the curves has the SE on
  the link scale (`SE`, from which the limits h(eta -+ crit SE) come) next
  to the delta-method SE of the means (`SE_mu`); before, only `SE_mu` was
  given. The bands of generalized fits and of their differences were
  checked against independent calculations (`dev/semireg-golden/
  check-bands-glmm.R`: M. P. Wand's binary example, a Poisson split plot)
  and the gaussian contrast band against Wand's growth-curve script
  (`check-bands.R`): equal to rounding error given the variance
  parameters.

* **Clearer colours in `semipred()` plots, set in one place, and
  `semipred(plot_args = )` to change them.** Each band is drawn first and
  its line on top in a darker shade of the same hue (1.x drew the band over
  the line, which dimmed it): one curve and its slopes in blue, the
  difference of two curves (`contr`) in green as in 1.x, the reference line
  at 0 (1 for ratios) dashed red. Under a difference, and under one slope,
  red bars mark the stretches where the band leaves out the reference.
  Curves by group take a palette that colour-blind readers can tell apart
  (Okabe and Ito; ggplot2's own beyond eight groups); the data are grey and
  half transparent. `plot_args` holds any of `line`, `band`, `diff_line`,
  `diff_band`, `ref`, `alpha`, `linewidth`, `points`, `point_alpha`,
  `point_size`, `palette`, `surface_points`. The difference plots of
  `covariatemeans(contr = )` use the same colours. `ci = FALSE` now also
  leaves out the band of a single curve.

* **`semiboundary()`: the region covered by two-dimensional data** (such as
  longitude and latitude), for `semipred(boundary = )`; and
  `semipred(boundary = TRUE)` uses it on the model's data. The region is
  the closing of the data with a distance d (every place within d of a
  point, shrunk back by d less a small margin), so it follows concave
  shapes, fills small gaps and keeps large empty areas as holes and
  distant clusters as separate pieces; d = 2 x the 95th percentile of the
  nearest-neighbour distances, so it follows the density of the data
  (`mult`, `d` to change it). `scale` puts the covariates on comparable
  scales: `"lonlat"` (longitude shrunk by the cosine of the latitude; found
  from the column names), `"same"` (map coordinates in one unit), `"sd"`
  (two different variables). Holes smaller than the area of 20 data points
  (`min_hole`) are filled. The outline is traced on a grid that gets finer
  with the density of the data: 0.1 s for 150 points, 4 s for 10,000, 9 s
  for 50,000. The help page has three examples with graphs (the scallop
  survey against its original boundary, copper drill holes with three
  values of `mult`, random points with a hole and a separate piece),
  surfaces drawn inside the boundary, and a section of suggestions. With twice the median
  spacing, uniform random designs of 50-1000 points gave 4 pieces and 12
  spurious holes on average; with the 95th percentile one piece and no
  holes, every point inside. On the `scallop` survey it covers 85% of the
  area of the boundary drawn with the original analysis (`scallopBdry`).
  `boundary` may now hold several rings separated by `NA` rows (separate
  pieces and holes); one polygon works as before. Only base R
  (`grDevices::contourLines()`), no new dependency.

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

* `covariatemeans()` and `semipred()` results have the same layout:
  `means`, `slopes` (`semipred(deriv = )`), `pairs` (with `contr`) and
  `settings`, with `tidy()`, `as.data.frame()`, `autoplot()` and `print()`
  methods (`?tidy.pm_curves`). In all of them `conf.low`/`conf.high` are
  pointwise and a simultaneous band is in `adj.conf.low`/`adj.conf.high`,
  as the adjusted limits of `predictmeans()$pairs`; a simultaneous
  `covariatemeans(contr = )` band also gives the response-scale
  `response.adj.conf.*`.

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

* Simultaneous bands of `semipred()` with t limits: the critical value
  took each point's chi-square by the Wilson-Hilferty approximation, which
  made it up to 0.8% too large at small df (df 8, against
  `mvtnorm::qmvt()`); it now uses the exact quantile transform. Checked
  against mvtnorm for df Inf, 30 and 8 (`dev/check-sim-crit.R`: within the
  Monte Carlo error, about 0.3% with the default 10,000 draws). The
  smoothing-parameter correction (`unconditional = TRUE`) was checked
  against mgcv's (`dev/check-uncond-mgcv.R`: SEs equal to 0.1% gaussian,
  0.5% Poisson), and the bands of `covariatemeans()` against emmeans
  (`dev/check-covariatemeans-bands.R`).

* **`semireg_tmb()` fits with zero inflation: bands and SEs were too
  narrow.** The covariance of the fixed and smooth coefficients used the
  working weights of the count part alone, which treat every zero as a
  count: SEs of the fixed effects were 5-22% too small (an example with 30%
  extra zeros), and so were the bands of `semipred()` and the SEs of
  `predictmeans()` on the fit. For zero-inflated fits the covariance is now
  the inverse of the block of TMB's joint precision for those coefficients
  (given the variance and zero-inflation parameters), the exact one; for
  other fits the two were already equal.

* **glmmTMB fits with zero inflation: the means were those of the
  conditional (count) part, mu, without saying so.** `predictmeans()`,
  `covariatemeans()` and `semipred()` (also on `semireg_tmb()` fits) give
  mu, as emmeans and `predict(type = "conditional")` do, not the overall
  mean (1 - pi) mu; the two can order treatments differently when the
  probability of an extra zero differs (a Poisson example: conditional
  means 4.80 and 5.94, overall means 4.00 and 3.83). A message now says
  which, `settings$component` is `"conditional"`, and the help pages have a
  section "Zero inflation (glmmTMB)". `predictmeans(mean_type =
  "overall")` gives the overall means (1 - pi) mu, with SEs from the joint
  covariance of the conditional and zero-inflation coefficients, averaged
  over the cells of the grid, and their ratios in `pairs`; equal to
  glmmTMB's `predict(type = "response")` and to emmeans with `component =
  "response"` (means and SEs). `mean_type = "conditional"` stays the
  default (the maintainer's decision of 2026-11-17); for models without
  zero inflation the argument has no effect. `covariatemeans()` and
  `semipred()` (also on `semireg_tmb()` fits) take `mean_type` too: curves
  of overall means with bands exp(log m -+ crit SE), and with `contr` their
  ratios and differences; equal to glmmTMB's `predict(type = "response")`
  and to emmeans with `component = "response"` (curves, SEs, differences
  and their SEs).

* `semipred()` and `predictmeans()` on `semireg()` fits stopped with "Can't
  find the variable(s) ..." when a random curve `(s(x, ...) | g)` or a
  `t2()` term had an object of the session in its arguments, such as
  `xt = list(range = rg)` or `k = kk`: those names were taken for data
  columns. Only the covariates, `by` and the group are now (as for `s()`
  already).

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

* A package logo (`man/figures/logo.png`: the README, the website and the
  package help page) and a workflow chart from the fitted model to the report
  (`man/figures/workflow.png`: the README, `?predictmeans-package` and
  `vignette("predictmeans")`).

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

* `R/defunct.R`: one table of the arguments and functions removed in 2.0 with their 2.0 spelling; `pm_defunct()` stops a call that names one (with the value given, e.g. `level = 0.05` -> `conf.level = 0.95`), `pm_defunct_fun()` is the body of a removed function. A test checks that no removed name can be partially matched to a remaining argument.
* `dev/review-options.md`: inventory of every argument of 1.1.1 and 2.0 (keep, rename, merge, remove) for the clean-up of the options before 2.0.
* Shorter empty-cell messages: the confounding warning (about 450
  characters) now names the two factors and the three ways to state the
  structure in one sentence; the empty-cell warning, the message on
  non-estimable means and the `anova()` heading point to
  `vignette("missing-cells")` / `?anova.lmerMod` for the explanation.
* One helper (`pm_result_table()`, `R/tidy.R`) gives the tables of every
  classed result to `tidy()` and `as.data.frame()` (`predictmeans`,
  `covariatemeans`/`semipred`, `pm_slopes`, `contrastmeans`); the error
  for a missing table reads "There is no 'pairs' table in this result".
  `contrastmeans()` builds its table once and reads the 1.x table from it
  (identical at tolerance 0 in 9 cases: permutations, equivalence, robust
  `vcov`, user `Df`, empty cells).
* The package's own numerical derivatives are in one file (`R/deriv.R`)
  with one rule: central differences with the step c max(|x|, scale),
  c = 1e-5 for first derivatives, 1e-4 for second and 1e-3 for third (the
  usual balance of truncation and rounding error), the scale saying what x
  is (1 for unconstrained variance parameters, 0.01 for lme4's theta, a
  covariate's spread with no location). Five copies with their own steps
  used it before (`semireg()` df, the covariance pieces of gls/lme,
  back-transformations, `pm_slopes()`, `semipred(deriv = )`).
  `semipred(deriv = )` and the second differences keep their steps; the
  first derivatives that stepped 1e-4 now step 1e-5 (more accurate):
  results move by at most 1e-9 relative (the 52 calls of
  `dev/snapshot-predictmeans.R` equal to 1e-12; the df of `semipred()`
  bands 3e-10). The derivative of a
  back-transformation (`trans`) for the response-scale SEs is the same
  function as for the bias adjustment and response-scale averaging (it was
  numDeriv's in one place, central differences in the other; the SEs
  moved by 1e-11). numDeriv
  stays where a derivative is computed once per model (lmer Satterthwaite
  df, equal to lmerTest's; Hessians of lme4 variance parameters;
  `varcomp()`'s delta method).

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
