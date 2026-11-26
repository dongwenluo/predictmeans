# predictmeans 2.0.0

predictmeans 2.0.0 is a major release. Defaults changed where simulations or the literature showed a more
reliable method; each such change names the argument that gives the 1.x result. Results are checked in the
tests against SAS PROC MIXED and Genstat runs, pbkrtest, nlme and emmeans. Every removed argument or function
stops with its 2.0 spelling, the value given filled in; nothing is silently ignored. `vignette("migrating")`
shows 1.1.1 and 2.0 code side by side.

## Breaking changes

### Removed functions

| Removed | Use instead |
|---|---|
| `smZ()`, `smoothZ =` | smooths in the formula of `semireg()` / `semireg_tmb()`, e.g. `s(x, bs = "os")` |
| `permglmer()`, `permglmmTMB()` | `permlmer()` (they were the same function) |
| `permanova.lmer()` | `permmodels(model, type = )` (Freedman-Lane; other p-values) |
| `R2_glmm()` | `R2_model(model)` (same values, as numbers) |
| `ci_mcp()` | `predictmeans(..., letter_ci = TRUE)` |

### Renamed arguments

* Dotted names are snake_case: `conf.level` -> `conf_level`, `letterCI` -> `letter_ci`, `as.is` -> `as_is`,
  `test.statistic` -> `test_statistic`, `max.degree` -> `max_degree`, `drop.unused.levels` ->
  `drop_unused_levels`. Result columns keep the broom names (`conf.low`, `std.error`, `p.value`).
* `level` (a significance level) -> `conf_level = 1 - level` in `predictmeans()`, `contrastmeans()`,
  `covariatemeans()`, `pm_slopes()`, `semipred()`, `permstats()`, `PMplot()`. In `varcomp()` `level` was a
  confidence level already: `conf_level` with the same value.
* `Df` -> `df_method = <number>`; `df_term(type = )` -> `df_term(df_method = )`.
* `covariateV` -> `at` in `covariatemeans()`, `semipred()`, `Kmatrix()`.
* `semipred(semireg = )` -> `semipred(model = )`; `simultaneous = TRUE` -> `band = "simultaneous"`.
* `permlmer(lmer0, lmer1)` -> `permlmer(model0, model1)`.
* `residplot(level = )` -> `re_level`; `CookD(idn = )` -> `id`.
* `semireg(ngbinomial = TRUE)` -> `semireg_tmb(family = glmmTMB::nbinom2)`.

### Removed arguments

* Back-transformation: `transOff` (put the offset inside `trans`, e.g. `trans = function(x) exp(x) - 1`);
  `count` (Poisson and negative binomial means are detected); `responsen` (read from the model).
* Printing: `prtnum`, `ndecimal`, `prt`. Nothing is printed during a call; use `print(x)` or
  `print(x, digits = )`.
* `predictmeans()` plot arguments: `plot`, `plot_args`, `which.plots`, `mplot`, `barplot`, `pplot`,
  `bkplot`, `prtplt`, `plotord`, `lineplot`, `lsd_bar`, `jitterv`, `plottitle`, `plotxlab`, `plotylab`,
  `basesz`. `plot = FALSE` is still accepted.
* `covariatemeans()` (`trellis`, `plotord`, `mtitle`, `ci`, `point`, `jitterv`) and `semipred()`
  (`trellis`, `scales`, `plotord`, `ci`, `point`, `jitterv`, `threeD`): use `autoplot()` on the result;
  `plot_args = list(...)` still works (superseded).
* `newwd` (`predictmeans()`, `covariatemeans()`, `residplot()`, `CookD()`, `PMplot()`), `residplot(ask)`,
  `PMplot(margin, legendx)`: ggplot2 draws in the current device.
* `semireg()`: `ncenter`, `nscale`, `resp_scale` (variables are used as they are); `predict_info` in
  `semireg()` and `semireg_tmb()` (always kept).
* `permmodels(fo = )`: the formula is read from the model.

### Other changes to scripts

* `predictmeans()` draws no plots. Use `autoplot(pm)` or `plot(pm)` with `type = "means"`, `"bar"`, `"ci"`
  (the 1.x back-transformed plot) or `"pvalue"`. The result has no `predictmeansPlot` or
  `predictmeansBarPlot`.
* `covariatemeans()` and `semipred()` draw only when asked: `plot = TRUE`, `autoplot()`, or
  `options(predictmeans.plot = TRUE)`.
* `pairwise` defaults to `NULL`: comparisons are made when `adj`, `atvar`, `permlist`, `letter_ci` or
  `pairwise = TRUE` asks for them. `pairwise` only sets what is printed.
* The "Pairwise LSDs" matrix is no longer printed; each pair's LSD is in `pairs$lsd`.
* `covariatemeans()$pltdf` is now `$pred_df` (as in `semipred()`); `$pltdf` still works with a message.
* `semipred(threeD)` / `autoplot(type = "3d")` gives one plotly widget, not a list.
* `contrastmeans()` returns class `"contrastmeans"` with `adj.conf.low`/`adj.conf.high` columns;
  `pm_slopes()` column `slope` is now `estimate`.
* `se_ranef()` returns a data frame (class `pm_ranef`); `as.list()` gives the 1.x list of SE matrices.
* `summary()` of a `semireg()` / `semireg_tmb()` fit returns a list with the mixed-model summary in `$model`.
* Results keep the data needed for plots, not the fitted model.
* `fit$knots_lst`, `range_lst`, `type_lst` and `smoothZ_call` are gone; refit models saved with `smZ()`.
* `library(predictmeans)` attaches lme4 and nlme only. glmmTMB moved to Imports and lmerTest to Suggests: call
  `library(glmmTMB)` or `library(lmerTest)` when needed.
* car, pbkrtest, HRW, lmeInfo, lmeSplines, plotly and splines2 moved to Suggests; plyr is no longer used.
  mgcv moved to Imports.
* `ncore` defaults to 1 (was 3) in `predictmeans()`, `permlmer()` and `permmodels()`.
* `residplot(id)` defaults to 3: the three largest residuals are labelled.

## Changes that affect results

* `anova(<lmer fit>, type = "I")` now uses Kenward-Roger df and F, like Types II and III, and `ddf =
  "Satterthwaite"` gives the previous table; both equal SAS PROC MIXED `HTYPE=1`. Type III tests a nested term
  (`A + A:B`) within each level of the outer factor, as SAS (it had 0 df), and the "empty cells" note no longer
  appears for nested factors or for factors that only have main effects.
* REML `lmer` fits use Kenward-Roger adjusted standard errors in `predictmeans()`, `contrastmeans()` and
  `covariatemeans()`, as SAS `DDFM=KR` and emmeans. Unbalanced designs change slightly. 1.x:
  `vcov_adjust = FALSE`.
* Kenward-Roger uses the observed information of the variance parameters, as SAS; it stays valid with
  missing values. pbkrtest's version: `df_method = "KR-expected"`.
* Degrees of freedom are chosen by model class (`df_method = "auto"`): Kenward-Roger for REML `lmer`, `lme`
  and `gls`; Satterthwaite for ML fits; inner-outer df for GLMMs; z tests for Poisson, binomial and negative
  binomial `glm`. Simulations showed the old df rejected too often (or too rarely for `glm`). 1.x:
  `df_method = "residual"` (or `"containment"` for `lme`/`gls`).
* `permmodels()` uses Freedman-Lane permutations (`method = "auto"`), so other terms keep their effects.
  1.x: `method = "raw"`.
* Permutation tests of mixed models decorrelate residuals with the symmetric square root of their covariance;
  the 1.x Cholesky factor depended on the order of the data. Raw p-values for mixed models also differ.
* Permutation p-values are (b + 1)/(B + 1) unless `perms` holds a complete enumeration.
* `permlmer()` with different fixed effects refits both models by ML, so the permuted statistic is the
  printed likelihood ratio.
* `se_ranef()` defaults to the prediction error SE (`type = "prediction"`, SAS "Std Err Pred"), larger
  than 1.x. 1.x: `type = "conditional"`.
* Type III tests do not depend on the model's contrasts (sum-to-zero hypotheses, as SAS Type 3).
* Means over empty cells keep the 1.x average of the observed cells, but their pairs are compared on the
  cells both means have; the 1.x comparison rejected too often in simulations.
* Means that depend on aliased coefficients are `NA`, with a message; estimable differences are shown.
* `aov()` with `Error()` is analysed by strata, as Genstat, with missing values estimated as Genstat's ANOVA.
  It was refitted as `lmer`. 1.x: fit the mixed model yourself.
* A response transformed in the formula (`log(y)`, `log(y + c)`, `sqrt(y)`, ...) is back-transformed by
  default. `trans = I` keeps the model scale. A given `trans` is now used also for models with a link.
* `bias_adjust = TRUE` gives the mean exactly by quadrature. The 1.x-style second-order formula:
  `bias_adjust = "taylor"`.
* `covariatemeans()` with a log link or `trans = exp` draws the back-transformed fitted value; 1.x added half
  the SE.
* The printed LSD goes with the adjusted p-values (Tukey HSD, Bonferroni, Sidak, Scheffe, mvt critical
  values); with `atvar` it uses the pairs within each level. `adj = "none"` is unchanged.
* `letter_ci = TRUE` reads letters from comparison intervals at the comparisons' level; drop a 1.x
  `level = 0.166`.
* `semipred()` limits are two-sided at the stated level (1.x defaults gave 90%). Gaussian `semireg()` fits
  get t limits with Satterthwaite df; generalized fits include smoothing-parameter uncertainty
  (`unconditional`).
* `semireg()` uses the `nloptwrap` optimizer. 1.x: `control = lmerControl(optimizer = "bobyqa")`.
* Gaussian `semireg_tmb()` fits are REML, as documented (they were ML).
* `glmer` Poisson back-transformed means are counts, as for `glm` and `glmmTMB`.

## New features

### Means and comparisons

* Stable result tables `means`, `pairs` and `settings` with fixed column names; `tidy()`, `glance()`,
  `as.data.frame()`, `coef()` and `vcov()` methods. The 1.x elements are kept.
* `weights = "equal"`, `"proportional"` (SAS `OM`), `"outer"`, `"cells"`, `"flat"` for averaging cells.
* Adjustments `"sidak"`, `"scheffe"` and `"mvt"`; adjusted confidence limits `adj.conf.low`/`adj.conf.high`.
* Comparison intervals: two means differ when their intervals do not overlap.
* `contrast_family()`, `specify_contrast()` (contrasts by level name), `pm_joint()` (joint F tests and
  slices), `pm_slopes()` (slopes per level), `contrastmeans(joint = TRUE)`.
* Equivalence tests (`equivalence =`); robust covariances (`vcov = "HC0"` to `"HC3"`, `~ cluster`, or a
  matrix or function).
* Response-scale means with delta-method SEs, ratios and odds ratios; `average = "response"`;
  `pm_boxcox()` and `pm_yeojohnson()` for formulas.
* `mean_type = "overall"` for zero-inflated glmmTMB models.
* `Kmatrix()` handles `poly()`, splines, transformations, `gam` smooths, covariates by name and
  combinations of factors that are not a model term; `cov_reduce` and `exclude`.
* Faster with many means and pairs; more than 100 means without comparisons skip the pair table.

### Tests and degrees of freedom

* `anova()` of `lmer` fits: `type1_basis = "fitted"` builds Type I hypotheses from the fitted model
  (X'V^-1X), as Genstat, instead of X'X (SAS, lmerTest; the default); `ddf = "KR-average"` gives Genstat's
  d.d.f. with an unadjusted F. Together they reproduce Genstat's Wald and F statistics; Genstat prints each
  term's Type II d.d.f., here each test has the d.d.f. of its own hypothesis.
* `df_method = "Kenward-Roger"`, `"KR-expected"`, `"KR-average"` (Genstat), `"Satterthwaite"`,
  `"containment"`, `"between-within"`, `"inner-outer"`, `"residual"` or a number, also for `lme` and `gls`.
* Inner-outer df of aggregated binomial GLMMs (`cbind(y, m - y)`, `weights` = trials) and of `clmm` with case
  weights count the cases in the residual stratum, so aggregated data and one row per case get the same df
  (`settings$df_cases`; weights that are not whole numbers count rows). Counting rows rejected 0-2.6% of true
  null hypotheses at 5% in simulations, counting trials 3.7-4.3%; `options(predictmeans.inner_outer_count =
  "rows")` counts rows. Trials assume binomial variation within clusters (with overdispersion the model leaves out
  they gave 9-12%, rows 3.5%): the advice `overdispersion_trials` (a caution when the fit's Pearson test rejects, a
  note when its ratio is above 1) says to add an observation-level random effect or to count rows.
* `anova()` of `lmer` fits with Kenward-Roger or Satterthwaite F tests, Types I, II, III, without lmerTest.
* Kenward-Roger and Satterthwaite computed without n x n matrices; much faster.

* Ordinal models: `MASS::polr()` and `ordinal::clm()` fits (location effects). Means on the latent scale (as
  emmeans' `mode = "latent"`), pairs as log odds ratios of a higher class, z tests; `mean_type = "prob"` adds the class
  probabilities at each mean and `mean_type = "mean_class"` the expected class, with delta-method SEs. Also
  `ordinal::clmm()` fits, for a typical subject (random effects at 0), with the inner-outer df of GLMMs (z tests
  rejected 11% of true nulls at 5% with 6 clusters). Checked against emmeans. Summaries say "z tests" when the df
  are infinite. `VGAM::vglm()` fits of the ordinal families `cumulative`, `propodds`, `acat`, `cratio` and `sratio`
  with `parallel = TRUE`: the same latent means, oriented so that a positive difference means higher classes,
  whatever the family and `reverse` (SEs from VGAM's expected information). vglm fits of other families
  (e.g. `poissonff`), which gave link-scale means with t tests before, now stop with a message.
  `nnet::multinom()`, `VGAM::vglm()` with `multinomial()` and ordinal vglm fits with non-parallel slopes: the class
  is a factor of the means, pairs within each class. `mean_type = "prob"` (default): class probabilities on the
  logit scale, back-transformed with limits inside 0-1, pairs as log odds ratios (with two classes the same as a
  binomial glm); `"prob_diff"`: differences of probabilities (as emmeans' `mode = "prob"`); `"link"`: the model's
  linear predictors; `"mean_class"` for ordinal families. Probabilities averaged over other factors on the
  probability scale; delta-method SEs, z tests; advice on separation (a class absent at some levels). Checked
  against emmeans; other functions stop for these fits. Where the cumulative curves of a non-parallel fit cross
  (negative probabilities), the means that use those points are NA with a caution.
* `ordinal::clm()` fits with nominal or scale effects (`nominal = `, `scale = `), refused before: class
  probabilities as for non-parallel vglm fits (`"link"` only without scale effects); equal to `predict.clm()`
  with its SEs (logit, probit, cloglog, loglog links).
* `ordinal::clmm()`: `bias_adjust = TRUE` with `mean_type = "prob"` or `"mean_class"` gives population-average
  class probabilities (integrated over the random intercepts by adaptive numerical integration), with delta-method SEs
  that include the random-effect variances; checked against numerical integration and the probit closed form.
  Random slopes stop.
* `bootmodels()` and `permmodels()` take `ordinal::clmm()` fits (tests with few clusters): a parametric bootstrap
  of the Type II Wald tests under the global null (random intercepts), whose `permlist` gives bootstrap p-values of
  the latent pairs in `predictmeans()`; and re-randomisation of each variable applied to whole clusters, within
  strata (fixed blocks, other cluster-level variables, coarser random effects), exact (all allocations enumerated)
  when there are at most `nperm`, with a note on the smallest attainable p-value (6 clusters split 3/3: 20
  allocations, p at least 0.1). In simulations with 6-12 clusters the default inner-outer t tests (3.4-5.4% at
  5%) and the bootstrap (3.7-5.1%) held the level, z tests did not (9-11%); the permutation was valid but weak.
* `weights = "proportional"` (and `"outer"`, `"cells"`, `"flat"`) count cases instead of rows where a row holds
  several: binomial and quasibinomial responses of any class (trials) and class responses (multinom, polr, clm, clmm,
  vglm with case weights). Aggregated data and one row per case then weight the factor combinations alike (as
  emmeans); covariates are still set from the rows. `settings$weights_by` says which. semireg means stop for
  weights other than `"equal"` (they were ignored).
* New article "Ordinal and categorical responses": latent means of cumulative link models, class probabilities
  of multinomial and non-parallel fits, the scales and the checks.

* New `corAR1xAR1()`: the separable AR1 x AR1 correlation of row x column field trials (as Genstat's
  `VSTRUCTURE AR,AR` and ASReml's `ar1(row):ar1(col)`), with an optional nugget, for `nlme::gls()` and `lme()`;
  every predictmeans tool works with it (Kenward-Roger df, means, SEDs, `residvariogram()`). On Genstat's Slate
  Hall trial phi and the variances agree to Genstat's printed digits, the 25 means and all 300 SEDs within 5e-6,
  and the variety test's d.d.f. with `df_method = "KR-average"`. Kenward-Roger and Satterthwaite df for gls and
  lme now leave out variance components (and a nugget) estimated at 0, and a failed Kenward-Roger falls back to
  Satterthwaite; `na.exclude` fits no longer break them. `residvariogram()` advice now points to it (or to a smooth `t2(row, col)` trend). New article "Spatial analysis of a field trial"
  (Slate Hall: variogram, AR1 x AR1 with and without a nugget, the `t2()` trend).
* `corAR1xAR1(by_group = TRUE)` with `form = ~ row + col | trial`: a phi pair (and nugget) for each trial; with
  `trial / variety` and `varIdent(~ 1 | trial)` each trial gets the estimates, SEDs and Kenward-Roger df of
  its own fit (checked against separate fits).
* `residplot()`, `pm_check()` and `residvariogram()` of `gls` and `lme` fits with a grouped correlation: the
  normalized residuals were misplaced when the data were not sorted by the groups (nlme's
  `residuals(type = "normalized")`, off by up to 2 standard deviations); they are now whitened group by group.
  For `corAR1xAR1()` fits the plots are taken in row and column order, so the residuals do not depend on the
  order of the data.

### Mixed models and variance components

* New `heritability()`: line-mean (entry-mean) heritability from an `lmer` fit with random genotypes, in the
  standard form and the mixed-model forms of Cullis et al. (2006), Piepho and Möhring (2007) and Oakey et
  al. (2006), for unbalanced trials; `method = "all"` shows them side by side. `method = "family"` gives the
  narrow-sense heritability of a single plant from half-sib or full-sib families (`relationship`), and
  `"repeatability"` that of single records; `"reliability"` (mean reliability) and `"simulated"` (Piepho and
  Möhring's simulated r2). Plot means: `within_variance`, `within_n` give the single-plant variance. SEs:
  delta method, or a parametric bootstrap (`nboot`). Works with weights, several genotype x environment terms, genotypes
  nested in a fixed factor (with a note) and boundary fits (genotypic variance 0). Reproduces every value of
  Schmidt et al. (2019, Figure 1) for the oats alpha design.
* `varcomp()` for every model class, with observed-information SEs and Satterthwaite intervals;
  `test = TRUE` tests each component (likelihood ratio, boundary mixture or parametric bootstrap).
* `se_ranef()` for `lmer`, `glmer`, `glmmTMB`, `lme` and `semireg` fits, with any covariance structure; checked
  against SAS PROC MIXED; `predict = TRUE` adds the fixed part to each BLUP.
* `R2_model()` with `method = "piepho"` or `"nakagawa"`, for every model class; reproduces Piepho (2023).

### Permutation and bootstrap

* Freedman-Lane permutations for `lm`, `glm`, `lmer`, `lme`, `gls`, `glmer` and `glmmTMB`, with or without
  refits (`refit`).
* `permlmer()` tests GLMMs on the working scale; `method = "quantile"` for 0/1 responses.
* Permutation confidence limits, permutation Tukey and `adj = "westfall"`; `robust = TRUE` (HC3, sign flips).
* `bootmodels()` for parametric bootstrap p-values, usable as `permlist`.
* `permindex()` and `perms` for restricted permutations; `permstats()`; `permcombine()`.

### Semiparametric models

* Smooths in the formula of `semireg()` and `semireg_tmb()`: mgcv bases, `by =`, `s(x, f)` and
  `s(x, f1*f2)` deviation curves, `bs = "os"`, `bs = "ztps"`, `bs = "user"`, `t2()` tensor products, random
  curves `(s(x) | g)`, glmmTMB covariance structures.
* `semipred()`: simultaneous bands, `deriv = 1, 2`, curve differences with a test, `bias_adjust`,
  `boundary`; `semiboundary()` for two-dimensional data.
* `anova()` of `semireg()` fits with tests of smooth variances; `predictmeans()` on `semireg()` fits.
* Fits print a short account of the smooth terms; `coef()`, `vcov()`, `fitted()` and other methods work.

### Plots

* `autoplot()` for every result: means, bars, confidence limits, p-value heat maps, curves, surfaces, 3-D,
  differences; `x`, `group`, `panel`, `interval`, `data`, `scales`.
* `residplot()`, `CookD()` and `PMplot()` are drawn with ggplot2 for every model class;
  `residplot(interactive = TRUE)` uses plotly. Many overlapping points are
  drawn as transparent open circles (or, with `dense = "bins"`, as shaded bins with the outstanding points on top); the Q-Q plot
  of many residuals keeps every tail point.

### Model checking and advice

* `pm_check()`: a table of checks of a fitted model (outliers, influence, normality, equal variances,
  dispersion, separation, collinearity, ...).
* New `residvariogram()`: the sample variogram of a model's residuals over a row x column field layout or
  coordinates (also latitude/longitude, by direction), one per trial with `by`, with a global envelope from data
  simulated from the fitted model and a test of spatial pattern; advice suggests a spatial model when neighbouring
  plots are more alike than the model allows. The estimator reproduces Genstat's `FVARIOGRAM` (Brooms Barn
  example). `pm_check(row = , col = )` adds the check. Rows or columns numbered within blocks stop with a message.
* Advice printed with results when the fit or data call for it; `advice()` returns it,
  `options(predictmeans.advice = FALSE)` hides it.

### Results and extension

* Exported generics `pm_coef()`, `pm_vcov()`, `pm_basis()`, `pm_data()`, `pm_df()`, `pm_linkinv()` for new
  model classes (`?pm_coef`). The columns of the design are matched to the coefficients by name; a name
  without a column stops with a message, never a silent mismatch.
* `pm_coef()` of an `aov` fit keeps its aliased coefficients as `NA`, as for `lm` and `glm`; every front end
  now reads a fit's coefficients and covariance only through `pm_coef()` and `pm_vcov()`.
* `?pm_coef` says what a new model class gets: means and pairs need `coef()`, `vcov()` and `terms()`;
  `residplot()` and `pm_check()` also `fitted()` and `residuals()` (else they stop and say so); the other
  functions list the classes they handle.
* New data sets, including `Caschool`, `DataExam5.2`, `copper` and M. P. Wand's examples.

## Bug fixes

* Estimability checks (means, their pairs, empty-cell contrasts, slopes) no longer depend on the units of a
  covariate: with a site mean in large units an aliased slope was reported as 0 with SE 0.
* `anova()` of `lmer` fits: a factor nested in another only in the data (`y ~ G + T`) is tested within it
  in Type III, as with `G/T`, with a caution to write `G/T`; balanced nested designs are no longer called
  unbalanced; ML fits get Satterthwaite df with a message instead of a Kenward-Roger label; a failed
  Kenward-Roger adjustment in Types II and III falls back to Satterthwaite.
* `pm_slopes()` gives no pairs when the slopes are equal by the model (no covariate x term interaction).
* `covariatemeans()` drew curves that the data do not identify (a covariate constant within the levels of a
  factor in the model, as a site mean with site; or a function of other factors). It now stops with a
  message pointing to `contr =` and `pm_slopes()`, which estimate the differences between curves; single
  points that cannot be estimated are NA, with advice.
* `pm_slopes()` gave a slope of 0 with SE 0 when the covariate is aliased with a factor (a site mean with
  site in the model). Such slopes are now NA, with advice; their differences are kept. The slopes table has
  an `estimable` column, as the means table.
* `autoplot()` of means by group with the data (`data = TRUE`) works with ggplot2 3.4, as `DESCRIPTION`
  promises (it stopped: "position_jitterdodge() requires at least one aesthetic to dodge by").
* Permutation tests, `bootmodels()` and variance-component tests of `lmer`/`glmer` fits with missing values
  work with older lme4 (1.1-35 dropped the missing rows twice in `refit()`).
* A `glmmTMB` fit with aliased coefficients (e.g. combinations of levels without data) no longer stops
  every means function with "logical subscript too long"; the means equal those of the same `lmer` ML fit.
* `semireg()` with a deviation smooth `s(x, f1*f2)` no longer gives the warning "number of items to replace".
* `permlmer()` reports lme4 convergence warnings of the permuted refits once, as a count, instead of one
  warning per refit.

* Tukey p-values within `atvar` levels used the number of pairs instead of means.
* `predictmeans()` failed for `lmer` models with a missing treatment combination, and for ML `lmer` fits.
* `predictmeans()` with `permlist` and `atvar` failed; a user `Df` for `lmer` gave NA p-values.
* Factor levels containing `":"` broke the pairwise tables.
* Within-`atvar` adjusted p-values could be matched to the wrong pairs.
* `contrastmeans()` failed with any missing treatment combination.
* Kenward-Roger ignored prior weights and offsets of `lmer` fits.
* `semireg()` 0/1 responses had SEs 15-25% too small; Gamma and inverse-Gaussian SEs were too small.
* `semireg()` gaussian fits had wrong `anova()` df and `varcomp()` SEs.
* `semireg_tmb()` negative binomial SEs were scaled by theta; zero-inflated bands were too narrow.
* `semipred(contr = )` bands of generalized models were wrong.
* `varcomp()` SEs for `glmer` and limits for `glmmTMB` were wrong.
* `R2_glmm()` was wrong for some glmmTMB families, offsets and correlated random effects.
* `permindex()` with `block` or `group` returned indices for sorted data.
* `permmodels()` permuted only one column of `cbind()` responses; type I `lmer` tables compared different
  statistics.
* `permmodels()` and `CookD()` could refit with a data set of the same name found elsewhere.
* `CookD()` failed for every `lmer` model.

## Documentation

* The help pages of `anova()` for `lmer` and `semireg` fits and of `se_ranef()` show their usage (CRAN
  pre-test: "Rd files without \usage").
* New vignettes: "Migrating from predictmeans 1.x", "Checking the model", "Missing cells and non-estimable
  means", "Smooth curves in mixed models", "Permutation and bootstrap tests", "Link and response scale",
  "Contrasts, joint tests and slopes", and "Matching SAS PROC MIXED" (the SAS and Genstat reference runs).
* `vignette("predictmeans")` rewritten, with a table of default df by model class.
* Help pages written with roxygen2; `?advice` lists every advice id.
* A pkgdown website with worked-example articles: an overview, tree improvement, and adjusted R² for GLMMs
  reproducing Piepho (2023).

The full development log of 2.0.0 is in the package source, file `dev/NEWS-2.0-development.md`.
