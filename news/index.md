# Changelog

## predictmeans 2.0.0

predictmeans 2.0.0 is a major release. Defaults changed where
simulations or the literature showed a more reliable method; each such
change names the argument that gives the 1.x result. Results are checked
in the tests against SAS PROC MIXED and Genstat runs, pbkrtest, nlme and
emmeans. Every removed argument or function stops with its 2.0 spelling,
the value given filled in; nothing is silently ignored.
[`vignette("migrating")`](https://dongwenluo.github.io/predictmeans/articles/migrating.md)
shows 1.1.1 and 2.0 code side by side.

### Breaking changes

#### Removed functions

| Removed | Use instead |
|----|----|
| [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md), `smoothZ =` | smooths in the formula of [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md) / [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md), e.g. `s(x, bs = "os")` |
| [`permglmer()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md), [`permglmmTMB()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) | [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md) (they were the same function) |
| [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) | `permmodels(model, type = )` (Freedman-Lane; other p-values) |
| [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) | `R2_model(model)` (same values, as numbers) |
| [`ci_mcp()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) | `predictmeans(..., letter_ci = TRUE)` |

#### Renamed arguments

- Dotted names are snake_case: `conf.level` -\> `conf_level`, `letterCI`
  -\> `letter_ci`, `as.is` -\> `as_is`, `test.statistic` -\>
  `test_statistic`, `max.degree` -\> `max_degree`, `drop.unused.levels`
  -\> `drop_unused_levels`. Result columns keep the broom names
  (`conf.low`, `std.error`, `p.value`).
- `level` (a significance level) -\> `conf_level = 1 - level` in
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md),
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md),
  [`permstats()`](https://dongwenluo.github.io/predictmeans/reference/permstats.md),
  [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md).
  In
  [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  `level` was a confidence level already: `conf_level` with the same
  value.
- `Df` -\> `df_method = <number>`; `df_term(type = )` -\>
  `df_term(df_method = )`.
- `covariateV` -\> `at` in
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md),
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).
- `semipred(semireg = )` -\> `semipred(model = )`; `simultaneous = TRUE`
  -\> `band = "simultaneous"`.
- `permlmer(lmer0, lmer1)` -\> `permlmer(model0, model1)`.
- `residplot(level = )` -\> `re_level`; `CookD(idn = )` -\> `id`.
- `semireg(ngbinomial = TRUE)` -\>
  `semireg_tmb(family = glmmTMB::nbinom2)`.

#### Removed arguments

- Back-transformation: `transOff` (put the offset inside `trans`,
  e.g. `trans = function(x) exp(x) - 1`); `count` (Poisson and negative
  binomial means are detected); `responsen` (read from the model).
- Printing: `prtnum`, `ndecimal`, `prt`. Nothing is printed during a
  call; use `print(x)` or `print(x, digits = )`.
- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  plot arguments: `plot`, `plot_args`, `which.plots`, `mplot`,
  `barplot`, `pplot`, `bkplot`, `prtplt`, `plotord`, `lineplot`,
  `lsd_bar`, `jitterv`, `plottitle`, `plotxlab`, `plotylab`, `basesz`.
  `plot = FALSE` is still accepted.
- [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  (`trellis`, `plotord`, `mtitle`, `ci`, `point`, `jitterv`) and
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  (`trellis`, `scales`, `plotord`, `ci`, `point`, `jitterv`, `threeD`):
  use
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  on the result; `plot_args = list(...)` still works (superseded).
- `newwd`
  ([`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
  [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md),
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md),
  [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md)),
  `residplot(ask)`, `PMplot(margin, legendx)`: ggplot2 draws in the
  current device.
- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md):
  `ncenter`, `nscale`, `resp_scale` (variables are used as they are);
  `predict_info` in
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  and
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  (always kept).
- `permmodels(fo = )`: the formula is read from the model.

#### Other changes to scripts

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  draws no plots. Use `autoplot(pm)` or `plot(pm)` with
  `type = "means"`, `"bar"`, `"ci"` (the 1.x back-transformed plot) or
  `"pvalue"`. The result has no `predictmeansPlot` or
  `predictmeansBarPlot`.
- [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  and
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  draw only when asked: `plot = TRUE`,
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html),
  or `options(predictmeans.plot = TRUE)`.
- `pairwise` defaults to `NULL`: comparisons are made when `adj`,
  `atvar`, `permlist`, `letter_ci` or `pairwise = TRUE` asks for them.
  `pairwise` only sets what is printed.
- The “Pairwise LSDs” matrix is no longer printed; each pair’s LSD is in
  `pairs$lsd`.
- `covariatemeans()$pltdf` is now `$pred_df` (as in
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md));
  `$pltdf` still works with a message.
- `semipred(threeD)` / `autoplot(type = "3d")` gives one plotly widget,
  not a list.
- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  returns class `"contrastmeans"` with `adj.conf.low`/`adj.conf.high`
  columns;
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md)
  column `slope` is now `estimate`.
- [`se_ranef()`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  returns a data frame (class `pm_ranef`);
  [`as.list()`](https://rdrr.io/r/base/list.html) gives the 1.x list of
  SE matrices.
- [`summary()`](https://rdrr.io/r/base/summary.html) of a
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  /
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  fit returns a list with the mixed-model summary in `$model`.
- Results keep the data needed for plots, not the fitted model.
- `fit$knots_lst`, `range_lst`, `type_lst` and `smoothZ_call` are gone;
  refit models saved with
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md).
- [`library(predictmeans)`](https://dongwenluo.github.io/predictmeans/)
  attaches lme4 and nlme only. glmmTMB moved to Imports and lmerTest to
  Suggests: call
  [`library(glmmTMB)`](https://github.com/glmmTMB/glmmTMB) or
  [`library(lmerTest)`](https://github.com/runehaubo/lmerTestR) when
  needed.
- car, pbkrtest, HRW, lmeInfo, lmeSplines, plotly and splines2 moved to
  Suggests; plyr is no longer used. mgcv moved to Imports.
- `ncore` defaults to 1 (was 3) in
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  and
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md).
- `residplot(id)` defaults to 3: the three largest residuals are
  labelled.

### Changes that affect results

- REML `lmer` fits use Kenward-Roger adjusted standard errors in
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  and
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
  as SAS `DDFM=KR` and emmeans. Unbalanced designs change slightly. 1.x:
  `vcov_adjust = FALSE`.
- Kenward-Roger uses the observed information of the variance
  parameters, as SAS; it stays valid with missing values. pbkrtest’s
  version: `df_method = "KR-expected"`.
- Degrees of freedom are chosen by model class (`df_method = "auto"`):
  Kenward-Roger for REML `lmer`, `lme` and `gls`; Satterthwaite for ML
  fits; inner-outer df for GLMMs; z tests for Poisson, binomial and
  negative binomial `glm`. Simulations showed the old df rejected too
  often (or too rarely for `glm`). 1.x: `df_method = "residual"` (or
  `"containment"` for `lme`/`gls`).
- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  uses Freedman-Lane permutations (`method = "auto"`), so other terms
  keep their effects. 1.x: `method = "raw"`.
- Permutation tests of mixed models decorrelate residuals with the
  symmetric square root of their covariance; the 1.x Cholesky factor
  depended on the order of the data. Raw p-values for mixed models also
  differ.
- Permutation p-values are (b + 1)/(B + 1) unless `perms` holds a
  complete enumeration.
- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  with different fixed effects refits both models by ML, so the permuted
  statistic is the printed likelihood ratio.
- [`se_ranef()`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  defaults to the prediction error SE (`type = "prediction"`, SAS “Std
  Err Pred”), larger than 1.x. 1.x: `type = "conditional"`.
- Type III tests do not depend on the model’s contrasts (sum-to-zero
  hypotheses, as SAS Type 3).
- Means over empty cells keep the 1.x average of the observed cells, but
  their pairs are compared on the cells both means have; the 1.x
  comparison rejected too often in simulations.
- Means that depend on aliased coefficients are `NA`, with a message;
  estimable differences are shown.
- [`aov()`](https://rdrr.io/r/stats/aov.html) with `Error()` is analysed
  by strata, as Genstat, with missing values estimated as Genstat’s
  ANOVA. It was refitted as `lmer`. 1.x: fit the mixed model yourself.
- A response transformed in the formula (`log(y)`, `log(y + c)`,
  `sqrt(y)`, …) is back-transformed by default. `trans = I` keeps the
  model scale. A given `trans` is now used also for models with a link.
- `bias_adjust = TRUE` gives the mean exactly by quadrature. The
  1.x-style second-order formula: `bias_adjust = "taylor"`.
- [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  with a log link or `trans = exp` draws the back-transformed fitted
  value; 1.x added half the SE.
- The printed LSD goes with the adjusted p-values (Tukey HSD,
  Bonferroni, Sidak, Scheffe, mvt critical values); with `atvar` it uses
  the pairs within each level. `adj = "none"` is unchanged.
- `letter_ci = TRUE` reads letters from comparison intervals at the
  comparisons’ level; drop a 1.x `level = 0.166`.
- [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  limits are two-sided at the stated level (1.x defaults gave 90%).
  Gaussian
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fits get t limits with Satterthwaite df; generalized fits include
  smoothing-parameter uncertainty (`unconditional`).
- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  uses the `nloptwrap` optimizer. 1.x:
  `control = lmerControl(optimizer = "bobyqa")`.
- Gaussian
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  fits are REML, as documented (they were ML).
- `glmer` Poisson back-transformed means are counts, as for `glm` and
  `glmmTMB`.

### New features

#### Means and comparisons

- Stable result tables `means`, `pairs` and `settings` with fixed column
  names; `tidy()`, `glance()`,
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) methods. The 1.x
  elements are kept.
- `weights = "equal"`, `"proportional"` (SAS `OM`), `"outer"`,
  `"cells"`, `"flat"` for averaging cells.
- Adjustments `"sidak"`, `"scheffe"` and `"mvt"`; adjusted confidence
  limits `adj.conf.low`/`adj.conf.high`.
- Comparison intervals: two means differ when their intervals do not
  overlap.
- [`contrast_family()`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md),
  [`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md)
  (contrasts by level name),
  [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
  (joint F tests and slices),
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md)
  (slopes per level), `contrastmeans(joint = TRUE)`.
- Equivalence tests (`equivalence =`); robust covariances
  (`vcov = "HC0"` to `"HC3"`, `~ cluster`, or a matrix or function).
- Response-scale means with delta-method SEs, ratios and odds ratios;
  `average = "response"`;
  [`pm_boxcox()`](https://dongwenluo.github.io/predictmeans/reference/pm_boxcox.md)
  and
  [`pm_yeojohnson()`](https://dongwenluo.github.io/predictmeans/reference/pm_boxcox.md)
  for formulas.
- `mean_type = "overall"` for zero-inflated glmmTMB models.
- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  handles [`poly()`](https://rdrr.io/r/stats/poly.html), splines,
  transformations, `gam` smooths, covariates by name and combinations of
  factors that are not a model term; `cov_reduce` and `exclude`.
- Faster with many means and pairs; more than 100 means without
  comparisons skip the pair table.

#### Tests and degrees of freedom

- `df_method = "Kenward-Roger"`, `"KR-expected"`, `"KR-average"`
  (Genstat), `"Satterthwaite"`, `"containment"`, `"between-within"`,
  `"inner-outer"`, `"residual"` or a number, also for `lme` and `gls`.
- [`anova()`](https://rdrr.io/r/stats/anova.html) of `lmer` fits with
  Kenward-Roger or Satterthwaite F tests, Types I, II, III, without
  lmerTest.
- Kenward-Roger and Satterthwaite computed without n x n matrices; much
  faster.

#### Mixed models and variance components

- [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  for every model class, with observed-information SEs and Satterthwaite
  intervals; `test = TRUE` tests each component (likelihood ratio,
  boundary mixture or parametric bootstrap).
- [`se_ranef()`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  for `lmer`, `glmer`, `glmmTMB`, `lme` and `semireg` fits, with any
  covariance structure; checked against SAS PROC MIXED; `predict = TRUE`
  adds the fixed part to each BLUP.
- [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
  with `method = "piepho"` or `"nakagawa"`, for every model class;
  reproduces Piepho (2023).

#### Permutation and bootstrap

- Freedman-Lane permutations for `lm`, `glm`, `lmer`, `lme`, `gls`,
  `glmer` and `glmmTMB`, with or without refits (`refit`).
- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  tests GLMMs on the working scale; `method = "quantile"` for 0/1
  responses.
- Permutation confidence limits, permutation Tukey and
  `adj = "westfall"`; `robust = TRUE` (HC3, sign flips).
- [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  for parametric bootstrap p-values, usable as `permlist`.
- [`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)
  and `perms` for restricted permutations;
  [`permstats()`](https://dongwenluo.github.io/predictmeans/reference/permstats.md);
  [`permcombine()`](https://dongwenluo.github.io/predictmeans/reference/permcombine.md).

#### Semiparametric models

- Smooths in the formula of
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  and
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md):
  mgcv bases, `by =`, `s(x, f)` and `s(x, f1*f2)` deviation curves,
  `bs = "os"`, `bs = "ztps"`, `bs = "user"`, `t2()` tensor products,
  random curves `(s(x) | g)`, glmmTMB covariance structures.
- [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md):
  simultaneous bands, `deriv = 1, 2`, curve differences with a test,
  `bias_adjust`, `boundary`;
  [`semiboundary()`](https://dongwenluo.github.io/predictmeans/reference/semiboundary.md)
  for two-dimensional data.
- [`anova()`](https://rdrr.io/r/stats/anova.html) of
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fits with tests of smooth variances;
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  on
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fits.
- Fits print a short account of the smooth terms;
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and other
  methods work.

#### Plots

- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  for every result: means, bars, confidence limits, p-value heat maps,
  curves, surfaces, 3-D, differences; `x`, `group`, `panel`, `interval`,
  `data`, `scales`.
- [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md),
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  and
  [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md)
  are drawn with ggplot2 for every model class;
  `residplot(interactive = TRUE)` uses plotly. Many overlapping points
  are drawn as transparent open circles (or, with `dense = "bins"`, as
  shaded bins with the outstanding points on top); the Q-Q plot of many
  residuals keeps every tail point.

#### Model checking and advice

- [`pm_check()`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md):
  a table of checks of a fitted model (outliers, influence, normality,
  equal variances, dispersion, separation, collinearity, …).
- Advice printed with results when the fit or data call for it;
  [`advice()`](https://dongwenluo.github.io/predictmeans/reference/advice.md)
  returns it, `options(predictmeans.advice = FALSE)` hides it.

#### Results and extension

- Exported generics
  [`pm_coef()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md),
  [`pm_vcov()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md),
  [`pm_basis()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md),
  [`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md),
  [`pm_df()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md),
  [`pm_linkinv()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  for new model classes
  ([`?pm_coef`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)).
  The columns of the design are matched to the coefficients by name; a
  name without a column stops with a message, never a silent mismatch.
- New data sets, including `Caschool`, `DataExam5.2`, `copper` and M. P.
  Wand’s examples.

### Bug fixes

- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  of means by group with the data (`data = TRUE`) works with ggplot2
  3.4, as `DESCRIPTION` promises (it stopped: “position_jitterdodge()
  requires at least one aesthetic to dodge by”).

- Permutation tests,
  [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  and variance-component tests of `lmer`/`glmer` fits with missing
  values work with older lme4 (1.1-35 dropped the missing rows twice in
  [`refit()`](https://rdrr.io/pkg/lme4/man/refit.html)).

- A `glmmTMB` fit with aliased coefficients (e.g. combinations of levels
  without data) no longer stops every means function with “logical
  subscript too long”; the means equal those of the same `lmer` ML fit.

- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  with a deviation smooth `s(x, f1*f2)` no longer gives the warning
  “number of items to replace”.

- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  reports lme4 convergence warnings of the permuted refits once, as a
  count, instead of one warning per refit.

- Tukey p-values within `atvar` levels used the number of pairs instead
  of means.

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  failed for `lmer` models with a missing treatment combination, and for
  ML `lmer` fits.

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  with `permlist` and `atvar` failed; a user `Df` for `lmer` gave NA
  p-values.

- Factor levels containing `":"` broke the pairwise tables.

- Within-`atvar` adjusted p-values could be matched to the wrong pairs.

- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  failed with any missing treatment combination.

- Kenward-Roger ignored prior weights and offsets of `lmer` fits.

- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  0/1 responses had SEs 15-25% too small; Gamma and inverse-Gaussian SEs
  were too small.

- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  gaussian fits had wrong
  [`anova()`](https://rdrr.io/r/stats/anova.html) df and
  [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  SEs.

- [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  negative binomial SEs were scaled by theta; zero-inflated bands were
  too narrow.

- `semipred(contr = )` bands of generalized models were wrong.

- [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  SEs for `glmer` and limits for `glmmTMB` were wrong.

- [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  was wrong for some glmmTMB families, offsets and correlated random
  effects.

- [`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)
  with `block` or `group` returned indices for sorted data.

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  permuted only one column of
  [`cbind()`](https://rdrr.io/r/base/cbind.html) responses; type I
  `lmer` tables compared different statistics.

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  and
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  could refit with a data set of the same name found elsewhere.

- [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  failed for every `lmer` model.

### Documentation

- New vignettes: “Migrating from predictmeans 1.x”, “Checking the
  model”, “Missing cells and non-estimable means”, “Smooth curves in
  mixed models”, “Permutation and bootstrap tests”, “Link and response
  scale”, “Contrasts, joint tests and slopes”, and “Matching SAS PROC
  MIXED” (the SAS and Genstat reference runs).
- [`vignette("predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/predictmeans.md)
  rewritten, with a table of default df by model class.
- Help pages written with roxygen2;
  [`?advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)
  lists every advice id.
- A pkgdown website with worked-example articles: an overview, tree
  improvement, and adjusted R² for GLMMs reproducing Piepho (2023).

The full development log of 2.0.0 is in the package source, file
`dev/NEWS-2.0-development.md`.
