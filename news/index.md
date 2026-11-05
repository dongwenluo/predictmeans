# Changelog

## predictmeans 1.1.1.9000 (development version)

### Changes that affect results or scripts

- **Kenward-Roger adjusted standard errors are now the default** for
  `lmer` models fitted by REML, in
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  and
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  (`vcov_adjust = TRUE`). This matches SAS `DDFM=KR` and emmeans.
  Standard errors, confidence limits, SEDs, LSDs and p-values can differ
  slightly from predictmeans 1.x (e.g. 7.896 instead of 7.887 for an
  unbalanced Oats example); balanced designs barely change. Use
  `vcov_adjust = FALSE` to reproduce 1.x results.
- **glmmTMB is no longer attached** by
  [`library(predictmeans)`](https://dongwenluo.github.io/predictmeans/);
  it moved from Depends to Imports. Scripts that fit models with
  `glmmTMB()` need
  [`library(glmmTMB)`](https://github.com/glmmTMB/glmmTMB). lme4, nlme
  and lmerTest are still attached.
- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  gives the same tables (including Satterthwaite df and p-values)
  whether an lmer model was fitted with
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) or
  [`lmerTest::lmer()`](https://rdrr.io/pkg/lmerTest/man/lmer.html).

### Stable result tables

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  results gain three slots with fixed names and unrounded values.
  `tidy()` (generics/broom) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
  them:
  - `means`: one row per mean, with `estimate`, `std.error`, `df`,
    `conf.low`, `conf.high`, back-transformed `response` columns,
    `group` letters, `n` and `estimable`;
  - `pairs`: one row per comparison, with `level1`, `level2`, `by` (for
    `atvar`), `estimate`, `std.error`, `df`, `statistic`, `p.value`,
    `p.adj`, `lsd` and `p.method`;
  - `settings`. They are checked against emmeans in the tests. The
    legacy slots and the printed output are unchanged, and the result
    has class `c("predictmeans", "pdmlist")`.

### Plots and summaries from the result tables

- [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  for
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  results gives ggplot objects that work in Shiny, reports and
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html):
  - `type = "means"`, with confidence limits, half-LSD bars or no bars;
  - `type = "pvalue"`, a heatmap of the adjusted pairwise p-values in
    the colour bands of
    [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md),
    one panel per `atvar` level;
  - `type = "response"`, the back-transformed means.
- [`summary()`](https://rdrr.io/r/base/summary.html) prints the means
  with a footer stating the model and df method, the SED and LSD ranges,
  the p-value method and adjustment, what the letters mean, and the
  covariate values used.
- Both read only the `means` and `pairs` tables: they compute nothing
  new, open no graphics device and print nothing unless printed.

### Supporting new model classes

- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  and
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  now get what they need from a model through six exported generics:
  - [`pm_coef()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md),
    [`pm_vcov()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md):
    estimates and their covariance;
  - [`pm_basis()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md):
    design matrix for new data;
  - [`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md):
    the fitting data;
  - [`pm_df()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md):
    df for each contrast;
  - [`pm_linkinv()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md):
    inverse link.

  The defaults work for any class with the usual `coef`, `vcov`,
  `terms`, `df.residual` and `family` methods. Supporting another class
  means writing methods for these generics only (see
  [`?pm_coef`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)),
  and model classes predictmeans doesn’t know are no longer rejected.
  Built-in classes behave exactly as before.
- Ordered factors keep their contrasts in the reference grid even when a
  model class stores no contrasts.

### Named contrasts

- New
  [`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md)
  writes contrasts with level names instead of counted columns, e.g.
  `specify_contrast(Kmatrix(fm, "A:B"), a = (`x:b1`+`y:b1`)/2 -`z:b1`)`.
  The result goes straight into `contrastmeans(ctrmatrix = )`. It checks
  that the weights sum to zero (unless `.zero_sum = FALSE`) and that no
  missing or non-estimable mean is used.
- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md):
  - accepts a contrast matrix whose columns are named by level, in any
    order and for a subset of levels;
  - keeps the contrast names in its result;
  - gains `weights`.

### New options

- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  and
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  gain `weights = c("equal", "proportional", "outer", "cells", "flat")`,
  which sets how the cells of the factors not in `modelterm` are
  weighted. `"proportional"` uses the observed frequencies, the same for
  every mean, like SAS `LSMEANS / OM`. The definitions are emmeans’, and
  the results match emmeans. The default `"equal"` is unchanged.

- `predictmeans(vcov_adjust = TRUE)` and
  `contrastmeans(vcov_adjust = TRUE)`: for `lmer` models fitted by REML,
  standard errors use the Kenward-Roger adjusted covariance (as SAS
  `DDFM=KR` and emmeans do) to go with the Kenward-Roger df. The default
  keeps the earlier, unadjusted standard errors.

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  and
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  pass `cov_reduce` and `exclude` on to
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).
  For example, `cov_reduce = "average"` gives `gam` means adjusted for a
  spatial smooth.

### Bug fixes

- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  failed with “non-conformable arguments” for every model with a missing
  treatment combination, even for contrasts that do not involve it. It
  now works, and stops with a clear message only when a contrast uses a
  mean that cannot be estimated.
- Tukey-adjusted p-values within `atvar` levels used the number of pairs
  instead of the number of means. With 4 means per level, Tukey was
  computed for 6 means, so p-values were too large (0.350 instead of
  0.197 in the Oats example). For `lmer` models they also used one
  averaged df instead of each pair’s df. Both now match emmeans.
- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  failed for every `lmer` model with a missing treatment combination
  (“replacement has 12 rows, data has 11”).
- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  failed for `lmer` models fitted by ML (`REML = FALSE`). Satterthwaite
  df are now used for such fits, with a message;
  [`df_term()`](https://dongwenluo.github.io/predictmeans/reference/df_term.md)
  falls back to Satterthwaite with a warning instead of stopping.
- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  and
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  refitted `lme`/`gls` models by re-evaluating the model call from
  inside the package. Data fitted inside a function were not found, and
  a data set with the same name elsewhere on the search path could be
  used silently. For example, a local `Oats` with `nitro` as a factor
  was replaced by
  [`nlme::Oats`](https://rdrr.io/pkg/nlme/man/Oats.html), so the
  permutation tests were for a different model. Refits now use the
  model’s own data and check that they reproduce its coefficients.
- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  and
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  now find the data of `gls` models fitted inside a function;
  [`nlme::getData()`](https://rdrr.io/pkg/nlme/man/getData.html) does
  not.
- [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  failed for every `lmer` model: it refitted with ML and then asked for
  Kenward-Roger df, which needs REML. It now works and agrees with
  `lme`.
- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  failed for `gls` models without a correlation structure.
- [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md):
  - for `glmmTMB` models it now finds data that are not in the global
    environment;
  - `over_disp = TRUE` gives a clear message when the model has no
    observation-level random effect.
- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md):
  - passing two or more of the deprecated `mplot`, `barplot`, `pplot`,
    `bkplot` arguments gave an error instead of a deprecation warning;
  - `trans = NULL` gave an error.
- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  with `permlist` counted permuted statistics strictly greater than the
  observed one. It now counts “greater or equal” (after rounding), as
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  does, so both give the same permutation p-value for the same contrast.
- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  and
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  no longer leave `options(scipen = 6)` set after they return; printed
  results look the same.
- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  without `ctrmatrix` only opens the contrast editor in an interactive
  session.

### Speed

- Satterthwaite df (used for `lmer` models fitted by ML, and by
  `df_term(type = "Satterthwaite")`) are computed for all contrasts at
  once from ingredients computed once per model. They equal lmerTest’s
  values, about 25x faster for many contrasts (780 pairs: 0.8 s -\> 0.03
  s), and no longer convert the model to lmerTest’s class.

- Kenward-Roger computations no longer form or invert the n x n
  covariance of the data. They use the model’s low-rank structure
  (Woodbury identity), with the same results as pbkrtest. n = 3,000
  takes under 0.05 s instead of 21.5 s; n = 20,000 takes 0.3 s
  (previously out of reach).

- Kenward-Roger df are computed for all contrasts at once (identical to
  [`pbkrtest::Lb_ddf()`](https://rdrr.io/pkg/pbkrtest/man/get_ddf_Lb.html)
  row by row), and the adjusted covariance is computed once per
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  call instead of once per model term, for the means, and again for the
  pairs.
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  on `lmer` models is about 10-50x faster: Oats split-plot 0.96 s -\>
  0.09 s; 40 means with 780 pairwise comparisons 5.3 s -\> 0.1 s.

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  no longer computes Kenward-Roger df for every permutation refit (the
  df were discarded). About 10x faster for `lmer` models (99
  permutations of the Oats split-plot: 21.8 s -\> 2.2 s), with identical
  permutation p-values.

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  and
  [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  use a sparse Cholesky factor and triangular solves instead of a dense
  n x n inverse.

- [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md)
  computes its traces without n x n matrices. For n = 4,000 it goes from
  4.4 s and 529 MB to 0.01 s with identical results; n = 100,000 now
  takes 0.02 s.

### Non-estimable means

- Means that cannot be estimated (they depend on aliased coefficients)
  are now NA in every table, with a message naming them; their pairwise
  differences are still shown when estimable. Before, a number that
  depended on the arbitrary choice of aliased coefficient was printed.
  Means over missing cells keep the 1.x rule (average of the observed
  cells, with the “Missing treatments’ combination” warning).

### Containment df

- New argument `df_method` in
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  and
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  (and `type = "containment"` in
  [`df_term()`](https://dongwenluo.github.io/predictmeans/reference/df_term.md))
  for lmer models: `"auto"` (the default, unchanged: Kenward-Roger for
  REML fits, Satterthwaite for ML fits), `"Kenward-Roger"`,
  `"Satterthwaite"` or `"containment"`. Containment gives the SAS
  `DDFM=CONTAIN` df (the smallest rank contribution to \[X Z\] of the
  random effects containing the term, else the residual df), with the
  unadjusted covariance unless `vcov_adjust = TRUE` is given. For the
  oats split plot this is 10 for Variety and 45 for Nitro and
  Variety:Nitro; with a whole plot missing, 9 and 42. The ranks use M’M,
  so a split plot with n = 6000 and 750 whole plots takes about 3 s.

### Bug fixes found by the refactor

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  with both `permlist` and `atvar` failed with “undefined columns
  selected”. It now gives permutation p-values within each `atvar`
  level.
- Factor levels containing `":"` (e.g. `"Golden:Rain"`) broke the legacy
  pairwise tables: `atvar` failed and the within-factor SEDs and LSDs
  were computed from misread labels (B4). Level labels are now taken
  from the label table, never split from row names.
- A user-supplied `Df` for an lmer model gave NA p-values for all pairs
  but the first (“missing value where TRUE/FALSE needed” when letters
  were computed). The given `Df` is now used for every pair.

### Documentation

- The deprecated plot arguments (`mplot`, `barplot`/`bar`, `pplot`,
  `bkplot`) turned off every plot not given, so `bar = TRUE` gave a NULL
  mean plot. Plots not given now keep their `which.plots` setting, and
  the deprecation warning shows the equivalent `which.plots` call.

- Fixed (regression in this development version): with nested treatments
  and one observation per block cell (e.g. `treat/contcomp` with
  `repl`), means were averaged over structurally empty cells and were
  wrong. Empty cells are now judged on the factors that interact only.

- `atvar` with several factors failed (“Names required”) when a group of
  `atvar` levels had only one or no observed level of the other factors,
  or (“replacement has 8 rows, data has 12”) when a group had one level.

- Within-`atvar` comparisons: adjusted p-values could be assigned to the
  wrong group or pair when the groups were not in row order (for example
  two `atvar` factors, or empty groups). This came from 1.x. P-values
  and letters are now matched to groups and levels by name.

- Three vignettes: “Getting started with predictmeans”, “Adjusted R² for
  generalized linear mixed models” and “Migrating from predictmeans
  1.x”.

- A pkgdown website (built by GitHub Actions and published to GitHub
  Pages), with two longer worked-example articles: an overview of the
  package and experimental design and analysis for tree improvement.

- [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  works with glmmTMB \>= 1.1.12, which no longer stores `Lambdat`,
  `Lind` and `lower` (they are rebuilt when missing).

- [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md)
  failed for gaussian glmmTMB fits (including
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md))
  with glmmTMB \>= 1.1.12, whose dispersion parameter is named
  `betadisp`; the fixed effects now come from `fixef(model)$cond`.

- New data sets `Body_length` and `Poisson_ind` (used in the R²
  vignette).

### Internal

- Help pages are now written with roxygen2: the sources are the
  `R/doc-*.R` files, and `man/` is generated with
  `roxygen2::roxygenise(roclets = "rd")`. NAMESPACE is still maintained
  by hand. The rendered pages are unchanged.

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  is split into stages (setup, estimate, compare, mean table,
  back-transform, plots, assemble). All statistics are computed before
  any plotting, so `plot = FALSE` skips all graphics work. The work of
  each stage is done by functions that take their inputs as arguments
  and return their results, so they can be tested and changed one at a
  time. Results are unchanged.

### Packaging

- `which.plots` is documented. DESCRIPTION no longer carries
  CRAN-generated fields, and `NeedsCompilation` is “no”.
- GitHub Actions R CMD check on Linux, macOS and Windows.

### Kmatrix() and semipred()

- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  now builds predicted means from a reference grid evaluated with the
  model’s own basis information. It handles covariates that:
  - enter through transformations:
    [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`, `bs()`,
    [`log()`](https://rdrr.io/r/base/Log.html),
    [`scale()`](https://rdrr.io/r/base/scale.html),
    [`I()`](https://rdrr.io/r/base/AsIs.html);
  - appear only inside interactions (`A:x` without `x`, `x:z`);
  - are mgcv smooth terms (`s()`, `te()`, `s(x, by = A)`). These cases
    previously stopped with “subscript out of bounds” or “object not
    found”.
- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  and
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  now support
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html)/`bam()` models,
  using the lpmatrix of the smooths. Random-effect smooths (`bs = "re"`,
  `bs = "fs"`) are excluded by default (`exclude` argument).
- New
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  argument `cov_reduce`: fix covariates with any summary function
  (default `mean`, e.g. `median`), or use `"average"` to average the
  predictions over the observed covariate values.
- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  now also handles:
  - character and logical predictors;
  - `factor(x)` terms in the formula (`modelterm` may be `"x"` or
    `"factor(x)"`);
  - models fitted without a `data` argument;
  - a transformed response together with `data =`.
- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  returns the reference grid (`refgrid`) and flags non-estimable rows in
  `attr(K, "estimable")`.
- A matrix `covariateV` for several covariates is now used as a set of
  joint prediction points, as documented; previously it was ignored.
- [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  uses
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  on the `semireg` object directly, without refitting an auxiliary
  [`lm()`](https://rdrr.io/r/stats/lm.html). The smooth basis is
  averaged over the same reference grid as the fixed part. As a result,
  factors not in `modelterm`, including `by` factors of the smooths, are
  averaged correctly instead of causing an error.
- The covariate summary printed by `prtnum = TRUE` now goes through
  [`message()`](https://rdrr.io/r/base/message.html), so it can be
  suppressed.
- Kmatrix() is faster on large data (about 4-15x on a 20,000-row,
  360-cell model).
- New tests check that results for every case supported by 1.1.1 are
  unchanged, and compare the new cases with emmeans and mgcv.

### Fewer dependencies

- `pbkrtest` moved from Imports to Suggests. Kenward-Roger adjusted
  covariances, df and F tests are computed by a port of pbkrtest’s
  `vcovAdj()` and `KRmodcomp()` that gives the same results in the
  tests. pbkrtest alone brought 11 tidyverse packages (dplyr, broom,
  tidyr, …). Installing predictmeans now pulls in 56 packages (112 for
  1.1.1).

- `car` moved from Imports to Suggests.
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  computes its Type II and III tests itself: lm F; glm LR, Wald and F;
  lme and gls Wald chi-square; lmer Wald chi-square and Kenward-Roger F.
  They are ported from car and equal
  [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) in the tests.
  Installing predictmeans now pulls in 68 packages (112 for 1.1.1).

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  with `type = "II"` now also works for `lm` models with aliased
  coefficients. 1.1.1 stopped with an error from
  [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html).

- `plyr` is no longer used.

- `HRW`, `lmeInfo`, `lmeSplines`, `plotly` and `splines2` moved from
  Imports to Suggests. They are only needed for, respectively:

  - HRW: example data;
  - lmeInfo:
    [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
    on `lme` models;
  - lmeSplines: `smZ(type = "smspline")`;
  - plotly: `semipred(threeD = TRUE)`;
  - splines2: the splines2 spline types in
    [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/smZ.md).
    A clear message asks for the package when such a feature is used.
    Installing predictmeans now pulls in 76 packages instead of 112.
    This also avoids current HRW’s requirement of R \>= 4.5.

- [`HRW::pointsInPoly()`](https://rdrr.io/pkg/HRW/man/pointsInPoly.html)
  is replaced by an internal copy. It also closes open polygons
  correctly; HRW’s version adds a spurious vertex when the last vertex
  does not repeat the first.
