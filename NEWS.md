# predictmeans (development version)

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

## New options

* `predictmeans(vcov_adjust = TRUE)` and `contrastmeans(vcov_adjust = TRUE)`:
  for `lmer` models fitted by REML, standard errors use the Kenward-Roger
  adjusted covariance (as SAS `DDFM=KR` and emmeans do) to go with the
  Kenward-Roger df. The default keeps the earlier, unadjusted standard errors.
* `predictmeans()` and `covariatemeans()` pass `cov_reduce` and `exclude` on
  to `Kmatrix()`. For example, `cov_reduce = "average"` gives `gam` means
  adjusted for a spatial smooth.

## Bug fixes

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

## Speed

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
* `R2_glmm()` computes its traces without n x n matrices. For n = 4,000 it
  goes from 4.4 s and 529 MB to 0.01 s with identical results; n = 100,000
  now takes 0.02 s.

## Packaging

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
