# predictmeans — notes for Claude Code

R package for predicted means, pairwise comparisons, permutation tests
and semiparametric models. Development is toward 2.0; read
`dev/predictmeans-2.0-plan.md` first (aims, decisions, bug list, status
table; the row marked “Next” is the current task) and `NEWS.md`.

## Aims (from the maintainer)

- Statistically right first. Each default is the most reliable method on
  the evidence (simulation, literature); matching SAS or Genstat is an
  option the user can choose, not the goal. SAS/Genstat runs check that
  a formula is computed as documented.
- Clear logic, not chaos: each statistical decision (df method,
  covariance, estimability, …) is made in one place and used by every
  function.
- Run fast; add emmeans/marginaleffects-style features; keep
  dependencies few; easy to extend to new model classes (adapter
  generics in `R/adapter.R`, see
  [`?pm_coef`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)).

## Decisions already made

- Kenward–Roger adjusted SEs are the default for REML lmer fits;
  `vcov_adjust = FALSE` reproduces 1.x numbers.
- glmmTMB and mgcv in Imports; lmerTest, car, pbkrtest and TMB in
  Suggests (own code in `R/anova.R`, `R/kr.R`, `R/satterthwaite.R`,
  `R/df-containment.R`; lmerTest dropped from Imports on 2026-11-17, the
  maintainer’s decision:
  [`anova.lmerMod()`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md)
  gives F tests with KR df for lmer fits,
  [`anova.semireg()`](https://dongwenluo.github.io/predictmeans/reference/anova.semireg.md)
  tests fixed effects and smooth terms; gaussian
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fits are class `semerMod`). Internal code that needs lme4’s own
  sequential table calls `pm_anova_seq()`, not
  [`anova()`](https://rdrr.io/r/stats/anova.html).
- Kenward-Roger leaves variance components estimated at zero out of the
  adjustment (as SAS); `R/kr.R`.
- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  draws plots only when asked (`plot = FALSE` default since 2.0;
  `which.plots` etc. turn plotting on); `autoplot()` draws from results.
- Means over empty cells keep the 1.x rule: the average of the observed
  cells, with a warning (the maintainer’s decision of 2026-11-13,
  replacing an earlier “NA by default”). Combinations absent by design
  (nested factors) are not empty cells. Truly non-estimable means are
  NA. The cells are those of each interaction term of the model. Factors
  that do not interact but are confounded in the data (`y ~ G + T` with
  T nested in G, the spray trial, each pasture on one farm) have no mean
  over all their combinations: the means then average over the
  combinations that exist, with a warning that names the two factors and
  the three ways to say it (`G/T` in the formula, leave out the aliased
  factor, `weights = "cells"`) (decision of 2026-11-16, replacing “NA
  with a note” of 2026-11-15, after the maintainer’s pasture data went
  from 0 to 60 NA pairs). The pair is found by trying
  ([`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
  `pm_confounded_pairs()`): pairs among the factors averaged over first.
- `aov` models with `Error()` (aovlist) are analysed by strata, as
  Genstat: standard errors from the stratum mean squares, never by
  refitting as lmer. Missing values in an orthogonal design are
  estimated as Genstat’s ANOVA does (exact SEDs), with a message
  recommending REML (decision of 2026-11-14).
- emmeans (and SAS) are not the reference for nested treatment
  structures or missing combinations (the EOD tree-improvement
  examples): there predictmeans averages with equal weights over the
  combinations that exist (`tests/testthat/test-nested.R`), where
  emmeans weights by the nesting hierarchy or gives NA.
- `R/inference.R` (`pm_inference()`) decides covariance, df and their
  labels for every front end; change a df or covariance rule there only.
- `R/p-adjust.R` decides multiplicity for every front end: families and
  their sizes (`pm_families()`), adjusted p-values (`pm_adjust_p()`) and
  the simultaneous limits that go with them (`pm_adjust_crit()`). In
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  every statistic of the pairs is computed once, in `pm_pair_table()`;
  the `pairs` table, the printed matrices and the letters read it.
- The maintainer can run SAS and Genstat on request: scripts go in
  `dev/sas/` and `dev/genstat/`. Reference values come only from real
  runs (uploaded results pages, output or CSV files). Do not use SAS
  tables typed into chat as references (see `dev/sas/README.md`).

## Layout

- `R/predictmeans.R` calls the stages in `R/predictmeans-stages.R`,
  ordinary functions that each take the results of those before and
  return a named list: `opt <- pm_setup()`, `est <- pm_estimate(opt)`,
  `cmp <- pm_compare(opt, est)`, `tab <- pm_tables(opt, est, cmp)`;
  `pm_flat()` joins them for `pm_make_plots()` and `pm_assemble()`. Plot
  settings go through `plot_args` (`pm_plot_args()`).
- `R/Kmatrix.R` + `R/refgrid.R`: reference grid and K matrix.
- `R/tidy.R`, `R/autoplot.R`: stable `means`/`pairs`/`settings` tables.
- `R/semireg-engine.R` (one fitting engine of
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)/[`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md))
  and `R/semireg-smooth.R` (smooths in the formula via mgcv; own
  `os`/`ztps` bases as mgcv smooth classes); plan in
  `dev/review-semireg.md`. Golden examples (the maintainer’s lme
  scripts, reference runs) in `dev/semireg-golden/`, data in
  `inst/extdata/semireg/`, checked by
  `tests/testthat/test-semireg-golden.R`.
- Help pages: edit the roxygen2 block above each function (data sets in
  `R/data.R`, the package page in `R/predictmeans-package.R`), then
  `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`. NAMESPACE is
  hand-written — edit it directly.

## Checks

- Tests: `NOT_CRAN=true Rscript -e 'devtools::test()'` (the 1.1.1
  reference tests only run with NOT_CRAN). One file:
  `testthat::test_file(<file>, package = "predictmeans", load_package = "source")`
  (with
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html) +
  `test_file()` tests can’t see unregistered internal methods).
  Mixed-model reference values are compared to 2% because lme4/nlme
  versions move them; everything else to 1e-6.
- Refactors must not change output: run
  `Rscript dev/snapshot-predictmeans.R <lib> before.rds` on the old
  build and `... after.rds` on the new one, then compare with
  `all.equal(readRDS("before.rds"), readRDS("after.rds"), tolerance = 0)`.
- Timings: `Rscript dev/benchmarks.R`. Pull requests are timed against
  their base branch in CI (`dev/benchmarks-ci.R`; fails when a key
  timing is over 1.5 times the base’s and at least 0.1 s longer).
- Run `R CMD check` once, just before opening a pull request (CI runs it
  on the PR as well); during the work run only the related tests.
- Add a NEWS entry and update the plan’s status table with each change.

## Working style the maintainer asked for

- One step per message; run only the tests related to the change, and R
  CMD check once just before a pull request.
- Keep explanations short and plain.
