# predictmeans — notes for Claude Code

R package for predicted means, pairwise comparisons, permutation tests and
semiparametric models. Development is toward 2.0; read
`dev/predictmeans-2.0-plan.md` first (aims, decisions, bug list, status
table; the row marked "Next" is the current task) and `NEWS.md`.

## Aims (from the maintainer)
- Statistically right first. Each default is the most reliable method on the
  evidence (simulation, literature); matching SAS or Genstat is an option the
  user can choose, not the goal. SAS/Genstat runs check that a formula is
  computed as documented.
- Clear logic, not chaos: each statistical decision (df method, covariance,
  estimability, ...) is made in one place and used by every function.
- Run fast; add emmeans/marginaleffects-style features; keep dependencies few;
  easy to extend to new model classes (adapter generics in `R/adapter.R`, see
  `?pm_coef`).

## Decisions already made
- Kenward–Roger adjusted SEs are the default for REML lmer fits;
  `vcov_adjust = FALSE` reproduces 1.x numbers.
- glmmTMB in Imports; car and pbkrtest in Suggests (own code in `R/anova.R`,
  `R/kr.R`, `R/satterthwaite.R`, `R/df-containment.R`).
- Means over empty cells keep the 1.x rule (average of observed cells, with a
  warning); truly non-estimable means are NA.
- The maintainer can run SAS and Genstat on request: scripts go in
  `dev/sas/` and `dev/genstat/`. Reference values come only from real runs
  (uploaded results pages, output or CSV files). Do not use SAS tables typed into chat as references
  (see `dev/sas/README.md`).

## Layout
- `R/predictmeans.R` calls the stages in `R/predictmeans-stages.R`; the work is
  done by pure functions there (`pm_estimate_means()`, `pm_compare_*()`,
  `pm_make_plots()`, ...). Only the setup stage edits the call's state.
- `R/Kmatrix.R` + `R/refgrid.R`: reference grid and K matrix.
- `R/tidy.R`, `R/autoplot.R`: stable `means`/`pairs`/`settings` tables.
- Help pages: edit the roxygen2 block above each function (data sets in
  `R/data.R`, the package page in `R/predictmeans-package.R`), then
  `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`. NAMESPACE is
  hand-written — edit it directly.

## Checks
- Tests: `NOT_CRAN=true Rscript -e 'devtools::test()'` (the 1.1.1 reference
  tests only run with NOT_CRAN). One file: `testthat::test_file(<file>, package =
  "predictmeans", load_package = "source")` (with `pkgload::load_all()` +
  `test_file()` tests can't see unregistered internal methods). Mixed-model reference values are compared to
  2% because lme4/nlme versions move them; everything else to 1e-6.
- Refactors must not change output: run
  `Rscript dev/snapshot-predictmeans.R <lib> before.rds` on the old build and
  `... after.rds` on the new one, then compare with
  `all.equal(readRDS("before.rds"), readRDS("after.rds"), tolerance = 0)`.
- Timings: `Rscript dev/benchmarks.R`.
- Run `R CMD check` once at the end of a piece of work, not after every edit.
- Add a NEWS entry and update the plan's status table with each change.

## Working style the maintainer asked for
- One step per message; run only the tests related to the change, and
  R CMD check once at the end.
- Keep explanations short and plain.
