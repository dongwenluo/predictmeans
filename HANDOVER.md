# Handover: predictmeans 2.0.0 (state at 2026-11-24)

Read first: `CLAUDE.md` (the maintainer's rules), then this file. Look up
details on demand in row **O16** of `dev/predictmeans-2.0-plan.md` (one long
row: every change of the last sessions, in order) and in
`dev/NEWS-2.0-development.md` (the full development log). `NEWS.md` is now the
short release note for users; don't let it grow back into a log.

## Repository

- R package **predictmeans 2.0.0** (DESCRIPTION), ready for CRAN submission.
  Windows 11, R 4.6.2, Positron, renv (`renv.lock` includes Suggests:
  `renv::snapshot(dev = TRUE)`).
- git, branch `main`, **no remote yet**. Last commit `8b3f7da`; working tree
  clean except `dev/conversations/` (git-ignored chat exports; do not publish).
- Last checks (2026-11-24): full suite 3,853 expectations, 0 failures,
  0 warnings, 2 skipped (all Suggests installed); `R CMD check --as-cran`
  with remote incoming checks and the PDF manual: 0 errors, 0 warnings,
  0 notes; `urlchecker::url_check()` clean; pkgdown site built to `docs/`
  (git-ignored), every page checked.
- `AGENTS.md`, `CLAUDE.md` and `HANDOVER.md` are in the package root (Positron
  reads `AGENTS.md` as project memory). They are excluded from the package by
  `.Rbuildignore`, but pkgdown publishes every top-level `.md`, so build the
  site with `dev/build-site.R` (it deletes their pages); the GitHub pkgdown
  workflow has the same step. Never add other private notes to the root.

## Maintainer's working rules (short; full text in CLAUDE.md)

- Rule no. 1: correct, accurate, robust, effective, explained.
- Every change: a NEWS.md entry (user-facing, short, in the right section of
  the 2.0.0 note) and an O16 update. Related tests while working; full suite,
  `R CMD check` and one commit at the end of a block (the maintainer asks for
  commits).
- Tests: `options(warnPartialMatchDollar = FALSE)`; `Sys.setenv(NOT_CRAN =
  "true")`; one file with `testthat::test_file(path, package = "predictmeans",
  load_package = "source")`; full suite `devtools::test(reporter = "silent")`.
  Keep the global environment clean (tests can see stray objects).
- Results that should not change: snapshot before/after at tolerance 0 (scripts
  `dev/snapshot-*.R`, or a quick list of results in the session).
- Advice: only when triggered by this fit/data; every id in the `?advice`
  table (a test checks this).
- Argument names snake_case; result columns broom names.
- Help pages from roxygen in `R/*.R`, then `devtools::document()`; never edit
  NAMESPACE by hand.
- SAS/Genstat references only from real runs: scripts in `dev/sas/` with data
  inline (no `infile`); the maintainer runs them (SAS OnDemand) and uploads
  the output (PDF/CSV) to the project.
- Plots: the 1.x look; clean (no unnecessary legends or colour bars).
- Python is not available; edit with the edit tools or R.
- eda4treeR pulls CRAN predictmeans 1.1.1 into the library: run
  `devtools::install(dependencies = FALSE, upgrade = FALSE)` before rendering
  articles or building the site (`dev/build-site.R`).

## Open: next steps, in order

1. **Builders.** win-builder (R-release, R-devel) submitted 2026-11-24, results
   go by email to the maintainer (`dongwen.luo@agresearch.co.nz`); mac-builder
   returned HTTP 502 twice: resubmit `devtools::check_mac_release()`. Fill the
   "(pending)" lines in `cran-comments.md` from the emails.
2. **GitHub.** Create `github.com/dongwenluo/predictmeans`, push `main`; then
   the workflows in `.github/workflows/` run: R-CMD-check, pkgdown (publishes
   the site at https://dongwenluo.github.io/predictmeans/), benchmarks, and
   revdep-check (reverse dependencies; results go into `cran-comments.md`).
3. **Submit 2.0.0 to CRAN** when builders and revdeps are clean.
4. **Piepho note.** `dev/sas/piepho-table6-note.md` (Table 6 of Piepho 2023 is
   wrong for 12 of 16 rows; SAS run `dev/sas/piepho-table6.sas`,
   `piepho-table6-results.csv`; cause traced in his own code,
   `dev/sas/Omega.sas`, `dev/sas/Tables 5 & 6 foot_mv in one go.sas`). The
   maintainer checks signature/tone and sends it.
5. Optional: a SAS check of `se_ranef()` for a GLMM or correlated random slopes
   (the lmer prediction SE is checked: `dev/sas/se-ranef.sas`); `spelling`
   package not installed, so no spelling check yet.

## What was done in the last two sessions (2026-11-23/24), by area

Details and file names in O16; user-facing summary in NEWS.md.

- **Documentation**: help pages, vignettes, README, DESCRIPTION reviewed against
  the 2.0 API, with the package's aims and their test evidence. Articles
  (`vignettes/articles/`): overview (reworded, 2.0 calls), tree-improvement,
  new R2-GLMM (Piepho 2023 examples; data in `vignettes/articles/data/`); code
  folded by default; `pkgdown/extra.css` (tab heights, folded code).
- **pkgdown**: Bootswatch flatly, primary navy, menus (Get started, Articles,
  Worked examples, News), reference sections with descriptions, favicons.
- **predictmeans() and friends**: `permmodels()` says when `test_statistic` is
  replaced or unused.
- **residplot()**: observed values from `pm_data()`, axis = response as modelled;
  labels = original row numbers; gam shows fitted vs observed; many points as
  open transparent circles (`dense = "bins"` optional, no colour bar); no
  legends; Q-Q of many residuals thinned (tails kept), line from all quartiles.
- **CookD()**: ggplot2 plot (`pm_cookd_plot()`).
- **covariatemeans()/semipred() plots**: `autoplot.pm_curves()`
  (`R/autoplot-curves.R`): lines, surface, 3d (plotly, one widget), difference;
  `group`/`panel` for any number of factors; uses the call's `plot_args` as
  defaults (`attr(x, "plot_args")`); one complete legend; marks where a band
  leaves out 0/1 drawn as a rug (worked on log axes). `$pltdf` is `$pred_df`.
- **semipred()**: `sm_term` by label, term, covariate or block; advice
  `smooth_other_covariate`; with `contr` a max-|t| test of equal curves
  (`settings$test`; level 3.0% in 200 null simulations); back-transformation
  shared with `covariatemeans()` (`pm_cov_back()`: trans, link, formula
  transform; `bias_adjust`, `sigma`; smooths excluded from the bias variance).
- **semireg fits**: `R/semireg-methods.R`: short print with a smooth table
  (smooth, term, covariate, basis, basis_dim, variance, linear) and how to pass
  them to `semipred(sm_term = )`; `summary()` (list with `$model`), `tidy()`,
  `coef`, `fixef`, `vcov`, `fitted`, `residuals`, `logLik`, `nobs`, `formula`
  (as written), `family`, `VarCorr`, `ranef`; the fit keeps `call`.
- **se_ranef()**: rewritten (`R/se_ranef.R`, `R/ranef-cov.R`): lmer, glmer,
  glmmTMB, lme, semireg; any covariance structure; `type = "prediction"`
  (default; Henderson C22, SAS Std Err Pred, checked against a real SAS run on
  unbalanced oats to 1e-4) or `"conditional"` (1.x); `predict = TRUE` adds the
  fixed part; data frame output, `as.list()` gives the 1.x list.
- **R2_model()**: reproduces Table 2 of Piepho (2023) (test with fixture
  `tests/testthat/fixtures/Binary_color.RData`).
- **Release**: Version 2.0.0; NEWS.md rewritten as a concise release note;
  `cran-comments.md` updated; `semipred` example under 5 s (slow parts in
  `\donttest{}`); `pm_coef` example uses `registerS3method()`; `s(x, f1*f2)`
  mgcv warning muffled; `permlmer()` reports refit convergence warnings once.

## Earlier work

See `dev/NEWS-2.0-development.md` and the plan's status table (rows O1-O16).
Earlier plan files are in `.posit/assistant/plans/` (workspace) and
`~/.posit/assistant/plans/` (outside the workspace).
