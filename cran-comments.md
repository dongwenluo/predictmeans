## Status

Version 2.0.0, not submitted yet. Submitted to win-builder (R-release and
R-devel) on 2026-11-24; results come by email to the maintainer. mac-builder
returned "HTTP 502 Bad Gateway" twice on 2026-11-24; resubmit with
`devtools::check_mac_release()`. Fill in the results below before submitting.

## Test environments

* Local: Windows 11 x64, R 4.6.2: `R CMD check --as-cran` with the remote
  CRAN incoming checks and the PDF manual, tests run as on CRAN
  (`NOT_CRAN=false`).
* win-builder: R-release and R-devel (results pending).
* mac-builder: R-release (to be resubmitted).
* GitHub Actions (`.github/workflows/R-CMD-check.yaml`, `--as-cran`): macOS,
  Windows and Ubuntu with R-release, Ubuntu with R-devel, once the repository
  is on GitHub.

## R CMD check results

Local (Windows, R 4.6.2): 0 errors | 0 warnings | 0 notes.

win-builder R-release: (pending)
win-builder R-devel (2026-11-17 r90926 ucrt): 1 NOTE (CRAN incoming
feasibility): the maintainer, and "Possibly misspelled words in DESCRIPTION:
Kenward, PROC". Both are names: Kenward of Kenward and Roger (1997), cited
with its DOI, and the 'SAS' procedure PROC MIXED.
mac-builder R-release: (pending)

The full test suite (3,853 expectations, `NOT_CRAN=true`) passes with 0
failures and 0 warnings. All URLs checked with `urlchecker::url_check()`.

## Reverse dependencies

One reverse dependency, eda4treeR 1.1.0, checked with
`.github/workflows/revdep-check.yaml` (`dev/revdep-check.R`) on 2026-11-24
(R 4.6.2, Ubuntu): `R CMD check` status OK with predictmeans 1.1.1 (CRAN) and
OK with 2.0.0. No new problems.

## Notes for CRAN

This is a major release that breaks some 1.x code on purpose, to give the
functions one vocabulary (snake_case arguments such as `conf_level`,
`df_method`, `at`). Every removed argument or function stops with a message
giving its 2.0 spelling; `vignette("migrating")` shows the 1.1.1 and 2.0 code
side by side.

* Removed functions: `smZ()`, `permglmer()`, `permglmmTMB()`,
  `permanova.lmer()`, `ci_mcp()`, `R2_glmm()` (each stops with its
  replacement, `?"predictmeans-defunct"`).
* Defaults that change, each with the argument that restores the 1.x result:
  Kenward-Roger adjusted standard errors for REML mixed models
  (`vcov_adjust = FALSE`), degrees of freedom chosen by model class
  (`df_method = "residual"`), Freedman-Lane permutations (`method = "raw"`),
  prediction-error standard errors in `se_ranef()` (`type = "conditional"`).
  `predictmeans()` no longer draws plots; `autoplot()` draws them from the
  result.
* lmerTest moved from Depends to Suggests (the package computes its own
  Kenward-Roger and Satterthwaite df); mgcv and glmmTMB are in Imports.
* Examples that take longer than a few seconds (simulated bands, refits) are in
  `\donttest{}`.
