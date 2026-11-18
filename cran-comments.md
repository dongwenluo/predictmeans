## Status

Draft for the 2.0.0 submission; not submitted yet. Update the results below
when the release is made (set `Version: 2.0.0` in DESCRIPTION and the NEWS
heading first).

## Test environments

* GitHub Actions (`.github/workflows/R-CMD-check.yaml`, `R CMD check
  --as-cran`): macOS, Windows and Ubuntu with R release, Ubuntu with R devel;
  an Ubuntu job without `NOT_CRAN` runs the tests as CRAN does.
* Local: Ubuntu, R 4.3.3.
* Before submitting: win-builder (release and devel) and mac-builder.

## R CMD check results

0 errors | 0 warnings | 0 notes (to be confirmed on the builders above).

## Reverse dependencies

Checked with `.github/workflows/revdep-check.yaml` (`dev/revdep-check.R`),
which checks every CRAN package that depends on, imports or suggests
predictmeans with the CRAN version and with this version. Results to be
added here.

## Notes for CRAN

This is a major release that breaks some 1.x code on purpose, to give the
functions one vocabulary (`conf.level`, `df_method`, `at`, `plot`). Every
removed argument or function stops with a message giving the 2.0 spelling;
`vignette("migrating")` shows the 1.1.1 and 2.0 code side by side.

* Removed functions: `smZ()`, `permglmer()`, `permglmmTMB()`,
  `permanova.lmer()`, `ci_mcp()`, `R2_glmm()` (each stops with its
  replacement, `?"predictmeans-defunct"`).
* Defaults that change, each with the argument that restores the 1.x result:
  Kenward-Roger adjusted standard errors for REML mixed models
  (`vcov_adjust = FALSE`), degrees of freedom chosen by model class
  (`df_method = "residual"` or `"containment"`), Freedman-Lane permutations
  (`method = "raw"`), and plots drawn only when asked for (`plot = TRUE`).
* lmerTest moved from Depends to Suggests (the package computes its own
  Kenward-Roger and Satterthwaite df); mgcv and glmmTMB are in Imports.
