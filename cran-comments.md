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

This is a major release. Some defaults change, each documented in NEWS with
the argument that restores the 1.x result: Kenward-Roger adjusted standard
errors for REML mixed models (`vcov_adjust = FALSE`), degrees of freedom
chosen by model class (`df_method = "residual"` or `"containment"`), and
Freedman-Lane permutations (`method = "raw"`). The 1.x functions and their
printed output are kept.
