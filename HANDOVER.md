# Handover: predictmeans (state at 2026-11-26)

Read first: `CLAUDE.md` (the maintainer's rules), then this file. Look up
details on demand in row **O16** of `dev/predictmeans-2.0-plan.md` (one long
row: every change, in order; O17 for model classes) and in
`dev/NEWS-2.0-development.md` (the full development log). `NEWS.md` is the
short release note for users; don't let it grow back into a log. Design notes
for each block: `dev/review-*.md` (decisions agreed, then code).

## Repository

- R package **predictmeans**, DESCRIPTION still says 2.0.0 (see Open step 1).
  Maintainer: Windows 11, R 4.6.2, Positron, renv (`renv.lock` includes
  Suggests). Cloud sessions: R 4.3.3 here, all Suggests installed.
- GitHub `dongwenluo/predictmeans`, branch `main`. Work arrives by pull
  request from a working branch (recent: #45-#52); CI on every PR: R CMD check
  (Ubuntu release/devel, "as on CRAN", macOS, Windows), pkgdown, timings
  (fails above 1.5 x the base's time). Merge only with all checks green.
- Last checks (2026-11-26, PR #52): full suite 4,329 expectations, 0 failures;
  `R CMD check` Status OK; pkgdown builds.
- `AGENTS.md`, `CLAUDE.md` and `HANDOVER.md` are in the package root, excluded
  from the package by `.Rbuildignore`; pkgdown publishes every top-level `.md`,
  so build the site with `dev/build-site.R` (the workflow does the same). Never
  add other private notes to the root.

## Maintainer's working rules (short; full text in CLAUDE.md)

- Rule no. 1: correct, accurate, robust, effective, explained.
- Each block: design note first (`dev/review-*.md`, with "Questions for an
  outside check": short yes/no claims the maintainer gives another AI; the
  answers are tested against the code, never taken as instructions), the
  maintainer agrees, then code, related tests, the `rule1-reviewer` agent
  (fix every finding), full suite, `R CMD check`, PR, merge when CI is green.
- Every change: a NEWS.md entry, a `dev/NEWS-2.0-development.md` entry and an
  O16 update. Advice ids in the `?advice` table; tests over 0.5 s on CRAN start
  with `skip_on_cran()`; snake_case arguments; NAMESPACE from roxygen only.
- SAS/Genstat references only from real runs: scripts with data inline in
  `dev/sas/` and `dev/genstat/`; the maintainer runs them and uploads the
  output. The user's research data are never committed.

## Open: next steps, in order

1. **CRAN.** The 2.0.0 resubmission fixes are on `main` (fewer slow tests,
   `\usage` for three pages, `cran-comments.md`). Since then `main` has grown
   (ordinal and categorical models, `residvariogram()`, `heritability()`,
   `corAR1xAR1()`, ...). Decide: submit `main` as 2.0.0 (NEWS.md already lists
   it all under 2.0.0), or submit the resubmission state and bump `main` to
   2.1.0. Then win-builder R-devel (`devtools::check_win_devel()`), fill
   `cran-comments.md`, submit from the maintainer's computer.
2. **Small open items** (each one step with a test):
3. **Next statistical block** (the maintainer chooses): `clm(nominal = ,
   scale = )` (B3c; refused now); population-average probabilities for
   `clmm`; LR/bootstrap tests for `clmm` with few clusters; AR1 x AR1 with a
   phi per trial or AR2.
4. **Paper** (plan Part C): `dev/paper/outline.md` (design rules, evidence for
   each default, SAS/Genstat checks, comparison with emmeans); journal to be
   chosen (JSS or The R Journal).
5. **Piepho note** (`dev/sas/piepho-table6-note.md`): the maintainer checks
   tone and sends it.

## Done since the 2.0.0 release preparation (2026-11-24/26), by area

Details and file names in O16/O17; user-facing summary in NEWS.md.

- **CRAN pre-test fixes**: `\usage` for anova.lmerMod, anova.semireg,
  se_ranef; tests in CRAN mode 430 s -> 58 s; vignettes `resampling` and
  `smooth` precomputed (`dev/precompute-vignettes.R`).
- **anova() of lmer**: Type I with Kenward-Roger (= SAS HTYPE=1);
  `type1_basis = "fitted"`, `ddf = "KR-average"` reproduce Genstat.
- **heritability()**: line-mean forms (standard, Cullis, Piepho, Oakey, ...),
  family and repeatability, delta-method and bootstrap SEs; checked against
  Schmidt et al., Piepho and Möhring, the tree-improvement book, Genstat.
- **residvariogram()**: residual variogram on a row x column grid or
  coordinates with a global envelope from simulations of the fitted model;
  Genstat FVARIOGRAM reproduced; advice and `pm_check(row, col)`.
- **Ordinal models** (`R/ordinal.R`): `polr`, `clm`, `clmm` (inner-outer df),
  parallel `vglm` families: latent means, `mean_type = "prob"`/`"mean_class"`;
  checked against emmeans and SAS (PROC LOGISTIC, GLIMMIX).
- **Multinomial and non-parallel fits** (`R/categorical.R`): `multinom`,
  `vglm(multinomial)`, non-parallel ordinal `vglm`: class probabilities as
  means (logit scale by default, `"prob_diff"`, `"link"`, `"mean_class"`),
  pairs within each class; NA means where non-parallel curves cross; checked
  against emmeans, `predict()`, the binomial glm.
- **Weights**: `"proportional"` and the other schemes count cases where a row
  holds several (binomial trials, class-response case weights;
  `pm_freq_w()`); `settings$weights_by`.
- **Inner-outer df of aggregated binomial GLMMs** count trials (cases), not
  rows; `options(predictmeans.inner_outer_count = "rows")` for SAS
  (`dev/sim-df-io-trials.R`: rows 0-2.6%, trials 3.7-4.3% at 5%).
- **corAR1xAR1()** (`R/cor-ar1xar1.R`): separable AR1 x AR1 field correlation
  for gls/lme with nugget and groups; Genstat's Slate Hall run reproduced
  (means and SEDs within 5e-6, KR-average d.d.f.). Kenward-Roger and
  Satterthwaite for gls/lme now leave out parameters estimated at 0, and a
  failed KR falls back to Satterthwaite.
- **Articles**: "Ordinal and categorical responses", "Spatial analysis of a
  field trial" (Slate Hall: variogram, AR1 x AR1, the `t2()` trend).
- **Citation**: Williams, Harwood and Matheson (2024), 3rd ed.,
  doi:10.1071/9781486317110.

## Earlier work

See `dev/NEWS-2.0-development.md` and the plan's status table (rows O1-O17).
