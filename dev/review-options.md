# Review of the options: what 2.0 keeps, renames, merges or removes

Step 0 of the clean-up (maintainer's request, 2026-11-18). No code is changed
by this note; the maintainer decides each row, then steps 1-7 carry it out.

Sources: the CRAN release 1.1.1 (2024-11-29, tarball supplied by the
maintainer; its arguments read by parsing `R/`, not by running it) and
`main` at 622c56a. 1.1.1 exported 20 functions; `main` exports 41.

## Rules for every removal (Rule no. 1: robust, never a silent wrong answer)

1. A removed argument or function **stops** with one line giving the 2.0
   spelling, e.g. `'Df' was removed in 2.0: use df_method = 10 (see
   vignette("migrating"))`. Never ignored through `...`, never a warning.
2. The messages live in one table, `R/defunct.R`; every front end calls one
   check (`pm_defunct(match.call())`).
3. Every removed spelling has a test: the 1.1.1 call, rewritten as the table
   says, gives the same numbers as before (tolerance 0, snapshot).
4. `vignette("migrating")` lists every row of the table below that is not
   "keep", with old and new code side by side.

## Decisions for the maintainer

| | Question | Proposal |
|---|---|---|
| D1 | `level` is a significance level (0.05) in `predictmeans()`, `contrastmeans()`, `covariatemeans()`, `semipred()`, `pm_slopes()`, `permstats()`, `PMplot()`, but a confidence level (0.95) in `varcomp()`; in `residplot()` it is the grouping level of random effects | `conf.level = 0.95` everywhere (R's `t.test()`, `confint()`, emmeans' `level`), `level` removed with the message "use conf.level = 1 - level". A new name, so an old `level = 0.05` cannot be read as 5% confidence. `residplot(level =)` becomes `re_level` |
| D2 | `Df` (a fixed df) beside `df_method` | `df_method` takes a number: `df_method = 10` |
| D3 | `vcov_adjust` beside `vcov` and `df_method` | keep `vcov_adjust` (one meaning, used to reproduce 1.x); or fold into `vcov = "unadjusted"`. Proposal: keep (it is a statistical choice, documented, tested) |
| D4 | `covariateV` (`Kmatrix`, `covariatemeans`, `semipred`) and `at` (`pm_slopes`) | `at` everywhere (a named list of covariate values) |
| D5 | `semireg()` and `semireg_tmb()`: one engine, two front ends with different arguments | keep both (maintainer, 2026-11-18, step 3b: a merge removes no code and leaves half the options for one engine only); remove `prt`, `predict_info`, `ngbinomial` |
| D6 | `permanova.lmer()` beside `permmodels()` (which tests lmer/glmer fixed effects on the shared engine) | remove `permanova.lmer()`: `permmodels(model, type =)` |
| D7 | `ci_mcp()`: letters from overlap of any intervals given by the user; overlapping 95% (or 83.4%) intervals do not answer the tests (Rule no. 1, `R/comparison.R`) | remove: `predictmeans(letterCI = TRUE)` gives letters from comparison intervals that agree with the tests |
| D8 | `R2_glmm()` beside `R2_model()` (same Piepho calculation, every class) | remove `R2_glmm()`: `R2_model(model)` |
| D9 | `smZ()` and `smoothZ =` (superseded 2026-11-17) | remove both: smooths in the formula (`s(x, bs = "os")` = 1.x `ZOSull`, `bs = "ztps"` = `Ztps`); the message names the formula term for each `type` |
| D10 | plot and print settings in the call of the analysis functions | analysis functions keep only `plot =` (TRUE draws the default plot); every look setting goes to `autoplot()` / `plot_args`; printing to `print(x, digits =)` |

## Function by function

Legend: **keep**; **rename** (old name stops with the new one); **merge**
(into another argument or function); **remove** (with the replacement in the
message); *new* = added in 2.0, not in 1.1.1. Default changes are listed at the
end.

### predictmeans()

| Argument | 1.1.1 | Proposal | 2.0 spelling |
|---|---|---|---|
| model, modelterm, data, pairwise, atvar, adj, covariate, letterCI, permlist, ncore | yes | keep | |
| Df | yes | merge (D2) | `df_method = 10` |
| level | 0.05 | rename (D1) | `conf.level = 0.95` |
| meandecr | yes | keep (step 2a: it sets which letter goes to the largest mean, so it cannot move to `print()`) | |
| trans | `I` | keep (only for a transformation outside the formula; detected from `log(y)` etc.) | |
| transOff | yes | remove | `trans = function(x) exp(x) - 1` |
| count | yes | remove: Poisson/neg. binomial already detected | `trans = ` for other cases |
| responsen | yes | already gone in 2.0; stop with message | not needed (response read from the model) |
| lsd_bar, plotord, lineplot, plottitle, plotxlab, plotylab, jitterv, basesz | yes | already moved to `plot_args`; now stop instead of warning | `plot_args = list(...)` or `autoplot()` |
| mplot, barplot, pplot, bkplot, prtplt | yes | already replaced by `which.plots`; stop | `which.plots = c("mean", "bar", ...)` |
| which.plots + plot | *new* / TRUE | merge into one | `plot = TRUE` (default set) or `plot = c("mean", "pval")` |
| prtnum, ndecimal | TRUE, 4 | remove | `print(x, digits = 4)`; nothing printed during the call |
| newwd | yes | remove | (ggplot2 draws in the current device) |
| cov_reduce, exclude, weights, df_method, vcov_adjust, vcov, equivalence, bias_adjust, sigma, mean_type | *new* | keep | |
| average | *new* | keep (how back-transformed means average: link or response) | |
| plot_args | *new* | keep | |

### contrastmeans(), pm_joint(), pm_slopes(), df_term()

| Function | Argument | Proposal | 2.0 spelling |
|---|---|---|---|
| contrastmeans | Df | merge (D2) | `df_method = 10` |
| contrastmeans, pm_slopes | level | rename (D1) | `conf.level` |
| pm_slopes | at | keep (D4: the name everywhere) | |
| df_term | type | rename to the package's name | `df_method` |
| df_term | vcov_adj | remove (step 2a: it was an internal covariance matrix, now `pm_df_term()`) | none |
| all | the rest | keep | |

### covariatemeans()

| Argument | Proposal | 2.0 spelling |
|---|---|---|
| model, modelterm, covariate, data, contr, contr_scale, band, cov_reduce, exclude, weights, df_method, vcov_adjust, bias_adjust, sigma, mean_type, plot | keep | |
| covariateV | rename (D4) | `at` |
| as.is (each line over its group's observed range) | keep (step 2b: `at` cannot do it; it changes the rows of the result) | |
| level | rename (D1) | `conf.level` |
| Df | merge (D2) | `df_method = 10` |
| trans, transOff, responsen | as `predictmeans()` | `trans` only |
| trellis, plotord, mtitle, ci, point, jitterv, newwd | remove (D10) | `autoplot(x, ...)` / `plot_args` |

### semipred()

| Argument | Proposal | 2.0 spelling |
|---|---|---|
| semireg | rename: the first argument is a model like everywhere else | `model` |
| modelterm, covariate, sm_term, contr, boundary, trans, contr_scale, band, unconditional, deriv, mean_type, plot, plot_args | keep | |
| covariateV | rename (D4) | `at` |
| level | rename (D1) | `conf.level` |
| simultaneous | remove (same as `band`) | `band = "simultaneous"` |
| prt | remove (D10) | `print(x)` |
| trellis, scales, plotord, ci, point, jitterv, threeD | remove (D10) | `plot_args` / `autoplot()` |

### semireg(), semireg_tmb(), smZ()

| Argument | Proposal | 2.0 spelling |
|---|---|---|
| semireg_tmb() | keep (D5) | |
| smZ(), smoothZ | removed (D9, step 3c; side-by-side code in `vignette("migrating")`) | `s(x, bs = "os")` etc. in the formula |
| ncenter, nscale, resp_scale | remove (scaling only served `smoothZ`; formula smooths are location-invariant) | none |
| ngbinomial | remove (not the ML fit) | `semireg_tmb(family = glmmTMB::nbinom2)` |
| prt | remove (D10) | `print(fit)` / `summary(fit)` |
| predict_info | remove (always stored; small) | none |
| formula, data, family, REML, control, start, verbose, drop.unused.levels, subset, weights, offset, contrasts | keep | |

### Permutation and resampling

| Function | Argument | Proposal | 2.0 spelling |
|---|---|---|---|
| permlmer | lmer0, lmer1 | rename (any class: lmer, glmer, glmmTMB) | `model0`, `model1` |
| permglmer(), permglmmTMB() | (aliases of `permlmer`) | remove | `permlmer()` |
| permanova.lmer() | | remove (D6) | `permmodels(model, type =)` |
| permmodels | prt | remove (D10) | `print()` |
| permmodels | fo | removed (step 4): a formula held in a variable is detected (`pm_formula_is_object()`) | none |
| permmodels | exact, data, method, perms, refit, robust, test.statistic, type, nperm, ncore, seed | keep | |
| bootmodels | prt | remove (D10) | `print()` |
| permstats | level | rename (D1) | `conf.level` |
| permindex, permcombine | all | keep | |

### Diagnostics, plots and small functions

| Function | Argument | Proposal | 2.0 spelling |
|---|---|---|---|
| residplot | level | rename (D1) | `re_level` |
| residplot | newwd, ask | remove | (ggplot2) |
| CookD | newwd | remove | |
| CookD | idn | rename (as `residplot`) | `id` |
| PMplot | level | rename (D1): the significance shading | `conf.level` (shading at p < 1 - conf.level) |
| PMplot | margin, legendx, newwd | remove (not used since the ggplot2 version) | |
| PMplot | mtitle, xylabel | keep, or move to `autoplot(type = "pvalue")` — keep | |
| ci_mcp() | | remove (D7) | `predictmeans(letterCI = TRUE)` |
| R2_glmm() | | remove (D8) | `R2_model(model)` |
| se_ranef, varcomp, pm_check, semiboundary, Kmatrix, contrast_family, specify_contrast, adapter generics (`pm_coef` ...), pm_boxcox, pm_yeojohnson | | keep |
| varcomp | level | rename (D1) | `conf.level` |
| Kmatrix | covariateV, prtnum | rename (D4) / remove | `at` / none |

## Defaults that change between 1.1.1 and 2.0 (for the migration guide)

| Function | Argument | 1.1.1 | 2.0 |
|---|---|---|---|
| predictmeans, covariatemeans, semipred | plot | TRUE | FALSE (`options(predictmeans.plot = TRUE)` for the 1.x behaviour) |
| predictmeans, permlmer, permmodels, permanova.lmer | ncore | 3 | 1 |
| predictmeans, contrastmeans | standard errors for REML mixed models | unadjusted | Kenward-Roger (`vcov_adjust = FALSE`) |
| several | df | per class | `df_method = "auto"` (table in `vignette("predictmeans")`) |
| semireg | control | `bobyqa` | `nloptwrap` |
| residplot | id | FALSE | 3 |

## Size of the change

About 45 arguments and 7 functions (`semireg_tmb`, `smZ`, `permglmer`,
`permglmmTMB`, `permanova.lmer`, `ci_mcp`, `R2_glmm`) go; the means
functions share one vocabulary (`at`, `conf.level`, `df_method`, `weights`,
`vcov`, `plot`). Every 1.1.1 call keeps working once rewritten as the table
says, with the same numbers, and every 1.1.1 call that is not rewritten
stops with the rewrite in its message.
