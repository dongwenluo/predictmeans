# Review of `predictmeans()` and `Kmatrix()` (2026-11-14)

Findings and proposals; the items marked **done** are in the code, the others
wait for the maintainer's decision. Scripts and timings: this file only
(the scratch scripts were not kept).

## predictmeans()

### Speed
| Case | Before | After |
|---|---|---|
| 12 means, lmer, pairwise | 0.03 s | 0.04 s |
| 60 means, lm, means only | 0.10 s | 0.02 s |
| 60 means, lm, `adj = "tukey"` | 1.82 s | 0.19 s |
| 96 means, lmer, means only | 1.53 s | 0.32 s |
| 96 means, lmer, pairwise | 1.52 s | 0.38 s |
| 200 lines, lmer, means/SED/LSD | 46.8 s | 0.7 s |
| 200 lines, lmer, `pairwise = TRUE` | minutes (letters) | 5 s |

- **done** The SEDs built the covariance matrix of all pairs to read its
  diagonal (n pairs x n pairs); now `rowSums((L V) * L)`.
- **done** One Tukey quantile per distinct df, not per pair.
- **done** More than 100 means, no comparison asked for: SED and LSD summaries
  from the covariance of the means (`pm_sed_summary()`); no pair table.
- **done** No letter groups above 100 means (the insert-absorb algorithm
  takes minutes for 200 means with unequal SEDs).
- Open: where the time goes now is the inference, computed again in every
  call on the same model: Kenward-Roger for `lme` 0.45 s and `gls` 0.2 s
  (numerical Hessian), `lmer` 0.05 s. A call for several terms repeats it.
  Options: (a) accept a vector of `modelterm`s; (b) keep the last model's
  inference in a small cache. `ptukey()` with many means is slow by itself
  (200 lines, Tukey: 18 s).

### Structure
- Each stage is a quoted block evaluated in the call's environment
  (`pm_code_setup`, ...), sharing about 90 variables; results move with
  `list2env()`, `pm_existing()` and name lists (`pm_plot_inputs`, 49 names;
  `pm_legacy_inputs`, 20). The pure functions are fine; the glue is hard to
  follow, and a typo in a name list fails silently (`ifnotfound = NULL`).
- Proposal: ordinary functions with explicit results: `opt <- pm_setup(...)`,
  `est <- pm_estimate(opt)`, `cmp <- pm_compare(opt, est)`,
  `tab <- pm_tables(opt, est, cmp)`, then plots and the result from those four
  lists. Checked with the snapshot (no output change). The plot code (400
  lines of quoted base/ggplot code) moves to functions of `opt`, `est`, `tab`.

### Arguments (45)
- 16 only concern plots (`lsd_bar`, `responsen`, `plotord`, `lineplot`,
  `plottitle`, `plotxlab`, `plotylab`, `which.plots`, `mplot`, `barplot`,
  `pplot`, `bkplot`, `jitterv`, `basesz`, `prtplt`, `newwd`). Proposal: remove
  the four deprecated ones (`mplot`, `barplot`, `pplot`, `bkplot`) in 2.0;
  keep `plot` and `which.plots`; pass the other plot settings through one
  `plot_args = list(...)` (1.x names still accepted with a note), and give
  `autoplot()` the same settings.
- `covariate`: one named list for values (see Kmatrix below).
- `prtnum`: a message already; default could be `FALSE`.
- `count`: set automatically for Poisson models; rarely needed.
- `ncore`: only used with `permlist`; fine.
- `Df`, `df_method`, `vcov_adjust`: three ways into one decision; keep (each
  documented as the 1.x switch), no change proposed.

### Other
- `Kmatrix()` called directly on `aov(... Error())` fails (`predictmeans()`
  converts the model first): convert inside `Kmatrix()`.
- `modelterm` must be a term of the model (`pm_check_modelterm()`): see below.

## Kmatrix()

### Efficiency
- The reference grid is the product of the levels of *all* factors of the
  fixed formula. Factors that do not interact with the term multiply it for
  nothing: `y ~ Line(300) + Block(90) + Site(4) + Year(3)` builds 324,000
  rows x 394 columns (3.6 s, 475 MB) where the average factorises: a column of
  the design belongs to one term, a term's variables always lie in one group
  of interacting variables, so each group can be averaged on its own small
  grid (here 300 + 90 + 4 + 3 rows). Valid for equal weights (and `"outer"`);
  `"proportional"`, `"cells"`, `"flat"` and the empty-cell rule need the joint
  cells of the groups they involve.
- `pm_model_info()` calls `model.matrix(model)` (n x p) only for the column
  names and contrasts: 40% of the time for n = 200,000; use the coefficient
  names and stored contrasts.
- The model frame is built twice (`model.frame()` for X and for the levels).
- `cov_reduce = "average"` with a continuous covariate crosses every observed
  value with every cell (200,000 values x 20 cells: 20 s, 1.4 GB).

### Covariates (bugs found)
- `covariate = c(z = 0.2, x = 0.5)`: the names are ignored; values go by
  position (x = 0.2, z = 0.5). Wrong values, no warning.
- `covariate = 0.5` with two covariates: recycled to both. `c(1, 2, 3)` with
  two covariates: an R warning, first two used.
- `covariateV = list(z = ..., x = ...)`: names ignored, by position.
- No way to vary one covariate and fix another at a chosen value.
- Proposal: one rule. `covariate` may be a named list (or named vector):
  `list(x = c(-1, 0, 1), z = 0.3)`; a variable with several values gives rows
  (all combinations, or a data frame for joint points), one value fixes it;
  variables not named use `cov_reduce`. Unnamed numeric values stay accepted
  only when there is exactly one value per covariate (else an error); the
  `covariate = "x", covariateV = ...` form keeps working.

### Any combination of factors and covariates
- `Kmatrix()` itself already accepts combinations that are not model terms
  (`"A:B"` in `y ~ A + B`); `predictmeans()` stops them. A covariate in
  `modelterm` (`"A:x"`) is refused unless `covariate` equals `modelterm`.
- Proposal: `modelterm` is any `:` combination of fixed-effect variables;
  factors give their levels, covariates their values (as above). Remove
  `pm_check_modelterm()` as a stop (keep its "did you mean" hint for unknown
  variables). Then the df rules that look the term up by name need a rule
  for a combination that is not a term: the smallest df of the terms whose
  coefficients the rows use (as the inner-outer rule already does).
- Open question: a combination with no data in an additive model (`a3:b2`
  in `y ~ A + B`) is estimable; `predictmeans()` drops rows without data.
  Show it (as emmeans) or keep dropping?

### New model classes
- `pm_model_info()` has one branch per class for the call, contrasts and
  coefficient names. Proposal: these three become adapter generics with
  defaults (`pm_contrasts()` exists), so a new class needs only `pm_coef()`,
  `pm_vcov()` and, if `model.matrix()` does not work, `pm_basis()`.

### Missing cells and nesting (checked)
Right as they are: an empty cell of interacting factors (warning, average of
the observed cells, comparisons on common cells); cells judged on the
interacting factors only (a block does not make cells look single; an empty
A:C cell does not matter when A and C do not interact); nested factors coded
`G/T`; nested factors crossed with a third factor with an accidental empty
cell.

Two cases that are not handled:
1. Nested factors coded additively, `y ~ G + T` (T within G): every mean is
   NA with the default weights, though `y ~ G/T` (the same model) gives them.
   The nesting rule is applied only when the formula has an interaction.
   Nesting is a fact of the data: combinations that cannot exist (each level
   of T occurs with one level of G) should be left out of the averages for
   any formula. (`weights = "cells"` or `"flat"` gives the right means now.)
2. The exemption "every cell has 0 or 1 observation" (kept from 1.x, where
   the cells included blocks): for unreplicated interacting factors with an
   empty cell, the mean is NA, while with two replicates it is the average
   of the observed cells with a warning. The cells no longer include blocks,
   so the exemption can go.
