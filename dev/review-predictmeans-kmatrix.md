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

## Done on 2026-11-15 (Kmatrix)

Covariate values by name; any combination of variables; grouped grids;
cells per interaction term; no 0/1 exemption; `pm_contrasts()`; see NEWS.
Checked with 592 `Kmatrix()` cases (before/after, and grouped grid against
the full grid) and the 52-call snapshot.

**Kept for decision: nesting that is only in the data.** `y ~ G + T` with T
nested in G (control/treated and treatments) and the spray trial
`y ~ droplet + conc + vol` with volume 400 only at fine droplets are the same
structure in the data. A rule from the data alone gives means in both; the
recorded decision for the spray trial is NA (the effects of droplet and of
volume 400 cannot be told apart). So nesting still has to be said in the
formula (`G/T`); without it the means are NA and a note says how to get them.

## Questions of 2026-11-15

### A combination that is not a model term in predictmeans()
Reasonable, with a note. The means are what the model implies: sums of the
factors' effects, so the differences between the levels of one factor are the
same at every level of the other and add nothing to the main-effect
comparisons. Useful for: a table of predictions for every combination,
including combinations without data (line x site from an additive model);
reports that want all cell predictions on one page. Not useful for testing.

### 1. How `adj` and the LSD are linked (as the code is now)
| `adj` | p-values | printed LSD | half-width that matches the p-values |
|---|---|---|---|
| none | t | t x SED | the same |
| bonferroni | x m | t at level/m x SED | the same (m = all pairs, even with `atvar`, where the p-values use the pairs within a level) |
| tukey | studentized range | t x SED (not adjusted) | qtukey/sqrt(2) x SED |
| sidak, scheffe, mvt | adjusted | t x SED (not adjusted) | their own critical value x SED |
| holm, hochberg, hommel, BH, BY | adjusted | t x SED (not adjusted) | none exists (step-wise methods) |
| permlist | permutation | 2 x SED | permutation quantile |

Oats, nitro, lm: printed LSD 10.65 with Tukey p-values; the matching HSD is
14.08. With `atvar = "Variety"` and Bonferroni the printed LSD is 32.8
(0.05/66) while the p-values use 0.05/6 (25.2). The `pairs` table already has
the matching limits (`adj.conf.low`, `adj.conf.high`).

Proposal: one rule. The LSD is the half-width of the limits that agree with
the p-values shown, with the same family: t (none), HSD (tukey), adjusted t
(bonferroni, sidak), Scheffe, mvt. For the step-wise methods there is no such
value: print the unadjusted LSD and say so in its note.

### 2. Significance and overlap of confidence intervals
Two means differ at level alpha when |d| > t x SED. Their intervals
mean +- c x SE do not overlap when |d| > c (SE1 + SE2). The two agree when

    c = t x SED / (SE1 + SE2),   SED^2 = SE1^2 + SE2^2 - 2 cov

- independent, equal SEs: c = t/sqrt(2): 83.4% intervals for alpha = 0.05
  (the spreadsheet; Austin & Hux 2002; Payton, Greenstone & Schenker 2003;
  Velleman & Hoaglin's 1.39);
- independent, SE ratio K: c = t sqrt(1 + K^2)/(1 + K), 83.4% at K = 1 up to
  95% as K grows (Maghsoodloo & Huang 2010, eq. 4; Velleman & Hoaglin's
  compromise 1.7);
- correlated means (blocks, repeated measures, mixed models): the covariance
  makes the SED smaller, so c is smaller still (Maghsoodloo & Huang eq. 20c;
  Afshartous & Preston 2010). Oats with random blocks: SE 7.34, SED 5.10:
  intervals that just touch at p = 0.05 are 51% intervals, and 95% intervals
  that do not overlap mean p < 0.00001. The 83% rule is wrong there.
- several means: one level for all intervals chosen so that the average
  error rate of the pairs is alpha (Goldstein & Healy 1995), or so that
  overlap agrees with as many of the tests as possible (Armstrong & Poirier
  2024); emmeans draws "comparison arrows" fitted the same way.

In the package now: `letterCI = TRUE` uses a fixed 83.4% (level 0.166) and its
note assumes uncorrelated means with equal SEs; the mean plots show 95% CIs or
half-LSD bars (`autoplot(interval = "lsd")`: exact when the SEDs are equal).

Proposal: (a) confidence intervals of means stay intervals for the means and
are not described as tests; (b) "comparison intervals": mean +- c x SE with c
from the formula above, one c for the plot (exact when all pairs have the
same c, as in balanced designs; else the average over the pairs, with the
number of pairs whose test disagrees reported), with the family's critical
value in place of t when `adj` has one; the half-LSD bars are the same thing
when SEDs are equal; (c) `letterCI` uses that level computed from the model
instead of 83.4%; (d) `covariatemeans()` and `semipred()` bands are bands for
each curve: to compare groups, a band for the difference between the curves,
not overlap; `contrastmeans()` intervals are already intervals of
differences.

### 3. Back-transformation
- GLM (no random effects): the inverse link of the linear predictor is the
  mean at those predictor values. No adjustment. The only choice is whether
  the other factors are averaged on the link scale (default) or on the
  response scale (`average = "response"`).
- GLMM: the inverse link with the random effects at zero is the response of
  a typical (median) cluster, not the mean of the population of clusters.
  `bias_adjust = TRUE` gives the population mean, h(eta) + sigma^2/2 h''(eta),
  sigma^2 the sum of the random-intercept variances (as emmeans). Default
  FALSE.
- Transformed response in a linear (mixed) model (`trans = exp` etc.): the
  back-transformed mean is the median on the original scale (geometric mean
  for log). `bias_adjust = TRUE` gives the mean, with the residual variance
  added to sigma^2. For log the exact value is exp(mu + sigma^2/2), for sqrt
  mu^2 + sigma^2; the code uses the second-order formula for every
  transformation (exact for sqrt, low for log when sigma is large).

Proposal: (a) find the transformation of the response in the formula
(log, log1p, log10, log2, sqrt, 1/y, y^p, logit) and back-transform by
default, as for links (now `trans` must be given); (b) Box-Cox and
Yeo-Johnson as `trans` objects with their inverse, so the same back-
transformation and adjustment work; (c) exact adjustment for log (and
quadrature for the others) in place of the second-order formula; (d) the
output says which it is: "back-transformed (median)" or "mean (bias
adjusted)"; (e) `pm_check()` gives the Box-Cox profile of an lm fit (lambda
and its interval, read as log/sqrt/none) from the residuals of the model, and
points to a GLM/GLMM for counts and proportions. Not proposed: choosing the
transformation automatically by a normality test (bestNormalize): it is a
modelling decision, and rank-based normalising has no back-transformation
for means and differences.

References: Austin & Hux (2002) J Vasc Surg 36:194; Maghsoodloo & Huang
(2010) J Stat Plan Inference 140:3295; Velleman & Hoaglin (1981) ABCs of EDA,
3.12; Goldstein & Healy (1995) JRSS A 158:175; Schenker & Gentleman (2001) Am
Stat 55:182; Payton, Greenstone & Schenker (2003) J Insect Sci 3:34;
Afshartous & Preston (2010) Comput Stat Data Anal 54:2296; Armstrong &
Poirier (2024) Political Analysis, "Decoupling visualization and testing
when presenting confidence intervals".
