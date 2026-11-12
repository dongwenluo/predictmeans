# Research log: permutation tests for linear and generalized linear mixed models

Started 2026-11-11. Kept so the work can become a paper. Plan and
background: `dev/permutation-plan.md`. Earlier study of Freedman–Lane in
`lm`/`lmer`: `dev/sim-freedman-lane.md`.

## 1. Questions

1. Does the random-effect test of `permlmer()` (Lee and Braun 2012: whiten
   the H0 marginal residuals with the Cholesky factor of V0, permute,
   re-correlate, REML LRT) keep its level in designs typical of field and
   animal trials, including edge cases?
2. Does the maintainer's fixed-effect scheme for mixed models (permute
   fitted values and whitened residuals together, y* = xb[π] + U e[π]; one
   set of refits serves all terms and all pairwise comparisons) keep its
   level, and how does it compare with Freedman–Lane (y* = xb0 + U0 e0[π])
   and the Satterthwaite / Kenward–Roger F tests?
3. Do both schemes work for GLMMs when applied on the working (link) scale
   (linearise the GLMM fit, V = Z G Z' + diag(1/w), statistic from weighted
   LMMs of the working response)?
4. Where do they fail (edge cases), and do they protect better than the
   usual asymptotic tests when the model is misspecified (overdispersion,
   unequal variances, heavy tails)?

## 2. Methods compared (names used in the result files)

| name | test | applies to |
|---|---|---|
| `perm` | `permlmer()` as in the package (Lee–Braun, REML LRT) | LMM random effects |
| `chibar` | REML LRT, 50:50 mixture of χ²(k−1) and χ²(k) | LMM random effects |
| `joint` | maintainer's scheme; LMM: F of the term from full-model refits; GLMM: same on the working scale, linearised at the full fit | fixed effects |
| `fl` | Freedman–Lane; LMM: reduced-model fitted values kept, its whitened residuals permuted, F from full-model refits; GLMM: same on the working scale, linearised at the reduced fit | fixed effects |
| `satt` | Satterthwaite F (lmerTest) | LMM fixed effects |
| `lrt_chisq` | LMM: ML LRT, χ²; GLMM: Laplace LRT, χ² | fixed effects |
| `joint`, `fl` in `lmm_fix_lev` | permlmer()/permanova.lmer() form: LRT statistic (ML criterion at REML estimates), xb and V from the reduced model; `joint` permutes xb with the residuals | LMM fixed effects |
| `work_perm` | Lee–Braun on the working scale (linearised at H0), REML LRT of weighted working LMMs | GLMM random effects |
| `work_chibar` | the working-scale LRT with the χ̄² reference | GLMM random effects |
| `glmm_chibar` | Laplace LRT of the GLMMs with the χ̄² reference | GLMM random effects |

All permutation p-values use 99 permutations and the +1 rule, so the test
is exact at the 5% and 1% levels when exchangeability holds.

## 3. Scripts, settings and records

| script | results | cells |
|---|---|---|
| `dev/sim-freedman-lane.R` | `dev/sim-freedman-lane*.rds` | lm and lmer, raw vs Freedman–Lane (see its own write-up) |
| `dev/sim-perm-ranef.R` | `dev/sim-perm-ranef/<cell>.rds` | main cases: split plot, random slope, GLMM pilots, GLMM fixed effects |
| `dev/sim-perm-edge.R` (+ `dev/sim-perm-lib.R`) | `dev/sim-perm-edge/<cell>.rds` | edge cases (each cell's description is stored in its file and listed in the script) |
| `dev/sim-study/summarise.R` | `dev/sim-study/results.csv`, `results.md` | all of the above, with Monte Carlo SE and flags |
| `dev/sim-perm-pairs.R` | `dev/sim-perm-pairs/<cell>.rds` | study 1b part 1: pairwise comparisons (joint permlist, Freedman–Lane per pair and per term with V fixed, parametric t / Tukey), 199 permutations |
| `dev/sim-perm-fixed.R` | `dev/sim-perm-fixed/<cell>.rds` | study 1b part 2: term tests with the package's own `permmodels()` raw and Freedman–Lane, `permlmer()` for GLMMs, against KR F / LRT, 199 permutations |
| `dev/sim-study/summarise-pairs.R` | `dev/sim-study/results-1b.md` | study 1b: per-pair and familywise error, power; term tests |

- Seeds: data set i of a cell uses `set.seed(sum(utf8ToInt(cell)) * 1000 + i)`,
  so results do not depend on cores or platform; any cell can be rerun.
- Every p-value of every data set is saved (`$p`), with failures
  (`$errors`), timings, date and package versions (edge-case files), so other
  levels, QQ plots of p-values and failure rates can be studied later.
- 1000 data sets per cell: Monte Carlo SE of a 5% rate 0.0069; rates between
  0.036 and 0.064 are consistent with 5%. Bradley's liberal criterion
  (0.025–0.075) is the usual looser bound.

### Run log

| date (UTC) | what | machine | software |
|---|---|---|---|
| 2026-11-10 | `sim-freedman-lane.R` (lm 2000 sets; lmer 300; lmer leverage confirm 1000) | sandbox, 4 cores Xeon 2.1 GHz | R 4.3.3, lme4 1.1.35.1, Matrix 1.6.5 |
| 2026-11-11 08:02 | `sim-perm-ranef.R 1000 99` (10 cells) | same | R 4.3.3, lme4 1.1.35.1, Matrix 1.6.5, lmerTest 3.1.3, predictmeans 1.1.1.9000 |
| 2026-11-11 09:26 | `sim-perm-ranef.R 1000 99` GLMM fixed-effect cells (4) | same | same |
| 2026-11-11 10:00–13:00 | `sim-perm-edge.R 1000 99 all 4` (23 cells, about 7.5 min each) | same | same |
| 2026-11-11 (later; the sandbox clock jumped) | `sim-perm-followup.R 99 … 4` (4 cells) and `sim-perm-followup2.R 99 all 4` (2 cells) | same | same |
| 2026-11-12 00:07– | `sim-perm-pairs.R 1000 199 … 3` (8 cells; `pr_rcbd_power` 500 sets) | same | same, predictmeans 1.1.1.9000 after step 2 (engine) |
| 2026-11-12 (after the pairs run) | `sim-perm-fixed.R 1000 199 … 3` (LMM and glm cells), `500 199 … 3` (GLMM, lme and power cells) | same | same, predictmeans after steps 6–8, glmmTMB 1.1.8 |

## 4. Results so far

See `results.md` (regenerate with `Rscript dev/sim-study/summarise.R`).
Observations, added as cells finish:

- LMM, test of a variance component (`permlmer()`, Lee–Braun): split plot
  normal 0.059, skewed errors and block effects 0.049, 20% missing 0.054,
  random slope 0.052; χ̄² reference 0.046–0.049. All within the Monte Carlo
  range. Power (whole-plot SD 0.5): permutation 0.383, χ̄² 0.366.
- LMM fixed term with an extreme-value covariate, `permanova.lmer()` form
  (LRT statistic, reduced-model xb and V): joint 0.053, Freedman–Lane 0.050,
  ML LRT χ² 0.052. Unlike the `permmodels()` form with an F statistic (6.3%
  in `sim-freedman-lane.md`), the joint scheme kept its level here.
- GLMM variance component on the working scale (`work_perm`), split plot,
  6 blocks: binary 0.065 (just above the range), binomial with 10 trials
  0.049, Poisson 0.059. The Laplace LRT with the χ̄² reference was very
  conservative: 0.019, 0.013, 0.008. Power (Poisson, whole-plot SD 0.4):
  working-scale permutation 0.637, working χ̄² 0.620, Laplace χ̄² 0.547.
  So for GLMM variance components the working-scale permutation test is
  close to nominal and more powerful than the usual test.
- GLMM fixed effects on the working scale (joint / Freedman–Lane /
  Laplace LRT χ²): binomial with extreme covariate 0.045 / 0.044 / 0.040;
  Poisson with extreme covariate 0.040 / 0.056 / 0.054; binary with an
  associated factor 0.038 / 0.050 / 0.048. Power (Poisson): 0.544 / 0.625 /
  0.652 — the joint scheme is conservative in GLMMs and loses power.

### Edge cases (`dev/sim-perm-edge/`)

- LMM random effects (`permlmer()` / χ̄²): 3 blocks 0.052 / 0.040; t3 data
  0.039 / 0.035; unequal error SDs 0.028 / 0.018; singular H0 fit 0.052 /
  0.036; plots of 1–4 units 0.056 / 0.045; crossed rows × columns 0.064 /
  0.046; correlated random slope 0.055 / 0.044; power with 3 blocks 0.381 /
  0.358. Lee–Braun never above 6.4%.
- LMM fixed effects (joint / FL / Satterthwaite): 3 blocks with extreme
  covariate 0.060 / 0.064 / 0.059; associated factor with skewed errors
  0.036 / 0.041 / 0.041; empty cells 0.055 / 0.050 / 0.046; **global null
  0.071** / 0.061 / 0.052 (unexpected for the joint scheme; to investigate);
  **group sizes 1/2/3 with SDs 3/1/1: 0.243 / 0.226 / 0.237** (all fail; the
  model's variance assumption is wrong); power 0.673 / 0.633 / 0.632.
- GLMM random effects on the working scale (perm / work χ̄² / Laplace χ̄²):
  binary whole plots of 2 **0.078** / 0.069 / 0.017; **rare binary events
  0.358 / 0.483** / 0.023 (137 of 1000 fits failed); Poisson mean 0.4 0.058 /
  0.049 / 0.011; **negative binomial fitted as Poisson 0.098** / 0.091 /
  **0.631**; 3 blocks 0.046 / 0.033 / 0.002.
- GLMM fixed effects (joint / FL / LRT): overdispersed counts 0.077 / 0.103 /
  **0.352**; rare binary events 0.040 / 0.070 / 0.071 (15 failed fits);
  3 blocks 0.052 / 0.052 / 0.040; global null 0.060 / 0.067 / 0.053.

### Main findings

1. Lee–Braun (`permlmer()`) is robust in LMMs: 13 null cases, never above
   6.4%, more power than the χ̄² test.
2. On the GLMM working scale the same idea keeps the level for counts and
   binomial data, where the usual Laplace LRT with χ̄² is very conservative
   (0.2–2.3%) and less powerful (0.547 vs 0.637). It fails for rare binary
   events (0.358): the working response turns each event into an outlier of
   size about 1/μ.
3. The maintainer's fixed-effect scheme is close to nominal in most cases
   but liberal under the LMM global null (0.071) and with an extreme
   covariate (0.063, F form). In GLMMs it is conservative (0.038–0.060) and
   the most robust to rare events and overdispersion, at a cost in power.
4. Under unmodelled overdispersion the asymptotic GLMM tests break down
   (0.35, 0.63) while the permutation tests stay near 0.08–0.10 — the
   strongest argument for the methods.
5. With unequal group variances nothing works (0.23–0.24).

### Follow-ups (step 1a; `dev/sim-perm-followup.R`, `dev/sim-perm-followup2.R`, results in `dev/sim-perm-followup/`)

- LMM global null, 8 blocks, 2000 data sets: joint (Cholesky, full model)
  0.062, joint (reduced model) 0.057, joint with the **true** V 0.065,
  Freedman–Lane 0.057, exact within-block permutation 0.048, Satterthwaite
  0.048. 16 blocks (1000 sets): 0.084, 0.071, 0.081, 0.070, 0.057, 0.052. So
  the excess is not chance and not from estimating V; it grows with more
  blocks.
- Diagnosis (30 and 20 data sets × 99 permutations, 16 blocks): the
  permuted F statistics average 0.93–0.95 under the joint scheme against
  1.04 under the exact within-block permutation; the residual variance
  estimates agree (0.937 vs 0.935), the sum of squares for A is about 8%
  smaller (1.70 vs 1.85). Cause: the Cholesky factor U is order-dependent;
  within a block U·1 rises along the rows (1.41 … 2.89 for 6 units with
  equal variances), and A's levels sit in fixed positions in each block.
- Symmetric square root of V instead (2000 new data sets, 8 blocks):
  joint 0.052, Freedman–Lane 0.053 (Cholesky on the same data: 0.067,
  0.060), within 0.053, Satterthwaite 0.049. **Fixed.**
- Extreme covariate related to A (1000 sets): joint Cholesky 0.063, FL
  Cholesky 0.054, joint symmetric 0.064, FL symmetric 0.046, Satterthwaite
  0.052. The symmetric root does not cure the joint scheme here (the
  raw-permutation weakness with extreme covariates); FL with the symmetric
  root is fine.
- Binary GLMMs, test of the whole-plot variance (500 data sets each; fits
  with `nAGQ = 0`, about 150 times faster than Laplace refits on rare-event
  data, which also often failed): rare events route B 0.062, parametric
  bootstrap 0.040, Laplace χ̄² 0.026 (4 failed fits of 500); whole plots of
  2 route B 0.056, bootstrap 0.056, Laplace χ̄² 0.026. Both fix the
  working-scale failures (0.358, 0.078).
- Tie bug: rejections of `permlmer()` in data sets whose LRT was 0 up to
  optimizer noise were 0.1–0.8% of null data sets (from the saved p-values;
  e.g. split plot 0.059 → 0.056 without them, crossed 0.064 → 0.060, tiny
  plots 0.056 → 0.048). Fixed in the package (LRTs rounded to 4 decimals).
  The working-scale tests in the scripts already truncated at 0 and were not
  affected.

### Package check of the new whitening (`dev/sim-perm-confirm.R`, `dev/sim-perm-confirm/`)

With the symmetric root in the package's own `permlmer()` and
`permmodels()` (1000 data sets each): `permlmer()` split plot 0.059 (χ̄²
0.048), crossed rows × columns 0.054 (χ̄² 0.036; Cholesky 0.064);
`permmodels()` global null with `lmer` 0.048 (its F test 0.057; Cholesky
0.062–0.071) and with `lme` 0.059 (F test 0.057). Timing of the root with
crossed random effects: n = 3000 0.11 s (Cholesky 3.9 s); n = 20000 21 s
(Cholesky over 30 min and 2.9 GB, stopped).

### Study 1b, part 1: pairwise comparisons (`dev/sim-perm-pairs/`, 199 permutations)

1000 data sets per null cell (500 for power); per-pair rates are the mean
over the true null pairs (largest in brackets); familywise error (FWER) is
over the true null pairs.

| cell | joint (permlist) | FL per pair, V fixed | FL per term, V fixed | t (KR / nlme) | FWER max-\|t\| joint / FL term / Tukey |
|---|---|---|---|---|---|
| RCBD, global null | 0.054 (0.068) | 0.058 (0.071) | 0.057 (0.072) | 0.056 (0.069) | 0.058 / 0.065 / 0.059 |
| RCBD, one mean 3 SD away | 0.044 (0.047) | 0.052 (0.057) | 0.051 (0.052) | 0.048 (0.050) | 0.032 / 0.042 / 0.037 |
| RCBD, skewed, partial null | 0.046 (0.052) | 0.051 (0.057) | 0.052 (0.058) | 0.048 (0.055) | 0.024 / 0.030 / 0.029 |
| split plot, whole-plot pairs | 0.055 (0.059) | **0.062 (0.065)** | **0.062 (0.068)** | 0.050 (0.055) | 0.034 / 0.047 / 0.028 |
| split plot, sub-plot pair | 0.048 | 0.054 | 0.049 | 0.052 | (one null pair) |
| covariate with extreme values | 0.057 (0.061) | 0.054 (0.060) | 0.053 (0.058) | 0.052 (0.055) | **0.067** / 0.056 / 0.051 |
| lme, varIdent (SD 1,1,2,3) | 0.059 (0.067) | **0.069 (0.075)** | **0.069 (0.078)** | 0.065 (0.073) | 0.065 / 0.086 / 0.079 |
| gls, AR(1) | 0.039 (0.046) | 0.061 (0.064) | 0.062 (0.067) | 0.059 (0.062) | 0.043 / 0.072 / 0.067 |
| power (means 0, 0.4, 0.8, 1.2) | 0.276 | 0.282 | 0.285 | 0.278 | |

Reading:

- One set of refits under the global null (the maintainer's `permlist`)
  works for pairs: 3.9–5.9% per pair, never above the Monte Carlo range;
  slightly conservative when other means differ (4.4–4.6%), and its max-|t|
  familywise error is right except with a covariate with extreme values
  (6.7%, the known weakness of permuting fitted values).
- Freedman–Lane with V fixed (no refits) is right in block designs (5.1–5.8%)
  but liberal where V rests on estimated variance parameters with little
  information: whole-plot pairs with 6 blocks (6.2%), varIdent (6.9%), AR(1)
  (6.1%) — the same as the parametric t, since both treat V as known. Hence
  Freedman–Lane **with refits** for linear mixed models (V estimated again
  for each permuted data set), checked in part 2 (`fx2_sp_pairs`,
  `fx2_lme_varident`).
- Freedman–Lane per term (one reduced model for all pairs of a factor) is as
  good as per pair, and gives a max-|t| familywise adjustment that matches
  Tukey's (a permutation Tukey test): this is the scheme used for pairs with
  refits.
- Power was the same for all methods.

### Study 1b, part 2: tests of a fixed term with the package code (`dev/sim-perm-fixed/`, 199 permutations)

`raw` = `permmodels(method = "raw")` (1.x scheme); `FL fixed` = Freedman–Lane
with V fixed (`refit = FALSE`); `FL refit` = Freedman–Lane with refits (the
2.0 default for mixed models; for GLMMs the working linear mixed model);
`permlmer` = `permlmer(H0, H1)` (working scale; quantile route for 0/1
data). 500 data sets per cell (range 0.031–0.069), 1000 for glm (0.036–0.064).

| cell | raw | FL fixed | FL refit | permlmer | usual test |
|---|---|---|---|---|---|
| LMM split plot, whole-plot factor (6 blocks) | 0.058 | **0.086** | **0.050** | | KR 0.046 |
| LMM split plot, whole-plot pairs (null pairs; Holm) | | (pairs study 0.062) | 0.052 (Holm 0.042) | | KR 0.051 (0.040) |
| LMM covariate with extreme values | 0.060 | 0.050 | 0.054 | | KR 0.048 |
| LMM unbalanced (20% missing) | 0.052 | 0.056 | 0.048 | | KR 0.054 |
| lme, varIdent | 0.036 | 0.042 | 0.036 | | nlme F 0.042 |
| glm Poisson, extreme covariate | **0.008** | 0.058 | | | LRT 0.055 |
| glm binomial, extreme covariate | **0.006** | 0.042 | | | LRT 0.050 |
| glm, negative binomial counts fitted as Poisson | 0.019 | 0.061 | | | LRT **0.329** |
| GLMM Poisson, extreme covariate | | 0.048 | 0.058 | 0.060 | LRT 0.052 |
| GLMM, negative binomial fitted as Poisson | | **0.100** | 0.074 | 0.086 | LRT **0.339** |
| GLMM binomial, 3 blocks | | 0.066 | 0.062 | 0.062 | LRT 0.064 |
| GLMM binary, rare events | | 0.031 | 0.041 | 0.051 | LRT 0.071 |
| glmmTMB nbinom2 (right family) | | 0.066 | 0.068 | 0.056 | LRT 0.070 |
| power: LMM, extreme covariate | 0.626 | 0.634 | 0.626 | | KR 0.626 |
| power: GLMM Poisson (A effects 0, 0.2, 0.4) | | 0.433 | 0.429 | 0.453 | LRT 0.475 |

Notes: in the first runs of the LMM and glm cells the "raw" column was
Freedman–Lane (the script relied on the default, which had become
`"auto"`); the glm cells were rerun and the LMM raw column recomputed on the
same data sets (`dev/sim-perm-fixed-raw.R`; the old column is kept as
`flr_all`, a second Freedman–Lane-with-refits run: 0.054, 0.052, 0.046).

Reading:

- Linear mixed models: Freedman–Lane with refits kept the level everywhere
  (3.6–5.4%), with the same power as the others; with V fixed it was
  liberal for the whole-plot factor (8.6%), because the permuted statistics
  then treat V as known. Hence the 2.0 default (refits).
- glm: permuting the response (1.x) was extremely conservative with a
  covariate (0.6–1.9%: it scrambles the covariate's effect), Freedman–Lane on
  the working scale kept the level (4.2–6.1%), and it was robust to
  overdispersion where the LRT failed (6.1% vs 32.9%).
- GLMMs: Freedman–Lane with refits of the working LMM 4.1–6.8% in
  well-specified models (nbinom2 with 48 observations: every test 5.6–7.0%),
  and much more robust to overdispersion than the LRT (7.4% vs 33.9%); V
  fixed was liberal there (10.0%). `permlmer()` (quantile route) is the most
  accurate for rare binary events (5.1%).

### Open questions

- ~~Why the joint scheme is liberal under the LMM global null (0.071).~~
  Answered: Cholesky whitening; the symmetric root fixes it.
- ~~Whether route B or a parametric bootstrap fixes rare binary events.~~
  Answered: both do.
- ~~Does the symmetric root also improve `permlmer()`?~~ Yes for the crossed
  design (0.064 → 0.054); unchanged for the split plot. Still open: the GLMM
  working scale, where the order-invariant root is not symmetric.
- ~~Pairwise comparisons by permutation~~: study 1b part 1 (`dev/sim-perm-pairs.R`),
  results above and in `results-1b.md`. Designs: RCBD (global null; partial null with
  one mean 3 SD away; skewed errors), split plot (whole-plot and sub-plot
  pairs, partial nulls), covariate with extreme values, `lme` with
  `varIdent`, `gls` with AR(1), power. Question: is one set of refits
  under the global null (the maintainer's `permlist`) good enough for each
  pair and for the family (max-|t|), or is a reduced model per pair
  (Freedman–Lane) needed?
- ~~Term tests with the package code~~ (`dev/sim-perm-fixed.R`), results above:
  Freedman–Lane with V fixed (no refits) for LMMs, on the working scale for
  glm and GLMMs, against raw permutation and KR F / LRT, including
  overdispersed counts and rare binary events.
- How a parametric bootstrap from the reduced (Poisson) model behaves under
  overdispersion (expected to fail like the LRT, since it trusts the model).

- Note on the run: the first batch ended with an error after "SIM DONE"
  because its script was edited while it ran (Rscript reads the file as it
  goes); all cells had been saved. Scripts are not edited during runs now.

## 5. Is it worth a paper?

What would be new, as far as the review in `dev/permutation-plan.md` found
(a proper literature search is still needed before claiming novelty):

- Lee–Braun for LMM random effects is published; our results only confirm
  it in agricultural designs and edge cases.
- The maintainer's fixed-effect scheme for mixed models (whitened residuals
  permuted together with fitted values, one set of refits for all terms and
  all pairwise comparisons) was not found in the literature in this form.
  Related: Lee and Braun (2012) (whitening), Kherad-Pajouh and Renaud (2015)
  (residual permutation in mixed designs), Manly's raw-data permutation
  (which it equals when V = I).
- For GLMMs, permutation on the working scale is published for a variance
  component (Zeng et al. 2015, PQL working response). For fixed effects in
  GLMMs there is no established permutation test; parametric bootstrap is
  the usual reference.

**Update after the follow-ups (2026-11-11):** the whitening finding adds a
clear methodological point for a paper: Cholesky-whitened residual
permutation (as in Lee and Braun 2012 and in our schemes) is order-dependent
and liberal for within-block fixed effects, increasingly so with more
blocks; an order-invariant square root removes the problem. To be checked
against the literature (we have not seen it reported).

**Update after the edge cases (2026-11-11):** yes, probably worth writing
up, centred on GLMMs: the working-scale permutation tests keep their level
where the usual Laplace tests are far too conservative (variance components)
or far too liberal (overdispersion), with clear limits (rare binary events,
unequal variances) that a paper should report. Before writing: the
parametric bootstrap comparison, route B for binary data, the global-null
puzzle, power curves, and a literature search on the joint scheme.

Original criterion: worth writing up if the edge cases show that (a) the working-scale GLMM
tests keep their level where the Laplace LRT / Wald tests do not (small
samples, overdispersion, rare events), and (b) they compare well with the
parametric bootstrap at lower cost. Possible outlets: *Computational
Statistics & Data Analysis*, *Journal of Statistical Computation and
Simulation*, *Biometrical Journal*, *Statistical Papers* (methods), or *The
R Journal* (software).

Still needed for a paper:

1. Parametric bootstrap from the reduced model as a reference in every
   cell (costly for GLMMs: two GLMM fits per bootstrap sample).
2. Power curves (several effect sizes), not just one alternative.
3. Levels 1% and 10%, and QQ plots of null p-values.
4. Timings per method.
5. Real examples: Oats split plot (lmer), Clinical trial (binomial GLMM),
   a count example (Poisson / negative binomial GLMM).
6. A short theory section: second-moment exchangeability of whitened
   residuals (Commenges 2003), why the joint scheme equals raw permutation
   when V = I, and studentised statistics (Chung and Romano 2013).
7. Literature search for prior use of the joint scheme (mixed-model
   permutation literature since 2000; Genstat, SAS, permuco documentation).

## 6. Reproducing

From the package root, with predictmeans installed:

```
Rscript dev/sim-perm-ranef.R 1000 99 <cells> <cores>
Rscript dev/sim-perm-edge.R 1000 99 all <cores>
Rscript dev/sim-study/summarise.R
```

`sim-perm-edge.R` skips cells that already have a result file, so it can be
restarted after an interruption; on Windows it uses a socket cluster.
Do not edit a script while it runs (Rscript reads it as it goes).
