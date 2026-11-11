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
| 2026-11-11 (queued) | `sim-perm-ranef.R 1000 99` GLMM fixed-effect cells (4) | same | same |
| 2026-11-11 (queued) | `sim-perm-edge.R 1000 99 all 4` (23 cells) | same | same |

## 4. Results so far

See `results.md` (regenerate with `Rscript dev/sim-study/summarise.R`).
Observations, added as cells finish:

- LMM split plot, test of the whole-plot variance (`permlmer()`): normal
  0.059, skewed errors and block effects 0.049, 20% missing 0.054; the χ̄²
  reference 0.046–0.049. All within the Monte Carlo range.

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

Worth writing up if the edge cases show that (a) the working-scale GLMM
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
