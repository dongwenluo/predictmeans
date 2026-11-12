# Permutation and resampling tests in predictmeans: review and plan

Written 2026-11-11. Goal (from the maintainer): keep in predictmeans the most
robust permutation test for each model class, and extend `permlmer()` to
`glmer` and `glmmTMB` models by permuting on the transformed (link) scale.
This file reviews what the package does now, checks the existing ideas
against theory and simulation, surveys the literature, and sets out the steps.

Contents: 1 Summary · 2 What the package has now · 3 Are the existing
schemes right? (random effects; the maintainer's fixed-effect scheme) · 4
Extending both to GLMMs on the working scale · 5 What "most robust" means · 6
Literature by model class · 7 Candidates and current recommendation · 8
Design: one engine · 9 Simulation study to choose the defaults · 10 Steps ·
11 Is there enough evidence to decide? · 12 Decisions for the maintainer ·
13 Risks · 14 References

## 1. Summary

Results of 37 simulated cases (1000 data sets each, 99 permutations; full
tables in `dev/sim-study/results.md`, research log in
`dev/sim-study/README.md`). Rates are type I error at the 5% level; 3.6–6.4%
is consistent with 5%.

- **Random effects in LMMs (`permlmer()`, Lee–Braun): keep as it is.** In 13
  null cases it stayed between 3.9% and 6.4%, except with unequal error
  variances (2.8%, conservative; the χ̄² test 1.8%). It never exceeded 6.4%
  and had slightly more power than the χ̄² test (38% vs 36–37%).
- **The maintainer's fixed-effect scheme in LMMs: sound; slightly liberal
  in two cases with Cholesky whitening** — no fixed effects at all (7.1%) and
  a covariate with extreme values with the F statistic (6.3%, earlier study)
  — fixed by symmetric-root whitening (follow-ups below). Elsewhere
  3.6–6.0%. Freedman–Lane was within or at the edge of the range in every
  case (4.1–6.4%), and so was the Satterthwaite F test. In LMMs the
  permutation tests do not beat Satterthwaite/Kenward–Roger.
- **Random effects in GLMMs on the working scale (the new `permglmer()`):
  works for counts and binomial data, fails for rare binary events.** Poisson
  5.9%, low counts 5.8%, binomial with 10 trials 4.9%, 3 blocks 4.6%; binary
  data 6.5% (whole plots of 5) and 7.8% (whole plots of 2); rare binary events
  (about 5%) 35.8% — unusable. The usual Laplace LRT with the χ̄² reference
  was very conservative (0.2–2.3%) and less powerful (55% vs 64%).
- **Fixed effects in GLMMs on the working scale: both schemes keep the level
  in well-specified models** (joint 3.8–6.0%, Freedman–Lane 4.4–6.7%, LRT
  4.0–5.4%). The joint scheme is more conservative and less powerful (54%
  vs 63% for Freedman–Lane, 65% for the LRT) but more robust with rare events
  (4.0% vs 7.0%).
- **Overdispersion not in the model** (negative binomial counts fitted as
  Poisson): the LRT rejected 35% (fixed effect) and 63% (whole-plot variance)
  of true null hypotheses; the working-scale permutation tests 7.7–10.3%.
  Still above 5%, but this is where permutation helps most.
- **Unequal variances by group** (small group with large variance): every
  test failed (23–24%). Permutation cannot fix a wrong variance model.

**Follow-ups (step 1a, `dev/sim-perm-followup*.R`):**

- **The cause of the liberal fixed-effect results is the Cholesky
  whitening, not the scheme.** Under the LMM global null the Cholesky
  versions rejected 6.2–6.7% (8 blocks) and 7.0–8.4% (16 blocks), even with
  the true V. Whitening with the symmetric square root of V gave 5.2% (the
  maintainer's scheme) and 5.3% (Freedman–Lane), the same as the exact
  within-block permutation (5.3%) and Satterthwaite (4.9%). With a covariate
  with extreme values the symmetric root does not help the maintainer's
  scheme (6.4%; Cholesky 6.3%) — the known weakness of permuting fitted
  values with extreme covariates, as for `lm` — while Freedman–Lane with the
  symmetric root holds (4.6%). The Cholesky factor treats the first unit of
  each block differently, so the results also depend on how the data are
  sorted; the symmetric root does not.
- **Binary GLMMs: route B (quantile back-transform) and the parametric
  bootstrap both fix the failures of the working scale**: rare events 6.2%
  and 4.0% (working scale 35.8%), whole plots of 2 both 5.6% (7.8%), on 500
  data sets each, with fast `nAGQ = 0` refits.
- **A small bug in `permlmer()` was found and fixed**: when the tested
  variance is estimated at 0, an observed LRT of 1e-9 could look larger than
  permuted LRTs of 0 (0.1–0.8% of null data sets). Statistics are now
  rounded before they are compared.

Recommendations: §7 (table, updated with these results), §10 (steps) and
§11 (what is ready to decide).

## 2. What the package has now

| Function | Models | Hypothesis | How the data are permuted | Statistic |
|---|---|---|---|---|
| `permlmer(m0, m1)` | lmer | random effects (same fixed effects) | whitened H0 marginal residuals permuted, re-correlated, added to H0 fitted values | REML LRT |
| `permlmer(m0, m1)` | lmer | fixed effects (fixed effects differ) | H0 fitted values **and** whitened residuals permuted together | ML log-likelihood at REML estimates |
| `permanova.lmer()` | lmer | each fixed term (type I/II/III) | calls `permlmer()` on nested pairs | as above |
| `permmodels(method = "raw")` | lm, aov, glm | terms, coefficients, `permlist` | response permuted | F (Deviance for glm type I), t |
| `permmodels(method = "raw")` | gls, lme, lmer | terms, coefficients, `permlist` | full-model fitted values and whitened residuals permuted together | F / Chisq, t |
| `permmodels(method = "freedman-lane")` | lm, aov | terms, coefficients | reduced-model residuals permuted (batched QR) | F, t |
| `bootmodels()` | lmer, glmer | terms, coefficients, `permlist` | parametric bootstrap from the *global* null (intercept + random effects) | F / Chisq, t |
| `permindex()` | any | design-based (randomisation) tests | permutations within blocks and/or of whole groups | user's |

Known weaknesses:

- `glm` with `method = "raw"` permutes the response, which also scrambles
  the effects of the other terms and the mean–variance relation (for a
  Poisson or binomial model the variance follows the mean). This is a test of
  the global null, not of one term when others have effects.
- `bootmodels()` simulates from the global null for every term. For one term
  the reference model should be the model without that term (as
  `pbkrtest::PBmodcomp()` does).
- `permanova.lmer()` (fixed effects) permutes the null model's fitted values
  together with the residuals. That scrambles the other terms' effects, like
  `method = "raw"` for `lm`, which our study showed is liberal with an
  extreme-value covariate (`dev/sim-freedman-lane.md`). Its statistic is the
  ML log-likelihood evaluated at REML estimates (consistent between observed
  and permuted data, so valid, but not the usual LRT; the printed table uses
  ML refits).
- No permutation test for GLMMs (`glmer`, `glmmTMB`) at all; GLMM fixed
  effects only through `bootmodels()` (global null), random effects not at all.
- Residuals are whitened with the Cholesky factor of V, which depends on the
  order of the data and makes the residual-permutation schemes liberal for
  within-block terms (found in step 1a; §3.3).
- `permlmer()` compared LRTs without rounding (fixed in step 1a).
- Defaults use 3 cores (`ncore = 3`); CRAN and CI allow 2.

## 3. Are the existing schemes right?

Two ideas are in the package: the test of random effects in `permlmer()`
(§3.1), and the maintainer's own scheme for fixed effects in mixed models,
used by `permmodels()` for lme/gls/lmer and by `permlmer()`/`permanova.lmer()`
when the fixed effects differ (§3.2). Both are checked here and both are
extended to GLMMs in §4.

### 3.1 Random effects: `permlmer()`

**What it does.** For H0 model `m0` and H1 model `m1` with the same fixed
effects, both fitted by REML:

1. V0 = Z Λ(θ0) Λ(θ0)' Z' + I, the covariance of y under H0 relative to σ²:
   H0's variance parameters θ0 are put into H1's structure and the tested
   components are set to 0.
2. Marginal residuals under H0: r = y − X β̂0.
3. Whiten: e = U⁻¹ r, with U the lower Cholesky factor of V0.
4. For each permutation π: y* = X β̂0 + U e[π] (re-correlate).
5. Refit both models to y*; statistic = REML LRT (negative values set to 0).
6. p = (#{LRT* ≥ LRT} + 1) / (B + 1).

**Why it works.** Under H0, Var(y) = σ² V0, so U⁻¹(y − Xβ) has covariance
σ² I: uncorrelated with equal variances. With normal errors and random
effects that makes the whitened residuals independent and identically
distributed, hence exchangeable, so permuting them leaves their joint
distribution unchanged. Re-correlating with U gives data with exactly the H0
covariance, in which the tested component is absent, so LRT* follows the
null distribution of the LRT. Under H1 the whitened residuals still carry the
tested structure (e.g. correlation within whole plots), which the observed
LRT picks up and the permutations destroy: that gives the power.

**It is a published method.** This is the restricted-LRT permutation test of
Lee and Braun (2012), which `permlmer()` already cites; Genstat's
VRPERMTEST implements the same idea. Lee and Braun build on Fitzmaurice,
Lipsitz and Ibrahim (2007), who permuted cluster labels in GLMMs with a single
random intercept, and Drikvandi et al. (2013) give a related test for several
variance components.

**Approximations (where it can go wrong).**

- Exchangeability is only "second-moment" (uncorrelated, equal variances)
  unless the data are normal (Commenges 2003). With skewed errors or skewed
  random effects each whitened residual mixes random effects and errors in
  different proportions, so they are not identically distributed.
- β̂0 and θ0 are estimated: the whitened residuals have covariance I − H
  (H the whitened hat matrix), not I. The same approximation as
  Freedman–Lane in `lm`; it matters most with few observations per parameter.
- If H0's fit is singular (θ0 = 0 for a remaining component) V0 is still
  valid.

### 3.2 Fixed effects: the maintainer's scheme

**What it does.** With fitted values xb = X β̂ and whitened residuals
e = U⁻¹(y − xb) (U U' = V, the fitted covariance), each permutation uses the
same index π for both:

  y* = xb[π] + U e[π]

and the model is refitted to y*. `permmodels()` (lme, gls, lmer) takes xb and
V from the full model and keeps every refit, so one set of refits gives
p-values for all terms, all coefficients, and (through `permlist`) all
pairwise comparisons and contrasts. `permlmer(m0, m1)` with different fixed
effects (and so `permanova.lmer()`) takes xb and V from the reduced model
and uses a likelihood-ratio statistic.

**What it amounts to.** With independent errors (V = I) it is exactly
permuting the response, y* = y[π]: the raw-data permutation of Manly, exact
under the global null. With correlated data it carries that idea over:
each unit's fixed part moves with its residual, while U keeps the
correlation structure (blocks, plots) in place. So it is the mixed-model
form of the "raw" scheme, and shares its properties:

- it tests the global null well (no fixed effect but the intercept);
- for one term when other terms have effects, the permuted data have the
  other terms' effects scrambled into noise, so the reference distribution
  is for a different nuisance configuration. With a studentised statistic
  (F, t) this usually matters little; it matters when a covariate has
  extreme values and is related to the tested term (Anderson and Legendre
  1999 for linear models);
- even under the global null it is approximate in mixed models (the
  re-correlated permuted residuals do not reproduce the data's joint
  distribution exactly unless V = I), which matters little in practice.

**Evidence so far.** In the Freedman–Lane study (`dev/sim-freedman-lane.md`)
this scheme kept 5% for `lmer` with an associated factor (5.7% normal, 6.0%
skewed; 300 data sets) and with skewed or heavy-tailed errors, but rejected
6.3% (1000 data sets) with a covariate with extreme values, where the
Freedman–Lane version (y* = xb0 + U0 e0[π], xb0 and U0 from the model
without the term) kept 5.0%. The new simulation adds the `permlmer()` /
`permanova.lmer()` form and edge cases (§3.3): the scheme was within the
range in most cases but 7.1% when no fixed effect exists (8 blocks, n = 48),
which theory did not predict and needs checking (§10, step 6).

**Follow-up.** The 7.1% was not chance (2000 new data sets: 6.2%; 16 blocks:
8.4%) and not caused by estimating V (true V: 6.5% and 8.1%). The permuted F
statistics were about 10% too small, in the numerator. The cause is the
Cholesky factor U: within a block, U applied to a constant gives a pattern
that rises along the rows (1.41, 1.93, 2.27, 2.52, 2.72, 2.89 for six units
with equal block and residual variances), and in these data A's levels sit
in fixed positions within each block, so the pattern lines up with A.
Whitening with the symmetric square root V^(1/2), which treats the units of
a block alike, gave 5.2% for the maintainer's scheme and 5.3% for
Freedman–Lane (2000 data sets), as the exact within-block permutation
(5.3%).

With a covariate with extreme values related to the tested term, the
symmetric root does not help the maintainer's scheme (6.4% vs 6.3%), while
Freedman–Lane with the symmetric root gives 4.6% (1000 data sets).

**Verdict.** Sound and useful — above all because one set of refits serves
every term and every pairwise comparison, which Freedman–Lane cannot do
without a separate null model per comparison. With symmetric-root whitening
it keeps its level under the global null; like raw permutation in `lm` it
stays slightly liberal (6.4%) with an extreme-value covariate related to the
tested term. The plan keeps it (as
`method = "raw"`, and for `permlist`), offers Freedman–Lane, and changes the
whitening of all residual-permutation schemes (including `permlmer()`) to
an order-invariant square root.

### 3.3 Simulation check

`dev/sim-perm-ranef.R` (results in `dev/sim-perm-ranef/`, 1000 data sets per
case, 99 permutations; results between 0.036 and 0.064 are consistent with
5%):

**Random effects, LMM** (`perm` = `permlmer()`; `χ̄²` = REML LRT with the
50:50 mixture):

| case | perm | χ̄² |
|---|---|---|
| split plot (6 blocks), normal | 0.059 | 0.046 |
| split plot, skewed errors and block effects | 0.049 | 0.046 |
| split plot, 20% missing | 0.054 | 0.049 |
| random slope (H0 random intercept) | 0.052 | 0.049 |
| 3 blocks | 0.052 | 0.040 |
| t3 errors and block effects | 0.039 | 0.035 |
| error SD 0.5/1/2 by subplot factor | 0.028 | 0.018 |
| block variance 0 (singular H0 fit) | 0.052 | 0.036 |
| plots of 1–4 units | 0.056 | 0.045 |
| crossed rows × columns | 0.064 | 0.046 |
| correlated random slope (2 parameters) | 0.055 | 0.044 |
| power: split plot, whole-plot SD 0.5 | 0.383 | 0.366 |
| power: 3 blocks, whole-plot SD 0.7 | 0.381 | 0.358 |

**Fixed effects, LMM** (`joint` = maintainer's scheme; `FL` =
Freedman–Lane; F of the term from refits of the full model; `Satt` =
Satterthwaite F):

| case | joint | FL | Satt |
|---|---|---|---|
| extreme covariate, 8 blocks (LRT form of `permanova.lmer()`) | 0.053 | 0.050 | 0.052 (χ² LRT) |
| extreme covariate, 8 blocks (F form; `sim-freedman-lane`) | 0.063 | 0.050 | 0.046 |
| extreme covariate, 3 blocks | 0.060 | 0.064 | 0.059 |
| associated factor, skewed errors | 0.036 | 0.041 | 0.041 |
| empty cells, large effect of the other factor | 0.055 | 0.050 | 0.046 |
| no fixed effects (global null) | **0.071** | 0.061 | 0.052 |
| group sizes 1/2/3 with SD 3/1/1 | **0.243** | **0.226** | **0.237** |
| power: extreme covariate, A effect | 0.673 | 0.633 | 0.632 |

**Follow-up: whitening** (`dev/sim-perm-followup*.R`; "within" = y permuted
within blocks, exact under the global null; 2000 data sets, range
0.040–0.060, except 16 blocks and the covariate case with 1000):

| case | joint, Cholesky | joint, reduced model | joint, true V | FL, Cholesky | joint, symmetric | FL, symmetric | within | Satt |
|---|---|---|---|---|---|---|---|---|
| global null, 8 blocks (`fu_lmm_global`) | 0.062 | 0.057 | 0.065 | 0.057 | — | — | 0.048 | 0.048 |
| global null, 16 blocks (`fu_lmm_global16`) | **0.084** | 0.071 | **0.081** | 0.070 | — | — | 0.057 | 0.052 |
| global null, 8 blocks (`fu2_lmm_global_sym`) | **0.067** | — | — | 0.060 | **0.052** | **0.053** | 0.053 | 0.049 |
| extreme covariate, 8 blocks (`fu2_lmm_lev_sym`) | 0.063 | — | — | 0.054 | **0.064** | **0.046** | (not exact here) 0.056 | 0.052 |

## 4. Extending both schemes to GLMMs on the working (link) scale

### 4.1 The working linear mixed model

A GLMM fitted by Laplace or adaptive quadrature can be linearised at a fit
(the basis of PQL, Breslow and Clayton 1993, and of SAS GLIMMIX
pseudo-likelihood, Wolfinger and O'Connell 1993). At the H0 fit, with
η0 = X β̂0 + Z b̂0 (without the offset) and μ0 = g⁻¹(η0 + offset):

- working response: z = η0 + (y − μ0) · dη/dμ
- working weights: w = a · (dμ/dη)² / V(μ0), with a the prior weights
  (number of trials for binomial) and V the family's variance function;
  for `glmer` these are `m@resp$sqrtWrkWt()^2`
- working model: z ≈ X β + Z b + ε, Var(b) = G, Var(ε) = φ · diag(1/w)
- so Var(z) ≈ V0 = Z G0 Z' + φ · diag(1/w)

This is exactly the V that `R2_glmm()` builds (`ZGZ' + diag(s²/w)`, with
`prior_weight` = the working weights). Points to get right, family by family:

| family | V(μ) | φ | note |
|---|---|---|---|
| binomial (lme4: y = proportion, a = trials) | μ(1 − μ) | 1 | binary data: z takes two values per mean |
| poisson | μ | 1 | |
| Gamma | μ² | 1/shape (σ² in lme4) | |
| gaussian, non-identity link | 1 | σ² | |
| nbinom2 (glmmTMB) | μ(1 + μ/θ) | 1 | θ = `sigma(model)`; `nbinom2()$variance` needs θ |
| nbinom1 (glmmTMB) | μ(1 + φ) | 1 | φ = `sigma(model)` |
| beta (glmmTMB) | μ(1 − μ)/(1 + φ) | 1 | `beta_family()$variance` omits 1/(1 + φ): check `R2_glmm()` |
| tweedie (glmmTMB) | φ μ^p | — | |

- G must be on the absolute scale: Λ Λ' for binomial/Poisson `glmer`,
  σ² Λ Λ' when the family has a dispersion; for glmmTMB the covariance
  matrices of `VarCorr()` (as `pm_re_terms()` and `pm_tmb_lambda()` give).
  Each family's V should be checked against `simulate()` moments in tests.
- Zero-inflated and `dispformula` glmmTMB models are out of scope at first.
- The whitening needs only the q × q factor (Woodbury: V0⁻¹ = D⁻¹ −
  D⁻¹ Z (G⁻¹ + Z' D⁻¹ Z)⁻¹ Z' D⁻¹ with D = φ diag(1/w)), as lme4 already
  computes; `permlmer()` would also gain from this for large n.

### 4.2 Permuting on the working scale

The same formulas as for linear mixed models, with z for y and V from §4.1:

- **Random effects (Lee–Braun):** linearise the H0 fit; marginal working
  residuals r = z − X β̂0; whiten e = U⁻¹ r (U U' = V0); permute;
  re-correlate: z* = X β̂0 + U e[π].
- **Fixed effects, maintainer's scheme:** linearise the full fit;
  z* = xb[π] + U e[π] with xb = X β̂ and e = U⁻¹(z − xb); one set of refits
  serves all terms, coefficients and pairs (a `permlist`).
- **Fixed effects, Freedman–Lane:** linearise the fit without the tested
  term; z* = xb0 + U0 e0[π].

One difference from linear models: in a GLMM the working variance of a unit
depends on its mean (w is computed from μ). Under the maintainer's scheme a
unit receives another unit's fixed part xb[π] but keeps its own working
variance, so mean and variance no longer match; Freedman–Lane keeps each
unit's mean and variance together. The pilot (§4.4) shows whether this
matters.

### 4.3 Can the permuted data be back-transformed to the response scale?

Not directly. g⁻¹(z*) is not a possible response: for binary data it lies
strictly between 0 and 1, for counts it is not a whole number, and its
variance no longer follows the family's mean–variance relation. A GLMM
refitted to it is a different model (with warnings for non-integer
binomial or Poisson data). Two routes are sound:

- **A. Stay on the working scale (recommended first).** Compute the
  statistic from working linear mixed models fitted to z* with weights w
  (`lmer(z* ~ ..., weights = w)` for H0 and H1), and the observed statistic
  in the same way from z. Observed and permuted statistics then come from the
  same procedure, which is what a permutation test needs. Each permutation
  costs two `lmer` fits instead of two GLMM fits. Zeng et al. (2015) proposed
  a permutation LRT computed from PQL working responses for a GLMM variance
  component and found it held its level; Chen et al. (2019, arXiv
  1906.03320; R package glmmVCtest) also test variance components on
  normalised PQL working responses. For fixed effects with V0 fixed the
  statistic is a GLS Wald F, so all permutations can be done at once by QR
  (as for `lm`), no refits at all. With V0 from the H0 fit this statistic is
  close to the score test, the statistic of the sign-flipping tests below.
- **B. Back-transform through the probability integral transform.**
  Randomised quantile residuals (Dunn and Smyth 1996) map each response to a
  standard normal value given its fitted mean: q = Φ⁻¹(u), u uniform
  between F(y − 1; μ0) and F(y; μ0). Permute q (with the random effects
  handled at cluster level: BLUPs rescaled to have covariance Ĝ0, Carpenter,
  Goldstein and Rasbash 2003, then permuted or re-drawn), and map back:
  y* = F⁻¹(Φ(q[π]); μ*) with μ* = g⁻¹(X β̂0 + Z b*). This gives valid
  responses of the same family, and the GLMM can be refitted to them. It is
  close to a semi-parametric bootstrap and costs two GLMM fits per
  permutation; it is a candidate for the simulation, not a first choice.

### 4.4 Pilot results (route A)

**Random effects, GLMM** (`perm` = Lee–Braun on the working scale;
`work χ̄²` = the working-scale LRT with the χ̄² reference, no permutation;
`Laplace χ̄²` = the usual GLMM LRT):

| case (split plot unless stated) | perm | work χ̄² | Laplace χ̄² | failed fits |
|---|---|---|---|---|
| binary, whole plots of 5 | 0.065 | 0.046 | 0.019 | 0 |
| binary, whole plots of 2 | **0.078** | 0.069 | 0.017 | 3 |
| binary, about 5% events | **0.358** | **0.483** | 0.023 | 137 |
| binomial, 10 trials | 0.049 | 0.042 | 0.013 | 0 |
| Poisson | 0.059 | 0.052 | 0.008 | 0 |
| Poisson, mean about 0.4 | 0.058 | 0.049 | 0.011 | 0 |
| Poisson, 3 blocks | 0.046 | 0.033 | 0.002 | 0 |
| negative binomial fitted as Poisson | **0.098** | 0.091 | **0.631** | 0 |
| power: Poisson, whole-plot SD 0.4 | 0.637 | 0.620 | 0.547 | 0 |

**Fixed effects, GLMM** (working scale; `joint` linearised at the full
fit, `FL` at the reduced fit; `LRT` = Laplace LRT, χ²):

| case | joint | FL | LRT | failed fits |
|---|---|---|---|---|
| binomial (10 trials), extreme covariate | 0.045 | 0.044 | 0.040 | 0 |
| Poisson, extreme covariate | 0.040 | 0.056 | 0.054 | 1 |
| binary, associated factor | 0.038 | 0.050 | 0.048 | 0 |
| binomial, 3 blocks, extreme covariate | 0.052 | 0.052 | 0.040 | 0 |
| Poisson, no fixed effects | 0.060 | 0.067 | 0.053 | 0 |
| binary, 5–10% events, associated factor | 0.040 | 0.070 | 0.071 | 15 |
| negative binomial fitted as Poisson, extreme covariate | 0.077 | 0.103 | **0.352** | 0 |
| power: Poisson, extreme covariate, A effect | 0.544 | 0.625 | 0.652 | 0 |

**Follow-up: binary data** (500 data sets, range 0.031–0.069; route B and
bootstrap refits with `nAGQ = 0`, observed and resampled alike):

| case | route B | parametric bootstrap | Laplace χ̄² | working-scale perm (main study) |
|---|---|---|---|---|
| binary, about 5% events | 0.062 | 0.040 | 0.026 | 0.358 |
| binary, whole plots of 2 | 0.056 | 0.056 | 0.026 | 0.078 |

Reading: on the working scale the Lee–Braun idea works for GLMM variance
components except for binary data with tiny clusters or rare events, where
the working response turns each event into a large outlier (about 1/μ) and
the linear approximation breaks down. There the quantile back-transform
(route B) or a parametric bootstrap is needed, and both work (follow-up
table). For fixed effects both
schemes behave, the maintainer's more conservatively. The mean–variance
mismatch of the joint scheme (§4.2) shows as lower power, not as a too-high
level.

## 5. What "most robust" means here

A scheme is robust if its type I error stays close to the nominal level
across: non-normal errors and random effects; unbalanced data and missing
cells; covariates with extreme values related to the tested term; few
clusters (4–8 blocks, as in field trials); for GLMs/GLMMs, overdispersion and
small binomial denominators; and moderate model misspecification. Among
robust schemes we prefer power, then speed (no refits where possible), then
simplicity. Two general results guide the choice:

- Permute under the reduced model (Freedman–Lane) rather than raw data when
  other terms have effects (Anderson and Legendre 1999; Anderson and
  ter Braak 2003; Winkler et al. 2014).
- Use a studentised (pivotal) statistic — F, t, Wald, LRT — not a raw
  coefficient or sum of squares: the test is then exact under
  exchangeability and still asymptotically valid when exchangeability fails
  (Chung and Romano 2013; DiCiccio and Romano 2017).

## 6. Literature by model class

- **Linear models.** Freedman–Lane is the usual recommendation (Anderson and
  Legendre 1999; Winkler et al. 2014); exact alternatives (Huh and Jhun 2001;
  Kherad-Pajouh and Renaud 2010) rotate the residuals into n − p
  exchangeable values. With heteroscedastic errors, sign flipping (wild
  bootstrap-like) with a studentised statistic is valid where permutation is
  not (Winkler et al. 2014).
- **Designed experiments, mixed ANOVA.** Restricted permutations within
  strata (whole plots for whole-plot terms, within whole plots for subplot
  terms) give exact tests for some terms of balanced designs (Anderson and
  ter Braak 2003); permuting reduced-model residuals usually has more power.
  Kherad-Pajouh and Renaud (2015) give an approximate residual permutation
  for any repeated-measures or mixed design (R package permuco, Frossard and
  Renaud 2021). `permindex()` already makes restricted permutations.
- **Linear mixed models, random effects.** Lee and Braun (2012) (=
  `permlmer()`); Drikvandi et al. (2013) for several components; Samuh et
  al. (2012) and El-Horbaty and Hanafy (2024) for robust variants with one
  component.
- **GLMs.** Potter (2005) permutes the residuals of the tested covariate
  regressed on the others and refits (R package glmperm, Werft and Benner
  2010); exact conditional logistic and Poisson regression (Mehta and Patel
  1995) for small data; sign flipping of effective score contributions
  (Hemerik, Goeman and Finos 2020; R package flipscores) is robust to
  overdispersion and heteroscedasticity and asymptotically exact, with
  standardised scores for small samples (De Santis et al., cited by Andreella
  2025).
- **GLMMs, random effects.** Fitzmaurice, Lipsitz and Ibrahim (2007)
  (permute cluster labels; single random intercept); Zeng et al. (2015)
  (permutation LRT on PQL working responses); Chen et al. (2019) (approximate
  RLRT on normalised working responses); parametric bootstrap of the LRT
  (Sinha 2009; `pbkrtest::PBmodcomp`, Halekoh and Højsgaard 2014).
- **GLMMs, fixed effects.** No established permutation method. Basso and
  Finos (2011) give exact multivariate permutation tests (fit per cluster;
  low power). Andreella (2025) applies score sign flipping at cluster level
  ("two-stage summary statistics"; an earlier arXiv version was withdrawn).
  Parametric bootstrap from the reduced model is the usual reference.
  Design-based randomisation tests (permuting treatment labels as they were
  randomised, refitting the model) are exact under the sharp null for any
  model class.

## 7. Candidates and current recommendation

"Prior" = before the study; "Now" = after the 37 cases of §3.3 and §4.4.

| Model | Hypothesis | Candidates | Prior | Now |
|---|---|---|---|---|
| lm, aov | terms, coefficients, contrasts | raw; **Freedman–Lane** (done); sign flip + HC studentised (heteroscedastic) | Freedman–Lane | Freedman–Lane (lm study) |
| gls, lme, lmer | fixed terms, contrasts | raw (maintainer's scheme); **FL whitened, refit**; FL whitened, V fixed (batched); design-restricted | FL whitened | **Freedman–Lane with symmetric-root whitening** for term tests (4.6–5.3%); the maintainer's scheme (symmetric root) for `permlist`, noting 6.4% with extreme covariates; help says KR/Satterthwaite F are as good for LMMs |
| lmer | random effects | **Lee–Braun** (`permlmer()`) | keep | keep (13 cases fine); tie fix done; move to symmetric-root whitening after a confirming run |
| glm | terms, coefficients | raw (current); FL on working scale (batched, score-like); Potter; **sign-flip effective score** | sign flip or FL working | not yet studied |
| glmer, glmmTMB | random effects | **working-scale LRT (route A)**; cluster-label permutation; route B; parametric bootstrap LRT | route A | route A for counts and binomial with trials; **route B for binary (0/1) responses** (bootstrap as an option); `nAGQ = 0` refits |
| glmer, glmmTMB | fixed terms | joint and FL on the working scale (weighted-LMM refits); V fixed (batched); route B; design-based; per-term parametric bootstrap | open | FL by default (power); joint as the robust option and for `permlist`; check V-fixed batching |

## 8. Design: one engine for all schemes

- **Hypotheses as reduced models or L matrices.** A term, a coefficient, a
  pairwise difference or a contrast is a linear hypothesis Lβ = 0; its
  reduced model is X N with N a basis of the null space of L. One engine then
  serves `permmodels()` tables, `predictmeans()` pairs and
  `contrastmeans()`.
- **A permutation object instead of a list of refits.** For Freedman–Lane
  the null distribution differs for each hypothesis, so a list of refits
  (`permlist`) cannot serve pairs. Instead `permmodels()` would return the
  permutation indices, the scheme and what is needed to recompute statistics;
  `predictmeans(permlist = )` would compute each pair's statistic under its
  own reduced model. Old `permlist` lists stay accepted.
- **Layers:** (1) `pm_working(model)`: z, w, φ, η, μ, offset, X, Z, G for
  lm/glm/lmer/glmer/glmmTMB/lme/gls (validated against `R2_glmm()` and
  `simulate()`); (2) `pm_whiten(V)`: an order-invariant square root W with W W' = V,
  not the Cholesky factor (§3.3). For V = D + A A' (D diagonal: I for LMMs,
  φ diag(1/w) on the GLMM working scale; A = Z Λ, n × q) take
  W = D^(1/2) M^(1/2) with M = I + D^(-1/2) A A' D^(-1/2); M^(1/2) =
  I + Q (diag(sqrt(1 + s²)) − I) Q' from the thin SVD D^(-1/2) A = Q S R',
  so the cost is O(n q²) and no n × n matrix is formed; W is unchanged by
  reordering the data. For LMMs (D = I) W is the symmetric root tested in
  §3.3; for GLMMs (D ≠ I) it is not symmetric, and step 1b checks it; (3)
  permutation indices, free or restricted (`permindex()`), with the seed
  handled by `pm_with_seed()`; (4) statistics: batched GLS F/t with V fixed,
  or refits (weighted `lmer`, GLMM) when variance parameters must be
  re-estimated; (5) p-values with failure counts and the `+1` rule.
- **API.** `permmodels(model, method = "auto")` picks the default for the
  model class (from §9); the named methods stay available. `permlmer()`
  accepts `glmerMod` and `glmmTMB` pairs, with `permglmer()` and
  `permglmmTMB()` as aliases. `permanova.lmer()` becomes a wrapper (and so
  works for GLMMs). `bootmodels()` gains a per-term reduced-model bootstrap.

## 9. Simulation study to choose the defaults

Extends `dev/sim-freedman-lane.R` and `dev/sim-perm-ranef.R`.

- **Designs**: completely randomised with a covariate (with and without
  extreme values); RCBD with 4 and 8 blocks; split plot (6 blocks); repeated
  measures with random slopes; unbalanced versions (20% missing); crossed
  random effects (one case).
- **Responses**: normal, skewed and heavy-tailed errors; skewed random
  effects; binary; binomial with 10 trials; Poisson; overdispersed Poisson
  (fitted as Poisson, and as nbinom2 in glmmTMB); beta (one case).
- **Hypotheses**: one fixed term (2–3 df) with another term having a large
  effect; one coefficient; one pairwise difference; one variance component
  (whole-plot, random slope).
- **Methods**: the candidates of §7 plus the parametric references (F with
  Kenward–Roger df, Wald, LRT with the 50:50 χ² mixture, parametric
  bootstrap from the reduced model).
- **Size**: 1000 data sets for type I error (Monte Carlo SE 0.007; range
  0.036–0.064), 500 for power; 99 permutations (exact level with the +1
  rule at 5%).
- **Rule**: a scheme qualifies if its type I error is within the Monte Carlo
  range in at least 90% of null cases and never above 0.075; among
  qualifiers choose by power, then time.
- **Compute**: about 1–2 s per data set for `lmer`-based methods and 5–20 s
  for GLMM refits. Run case by case, results saved per case (resumable);
  the large GLMM-refit cases are best run on a machine with more cores.

## 10. Steps

Each step is one message / one PR-sized piece, with tests, NEWS and a plan
update; R CMD check at the end of each.

1. **Check `permlmer()`, the maintainer's scheme, and both on the GLMM
   working scale** — **done** (§3.3, §4.4; 37 cases).
1a. **Follow-ups** — **done**: (i) the 7.1% came from Cholesky whitening;
   the symmetric root fixes it; (ii) route B and the parametric bootstrap
   fix binary GLMMs; (iii) rule: route B for 0/1 responses, the working
   scale for counts and binomial with trials (simpler than a threshold on
   events per cluster, and covers every failing case); (iv) the bootstrap
   reference in the main cases is left for the paper; (v) the tie bug in
   `permlmer()` fixed.
1b. **Second, smaller study** (can run while steps 2–5 are built): pairwise
   comparisons and contrasts (per-pair and familywise error for `permlist`
   from the joint scheme vs Freedman–Lane per pair vs t with KR df);
   `permlmer()` and the GLMM working scale with symmetric-root whitening
   (rerun 4–5 cases); glm without random effects; glmmTMB families (nbinom2,
   beta); route B vs bootstrap under overdispersion; parametric bootstrap
   reference for the paper.
2. **Engine, part 1 (no change in results)** — **done** (`R/perm-engine.R`):
   `pm_perm_stream()` (indices with a private seed, drawable in chunks),
   `pm_perm_map()` (refits, parallel, errors → failures), `pm_perm_drop()`,
   `pm_perm_p()` (rounded comparison, +1 rule, NAs left out),
   `pm_perm_tstats()` (|t| of contrasts from a permlist). `permmodels()`
   (raw and Freedman–Lane), `permlmer()`, `bootmodels()` and the pair and
   contrast p-values use it. `dev/snapshot-perm.R` (33 calls) and
   `dev/snapshot-predictmeans.R`: identical, except that a seed no longer
   moves the user's random number stream. Found and fixed: `cbind(R, F)`
   responses permuted only `F`.
2b. **Order-invariant whitening** (`pm_whiten()`, §8) in `permmodels()`
   (lme, gls, lmer) and `permlmer()` — **done** (`R/whiten.R`): symmetric
   root from q × q eigen decompositions (block roots for nlme). Confirmed
   with the package functions (`dev/sim-perm-confirm.R`, 1000 data sets):
   `permlmer()` split plot 0.059, crossed 0.054 (Cholesky 0.064);
   `permmodels()` global null `lmer` 0.048 (Cholesky 0.062–0.071), `lme`
   0.059. Much faster with crossed random effects (n = 3000: 0.11 s vs
   3.9 s; n = 20000: 21 s vs > 30 min). The 1.1.1 reference test compares
   `perm_lmer`/`perm_lmer_III` without their permutation parts.
3. **Freedman–Lane for any linear hypothesis in lm/aov** — **done**:
   `permmodels(method = "freedman-lane")` returns a permutation object
   (`perm`: X, z, the random number state of the permutations, coefficient
   names); `pm_perm_fl_t()` tests each row of an L matrix with its own
   reduced model (X N, N the null space of the row), regenerating the same
   permutations in chunks; `predictmeans()`/`contrastmeans()` use it.
   Checked against direct refits (pair 2 v 3 by merging levels). Also fixed:
   with a `permlist`, observed t used the KR-adjusted SE while the refits
   use the unadjusted one; now `vcov_adjust` is not used with a `permlist`.
   The default for lm/aov is decided in step 9 (`method = "auto"`).
3b. **`permindex()`** (maintainer's question: still used? worth keeping?)
   — **done**: it was not used inside the package, but it is exported and
   the natural source of restricted (design-based) permutations, so kept and
   fixed: indices for the data's own order (1.x: for a sorted copy), a
   separator in the duplicate check, an error for whole-group swaps of
   unequal sizes, and a seed that leaves the user's stream alone (same
   draws for sorted data). `permmodels()` (raw and Freedman–Lane) and
   `permlmer()` take its matrix as `perms` (`pm_perm_source()`); the lme/gls
   path maps it through its internal reordering (tested: same p-values for
   shuffled data).
4. **Working-scale layer `pm_working()`** — **done** (`R/working.R`) for
   glm, glmer, glmmTMB (no zero inflation, dispformula ~1; families poisson,
   binomial, nbinom1, nbinom2, Gamma, gaussian, beta): z, working variances
   on the absolute scale (offsets handled), A = Z L with ZGZ' = AA'. Checked:
   glmer equals the pilot's `linearise()`; glmmTMB equals glmer for the same
   model; each glmmTMB family's variance against `simulate()` (ratio
   1.00 ± 0.03). `R2_glmm()` was wrong for nbinom1, nbinom2 and beta (and
   ignored offsets): fixed with these variances.
5. **`permlmer()` for GLMMs** — **done** (`R/permglmm.R`), aliases
   `permglmer()` and `permglmmTMB()`: route A (working scale, weighted
   `lmer` refits, REML for random effects, ML for fixed effects;
   Lee–Braun / Freedman–Lane from the H0 linearisation) by default; route B
   (quantile, GLMM refits with `nAGQ = 0`) for 0/1 responses
   (`method = "auto"`). Observation-level random effects are refused on the
   working scale (not separable from the working residual; `bootmodels()`).
   Cluster-label permutation not added (restricted permutations via
   `perms` cover it).
6. **Mixed-model fixed effects** — **done** (`R/perm-fl.R`):
   `permmodels(method = "freedman-lane")` for lmer, lme, gls with V fixed at
   the fit (whitened by the symmetric root, then the batched `lm` code; no
   refits; statistics = the model's Wald F and t, unadjusted covariance;
   type II/III car-style, type I sequential in the GLS metric). The
   refit variant (V re-estimated per permutation) was 4.6–5.3% in §3.3;
   the V-fixed variant is checked by study 1b (pairs: `flpair`, `flterm`).
   Found and fixed: raw `permmodels()` for lmer, type I, compared lmerTest's
   type III observed F with lme4's sequential permuted F; now both are
   lmerTest's type I (Wald with Doolittle/Cholesky hypotheses of X'X).
   `permanova.lmer()` is unchanged (term by term with `permlmer()`); GLMM
   terms are tested with `permmodels(method = "freedman-lane")` or
   `permlmer()`. Also fixed: refits of models fitted inside a function with
   the formula or family given as a variable (`pm_refit()`).
7–8. **GLM and GLMM fixed effects** — **done** in the same code: the
   working-scale model of `pm_working()`, whitened with V = AA' +
   diag(working variances), then batched Freedman–Lane; the dispersion is
   estimated from the working residuals (robust to overdispersion: pilot
   7.7–10% against LRT 35%). Sign-flip score tests not added (not needed
   if study 1b part 2 confirms the level; then for a later version).
   `permmodels()` raw is refused for glmer/glmmTMB (Freedman–Lane only).
9. **`method = "auto"`**, documentation (help pages, the resampling
   vignette with a "which test for which model" table), NEWS, R CMD check, PR.

## 11. Is there enough evidence to decide the structure?

Yes for the overall structure and for most defaults; not yet for four
parts. The table says what each decision rests on.

| Decision | Evidence | Ready? | Missing |
|---|---|---|---|
| One engine: hypotheses as L matrices, a permutation object instead of a list of refits, shared indices and seeds, failure counts, rounded statistics | design; the tie bug found in `permlmer()` shows why statistics must be compared after rounding | **yes** | — |
| lm/aov terms: Freedman–Lane by default | Freedman–Lane study (2000 sets per case), literature | **yes** | pairwise comparisons (below) |
| LMM random effects: keep Lee–Braun (`permlmer()`), with the tie fix | 13 null cases, 2 power cases | **yes** | — |
| Whitening: an order-invariant (symmetric) square root of V instead of the Cholesky factor, in every residual-permutation scheme | 2 follow-up cases: 6.0–6.7% → 5.2–5.3%; theory (order dependence) | **yes** | a confirming rerun of `permlmer()` cases (step 1b) |
| LMM fixed terms: Freedman–Lane with symmetric whitening for term tests; the maintainer's scheme for `permlist` | 8 null cases + 4 follow-up cases (FL symmetric 4.6–5.3%; joint symmetric 5.2–6.4%) | **yes** | — |
| GLMM random effects: working scale for counts and binomial with trials; route B for 0/1 responses | 9 + 2 cases | **yes** | route B vs bootstrap under overdispersion (option, not structure) |
| GLMM fixed terms: which default (joint or FL) | 7 null cases, 1 power case | **no** | pairwise comparisons; parametric bootstrap reference; overdispersion fitted as nbinom (glmmTMB) |
| **Pairwise comparisons and contrasts by permutation** (`predictmeans(permlist = )`), the main use in predictmeans | none: only whole-term tests were simulated | **no** | a study of per-pair and familywise error (none, Holm, Tukey-type) for the joint scheme's `permlist` vs Freedman–Lane per pair vs t tests with KR df |
| glm (no random effects) | none | **no** | raw vs Freedman–Lane on the working scale vs sign flipping |
| glmmTMB families (nbinom1/2, beta, Gamma) | none; variance table in §4.1 unchecked | **no** | family-by-family check of V, then a few null cases |
| gls/lme (nlme) | none directly; same whitening as lmer | assumed | one or two cases to confirm |

So the build can start with the parts that are ready — the engine (step
2), order-invariant whitening (2b), Freedman–Lane for lm/aov (step 3),
`permlmer()` for GLMMs (step 5) — while the second, smaller study (step 1b)
covers pairwise comparisons, glm, glmmTMB families and the GLMM fixed-effect
default. The pairwise study matters
most: predictmeans is mainly used for comparisons of means, and the choice
between one set of refits (the maintainer's `permlist`) and a null model per
comparison (Freedman–Lane) is the central structural decision.

## 12. Decisions for the maintainer

1. Make Freedman–Lane the default for lm/aov once it also gives pairs
   (step 3)? Recommended: yes.
2. Extend `permlmer()` to GLMMs with `permglmer()`/`permglmmTMB()` aliases,
   rather than separate functions? Recommended: yes.
3. Own code for sign flipping and the working scale (no flipscores, permuco,
   RLRsim dependencies)? Recommended: yes, as with car and pbkrtest.
4. Default `ncore = 1` for all resampling functions? Recommended: yes.
5. glmmTMB zero-inflation and `dispformula` models out of scope at first?
   Recommended: yes.
6. Where to run the large GLMM simulation cases (here they take hours on
   4 cores)?

## 13. Risks

- The working-scale approximation fails for binary data with rare events
  (35.8% instead of 5% for a variance component) and is liberal with whole
  plots of 2 (7.8%). Route B (6.2%, 5.6%) and the bootstrap (4.0%, 5.6%) fix
  both; route B is the plan for 0/1 responses.
- Changing the whitening changes `permmodels()` and `permlmer()` results for
  mixed models (not only in the hard cases), so it is a visible change with
  a NEWS entry, not a refactor.
- Unequal variances between groups break every test studied (23–24%);
  permutation is no cure for a wrong variance model. The help should say so.
- Second-moment exchangeability (Commenges 2003) is all residual-based
  schemes can offer outside normal models; small-sample type I error must be
  checked case by case.
- GLMM refits fail or warn more often; failures must be counted and reported,
  and the p-value based on successful refits only.
- Run time: GLMM-refit schemes are slow; batched V-fixed schemes are fast but
  ignore variance re-estimation.

## 14. References

Checked against publisher or index records (title, journal, year) during
this review; full texts were not all available, so method details for some
come from abstracts:

- Anderson, M. J. and Legendre, P. (1999). An empirical comparison of
  permutation methods for tests of partial regression coefficients in a
  linear model. *J. Statist. Comput. Simul.* 62, 271–303.
- Anderson, M. J. and ter Braak, C. J. F. (2003). Permutation tests for
  multi-factorial analysis of variance. *J. Statist. Comput. Simul.* 73,
  85–113.
- Andreella, A. (2025). Robust inference for generalized linear mixed models:
  a "two-stage summary statistics" approach based on score sign flipping.
  *Psychometrika*, doi 10.1017/psy.2024.22.
- Basso, D. and Finos, L. (2011). Exact multivariate permutation tests for
  fixed effects in mixed-models. *Commun. Statist. Theory Methods*, doi
  10.1080/03610926.2011.627103.
- Breslow, N. E. and Clayton, D. G. (1993). Approximate inference in
  generalized linear mixed models. *JASA* 88, 9–25.
- Chen, S. T. et al. (2019). An approximate restricted likelihood ratio test
  for variance components in generalized linear mixed models. arXiv
  1906.03320.
- Commenges, D. (2003). Transformations which preserve exchangeability and
  application to permutation tests. *J. Nonparametric Statist.* 15, doi
  10.1080/1048525031000089310.
- Drikvandi, R., Verbeke, G., Khodadadi, A. and Partovi Nia, V. (2013).
  Testing multiple variance components in linear mixed-effects models.
  *Biostatistics* 14, 144–159.
- Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
  *J. Comput. Graph. Statist.* 5, 236–244.
- El-Horbaty, Y. S. and Hanafy, E. M. (2024). A Monte Carlo permutation
  procedure for testing variance components using robust estimation methods.
  *Statistical Papers* 65, 335–356.
- Fitzmaurice, G. M., Lipsitz, S. R. and Ibrahim, J. G. (2007). A note on
  permutation tests for variance components in multilevel generalized linear
  mixed models. *Biometrics* 63, 942–946.
- Frossard, J. and Renaud, O. (2021). Permutation tests for regression,
  ANOVA, and comparison of signals: the permuco package. *J. Statist.
  Software* 99(15).
- Hemerik, J., Goeman, J. J. and Finos, L. (2020). Robust testing in
  generalized linear models by sign flipping score contributions. *JRSS B*
  82, 841–864.
- Kherad-Pajouh, S. and Renaud, O. (2015). A general permutation approach for
  analyzing repeated measures ANOVA and mixed-model designs. *Statistical
  Papers* 56, 947–967.
- Lee, O. E. and Braun, T. M. (2012). Permutation tests for random effects in
  linear mixed models. *Biometrics* 68, 486–493.
- Potter, D. M. (2005). A permutation test for inference in logistic
  regression with small- and moderate-sized data sets. *Statist. Med.* 24,
  693–708.
- Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M. and
  Nichols, T. E. (2014). Permutation inference for the general linear model.
  *NeuroImage* 92, 381–397.
- Wolfinger, R. and O'Connell, M. (1993). Generalized linear mixed models: a
  pseudo-likelihood approach. *J. Statist. Comput. Simul.* 48, 233–243.
- Zeng, P. et al. (2015). Permutation-based variance component test in
  generalized linear mixed model with application to multilocus genetic
  association study. *BMC Med. Res. Methodol.* 15, 37.

From memory (to be checked before they are cited in help pages):
Carpenter, Goldstein and Rasbash (2003), *JRSS C* 52, 431–443; Chung and
Romano (2013), *Ann. Statist.* 41, 484–507; DiCiccio and Romano (2017), *JASA*
112, 1211–1220; Freedman and Lane (1983), *J. Bus. Econ. Statist.* 1,
292–298; Halekoh and Højsgaard (2014), *J. Statist. Software* 59(9); Huh and
Jhun (2001), *Commun. Statist. Theory Methods* 30, 2023–2032; Kherad-Pajouh
and Renaud (2010), *Comput. Statist. Data Anal.* 54, 1881–1893; Mehta and
Patel (1995), *Statist. Med.* 14, 2143–2160; Samuh et al. (2012), *Commun.
Statist. Simul. Comput.*; Sinha (2009), *Canad. J. Statist.*; Werft and
Benner (2010), *R Journal* 2(1), 39–43.
