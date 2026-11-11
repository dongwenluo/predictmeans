# Permutation and resampling tests in predictmeans: review and plan

Written 2026-11-11. Goal (from the maintainer): keep in predictmeans the most
robust permutation test for each model class, and extend `permlmer()` to
`glmer` and `glmmTMB` models by permuting on the transformed (link) scale.
This file reviews what the package does now, checks the existing ideas
against theory and simulation, surveys the literature, and sets out the steps.

Contents: 1 Summary · 2 What the package has now · 3 Is `permlmer()` right?
· 4 Permuting GLMMs on the working scale · 5 What "most robust" means · 6
Literature by model class · 7 Candidates and current recommendation · 8
Design: one engine · 9 Simulation study to choose the defaults · 10 Steps ·
11 Decisions for the maintainer · 12 Risks · 13 References

## 1. Summary

(To be completed when the simulation in §3 and §4.4 has finished.)

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
- Defaults use 3 cores (`ncore = 3`); CRAN and CI allow 2.

## 3. Is `permlmer()` right? (random effects)

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

**Fixed effects (the `ref = FALSE` branch, used by `permanova.lmer()`).**
This part is not in Lee and Braun. Permuting H0's fitted values together
with the residuals is the "raw" scheme; the Freedman–Lane version keeps the
fitted values in place, y* = X β̂0 + U e[π], i.e. the same formula as for
random effects.

**Simulation check.** `dev/sim-perm-ranef.R` (results in
`dev/sim-perm-ranef/`, 1000 data sets per case, 99 permutations; results
between 0.036 and 0.064 are consistent with 5%):

(Running; results to be added.)

## 4. Permuting GLMMs on the working (link) scale

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

### 4.2 Permuting

As in `permlmer()`: marginal working residuals r = z − X β̂0, whiten
e = U⁻¹ r (U U' = V0), permute, re-correlate: z* = X β̂0 + U e[π].

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

### 4.4 Pilot results (route A, variance components)

(Running; results to be added.)

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

Recommendations before our study ("prior"); the study (§9) decides the
defaults.

| Model | Hypothesis | Candidates | Prior |
|---|---|---|---|
| lm, aov | terms, coefficients, contrasts | raw; **Freedman–Lane** (done); sign flip + HC studentised (heteroscedastic) | Freedman–Lane |
| gls, lme, lmer | fixed terms, contrasts | raw (current); **FL whitened, refit**; FL whitened, V fixed (batched); design-restricted | FL whitened |
| lmer | random effects | **Lee–Braun** (`permlmer()`) | keep |
| glm | terms, coefficients | raw (current); FL on working scale (batched, score-like); Potter; **sign-flip effective score** | sign flip or FL working, by simulation |
| glmer, glmmTMB | random effects | **working-scale LRT (route A)**; cluster-label permutation (single intercept); parametric bootstrap LRT | route A if the pilot holds |
| glmer, glmmTMB | fixed terms | FL working scale (V fixed, batched); FL working scale with weighted-LMM refits; route B; design-based; per-term parametric bootstrap | open: by simulation |

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
  `simulate()`); (2) `pm_whiten(V)`: sparse or Woodbury Cholesky; (3)
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

1. **Check `permlmer()` and pilot route A** (this step): theory (§3),
   simulation (`dev/sim-perm-ranef.R`).
2. **Engine, part 1 (no change in results):** `pm_whiten()`, permutation
   indices, failure counting; `permmodels()`, `permlmer()` and
   `bootmodels()` moved onto it; snapshot test shows identical output.
3. **Freedman–Lane for any linear hypothesis in lm/aov:** pairs and
   contrasts in `predictmeans()`/`contrastmeans()` from a permutation object;
   then the decision on making Freedman–Lane the default for lm/aov.
4. **Working-scale layer `pm_working()`** for glm, glmer, glmmTMB (no zero
   inflation): weights, dispersion, G on the absolute scale, checked family by
   family against `simulate()` and `R2_glmm()`; fix `R2_glmm()` if the beta
   or negative binomial variance is wrong.
5. **`permlmer()` for GLMMs (route A)**, aliases `permglmer()` and
   `permglmmTMB()`; cluster-label permutation for a single random
   intercept as an option; tests against the pilot.
6. **Mixed-model fixed effects:** FL whitened (refit and V-fixed) for
   lme/gls/lmer; `permanova.lmer()` on the engine; simulation; default.
7. **GLM fixed effects:** sign-flip effective score (own code, no new
   dependency), FL on the working scale, Potter; simulation; default.
8. **GLMM fixed effects:** FL working scale (V fixed; weighted-LMM refits),
   route B, design-based, per-term parametric bootstrap; simulation; default.
9. **`method = "auto"`**, documentation (help pages, the resampling
   vignette with a "which test for which model" table), NEWS, R CMD check, PR.

## 11. Decisions for the maintainer

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

## 12. Risks

- The working-scale approximation is poorest for binary data with few
  observations per cluster (PQL's known weakness). The permutation test only
  needs the statistic to rank data sets sensibly and the permuted data to
  mimic the null, so it may still hold its level — the pilot and study
  decide.
- Second-moment exchangeability (Commenges 2003) is all residual-based
  schemes can offer outside normal models; small-sample type I error must be
  checked case by case.
- GLMM refits fail or warn more often; failures must be counted and reported,
  and the p-value based on successful refits only.
- Run time: GLMM-refit schemes are slow; batched V-fixed schemes are fast but
  ignore variance re-estimation.

## 13. References

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
