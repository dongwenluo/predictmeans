# Part B0: more model classes, ordinal models first (review note, 2026-11-25)

**Decided 2026-11-25:** the maintainer agreed D1-D6 as proposed. B1 (`polr`, `clm`) done: `R/ordinal.R`,
`test-ordinal.R`; `mean_type = "prob"` / `"mean_class"` instead of a new `type` argument (one argument for
what the means are).

For the maintainer's decisions before any code. Packages here: MASS 7.3 (`polr`),
ordinal 2025 (`clm`, `clmm`), VGAM, nnet, emmeans 1.11 (reference).

## What an ordinal model estimates

`polr` and `clm` fit P(Y <= j) = F(theta_j - eta), with eta = x'beta (no intercept)
and cut points theta_1 < ... < theta_{J-1}; `clmm` adds random effects to eta. On
the wine data (`rating ~ temp + contact`, 5 classes) emmeans gives for temp:

| mode | cold | warm | what it is |
|---|---|---|---|
| `latent` | -1.33 (0.36) | 1.17 (0.35) | eta: mean of the latent variable, logistic scale |
| `mean.class` | 2.36 (0.13) | 3.48 (0.15) | expected class number, sum_j j P(Y = j) |
| `prob` (per class) | 5 probabilities each | | P(Y = j) |
| `exc.prob` / `cum.prob` | 0.34 / 0.66 | 0.62 / 0.38 | P(Y > j) / P(Y <= j), averaged over j |

Differences of latent means are log cumulative odds ratios: the one scale on which
the model's effects are linear and the same for every cut point.

## Decisions

**D1. Default "mean".** Proposal:
- `predictmeans()` means and pairs on the latent scale (eta), as for a GLM on the
  link scale. Every part of the package works there unchanged (pairs, letters,
  adjustments, contrasts, `pm_joint()`, permutation tests).
- The response scale: `type = "prob"` gives the table of class probabilities per
  mean (with delta-method SEs), and `type = "mean_class"` the expected class.
  The tests stay those of the latent scale, as for GLMs (back-transformed means
  are descriptions; NEWS of 2026-11-15).
- `mean_class` assumes equally spaced classes, so it is never the default.

**D2. Averaging over other factors.** Proposal: average eta (the link scale), then
transform, the rule predictmeans uses for GLMs. Averaging probabilities instead is
what `weights`/`bias_adjust` style options do; leave for later.

**D3. Degrees of freedom.** `polr`/`clm`: asymptotic (z), as emmeans and the
fixed-dispersion GLM rule (simulation `dev/sim-df-glm/`). `clmm`: z first; the
inner-outer df of GLMMs only after a small simulation shows it is needed.

**D4. `clmm` probabilities.** Conditional on random effects at 0 (a typical
subject), as emmeans; marginal (integrated over the random effects) with
`bias_adjust = TRUE`, as for GLMM means, with the same advice line.

**D5. Order.** B1 `polr` + `clm` (adapter methods only: `pm_coef`/`pm_vcov`
without the cut points, `pm_basis` without the intercept, sign eta = x'beta,
`pm_linkinv`; probabilities through `R/response-scale.R`). B2 `clmm`. B3 `vglm`
(cumulative and multinomial; several linear predictors, a larger change).
`nnet::multinom` later.

**D6. References.** emmeans to 1e-6 (`latent`, `mean.class`, `prob`) for
`polr`/`clm`, 2% for `clmm`. One SAS run (PROC LOGISTIC / GENMOD cumulative logit,
PROC GLIMMIX for clmm) if the maintainer wants a second reference.

**What other functions do** (B1 to B3): `anova()` for these classes: Wald
chi-square (Types II/III on eta); `residplot()`/`pm_check()`: surrogate or
randomized quantile residuals, or a clear stop; `R2_model()`, `varcomp()` (clmm),
`permmodels()`: stop with the class list until each is added and checked.
