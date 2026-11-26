# Part B4: clm nominal/scale, population-average clmm probabilities, few-cluster tests, AR1 x AR1 per trial (review note, 2026-11-26)

For the maintainer's decisions before any code. Four steps, each its own PR, in the
order proposed below. Packages here: ordinal 2026.7.26, VGAM 1.1-14, nlme.

## G1. clm with nominal or scale effects (B3c)

**What the fits are.** `clm(y ~ x, nominal = ~ z)`: P(Y <= j) = F(theta_j(z) - x'b), cut
points that depend on z: the non-parallel cumulative model (checked here: the same
log-likelihood as `vglm(cumulative(parallel = FALSE ~ z))`, -86.20855 on wine).
`clm(y ~ x, scale = ~ s)`: P(Y <= j) = F((theta_j - x'b) / exp(s'zeta)): a different spread
for each level of s, so the latent means of two scale levels are not on one scale.

**Proposal.** Both go to the class-probability path of B3b (`R/categorical.R`):
`mean_type = "prob"` (logit, default), `"prob_diff"`, `"mean_class"`, `"link"` (for nominal
only: the J - 1 predictors; refused with scale effects, where a predictor is not linear
in the coefficients). One new piece: B3b writes each predictor as linear in the
coefficients, eta(x) = H (theta * x[col]); a scale effect divides by exp(s'zeta), so the
parts get a general form, eta(x, theta) with its Jacobian (analytic for clm: d eta /
d theta, d eta / d b, d eta / d zeta = -eta s'), used by the same averaging and delta
method. Crossing curves (nominal) and separation are handled as in B3b.

**References.** nominal: the same model by `vglm(cumulative(parallel = FALSE))`, whose
results B3b already checks against `predict()` and a numerical delta method (to 1e-6);
scale: `predict.clm(type = "prob", se.fit = TRUE)` (clm's own delta method) at the cells
of the model; both: numDeriv delta method on `coef(clm)`.

## G2. Population-average class probabilities for clmm

**Now.** clmm probabilities are for a typical subject (random effects at 0), with the
advice `ordinal_typical_subject` saying they are more extreme than the average over
subjects.

**The quantity.** P_pop(Y = j | x) = integral of P(Y = j | x, u) over u ~ N(0, Var(z'u)).
The random part enters the linear predictor as one normal number z'u, with variance
z' Sigma z (sum of the intercept variances for crossed or nested intercepts; with
random slopes it depends on x), so a one-dimensional Gauss-Hermite rule is exact for
any random-effect structure (30 nodes: error below 1e-8 for sd up to 3 on the logit
scale; checked against 1e6 Monte Carlo draws in the tests).

**Interface. Proposal:** `bias_adjust = TRUE` with `mean_type = "prob"` or
`"mean_class"`, as for GLMM means: in predictmeans `bias_adjust = TRUE` already means "the
mean over the clusters instead of the typical cluster" (`pm_back_label()`). The pairs
stay on the latent scale (subject-specific comparisons, valid for every subject); the
class table gives population-average probabilities with delta-method SEs that include
the uncertainty of the random-effect variance (`vcov(clmm)` has `ST1`). The advice
`ordinal_typical_subject` then says the other way round (how to get the typical subject).

**References.** Monte Carlo integration (1e6 draws) for the probabilities; numerical
delta method (numDeriv over the fixed and ST parameters) for the SEs; the limit sd -> 0
equals the typical subject; for a probit link the integral is closed form,
Phi(eta / sqrt(1 + sigma^2)) (an exact check).

## G3. Tests for clmm with few clusters

**Now.** Inner-outer df (B2) for a between-cluster treatment: 3.7% at the 5% level in one
setting with 6 clusters (z: 11%), but 7.3% and 7.75% in the other two settings of that
simulation, so with few clusters it can be too liberal as well as too conservative; an
exact or resampling test is the robust answer.

**Options.**
- (a) **Cluster permutation** for between-cluster terms: permute the treatment labels
  among clusters (restricted permutations, `permindex()`), refit, compare the Wald or LR
  statistic. Exact under exchangeability of clusters; with 6 clusters split 3/3 only 20
  permutations (smallest p 0.05), so it is useful from about 8 clusters.
- (b) **Parametric bootstrap LRT**: simulate from the null fit, refit both models. Any
  term; about the nominal size when the model is right; costs 2 x B clmm fits.
- **Proposal:** (b) in `bootmodels()` for clmm (any term), and (a) in `permmodels()` for
  terms constant within clusters, with the smallest attainable p in the notes; both
  through the existing engines (`R/perm-engine.R`, `R/boot-parts.R`), so `ncore` and
  identical results for any `ncore` come with them. Evidence: a simulation of size and
  power at 6, 8 and 12 clusters against inner-outer t and z (`dev/sim-clmm-few/`).

## G4. AR1 x AR1 with a phi per trial

**Proposal.** `corAR1xAR1(form = ~ row + col | trial, by_group = TRUE)`: a pair (phi_r,
phi_c) (and a nugget share) for each trial, as fitting each trial alone but with shared
treatment effects; for trial-specific residual variances combine with `weights =
varIdent(~ 1 | trial)` (said on the help page). Parameters 2 (or 3) x number of trials,
named by trial. AR2 and other orders stay out unless asked for.

**References.** Exact: the joint fit with per-trial phi, `varIdent(~ 1 | trial)` and
trial-specific fixed effects (`yield ~ trial / variety`) has the log-likelihood and phi
of the separate per-trial fits (sum of their REML log-likelihoods up to the constant).
Genstat: a two-trial run (Slate Hall twice, or Slate Hall and another trial) if the
maintainer wants a second reference.

## Order and effort

G1 (half a day: reuses B3b) -> G2 (half a day) -> G4 (a few hours) -> G3 (a day, mostly
the simulation). Each: tests, the Rule no. 1 review, R CMD check, PR.

Also open from before (the maintainer decides): inner-outer df of aggregated binomial
GLMMs, rows or trials; a simulation of their size can settle it (half an hour).

## Questions for an outside check

1. Check: is `clm(y ~ x, nominal = ~ z)` the same model as
   `VGAM::vglm(y ~ x + z, cumulative(parallel = FALSE ~ z))`, so their class
   probabilities agree?
2. Check: for a cumulative link mixed model, is the population-average P(Y <= j | x)
   the integral of F(theta_j - x'b - v) over v ~ N(0, z' Sigma z), so that a
   one-dimensional Gauss-Hermite rule suffices for any random-effect structure?
3. Check: for a cumulative probit mixed model with random intercept variance sigma^2,
   is the population-average P(Y <= j | x) = Phi((theta_j - x'b) / sqrt(1 + sigma^2))?
4. Check: with a cumulative logit mixed model and 6 clusters (3 per treatment), is a
   permutation test of a between-cluster treatment limited to a smallest p of 1/20?
5. Check: does `ordinal::clmm`'s `vcov()` give the random-effect parameters on the
   log standard deviation scale?
