# Part B3: VGAM::vglm (and nnet::multinom) (review note, 2026-11-25)

**Decided 2026-11-25:** E1-E6 agreed as proposed. B3a done (`R/ordinal.R`, `test-ordinal.R`). B3b done
2026-11-26 (`R/categorical.R`, `test-categorical.R`), with one change to E2 chosen by the maintainer: the
default `mean_type = "prob"` compares each class probability on the logit scale (log odds ratios of the class
against the others, limits inside 0-1, as a binomial glm); `"prob_diff"` gives E2's differences of
probabilities; `"link"` E2's option (a).

For the maintainer's decisions before any code. Packages here: VGAM 1.1-14, nnet,
ordinal, emmeans 2.0.4 (which has no vglm method; it supports nnet::multinom).

## What vglm fits

A vglm has M linear predictors eta_1, ..., eta_M (M = J - 1 for J classes), each
x'b_j, with constraint matrices that tie coefficients across predictors. On the
wine data (`rating`, 5 classes), checked here:

| family | predictors | coefficients | same fit as |
|---|---|---|---|
| `cumulative(parallel = TRUE)` | logit P(Y <= j) = theta_j + x'b | one b, J - 1 intercepts | `clm`, `polr`, but **sign flipped** (b = -2.503 for warm; clm +2.503); `reverse = TRUE` flips it back |
| `cumulative(parallel = FALSE ~ x)` | logit P(Y <= j) = theta_j + x'b_j | b varies with j for x (partial proportional odds) | `clm(nominal = ~ x)` (refused by predictmeans now) |
| `propodds()` | = cumulative(parallel = TRUE, reverse = TRUE) | | `clm` |
| `acat()`, `cratio()`, `sratio()` | adjacent-category / continuation ratio logits | parallel or not | none in predictmeans |
| `multinomial(refLevel = r)` | log P(Y = j) / P(Y = r) = x'b_j | b for each non-reference class | `nnet::multinom` (emmeans supports it) |

Only the parallel models have one latent scale. For the others each mean is a
vector (one value per predictor, or per class).

## Decisions

**E1. Scope and order.** Proposal:
- B3a: parallel `cumulative()` / `propodds()` (and `acat`, `cratio`, `sratio` with
  `parallel = TRUE`): the B1 machinery, with the sign and the link of each family
  read from the fit, so a positive latent difference always means higher classes
  more likely. Reference: the same model by `clm()` (identical MLE) to 1e-6.
- B3b: non-parallel models and `multinomial()`, through the class probabilities
  (E2). Together with `nnet::multinom` (the plan's "later"), since emmeans covers
  it and it is the same model as `vglm(multinomial)`.
- Not: `vgam()` (smooths), `vglmm`/`rrvglm`, families that are not for a
  categorical response (they stop with a message).

**E2. What a "mean" is when there is no single latent scale (B3b).** Options:
- (a) Means of each linear predictor eta_j, and pairs within each predictor (a
  `by` column, as `atvar`): the model's own scale, every test exact on it, but
  J - 1 tables, and for multinomial each depends on the reference class.
- (b) Class probabilities as the means (mean_type = "prob" the only type), pairs
  as differences of probabilities with delta-method SEs, letters per class.
- Proposal: (b) as the default for B3b: it does not depend on the reference class
  or the coding, it is what users report, and emmeans does the same for multinom
  (`mode = "prob"`); (a) on request (`mean_type = "link"`) for the model-scale
  tests. Pairs of probabilities are Wald tests on the probability scale: one
  advice line says so (they are less exact near 0 and 1).

**E3. Degrees of freedom.** z (no dispersion), as polr/clm (pm_df_rule). vglm
with a dispersion family is not categorical and stops (E1).

**E4. Sign and labels.** The summary line of B1 ("log odds ratios of a higher
class") for parallel cumulative fits whatever `reverse =` was; for acat/cratio/
sratio the line names their scale (adjacent-category / continuation-ratio log
odds). Proposal: always report so that a positive difference means higher classes
more likely, and say which logit it is.

**E5. References.** B3a: `clm` on the same data (identical model) for cumulative,
and `clm(link = )` for other links; for acat/cratio/sratio, the probabilities of
`predict(type = "response")` and SEs against a numerical delta method
(`numDeriv` on the coefficients). B3b: `emmeans` on the equivalent `nnet::multinom`
fit (same MLE, different optimiser: 1e-4) and `predict.vglm()`; non-parallel
cumulative against `predict()` and the numerical delta method.

**E6. Other functions.** `residplot()`, `pm_check()`, `R2_model()`, `varcomp()`,
`permmodels()`, `bootmodels()`: stop with the class list (as B1), until each is
added and checked.

## Effort

B3a: small (adapter methods, sign, families; about B1's size). B3b: larger (a
mean is a vector; the pairs, letters and plots work per class), one or two days
of the B1 kind with its review.
