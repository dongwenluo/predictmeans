# Review of the defaults (2026-11-13)

The maintainer's rule: each default is the most reliable method on the
evidence; matching SAS or Genstat is an option. This note puts each default
next to the evidence we have (our simulations, the SAS and Genstat runs, the
literature) and lists what to change.

## 1. The defaults and their evidence

| choice | default now | evidence | verdict |
|---|---|---|---|
| df, REML `lmer`, `lme`, `gls` | Kenward–Roger, observed information | `dev/sim-df/`: 4.2–6.6% everywhere (split plots, AR(1), UN, missing values); the expected information 6.2–8.7% with missing values. Kenward and Molenberghs (1998): with missing data only the observed information is valid. SAS `DDFM=KR` = ours in every value (16 + 5 runs). Luke (2017): KR and Satterthwaite hold the level for `lmer`. | keep |
| covariance with KR df | KR-adjusted (`vcov_adjust = TRUE`) | Kenward and Roger (1997); the t tests of the simulations used it (4.2–5.7%). SAS's `KR` adds second-derivative terms for AR(1); ours (`KR(FIRSTORDER)`) does not depend on the parameterisation (Kenward and Roger 2009). | keep |
| df, ML `lmer` | Satterthwaite | lmerTest; KR needs REML. | keep |
| df, ML `lme`/`gls` | Satterthwaite, ML (was nlme's df / residual df) | `dev/sim-df-ml/`: unstructured, residual df t 6.9%, F 9.3–11.7%; Satterthwaite 5.5%, 6.3–8.8%. Now as ML `lmer`. | changed (step 3) |
| Satterthwaite F tests | Fai–Cornelius, SAS's basis | liberal with several covariance parameters (UN 9.7–11.9%) | keep as an option; KR is the default |
| df, GLMMs | inner-outer (by grouping level) | `dev/sim-df-glmm/`: 3.6–8.1% against 8.4–15.2% (residual) and 7.3–12.7% (containment). Li and Redden (2015, binary cluster trials): between-within holds the level with 10 clusters, residual and containment inflated; inner-outer equals between-within for between-cluster effects and is less conservative within. The 6–8% left with 4 clusters per group is the Wald SE of the Laplace fit: `bootmodels()`. | keep |
| df, `aov` with `Error()` | strata (Genstat) | Genstat runs (Exams 3.1, 4.4, 8.1) | keep |
| df, `lm` | residual | exact | keep |
| df, `glm` with fixed dispersion (Poisson, binomial) | residual df (t) | no evidence either way; emmeans and SAS GLIMMIX use z | open (step 5: simulate first) |
| weights over other factors | equal (population marginal means) | Searle, Speed and Milliken (1980); SAS, emmeans, Genstat | keep |
| nested factors | equal over the combinations that exist | Genstat (Exam 4.3) | keep |
| **means over an empty cell** | average of the observed cells, with a warning (1.x; maintainer's decision) | an estimable function (ASReml/Genstat `!PRESENT`: "averaging only the estimable cells"); SAS and emmeans give `Non-est`/NA | keep the means |
| **comparisons of such means** | the difference of the two averages | **biased** — see section 2 | **change** (step 1) |

## 2. Empty cells

A mean that averages over its own observed cells is a fair summary of those
cells. The difference of two such means is not a fair comparison when they
average over different cells: it carries the effects of the other factor.
With A1 missing at B4, A1 − A2 estimates mean(B1..B3) − mean(B1..B4) of the B
effects as well as the A difference.

`dev/sim-empty-cells.R` (3 × 4 factorial, 3 replicates, `y ~ A * B`, no A
effect, 2000 data sets; a true 5% lies in 3.9–6.1%), % of true nulls rejected
for A1 vs A2 (`dev/sim-empty-cells/results.md`):

| empty cells | B4 differs by | present (default now) | common cells | additive model |
|---|---|---|---|---|
| A1:B4 | 0 / 1 / 2 / 4 σ | 5.1 / 8.3 / 19.0 / 59.2 | 5.4 / 5.1 / 5.6 / 4.1 | 4.6 / 5.6 / 5.7 / 4.1 |
| A1:B4, A2:B3 | 0 / 1 / 2 / 4 σ | 5.3 / 10.0 / 27.6 / 78.0 | 4.7 / 5.7 / 4.4 / 5.4 | 4.4 / 6.4 / 4.3 / 4.9 |
| A1:B4, interaction at A2:B4 | 1 / 2 / 4 σ | 8.2 / 17.6 / 57.0 | 5.2 / 4.7 / 5.3 | 5.7 / 6.2 / 7.5 |
| two empty, interaction | 1 / 2 / 4 σ | 9.8 / 26.2 / 77.1 | 4.8 / 4.3 / 5.4 | 5.7 / 7.3 / 13.8 |

- **Present (now):** right only when the other factor has no effect; with a
  B effect of 2σ it rejects 19–28% of true nulls.
- **Common cells** (each pair compared over the levels of the other factors
  observed for both; a Type IV-like contrast): 4.1–5.7% everywhere, with or
  without interaction. It answers a clear question: the difference averaged
  over the cells both levels have.
- **Additive model:** right without interaction, biased with it (up to 13.8%).

Literature: the population marginal mean over all cells is not estimable
(Searle, Speed and Milliken 1980), which is why SAS and emmeans give NA;
Milliken and Johnson (1984, p. 185) call SAS's Type III hypotheses with empty
cells "the worst hypotheses to consider ... no reasonable way to interpret
them", and analyse missing cells with the cell-means model, comparing
observed cells; ASReml and Genstat average over the present cells only on
request (`!PRESENT`, Welham et al. 2004), and their guide says predictions
over empty cells "cannot be averaged in a meaningful way".

**Proposal (step 1):** keep the means (the maintainer's decision); compare
two means that average over different cells on the cells both have, and say
so (a message naming the levels used; a pair with no common cells is NA). The
joint test of the term (`pm_joint()`) uses the levels common to all the means
it compares. The warning is given only when a mean actually averages over an
empty cell, and names those means: now it is also given for the cell means
themselves (`"A:B"`), where nothing is averaged.

## 3. One decision in one place

- `df_term()` (`R/misc.R`) maps `type` to the information and checks REML
  itself, a second copy of `pm_df_method()`: make it call `pm_df_method()`.
- `pm_anova()` (the `permmodels()` tables) uses Kenward–Roger with the
  expected information, for `car::Anova()` parity: by the maintainer's rule
  its default should be the package's (observed), with the expected one as an
  option.
- The empty-cell rule (which cells count, nesting, the warning) sits inside
  `Kmatrix()`; move it into one function used by `Kmatrix()` and the
  comparisons of step 1.

## 4. Steps

1. Empty cells: comparisons over common cells, warning only where a mean
   averages over an empty cell (changes `pairs` for such means) — **done**
   (maintainer's yes, 2026-11-13; `R/empty-cells.R`).
2. `df_term()` through `pm_df_method()`; empty-cell rule in one function (no
   change in output; snapshot check) — **done**.
3. ML `lme`/`gls`: Satterthwaite from the ML Hessian, as ML `lmer` (changes
   df of ML nlme fits) — **done**. `dev/sim-df-ml/results.md`: unstructured
   covariance, t 5.5% (residual df 6.9%), F 6.3–8.8% (9.3–11.7%); AR(1) the
   same either way. F tests of ML fits stay liberal with several covariance
   parameters: REML and Kenward–Roger are the remedy.
4. `pm_anova()` tables: observed information by default (changes the KR F
   tests of the `permmodels()` tables slightly where data are unbalanced).
5. `glm` with fixed dispersion: simulate t (residual df) against z, then
   decide.

## References

- Gilmour, A. R., Thompson, R. and Cullis, B. R. (1995). Average information
  REML. *Biometrics* 51, 1440–1450.
- Kenward, M. G. and Molenberghs, G. (1998). Likelihood based frequentist
  inference when data are missing at random. *Statistical Science* 13,
  236–247.
- Kenward, M. G. and Roger, J. H. (1997). Small sample inference for fixed
  effects from restricted maximum likelihood. *Biometrics* 53, 983–997.
- Kenward, M. G. and Roger, J. H. (2009). An improved approximation to the
  precision of fixed effects from restricted maximum likelihood.
  *Computational Statistics and Data Analysis* 53, 2583–2595.
- Li, P. and Redden, D. T. (2015). Comparing denominator degrees of freedom
  approximations for the generalized linear mixed model in analyzing binary
  outcome in small sample cluster-randomized trials. *BMC Medical Research
  Methodology* 15, 38.
- Luke, S. G. (2017). Evaluating significance in linear mixed-effects models
  in R. *Behavior Research Methods* 49, 1494–1502.
- Milliken, G. A. and Johnson, D. E. (1984). *Analysis of Messy Data, Vol. 1:
  Designed Experiments.* Van Nostrand Reinhold.
- Searle, S. R., Speed, F. M. and Milliken, G. A. (1980). Population marginal
  means in the linear model: an alternative to least squares means. *The
  American Statistician* 34, 216–221.
- Welham, S., Cullis, B., Gogel, B., Gilmour, A. and Thompson, R. (2004).
  Prediction in linear mixed models. *Australian and New Zealand Journal of
  Statistics* 46, 325–347.
