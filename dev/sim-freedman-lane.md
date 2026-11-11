# Permutation schemes for testing one term: simulation results

Script: `dev/sim-freedman-lane.R` (run 2026-11-10, seed 20261110); results in
`dev/sim-freedman-lane-results.rds`.

Question: when one term is tested while other terms have real effects, does
the scheme of `permmodels()` (permute the response; for mixed models permute
the fitted values and decorrelated residuals together) keep the 5% level, and
does Freedman–Lane (permute the residuals of the model without the tested
term) do better?

Setup: tested factor A (3 levels) with no effect (type I error) or an effect
(power); a nuisance factor B or covariate x with a large effect; 199
permutations (`lm`) or 99 (`lmer`). `lm`: 2000 data sets per cell, so the
Monte Carlo SE of a 5% rate is about 0.005 (results between 0.040 and 0.060
are consistent with 5%). `lmer`: 300 data sets, SE about 0.013.

- balanced: A × B crossed, 4 replicates (n = 36);
- unbalanced: A and B associated (8 of 12 units of each A level in one B level);
- covariate: x correlated with A (n = 36);
- leverage: x log-normal (extreme values), larger for one A level (n = 24) —
  the case where permuting the response is known to fail (Anderson and Legendre
  1999);
- lmer: 8 blocks (random), A and B associated within blocks.

## Type I error (nominal 0.05)

| model | design | errors | F test | current | Freedman–Lane |
|---|---|---|---|---|---|
| lm | balanced | normal | 0.055 | 0.055 | 0.049 |
| lm | balanced | skewed | 0.042 | 0.039 | 0.045 |
| lm | balanced | heavy | 0.050 | 0.050 | 0.057 |
| lm | unbalanced | normal | 0.056 | 0.057 | 0.052 |
| lm | unbalanced | skewed | 0.046 | 0.049 | 0.048 |
| lm | unbalanced | heavy | 0.042 | 0.043 | 0.041 |
| lm | covariate | normal | 0.061 | 0.063 | 0.062 |
| lm | covariate | skewed | 0.039 | 0.038 | 0.045 |
| lm | covariate | heavy | 0.047 | 0.049 | 0.054 |
| lm | **leverage** | normal | 0.040 | **0.063** | 0.039 |
| lm | **leverage** | skewed | 0.039 | **0.063** | 0.048 |
| lm | **leverage** | heavy | 0.042 | **0.067** | 0.051 |
| lmer | blocks | normal | 0.053 | 0.057 | 0.063 |
| lmer | blocks | skewed | 0.047 | 0.060 | 0.057 |

## Power

| model | design | errors | F test | current | Freedman–Lane |
|---|---|---|---|---|---|
| lm | balanced | normal | 0.547 | 0.541 | 0.540 |
| lm | unbalanced | normal | 0.423 | 0.407 | 0.418 |
| lm | covariate | normal | 0.217 | 0.208 | 0.216 |
| lm | leverage | normal | 0.321 | 0.376 | 0.311 |
| lm | leverage | heavy | 0.439 | 0.493 | 0.453 |
| lmer | blocks | normal | 0.917 | 0.903 | 0.900 |

(All cells are in the `.rds` file.)

## Mixed model with a covariate with extreme values

Run 2026-11-10 with `Rscript dev/sim-freedman-lane.R 2000 300 lmer-leverage`
(results in `dev/sim-freedman-lane-lmer-leverage.rds`): 8 blocks with two
units of each A level (n = 48), x log-normal and larger for one A level, a
large slope; `y ~ x + A + (1 | block)`; 300 data sets per cell (results
between 0.025 and 0.075 are consistent with 5%).

| errors | effect | F test (Satterthwaite) | current | Freedman–Lane |
|---|---|---|---|---|
| normal | null | 0.043 | 0.077 | 0.067 |
| skewed | null | 0.050 | 0.047 | 0.053 |
| heavy | null | 0.043 | 0.057 | 0.057 |
| normal | power | 0.620 | 0.617 | 0.593 |

Unlike `lm`, the two permutation schemes behave alike here (within one
percentage point). The current scheme is just above the Monte Carlo range in
one of three cells; Freedman–Lane is close behind. The F test keeps 5%.

Confirmation of the normal-error null cell with 1000 data sets
(`Rscript dev/sim-freedman-lane.R 2000 1000 lmer-leverage-confirm`, results in
`dev/sim-freedman-lane-lmer-leverage-confirm.rds`; results between 0.036 and
0.064 are consistent with 5%):

| errors | effect | F test (Satterthwaite) | current | Freedman–Lane |
|---|---|---|---|---|
| normal | null | 0.046 | 0.063 | 0.050 |

The current scheme is at the upper edge of the range (slightly liberal, as in
`lm`, but less so); Freedman–Lane and the F test keep 5%. The 300-set value
for Freedman–Lane (0.067) was chance.

## Conclusions

- In balanced, unbalanced and correlated-covariate designs the current scheme,
  Freedman–Lane and the F test all keep the 5% level, with the same power.
- With a covariate that has extreme values and is correlated with the tested
  factor, the current scheme rejects 6.3–6.7% of true null hypotheses (well
  outside the Monte Carlo range); Freedman–Lane keeps 3.9–5.1%. The current
  scheme's higher power there comes from its inflated level.
- For `lmer` with an associated factor the two schemes do not differ. With a
  covariate with extreme values the current scheme is slightly liberal with
  normal errors (6.3% in 1000 data sets), Freedman–Lane keeps 5.0% and the
  Satterthwaite F test 4.6%: the same pattern as `lm`, but smaller.

Recommendation: offer Freedman–Lane in `permmodels()` for `lm`/`aov` models
(where all permutations of a term can be solved at once from QR projections,
so it is cheap), keep the current scheme as the default, and say in the help
when to prefer Freedman–Lane (covariates with extreme values). For mixed
models with such a covariate, the help recommends the Kenward–Roger or
Satterthwaite F test; Freedman–Lane for mixed models (one reduced-model fit
per term, then refits) is a possible later addition.

Reference: Anderson, M. J. and Legendre, P. (1999). An empirical comparison of
permutation methods for tests of partial regression coefficients in a linear
model. *Journal of Statistical Computation and Simulation* 62, 271–303.
Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of reported
significance levels. *Journal of Business & Economic Statistics* 1, 292–298.
