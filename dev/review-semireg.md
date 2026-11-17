# Semi-parametric modelling in predictmeans: review and upgrade plan

Step 1 of the upgrade asked for on 2026-11-17 (plan:
`dev/predictmeans-2.0-plan.md`, row "Semi-parametric upgrade"). Aim, in the
maintainer's words: make semi-parametric regression robust, flexible, fast,
easy to write and powerful, and show what only predictmeans can do.

Decisions already given:

1. Smooths are written **in the model formula**
   (`y ~ A*B + s(x, by = A) + (1|Block)`). The unpenalised (fixed) part of
   each smooth is added to the fixed effects automatically.
   `smoothZ = list(smZ(...))` keeps working.
2. `s(x, f1*f2)` is an **ANOVA-like split**: a common curve plus deviation
   curves for `f1`, `f2` and `f1:f2`, each with its own smoothing variance.
   `s(x, f1:f2)` is one curve per cell.
3. Bases come from **mgcv's constructors plus the package's own O'Sullivan
   and thin-plate (`Ztps`) bases**. mgcv moves from Suggests to Imports (it
   is a recommended package that ships with R); fits of 1.x stay reproducible.

## 1. The present code

| Where | What it does | Problems |
|---|---|---|
| `R/semireg.R` (403 lines) | Adds a placeholder `(1 | sm)` term per smooth block, fits `lmer`/`glmer` once to get the structure, replaces `Zt`, fits again with `mkLmerDevfun()`/`optimizeLmer()`, then converts to an lmerTest class | Two fits; the data are put into the calling frame by name (`assign(as.character(lmerc$data), data)`), which breaks when `data` is an expression; numeric variables are centred, all but the first response name, so `cbind(succ, fail)` responses fail (open question O4); one iid variance per Z block (no shared or structured penalties); negative binomial by a two-stage theta; needs lmerTest |
| `R/semireg_tmb.R` (279) | The same for glmmTMB (`doFit = FALSE`, replace `Z`, `fitTMB()`) | A near copy of `semireg()`; the parsing and the covariance of the coefficients are kept in step by hand |
| `R/smZ.R` (212) | Eleven basis types (`ZOSull`, `Ztps`, `ns`, `bs`, splines2 types, `smspline`); `by = factor`, `group = TRUE` for a variance per level | The `by`/`group`/missing-value code is repeated five times; `by` takes one factor only; no numeric `by` (varying coefficient), no tensor product, no shared-variance factor smooth; knots and range are kept as attributes and the stored call is evaluated again for prediction |
| Fixed part | Written by the user | Nothing checks it. If `x:f` is left out of the formula with `smZ(x, by = f)`, the group slopes are shrunk to a common one without any message |
| `R/semipred.R` (583), `pm_semi_design()` | Predicts by evaluating the stored `smZ()` call on the grid; curves, bands, differences (`pm_diff_table()`) | Bands are conditional on the smoothing parameters: in `dev/sim-curve-difference.R` the band of a difference between two curves covered 0.96 for a linear model but 0.92-0.94 pointwise and 0.87 simultaneous for Poisson and binary models (open question O5); long plotting branches |
| `R/semireg-means.R` | `predictmeans()` on a fit, via `pm_semi_design()` | Normal-based tests only (no KR/Satterthwaite df for the parametric terms) |

## 2. gamm4, and what predictmeans can do beyond it

gamm4 is not installed in this sandbox (CRAN not reachable), so this part
follows its documented method (Wood and Scheipl, gamm4 manual; Wood 2017,
section 7.7) and is to be checked against its source when it can be installed
(`dev/` only, never a dependency). The mgcv functions it relies on were checked
here (mgcv 1.9-1):

| Smooth | `smoothCon()` | `smooth2random()`: fixed columns `Xf` | random blocks (Z columns) |
|---|---|---|---|
| `s(x)` | 1 smooth | 1 (the linear term) | 1 (8) |
| `s(x, by = f)` | 2 smooths | 1 each (`x` within the level) | 1 each (8) |
| `s(x, f, bs = "fs")` | 1 | 0 | 3 (16, 2, 2): curves, and per-level intercepts and slopes as random |
| `t2(x, f, bs = c("cr", "re"))` | 1 | 0 | 2 (6, 4) |

### Method to adopt

- `mgcv::interpret.gam()` splits the formula into parametric terms and smooth
  specifications.
- `mgcv::smoothCon(absorb.cons = TRUE)` builds each smooth, and
  `mgcv::smooth2random()` re-parameterises it into fixed columns `Xf` (the
  penalty's null space) and random blocks with iid coefficients.
  - This is what guarantees the fixed terms of each smooth.
  - The level means of a `by` factor are not in `Xf`. mgcv and gamm4 expect the
    user to write `f` as a main effect; predictmeans will add it.
- `lme4::lFormula()` on a formula with one placeholder factor per random
  block, `Zt` replaced, then `mkLmerDevfun()` / `optimizeLmer()`. This is
  **one fit**, and no data are written into the caller's frame.
- `mgcv::PredictMat()` predicts from the stored smooth object: no stored call
  is evaluated again.

### gamm4's limits that predictmeans can lift

| gamm4 | predictmeans after the upgrade |
|---|---|
| Two objects (`$mer`, `$gam`); `predict()` only | Marginal means of factors with smooths in the model, pairwise comparisons, letters, contrasts (`predictmeans()`, `contrastmeans()`) |
| No difference between two curves with its band | `semipred(contr = )` / `pm_diff_table()`: model scale, ratio or odds ratio, simultaneous band, ranges where curves differ |
| Normal-based tests for parametric terms | Kenward-Roger / Satterthwaite df through `pm_inference()` |
| `te()` not allowed | `te()` rewritten as `t2()` with a message |
| lme4 families only | glmmTMB engine: negative binomial, beta, zero inflation, dispersion models |
| Bands conditional on smoothing parameters | Correction for their uncertainty (`unconditional = TRUE`), checked by simulation |
| No test that a curve is needed beyond a line, or that curves differ by a factor | Variance-component tests by permutation (`permlmer()`) or parametric bootstrap (`bootmodels()`) |
| No crossed-factor smooth syntax | `s(x, f)`, `s(x, f1:f2)`, `s(x, f1*f2)`, `s(x1, x2, f)` |
| No R² or variance shares | `R2_model()`, `varcomp()` per smooth |

## 3. Smooth syntax

| Written | Meaning | Fixed part added | Built as |
|---|---|---|---|
| `s(x)`, `s(x, k =, bs =)` | one curve | `x` | mgcv; `bs = "os"` the own O'Sullivan basis, `bs = "ztps"` the own thin plate |
| `s(x, by = f)` | a curve per level, a variance each (`id =` shares one) | `f`, `f:x` | mgcv `by` |
| `s(x, by = z)`, `z` numeric | varying coefficient | `z`, `z:x` | mgcv |
| `s(x, f, bs = "fs")` | random curves per level, one variance | none (all random) | mgcv |
| `s(x1, x2)`, `t2(x1, x2)` | surface | `x1`, `x2` (and `x1:x2` for `t2`) | mgcv |
| `te(x1, x2)` | surface | as `t2` | rewritten as `t2()`, with a message |
| `s(x, f)` | common curve + a deviation curve per level | `x`, `f`, `f:x` | expanded by predictmeans |
| `s(x, f1:f2)` | one curve per cell | `f1:f2`, `f1:f2:x` | `by = interaction(f1, f2)` |
| `s(x, f1*f2)` | common curve + deviations for `f1`, `f2`, `f1:f2`, a variance each | `x`, `f1*f2`, `(f1*f2):x` | expanded term by term |
| `s(x1, x2, f)` | common surface + deviations per level | as `s(x, f)` | a 2-D basis per piece |

How the fixed part and the deviations are built:
- The expansion is done before mgcv sees the terms.
- Each piece's `Xf` gives its fixed columns, so the fixed part is always
  complete.
- Fixed terms the user also wrote are not doubled: a column that is already
  in the span of the others is dropped.
- Deviation curves are identified as mgcv recommends for a `by` smooth next to
  a main smooth: fixed part with sum-to-zero contrasts, penalised deviations.
  A test checks that the common curve plus the deviations equals the
  per-level curves.

`smZ()` keeps its interface. Its repeated code becomes one helper, and it
gains `by = ~ f1:f2` and a numeric `by`.

## 4. HRW (Harezlak, Ruppert and Wand 2018): what there is to reproduce

The package's 45 demo scripts do not mostly use `lme`:

| Fitting function | Demos |
|---|---|
| `lme` (nlme) | 1: `maleGrowthIndiananlme` (group-specific curves: mean curve + black-vs-white contrast curve, `pdIdent` blocks for O'Sullivan Z, per-subject `pdSymm(~x)` + `pdIdent` curve) |
| `mgcv::gam` | 15 (`WarsawAptsSimpSemi`, `WarsawAptsSimpSemiInt`, `WarsawAptsGAMfit`, `BostMortGAMfit`, `BostMortFacByCurv`, `CaSchoolGAMfit`, `CaSchoolStepFit`, `OFPGAMfit`, `SydneyDisplaySophis`, `ozoneDisplayConvHull`, `ozone3Dspin`, `returnsDrvFig`, ...) |
| `rstan` (Bayesian mixed-model form of the same penalised splines, O'Sullivan Z) | 20 (`femSBMDbayes`, `indonRespirBayes`, `maleGrowthIndianaBayes`, `margAddMod`, `proteinMargAddInt`, `WarsawAptsSimpFacByCurv`, `coralAna`, `npReg*`, ...) |
| other (`vgam`, `rqss`, `svm`, `fda`) | 9: outside the scope (quantile, SVM, functional regression) |

The book's chapter 4 writes the penalised spline as a linear mixed model
fitted by `lme`, but the shipped code does that once. The reproduction set
(step 6) is therefore:

| HRW example | Their fit | predictmeans fit | What is compared |
|---|---|---|---|
| `maleGrowthIndiananlme` | `lme` with `pdIdent`/`pdSymm` blocks | `semireg(y ~ black*x + s(x, by = black) + (x | id) + s(x, id, bs = "fs"))` | variance components, mean and contrast curves, bands; lme vs lmer time |
| `femSBMDbayes` | Stan (group-specific curves by ethnicity) | `s(age, ethnicity)` + random curves per subject | posterior means vs REML curves |
| `indonRespirBayes` | Stan (binary, smooth of age, random intercepts) | `semireg(..., family = binomial)`; glmmTMB engine | curve and odds ratios |
| `WarsawAptsSimpSemi`, `WarsawAptsSimpFacByCurv`, `WarsawAptsGAMfit` | `gam`, Stan | `s(construction.date, by = district)` | curves vs `gam(method = "REML")` |
| `margAddMod`, `proteinMargAddInt` | Stan (additive mixed model; factor-by-curve interaction) | `s(x1) + s(x2, f)` with random effects | curves, interaction deviation |
| `ozoneDisplayConvHull`, `SydneyDisplaySophis` | `gam` with a bivariate smooth | `s(longitude, latitude)` / `bs = "ztps"` | surface vs `gam` |
| `BostMortFacByCurv` | `gam` binomial factor-by-curve | `s(dir, by = self)`, binomial | curves vs `gam` |

### 4.1 Results (step 6)

The maintainer supplied golden examples: Wand's course scripts fitting
penalised splines with `nlme::lme` (and `glmmPQL`, `gamm`, `gam`), in
`dev/semireg-golden/original/`; data in `inst/extdata/semireg/`.
`make-reference.R` reruns each original model (its own `ZOSull.r` and
`Ztps.r`) and stores the results in
`tests/testthat/fixtures/semireg-golden.rds`; `compare.R` refits them with
`semireg()`, smooths written in the formula. The bundled `ZOSull` differs
from the package's by a rotation of the columns only (Z Z' equal to 1e-10),
so the fits are the same.

| Example | Model | semireg formula | logLik semireg / reference | max fitted diff (SD of y) | time semireg (nloptwrap / bobyqa) vs reference |
|---|---|---|---|---|---|
| addMix | spline + binary x + random intercept (lme) | `y ~ x + s(t, bs = "os", k = 15, xt = list(range = c(0, 1))) + (1 \| id)` | -234.098 / -234.098 | 1e-7 | 0.3 / 0.6 s vs 0.08 s |
| copper | thin-plate (x, y) with given knots + O'Sullivan z (lme) | `s(xcoord, ycoord, bs = "ztps", xt = list(knots = kn)) + s(zcoord, bs = "os", xt = list(knots = zk, range = rz))` | -207.242 / -207.242 | 1e-6 | 0.3 / 0.3 s vs 0.07 s |
| growthIndiana | black boys: mean curve + subject curves, `pdSymm(~age)` + `pdIdent` spline (lme) | `s(age, bs = "os", k = 20, xt = list(range = c(5.5, 20))) + (s(age, k = 10, xt = ...) \| idnum)` | -935.745 / -935.745 | 2e-6 | 1.4 / 2.9 s vs 2.1 s |
| growthIndianaCont | all boys: mean + black-white contrast curve + subject curves (lme; = HRW `maleGrowthIndiananlme`) | `s(age01, bs = "os", k = 20) + s(age01, by = black, bs = "os", k = 20) + (s(age01, k = 10) \| idnum)` | 6386.102 / 6386.102 | 2e-7 | 8.1 / 11.2 s vs 31 s |
| ozoneCalifLMM | three O'Sullivan smooths (lme) | `s(x1, bs = "os", k = 15, xt = list(range = r1)) + ...` | -1017.000 / -1017.000 | 2e-5 | 0.5 / 0.6 s vs 0.10 s |
| scallop | thin-plate surface of asinh(catch), given knots (lme) | `s(longitude, latitude, bs = "ztps", xt = list(knots = sk))` | -308.020 / -308.020 | 5e-6 | 0.07 / 0.07 s vs 0.02 s |
| wageTUnion | binary, O'Sullivan (glmmPQL) | `s(wage, bs = "os", k = 15, xt = ...)`, `family = binomial` | Laplace vs PQL | 0.013 (probability) | 0.7 s vs 0.1 s |
| SBMDfem | `gamm(s(age) + ethnicity, random = ~1 \| id)` | `s(age) + black + hispanic + white + (1 \| idnum)` | equals `gam(REML)` to 1e-7; gamm's `lme` 2.4e-3 SD away (stopped short) | | 0.1 s vs 0.1 s |
| ozoneCalifGAM | `gam`, GCV | `s(x1) + s(x2) + s(x3)` | equals `gam(REML)` to 2e-7; GCV fit 0.11 SD away (other criterion) | | 0.1 s vs 0.02 s |

Findings:

- Every `lme` example is reproduced exactly (the REML optimum: lme and lmer
  agree to 1e-7 in log-likelihood). The addMix script's streamlined SE of the
  binary effect (0.0678) equals `vcov(fit$semer)`.
- HRW's group-specific curves (`pdSymm(~ x)` per subject + a spline
  deviation per subject with one variance) had no formula form: `bs = "fs"`
  gives uncorrelated intercepts and slopes (log-likelihood 2731.8 vs
  2747.1 with `(x | id)` added on top, which duplicates them). Added
  `(s(x, k) | g)`: `(x | g)` + the wiggly part of the smooth per level, one
  variance (`R/semireg-smooth.R`); population curves leave it out.
- The map from a smooth's basis to its random part (for prediction) was
  found by regression on the data, which fails when the covariate has fewer
  distinct values than basis columns (copper depth every 5 m) or an `fs`
  subject has few observations. It is now `smooth2random()`'s own
  transformation (exact), with the regression and an `fs` per-level map as
  fallbacks.
- The script's ethnicity test (REML likelihood ratio between models with
  different fixed effects) is not valid in general; `anova()` gives the
  Kenward-Roger F test (p = 6.1e-7 here; the LRT gave 4.6e-7).
- Speed: lmer's time is the optimizer. `bobyqa` (semireg's default since
  1.x) needs 500-3000 deviance evaluations here because the smooth blocks'
  theta lie far from the start (theta = 50-70 after mgcv's penalty scaling);
  `nloptwrap` needs a third of them, with the same optimum. O10, decided by
  the maintainer: `nloptwrap` is `semireg()`'s default (also for GLMMs, where
  binomial and Poisson fits reached the same optimum 1.3-4 times faster).
- Step 8: each smooth block is scaled as a random intercept (sum of squares
  n) before the fit. The model is the same (the block's variance is free),
  but lme4's start (theta = 1) is then near the optimum whatever the
  covariate's units. mgcv's own penalty scaling had left
  `s(x, by = f, id = 1)` on Dialyzer at a worse local optimum (-362.5 vs
  -355.6, the `smZ(by = f)` fit), and the fits got faster (times above,
  after the change; growth curves 18 s -> 8 s).
- PQL (the script) and Laplace (`semireg()`, glmer) differ by up to 0.013 in
  the fitted probabilities of the binary example: different
  approximations, not an error.

## 5. Cases that only predictmeans handles (step 7; simulated data and HRW data)

1. Marginal means and pairwise comparisons of factor levels at chosen
   covariate values, with a smooth in the model and Kenward-Roger df.
2. The difference of two curves with a simultaneous band and the covariate
   ranges where they differ, on the model scale or as a ratio (counts).
3. `s(x, f1*f2)`: which variance component is needed (does the `f1:f2`
   deviation curve matter?), by permutation.
4. A split-plot design with curves over time: treatments on whole plots,
   smooth over time per treatment, KR df for treatment comparisons.
5. Counts with negative binomial / zero inflation and a smooth (glmmTMB
   engine), with rate ratios between curves.
6. `R2_model()` and `varcomp()`: the share of each smooth.

### 5.1 Results (step 7)

All six cases run on the current code (`tests/testthat/test-semireg-cases.R`,
vignette `smooth.Rmd`, section "What predictmeans adds"). One gap was found
and closed:

- **Means at a covariate value had normal-based tests.** `predictmeans()` on
  a `semireg()` fit gave z tests (df = Inf) even in a split-plot design,
  where `anova()` gives Kenward-Roger df for the fixed effects. The means
  include the smooths' random coefficients, so the df must be those of a
  combination of fixed and random coefficients: `pm_sr_satt_ddf()`
  (`R/semireg-engine.R`), Satterthwaite with the prediction error variance
  L C(psi) L' as a function of psi = (theta, sigma), as SAS does for
  estimates that involve random effects. For fixed effects alone it equals
  lmer's Satterthwaite df. Gaussian lme4-engine fits only; generalized fits
  keep z.
- Type I error (`dev/sim-semireg-df.R`; 6 blocks, 3 treatments on whole
  plots, 11 times, no treatment effect, a1 vs a2 at t = 5, 500 data sets,
  Monte Carlo SE 0.010):

  | plot SD | z (before) | Satterthwaite t | median df |
  |---|---|---|---|
  | 0.2 | 0.052 | 0.046 | 27 |
  | 0.4 | 0.068 | 0.046 | 15 |

- `semipred()` bands of gaussian fits were still normal-based; since
  2026-11-17 they are t limits with these df at each point (one gradient
  for all points: `pm_sr_satt_ddf()` takes many rows at once) and, for
  simultaneous bands, a multivariate t (`pm_sim_crit(df = )`), the
  smoothing parameters held at their estimates (design variance
  components and residual only). Coverage of a treatment's curve in the
  split plot: normal 0.91, t 0.94-0.95 pointwise, 0.96-0.97 simultaneous;
  `dev/sim-semipred-df.R`.
- `anova()`'s likelihood ratio tests of the smooth variances (halved
  p-values, chi-square 0:1) are conservative (`dev/sim-semireg-lrt.R`;
  `s(x, f1*f2)`, 200 data sets per scenario, Monte Carlo SE 0.015):

  | true shapes | reject f1 | reject f2 | reject f1:f2 |
  |---|---|---|---|
  | one shape | 0.010 | 0.015 | 0.010 |
  | differ by f1 | 1.000 | 0.025 | 0.030 |

  The variance was estimated at zero in about 40% of the null data sets.
  The exact null distribution of the restricted LRT (Crainiceanu and
  Ruppert, 2004) or a permutation test (`permlmer()`) would bring the level
  closer to 5%. Closed by the tests of variance components (plan row O12):
  `anova()` and `varcomp(test = TRUE)` use the exact restricted LRT where
  it is exact (the only random effect) and the mixture otherwise;
  `test_method = "bootstrap"` gives a parametric bootstrap. For these
  deviation curves every reference stays conservative (0% at 5% in
  `dev/sim-vc-defaults.R`, scenario G), with full power.

## 6. Steps

1. This note.
2. **One engine** (`R/semireg-engine.R`) — **done**. Result check: 12 example fits (gaussian, ML, binomial, Poisson, `bs` basis, by-factor and group smooths, Dialyzer, indonRespir, glmmTMB Poisson and negative binomial) give the same coefficients, covariance matrices and `semipred()` output as before (log-likelihoods within 1e-7); total time 19 s -> 14 s. Found and fixed on the way: lmerTest's Satterthwaite df and `varcomp()` used the deviance of the placeholder model (now equal to `lme`); a zero smoothing variance stopped the fit; `cbind()` responses; data written into the caller's frame. Was: `semireg()` and `semireg_tmb()` run on
   it with the existing `smoothZ` form only, and must reproduce current
   results: the `semi_*` reference cases, `test-semipred.R` and
   `test-glmmtmb.R` pass unchanged. One fit instead of two (timed).
3. **Formula smooths**: `s()`, `by`, `fs`, `t2`/`te`, `os`, `ztps`, with
   automatic fixed terms; mgcv to Imports. Checked against
   `mgcv::gam(method = "REML")`.
   **Done.** `R/semireg-smooth.R`. Set-up by `gam(fit = FALSE)` (bases,
   knots, identifiability constraints as a gam), each smooth given the design
   `PredictMat()` produces (so that prediction uses the same columns), then
   `smooth2random()`; each random block's map from `PredictMat()` columns is
   kept for prediction. Fixed terms are found by checking which candidate
   terms span the smooth's unpenalised columns (the variables, then their
   interaction, then squares; times the `by` variable). Results: fitted values
   equal `gam(method = "REML")` in 11 cases (`s(x)`, `cr`, `ps`, `cs`,
   `by` factor, `id`, numeric `by`, `fs`, `s(x) + s(z)`, `s(x, z)`, with a
   random intercept); `semipred()`/`predictmeans()` equal `predict.gam()`;
   `bs = "os"`/`"ztps"` equal the `smZ()` fits exactly. Found on the way:
   variables of a smooth must not be centred (a numeric `by` multiplies the
   smooth); `semireg_tmb()` scaled the covariance by `sigma()`, which is theta
   for a negative binomial (fixed). For `t2()` the design of
   `gam(fit = FALSE)` and `smoothCon()` is not in the span of `PredictMat()`
   (mgcv 1.9-1; a fitted gam is consistent): the sum-to-zero constraint is
   not absorbed in that design. Closed 2026-11-17: `t2()` is set up without
   the constraint (`smoothCon(absorb.cons = FALSE)`), so its unpenalised
   part is exactly `1, x1, x2, x1:x2` (with mgcv's constraint on the whole
   product 5% of the penalised space leaks into it, and a random intercept
   then shifts the population surface: 0.04 on a test surface); the fit
   equals `gam()` with that basis and penalties (`paraPen`), and `gam()`'s
   `t2()` without other random effects (`test-semireg-tensor.R`). `te()` is
   fitted as `t2()` with a message; `ti()` is refused. `smZ()`'s own tidy-up is left:
   the formula form covers what it lacked (numeric `by`, crossed factors via
   `by = interaction(f1, f2)`).
4. **Crossed-factor syntax**.
   **Step 4 done.** `pm_sm_parse()` expands a factor term among the arguments
   of `s()`: a common piece `s(x)` and, per term of the factor formula, a
   piece `s(x, by = <term>)` whose level curves share one variance (merged
   into one block, named `s_x_f1`, `s_x_f1_f2`, ...); cells of several
   factors get a derived `by` factor, rebuilt for new data. Each piece is set
   up in its own `gam(fit = FALSE)` call, so mgcv puts no side constraints
   between them. Fixed terms are written with the original factors
   (`x + f1*f2 + f1:x + f2:x + f1:f2:x`) and checked to span each piece's
   unpenalised columns. Checks: `s(x, f)` equals the HRW-style `lme` with one
   `pdIdent` block for the common curve and one for the deviations (same
   log-likelihood); `s(x, f1:f2)` equals `s(x, cell)`; with data whose curve
   shape depends on `f1` only, the `f1:f2` deviation variance is 0 and the
   `f2` one small.
5. **Prediction and inference**:
   - prediction from the stored smooth objects;
   - the `unconditional` correction, with the coverage simulation;
   - KR/Satterthwaite df without lmerTest, then lmerTest out of Imports.
   **Step 5 done (part).** Smoothing-parameter correction (`pm_sr_uncond()`):
   first-order J V_theta J' (lmer: J by numerical derivative of the
   profiled estimates from the fit's deviance function, V_theta from
   `pm_lmer_varpar()`; glmer: the stage-1 deviance (theta only, beta and u by
   PIRLS) for both; glmmTMB: TMB's joint precision). Within 1% of mgcv's
   `unconditional = TRUE`; glmmTMB and lme4 agree. Coverage simulation
   (`dev/sim-curve-difference.R`): Poisson 0.933 -> 0.957 pointwise, 0.907 ->
   0.940 whole curve; binomial 0.936 -> 0.955, 0.893 -> 0.947; gaussian 0.968
   -> 0.976, 0.973 -> 0.980. Default: on for generalized fits only. The
   package's own Kenward-Roger df on a semireg fit equal those of the same
   model fitted by `lme` (281.495 both), so lmerTest is not needed for the
   tests; whether to drop it (it gives `anova()`, `summary()` p-values and
   `ranova()` on `fit$semer`) is the maintainer's decision. The `smoothZ`
   form still predicts by re-evaluating its `smZ()` call; it drops `k` and
   uses the stored knots, so it does not depend on the caller's objects in
   the usual case.
6. **HRW reproduction** — done with the maintainer's golden lme examples (`dev/semireg-golden/`); results in section 4.1.
7. **Only-predictmeans cases**, as tests and as a vignette section.
8. **Documentation**: `vignettes/smooth.Rmd`, help pages, `inst/llms-context.md`.
9. R CMD check, then a pull request.
