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

## 6. Steps

1. This note.
2. **One engine** (`R/semireg-engine.R`) — **done**. Result check: 12 example fits (gaussian, ML, binomial, Poisson, `bs` basis, by-factor and group smooths, Dialyzer, indonRespir, glmmTMB Poisson and negative binomial) give the same coefficients, covariance matrices and `semipred()` output as before (log-likelihoods within 1e-7); total time 19 s -> 14 s. Found and fixed on the way: lmerTest's Satterthwaite df and `varcomp()` used the deviance of the placeholder model (now equal to `lme`); a zero smoothing variance stopped the fit; `cbind()` responses; data written into the caller's frame. Was: `semireg()` and `semireg_tmb()` run on
   it with the existing `smoothZ` form only, and must reproduce current
   results: the `semi_*` reference cases, `test-semipred.R` and
   `test-glmmtmb.R` pass unchanged. One fit instead of two (timed).
3. **Formula smooths**: `s()`, `by`, `fs`, `t2`/`te`, `os`, `ztps`, with
   automatic fixed terms; mgcv to Imports. Checked against
   `mgcv::gam(method = "REML")`.
4. **Crossed-factor syntax**.
5. **Prediction and inference**:
   - prediction from the stored smooth objects;
   - the `unconditional` correction, with the coverage simulation;
   - KR/Satterthwaite df without lmerTest, then lmerTest out of Imports.
6. **HRW reproduction** (`dev/hrw/`), with a results table added to this note.
7. **Only-predictmeans cases**, as tests and as a vignette section.
8. **Documentation**: `vignettes/smooth.Rmd`, help pages, `inst/llms-context.md`.
9. R CMD check, then a pull request.
