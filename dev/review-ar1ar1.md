# Separable AR1 x AR1 residuals for field trials (review note, 2026-11-26)

For the maintainer's decisions before any code. **Decided 2026-11-26:** S1-S6 agreed; design D1-D7 below, checked by two outside AIs.

## Why

The standard spatial model for a row x column field trial is a separable
first-order autoregression in both directions (Gilmour, Cullis and Verbyla
1997): the residuals e = sigma^2 (Sigma_r (x) Sigma_c), with
corr(plot i, plot j) = phi_r^|dr| phi_c^|dc|. Genstat (`VSTRUCTURE ... AR,AR`)
and ASReml (`ar1(row):ar1(col)`) fit it; in R only ASReml-R (commercial), sommer
(the two phi fixed on a grid) and INLA do. nlme has `corAR1` (one direction) and
isotropic `corExp`/`corGaus`; glmmTMB 1.1.15 has `ar1()` in one direction and
isotropic `exp()`, `gau()`, `mat()`, but no product of two AR1s with their own phi.
So `residvariogram()` can show the need for the model, and its advice can only
point to isotropic stand-ins or a smooth surface.

## Proposal

A correlation structure for nlme written in R (no compiled code):
`corAR1xAR1(value = c(0.5, 0.5), form = ~ row + col | trial, nugget = FALSE)`,
used as `gls(yield ~ variety, correlation = corAR1xAR1(form = ~ row + col),
data = d)` or in `lme()` with random replicates or blocks. nlme lets a package
define its own `corStruct` class (methods `Initialize`, `corMatrix`, `coef`,
`coef<-`; nlme's `corFactor`/`recalc` then work from `corMatrix`), so every
nlme tool works with it, and so does everything predictmeans does for gls/lme:
Kenward-Roger and Satterthwaite df (`R/kr.R`, `R/satt-nlme.R`: they differentiate
the covariance numerically, whatever the structure), means, SEDs, `varcomp()`,
`residplot()` (normalized residuals) and `residvariogram()` (a flat variogram
when the model fits).

- Positions: integer row and column over the whole field (the check of
  `residvariogram()`); the correlation of two plots is phi_r^|dr| phi_c^|dc|, so
  missing plots need no special handling (the correlation matrix of the plots
  present is the sub-matrix).
- Parameters: phi_r, phi_c in (-1, 1) on the scale atanh(phi) for the
  optimiser (negative phi allowed, as Genstat).
- `nugget = TRUE`: an independent plot error besides the AR1 x AR1 part (Genstat's
  "units" term, the measurement error), as nlme's nugget of `corExp`.
- Several trials: `| trial`, the same phi in each (separate phi per trial: fit
  each trial alone, or later).

## Decisions

**S1. Name and place.** Proposal: exported `corAR1xAR1()` in a new
`R/cor-ar1xar1.R`, documented with the Slate Hall example; its methods internal.

**S2. Scope of the first step.** Proposal: gls and lme, with and without the
nugget, one or several trials; AR2 and other orders later (Genstat's ORDER = 2),
and glmmTMB not at all (no way to add a structure there).

**S3. Size.** The correlation matrix of a trial is dense (n x n per trial):
fine to about 2000 plots per trial (a few seconds per likelihood); larger trials
get a message pointing to the smooth surface `t2(row, col)` in `semireg()`.

**S4. References.** Exact values in Genstat's REML Guide, chapter 3, for the
Slate Hall trial (the data are in `tests/testthat/_data/variogram/slatehall.csv`):
- variety fixed, AR1 x AR1: sigma^2 3.876, phi 0.4586 (rows) and 0.6838 (columns),
  deviance 249.35 on 122 d.f.; variety Wald 313.04 on 24, F 13.04 on 80.0 d.d.f.;
- with a nugget (units 0.500), replicates: AR part 0.500 x 9.240 = 4.620, phi 0.6922 and 0.8442,
  deviance 242.20, average SED of variety means 0.6072, Wald 241.61, d.d.f. 74.7;
- the variogram of the fitted model's residuals (`residvariogram()`, Genstat's
  `FVARIOGRAM` of the AR1 x AR1 residuals already in `test-variogram.R`).
Proposal: a Genstat run by the maintainer printing the variety means and SEDs
(and the d.d.f. per comparison) of both models to 6 digits, for tests to 1e-4;
the deviance is compared as a difference between the two models (Genstat's
omits constants). ASReml-R, if the maintainer has it, as a second reference.

**S5. Use elsewhere.** Proposal: `residvariogram()` advice names
`corAR1xAR1()` for row x column layouts; `pm_check()` likewise; the
tree-improvement / field-trial article gets a Slate Hall section (variogram of the
RCB fit, the AR1 x AR1 fit, its flat variogram, variety means).

**S6. Checks of its own.** nlme's `corAR1` is the special case of one row (or
one column): equal to 1e-8; a simulated field with known phi recovered; KR df
and SEDs against Genstat (S4); the variogram of normalized residuals flat in
simulations with the right model (the envelope of `residvariogram()` at its
nominal size).

**Documentation (maintainer, 2026-11-26).** Wherever `corAR1xAR1()` is documented (its help page, the
Slate Hall article section, the `residvariogram()` advice), also show the smooth field trend, with its comment:

| model | call | comment |
|---|---|---|
| Smooth field trend | `semireg(y ~ variety + t2(row, col) + (1 \| rowf) + (1 \| colf), data = d)` | a smooth surface over the field, as SpATS fits; often as good as AR1 x AR1 in practice. |

(`row`, `col` numeric positions; `rowf`, `colf` the same as factors.) Slate Hall, first look
(2026-11-26): it fits in under a second; the average SED of the 25 variety means is 0.668 against
Genstat's 0.590 (AR1 x AR1) and 0.605 (with the nugget). So on this trial the AR1 x AR1 fit is the more
precise; the article shows both with their variograms and states this, rather than "as good" alone.

## Genstat run (2026-11-26)

Output: `dev/genstat/genstat-results-2026-11-26-ar1ar1-slatehall.txt` (script
`ar1ar1-slatehall.gen`); 25 means and 300 SEDs of each model to 8 decimals.

| model | variance parameters | phi row, col | deviance (d.f.) | variety Wald, F, d.d.f. |
|---|---|---|---|---|
| A AR1 x AR1 | sigma^2 3.876 | 0.4586, 0.6838 | 249.35 (122) | 313.04, 13.04, 80.0 |
| B + nugget | units 0.486, AR part 0.486 x 9.420 = 4.578 | 0.6827, 0.8438 | 242.35 (121) | 245.39, 10.21, 75.7 |
| C replicates + A | replicates **-0.024**, sigma^2 3.880 | 0.4602, 0.6827 | 249.34 (121) | 312.41, 13.01, 79.3 |
| D replicates + B | replicates **-0.0777**, units 0.500, AR 0.500 x 9.240 | 0.6922, 0.8442 | 242.20 (120) | 241.61, 10.05, 74.7 |

- A and B are the references for `gls`. Checked here at Genstat's rounded
  estimates: GLS means to 2e-5, sigma^2 and the REML deviance (without the
  2 pi constant) to the printed digits, so the data file and the
  parametrisation match.
- With a nugget Genstat makes `units` the residual and scales the AR1 x AR1 term
  by a ratio (its warning "More than one residual term" is about this); the
  same model as nlme's nugget form, nugget share 1 / (1 + 9.420) = 0.096.
- C and D have **negative** replicate components (Genstat leaves them
  unconstrained). `lme()` keeps variances positive, so it reaches the boundary
  (replicates about 0), where C is A and D is B. C and D are therefore not
  references for `lme()`; the lme tests check C = gls A and D = gls B at the
  boundary, and the lme code path against A/B values with replicates fixed
  near 0. (S4 first said 4.662 for D's AR part; the run gives 4.620.)

## Design for the code (2026-11-26, before any code)

S1-S6 stand. This section fixes how they are built, from nlme's own machinery
(read from nlme 3.1: `corFactor.corStruct()`, `recalc.corStruct()`,
`logDet.corStruct()` and `Initialize.corStruct()` all work from `corMatrix()`).

**D1. The class.** `corAR1xAR1(value = c(0.5, 0.5), form = ~ row + col | group,
nugget = FALSE, fixed = FALSE)` returns an object of class
`c("corAR1xAR1", "corStruct")`: the unconstrained parameters
`atanh(phi_r), atanh(phi_c)` (and `qlogis(nugget)` when `nugget = TRUE`), with
the formula in an attribute, as nlme's own classes. Methods (internal,
`@exportS3Method nlme::...`):
- `Initialize.corAR1xAR1()`: `Initialize.corStruct()` (groups, `Dim`,
  covariate = the row and column of each plot, integers checked), then the
  factor and log-determinant, as `Initialize.corSpatial()`;
- `corMatrix.corAR1xAR1(object, covariate, corr = TRUE)` handles both
  branches itself (nlme's `corMatrix.corStruct(corr = FALSE)` counts
  `length()` of the covariate, wrong for a row-column matrix); `corr = FALSE`
  returns the transposed inverse Cholesky factor of each group's matrix with
  `attr(, "logDet")` = -sum(log(diag(chol(R)))) (summed over groups), as
  nlme's classes; per group
  `phi_r^|dr| * phi_c^|dc|` (outer products of the two AR1 matrices indexed by
  the plots present: missing plots need nothing), times `1 - nugget` off the
  diagonal; `corr = FALSE` gives the transpose of the inverse Cholesky factor
  with attribute `logDet`, as nlme's classes;
- `coef.corAR1xAR1(object, unconstrained = TRUE)` (FALSE: phi and the nugget
  share; `numeric(0)` when `fixed = TRUE`, so the optimiser leaves them), `coef<-` (recompute the factor and log-determinant: a stale cached
  factor would be used by the next likelihood), `print`/`summary` via the corStruct methods, `Dim`.
Starting values: 0.5, 0.5 (and nugget 0.1); `fixed = TRUE` keeps them.

**D2. gls.** `gls(yield ~ variety, correlation = corAR1xAR1(form = ~ row + col),
data = d)` is model A; with `nugget = TRUE` model B. nlme's nugget is the share
of the independent plot variance in the total: n = units / (units + spatial).
In the run Genstat made `units` the residual and gave the AR1 x AR1 term the
ratio g = spatial / units = 9.420, so n = 1 / (1 + g) = 0.0960 (with AR1 x AR1
as the residual and units as a random term, Genstat's ratio is the inverse and
n = g / (1 + g)); the tests compare the variances themselves (units 0.486,
spatial 0.486 x 9.420). References, in two tiers (Genstat stops on a deviance
change below 1e-4, which does not make phi exact to 1e-4): (i) phi and the
variances fixed at Genstat's printed estimates (`fixed = TRUE`): means and SEDs
to the rounding of those estimates (2e-5 when checked by hand here); (ii) free
fits: REML deviance A - B = 7.00 and the deviance at our estimates no larger
than at Genstat's, phi to about 1e-3 for model A and about 1e-2 for model B
(with a nugget the likelihood surface in phi and the nugget can be a flat ridge).

**D3. lme: the correlation must span the field.** In `lme()` a correlation
structure acts within the groups of the innermost random effect: with
`random = ~ 1 | rep` the AR1 x AR1 would be within each replicate, with no
correlation across replicate boundaries. That is a different model from
Genstat's C and D (AR1 x AR1 over the whole field plus replicate effects).
To get Genstat's model in nlme, the replicate effects go in a single group:
`lme(yield ~ variety, random = list(all = pdIdent(~ rep - 1)),
correlation = corAR1xAR1(form = ~ row + col | all), data = d)` with
`d$all <- factor(1)`. `corAR1xAR1()` cannot see the random part, so the help
page shows both forms and says which model each is; `pm_check()`/advice could
flag `| rep` nested within a random factor that is a field block (later).
References for lme: Genstat C and D have negative replicate components, so
lme (variances >= 0) reaches the boundary; tested: lme C = gls A and lme D =
gls B at the boundary (to the optimiser's tolerance). Genstat's REML Guide
confirms it: with the replicate variance constrained to 0 it gives phi 0.6827
and 0.8438, the values of model B.

**D4. predictmeans' df and SEs.** `R/kr.R` (`pm_vcov_kr_nlme()`) and
`R/satt-nlme.R` differentiate the covariance of gls/lme fits numerically in
the variance parameters, through `corMatrix()`; they must be checked with the
new class: KR df of a variety difference against Genstat's d.d.f. of the
variety F test (80.0 for A, 75.7 for B: the per-pair df are not printed, so the
check is the Wald F d.d.f. from `anova()` with KR, to 1%), and SEDs to 1e-4.

**D5. One trial or several.** `| trial` gives a block-diagonal correlation
(one AR1 x AR1 per trial, the same phi); separate phi per trial is later.

**D6. Size.** Dense n x n per group; above 2000 plots per group a message
points to `semireg(y ~ variety + t2(row, col) + (1 | rowf) + (1 | colf))`.

**D7. Documentation (maintainer).** The help page and the Slate Hall article
section show AR1 x AR1 next to the smooth field trend (`t2(row, col)`, as
SpATS), with their variograms and average SEDs (Slate Hall: 0.590 AR1 x AR1,
0.605 with nugget, 0.668 t2 trend).

Tests (`test-cor-ar1xar1.R`): corAR1 (one row) equal to 1e-8 in logLik and
estimates; Genstat A and B (means, SEDs, phi, sigma^2, deviance difference);
lme single-group form = Genstat model at a positive replicate variance (when
fixed) and = gls at the boundary; missing plots (drop 10 plots: no error,
the correlation of the rest is the sub-matrix); simulated field with known phi
recovered (skip_on_cran); `residvariogram()` of the AR1 x AR1 fit within its
envelope.

## Outside check (2026-11-26, answers tested)

1. lme correlation within the innermost random group: holds (also finer groups); D3 stands.
2. Kronecker product, sub-matrix for missing plots: holds; built from the original coordinates.
3. nugget n = 1/(1+g): holds for this run (units the residual); n = g/(1+g) when AR1 x AR1 is the residual; D2 says which.
4. four methods: holds; corMatrix(corr = FALSE) handled by the class, coef<- refreshes the factor (D1).
5. phi to 3-4 decimals: not a rule (Genstat's tolerance is on the deviance); tests in two tiers (D2).
   (Two AIs asked, OpenAI and Gemini: they agree on 1-4; on 5 Gemini says "usually 3-5 decimals"
   without a nugget and 1-2 with one, OpenAI that it is no established rule. Gemini added the exact
   corr = FALSE contract and coef() = numeric(0) when fixed: both in D1.)

## Questions asked

1. Check: in `nlme::lme()`, does a `correlation =` structure act only within
   the groups of the innermost random effect, so that `random = ~ 1 | rep`
   with an AR1 x AR1 over row and column gives no correlation between plots in
   different replicates?
2. Check: is a separable AR1 x AR1 correlation on a complete rows x columns
   grid the Kronecker product of the two AR1 matrices, and for plots missing
   from the grid simply the sub-matrix of the plots present?
3. Check: is nlme's `nugget` (correlation (1 - n) rho off the diagonal) the same
   model as Genstat's "units + AR1 x AR1" with ratio g = (AR variance) / (units
   variance), with n = 1 / (1 + g)?
4. Check: can a user-defined nlme `corStruct` that supplies only
   `Initialize`, `corMatrix`, `coef` and `coef<-` be used in `gls()` and
   `lme()` with REML, `anova()`, `intervals()` and `predict()`?
5. Check: for Slate Hall-sized trials (150 plots), do REML estimates of
   phi_row and phi_col from nlme and from Genstat's AI algorithm usually agree
   to 3-4 decimals?

## Effort

About a day: the class (small), tests against corAR1 and Genstat, the advice
and help; the article section with the B3 articles.

## References

Gilmour, A. R., Cullis, B. R. and Verbyla, A. P. (1997). Accounting for natural
and extraneous variation in the analysis of field experiments. *Journal of
Agricultural, Biological, and Environmental Statistics*, 2, 269-293.

Pinheiro, J. C. and Bates, D. M. (2000). *Mixed-Effects Models in S and S-PLUS*.
Springer (section 5.3 and 6.3: user-defined correlation structures).
