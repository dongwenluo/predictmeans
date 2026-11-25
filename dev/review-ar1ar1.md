# Separable AR1 x AR1 residuals for field trials (review note, 2026-11-26)

For the maintainer's decisions before any code.

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

## Effort

About a day: the class (small), tests against corAR1 and Genstat, the advice
and help; the article section with the B3 articles.

## References

Gilmour, A. R., Cullis, B. R. and Verbyla, A. P. (1997). Accounting for natural
and extraneous variation in the analysis of field experiments. *Journal of
Agricultural, Biological, and Environmental Statistics*, 2, 269-293.

Pinheiro, J. C. and Bates, D. M. (2000). *Mixed-Effects Models in S and S-PLUS*.
Springer (section 5.3 and 6.3: user-defined correlation structures).
