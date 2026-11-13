# predictmeans: context for AI assistants and scripted use

Package version 2.0 (development). Compact reference to the result tables, defaults and
common recipes. Full documentation: the help pages and vignettes
(`browseVignettes("predictmeans")`) and https://dongwenluo.github.io/predictmeans/.

## What it does

Predicted (least-squares, "marginal") means of model terms, their standard errors,
confidence limits, pairwise comparisons, LSDs, contrasts, joint tests and slopes, for
`aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer`, `glmmTMB`, `gam` (mgcv) and
`semireg()` fits. Means average over the other factors of the model (equal weights by
default) with covariates at their means. Results match SAS `LSMEANS` and emmeans.

## Core call

```r
pm <- predictmeans(model, modelterm, atvar = NULL, adj = "none", level = NULL,
                   pairwise = FALSE, trans = I, transOff = 0, plot = FALSE, prtnum = FALSE)
```

- `modelterm`: a fixed-effect term written as in the model, e.g. `"A"` or `"A:B"`; its
  variables must be factors. A wrong order or case gives an error naming the right term.
- `atvar`: compare the levels of the other factor(s) within each level of `atvar`.
- `adj`: `"none"`, `"tukey"`, `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`,
  `"BY"`, `"fdr"`, `"sidak"`, `"scheffe"`, `"mvt"` (exact single-step, = SAS
  `ADJUST=SIMULATE`; needs mvtnorm). Any letter case. Giving `adj` or `atvar` turns on
  pairwise comparisons.
- `level`: significance level for confidence limits and LSDs, default 0.05.
- `plot = FALSE` in scripts (the default `TRUE` draws plots), or `options(predictmeans.plot = FALSE)`
  once; `prtnum = FALSE` suppresses printing of the reference-grid note.
- `generics::glance(pm)`: one row describing the result (term, numbers of means and pairs, df
  method, adjustment, SED range, mean LSD).
- `df_method`: `"auto"` (Kenward-Roger for REML `lmer`, `lme`, `gls`; Satterthwaite for ML fits;
  `"inner-outer"` for GLMMs: df of the grouping level at which each effect varies; `"residual"` = 1.x),
  `"Kenward-Roger"`, `"Satterthwaite"`, `"containment"` (SAS `DDFM=CONTAIN`), `"between-within"`
  (SAS `DDFM=BETWITHIN`; also for `lme` and `gls`), `"KR-expected"` (Kenward-Roger with the expected
  information, as pbkrtest/emmeans; `"Kenward-Roger"` uses the observed information, as SAS `DDFM=KR`),
  `"KR-average"` (the average information: Genstat's df).
  For `lme` and `gls`, `"Satterthwaite"`, `"Kenward-Roger"`, `"KR-expected"` and `"KR-average"` work with any
  correlation/variance structure (defaults: Kenward-Roger for REML fits, Satterthwaite for ML fits;
  `"containment"`: nlme's or the residual df, as 1.x).
- `vcov_adjust = TRUE` (default): Kenward-Roger adjusted covariance for REML `lmer`, as SAS
  `DDFM=KR` and emmeans; `FALSE` reproduces predictmeans 1.x standard errors.
- `weights`: `"equal"` (default), `"proportional"` (= SAS `OM`), `"outer"`, `"cells"`,
  `"flat"`, as in emmeans.
- `covariate`, `cov_reduce`, `exclude`: covariate values, how covariates are reduced
  (function or `"average"`), smooth terms to leave out (`gam`).

## Result: stable tables

`pm` has class `c("predictmeans", "pdmlist")`. Use these slots (fixed names, unrounded):

`pm$means`, one row per mean:

| column | meaning |
|---|---|
| factor (and covariate) columns | the level of each variable of `modelterm` |
| `estimate`, `std.error`, `df` | mean on the model (link) scale, SE, denominator df |
| `conf.low`, `conf.high` | confidence limits at `1 - level` |
| `response`, `response.std.error`, `response.conf.low`, `response.conf.high` | only when the means are back-transformed (non-identity link or `trans`); SE by the delta method |
| `group` | letters (only with pairwise comparisons); means sharing a letter do not differ |
| `n` | number of observations in the cell |
| `estimable` | FALSE when the mean cannot be estimated (then `estimate` is NA) |

`pm$pairs`, one row per pair (NULL when there is only one mean):

| column | meaning |
|---|---|
| `level1`, `level2` | the two means, labels `"a:b"` for several factors |
| `by` | the `atvar` level (only with `atvar`; NA for pairs across levels) |
| `over` | only when a mean averages over empty cells: the cells of the other factors the pair is compared on (those both means have; `"none"`: NA pair); NA for the other pairs |
| `estimate` | mean of `level1` minus mean of `level2` (link scale; with `over`: on the common cells) |
| `std.error`, `df`, `statistic`, `p.value` | t test of the difference |
| `p.adj` | p-value adjusted by `adj` (within each `by` level with `atvar`) |
| `p.equiv` | with `equivalence = d` or `c(lower, upper)`: TOST p-value (difference within the margins) |
| `lsd` | least significant difference at `level` |
| `p.method` | `"t"`, `"tukey"`, `"permutation"` or `"bootstrap"` |
| `ratio`, `ratio.conf.low`, `ratio.conf.high` | log link or `trans = exp`: ratio of the means |
| `odds.ratio`, `odds.ratio.conf.low`, `odds.ratio.conf.high` | logit link: odds ratio |

`pm$settings`: list with `modelterm`, `level`, `adjust`, `atvar`, `df_method`,
`vcov_adjust`, `vcov` (`"model"`, or the covariance given by `vcov =`: `"HC0"`-`"HC3"`, `"CR1, clusters: ..."`,
`"user matrix"`, `"user function"`), `equivalence` (the TOST margins, or NULL), `bias_adjust` (the SD
used with `bias_adjust = TRUE`, else NA), `average` (`"link"` or `"response"`), `singular` (TRUE when
the lmer/glmer fit is singular; NA for other classes), `nperm`, `covariates`,
`model_class`, `version`.

Also: `as.data.frame(pm)` / `as.data.frame(pm, what = "pairs")`, `generics::tidy(pm)`,
`summary(pm)`, `ggplot2::autoplot(pm, type = c("means", "pvalue", "response"))`,
`coef(pm)` and `vcov(pm)` (means and their covariance, named `"a:b"`).

The 1.x slots (`"Predicted Means"`, `"Standard Error of Means"`, `"LSD"`, `mean_table`,
`"Pairwise p-value"`, ...) are still there, with names that vary with the options;
prefer the stable tables.

## Recipes

```r
library(predictmeans)
Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)

# means of one factor, all pairs, Tukey; letters in pm$means$group
pm <- predictmeans(fm, "nitro", adj = "tukey", plot = FALSE, prtnum = FALSE)

# nitrogen compared within each variety
pm2 <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)

# contrasts by level names; standard families (Dunnett, polynomial, ...)
K <- Kmatrix(fm, "nitro:Variety")
cm <- specify_contrast(K, N_Victory = `0.6:Victory` - `0:Victory`)
contrastmeans(fm, "nitro:Variety", ctrmatrix = cm)[[1]]
contrastmeans(fm, "nitro", ctrmatrix = contrast_family(Kmatrix(fm, "nitro"), "trt.vs.ctrl"),
              adj = "mvt")[[1]]

# joint F tests: nitrogen within each variety (SAS SLICE)
pm_joint(fm, "nitro:Variety", by = "Variety")

# slopes of a numeric covariate by group (fit with nitro numeric)
Oats2 <- as.data.frame(nlme::Oats)
fm2 <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats2)
pm_slopes(fm2, "nitro", "Variety")

# GLM / GLMM: response-scale means, ratios or odds ratios in pm$pairs
data(Drug)
g <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
predictmeans(g, "Drug", pairwise = TRUE, plot = FALSE)$pairs

# transformed response
# predictmeans(lmer(log(y) ~ ...), "A", trans = exp, plot = FALSE)   # ratios of geometric means

# permutation p-values for terms and pairs (Freedman-Lane by default, any model class:
# lm, aov, glm, lme, gls, lmer, glmer, glmmTMB; mixed models refitted, refit = FALSE is faster)
perm <- permmodels(fm, nperm = 999, seed = 1, prt = FALSE)
predictmeans(fm, "nitro", permlist = perm, plot = FALSE, ncore = 1)$pairs
# unequal error variances in lm/glm: HC3 statistics with sign flips (wild bootstrap)
# permmodels(lm_fit, nperm = 999, seed = 1, robust = TRUE)
# effect sizes and a check of the parametric test; several responses combined
permstats(perm)
# permcombine(y1 = perm1, y2 = perm2)   # same rows, same seed and nperm
# the 1.x scheme (response / fitted values permuted, refits kept in permlist)
perm_raw <- permmodels(fm, nperm = 999, ncore = 1, seed = 1, prt = FALSE, method = "raw")
# random effects or any two nested mixed models (lmer, glmer, glmmTMB)
# permlmer(fm0, fm1, nperm = 999, ncore = 1, seed = 1)    # permglmer(), permglmmTMB() too
# randomisation test: permutations within blocks
idx <- permindex(Oats, block = "Block", nsim = 999, seed = 1)
permmodels(lm(yield ~ Block + nitro, data = Oats), perms = idx, prt = FALSE)$ANOVA

# R2 as numbers; variance components
R2_model(fm); R2_model(fm, method = "nakagawa"); varcomp(fm)
```

## Things to know

- Kenward-Roger standard errors are the default for REML `lmer` fits, so SEs differ
  slightly from predictmeans 1.x; use `vcov_adjust = FALSE` for the 1.x values.
- Tests and letters are on the link scale for GLMs/GLMMs; back-transformed GLMM means are
  those of a unit with all random effects zero (not population averages); ratios do not
  depend on this.
- Empty cells: cell means of unobserved combinations are not shown; a mean averaging over
  an empty cell follows the 1.x rule (average of the observed cells, with a warning naming
  it), and its pairs are compared on the cells both means have (column `over`), not as
  the difference of the two averages (since 2.0). Means that cannot be estimated (confounded
  factors) are NA, with their estimable differences kept. Check `estimable` and `n`.
- `letterCI = TRUE`: letters from overlap of the confidence intervals; `level` then sets
  the CIs (default 0.166, i.e. 83.4% CIs, about a 0.05-level test).
- `permmodels()`, `permlmer()`, `permanova.lmer()` and `predictmeans(permlist = )` use
  `ncore = 1` by default since 2.0 (3 before); set `ncore` for parallel refits. With `permlist`,
  `adj = "tukey"`/`"mvt"` is the single-step and `adj = "westfall"` the step-down
  permutation max-|t| (Westfall-Young); `"scheffe"` is not available.
- `permmodels()` defaults to Freedman-Lane since 2.0 (`method = "auto"`); `method =
  "raw"` reproduces 1.x. For GLMs/GLMMs it works on the link (working) scale.
  `permlmer()` for GLMMs uses the working scale, and randomised quantile residuals for
  0/1 responses. With a `permlist`, `vcov_adjust` (Kenward-Roger) is not used.
- `vcov = "HC3"` (or `"HC0"`-`"HC2"`, `~ cluster`) gives robust SEs for lm/glm in `predictmeans()`,
  `contrastmeans()`, `pm_slopes()` and `pm_joint()`; a matrix or a function works for any model (e.g.
  clubSandwich for mixed models). The df stay the model's.
- Back-transformed means (`response` columns) are those of a unit with random effects zero, averaged on
  the link scale. `bias_adjust = TRUE` (with `sigma`) gives population means (emmeans' `bias.adjust`);
  `average = "response"` averages the back-transformed cells (emmeans' `regrid()`).
- `permstats()` gives standardised effect sizes and `param.size`, the rejection rate of the
  parametric test in the permuted data sets (near `level` when the parametric test suits
  the data). `permcombine()` combines tests of several responses made with the same
  permutations (Fisher or Tippett) and gives min-p adjusted p-values over the responses.
- `semipred()` limits are two-sided at `level` (95% by default) since 2.0; 1.x gave 90%.
  `simultaneous = TRUE` gives simultaneous bands. `predictmeans()` on a `semireg()` fit
  gives means with the smooth included and normal-based tests.
- `marginaleffects::hypotheses()` works on `predictmeans()` results after
  `options(marginaleffects_model_classes = "predictmeans")`.
- Data must be available when `predictmeans()` runs: it re-reads the model's `data`
  argument, so fit models with a named data frame (not an expression).
