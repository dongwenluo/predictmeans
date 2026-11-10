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
- `plot = FALSE` in scripts (the default `TRUE` draws plots); `prtnum = FALSE` suppresses
  printing of the reference-grid note.
- `df_method`: `"auto"` (Kenward-Roger for REML `lmer`, Satterthwaite for ML),
  `"Kenward-Roger"`, `"Satterthwaite"`, `"containment"`.
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
| `estimate` | mean of `level1` minus mean of `level2` (link scale) |
| `std.error`, `df`, `statistic`, `p.value` | t test of the difference |
| `p.adj` | p-value adjusted by `adj` (within each `by` level with `atvar`) |
| `lsd` | least significant difference at `level` |
| `p.method` | `"t"`, `"tukey"`, `"permutation"` or `"bootstrap"` |
| `ratio`, `ratio.conf.low`, `ratio.conf.high` | log link or `trans = exp`: ratio of the means |
| `odds.ratio`, `odds.ratio.conf.low`, `odds.ratio.conf.high` | logit link: odds ratio |

`pm$settings`: list with `modelterm`, `level`, `adjust`, `atvar`, `df_method`,
`vcov_adjust`, `nperm`, `covariates`, `model_class`, `version`.

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

# permutation / bootstrap p-values for pairs
perm <- permmodels(fm, nperm = 999, ncore = 1, seed = 1, prt = FALSE)
predictmeans(fm, "nitro", permlist = perm, plot = FALSE, ncore = 1)$pairs

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
  an empty cell follows the 1.x rule (average of the observed cells, with a warning) and
  is not comparable with complete means. Means that cannot be estimated (confounded
  factors) are NA, with their estimable differences kept. Check `estimable` and `n`.
- `letterCI = TRUE`: letters from overlap of the confidence intervals; `level` then sets
  the CIs (default 0.166, i.e. 83.4% CIs, about a 0.05-level test).
- `permmodels()`, `permlmer()` and `predictmeans(permlist = )` use `ncore = 3` by
  default; use `ncore = 1` or 2 in package checks and on shared machines. `adj` options
  `"tukey"`, `"scheffe"` and `"mvt"` are not available with `permlist`.
- `semipred()` limits are two-sided at `level` (95% by default) since 2.0; 1.x gave 90%.
  `simultaneous = TRUE` gives simultaneous bands. `predictmeans()` on a `semireg()` fit
  gives means with the smooth included and normal-based tests.
- `marginaleffects::hypotheses()` works on `predictmeans()` results after
  `options(marginaleffects_model_classes = "predictmeans")`.
- Data must be available when `predictmeans()` runs: it re-reads the model's `data`
  argument, so fit models with a named data frame (not an expression).
