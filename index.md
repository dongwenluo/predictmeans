# predictmeans

Predicted (least-squares) means, standard errors, LSDs, pairwise
comparisons and contrasts for linear, generalized linear and mixed
models — `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer`, `glmmTMB`,
`gam` and semiparametric `semireg` models — with Kenward–Roger,
Satterthwaite or containment degrees of freedom, permutation and
bootstrap tests, R² for mixed models and diagnostics. Results match SAS
`LSMEANS` and emmeans.

## Installation

``` r

install.packages("predictmeans")                       # CRAN
# remotes::install_github("dongwenluo/predictmeans")    # development version
```

## Example

``` r

library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)

pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
pm$means     # one row per mean: estimate, std.error, df, conf.low, conf.high, group, n, ...
pm$pairs     # one row per comparison: estimate, std.error, df, statistic, p.value, p.adj, ...
ggplot2::autoplot(pm)
```

Every
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
result has three stable data frames, `means`, `pairs` and `settings`,
with fixed column names;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
`tidy()`, [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) work on it.

## Main functions

| Task | Function |
|----|----|
| Means and pairwise comparisons | [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) |
| Contrasts | [`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md), [`contrast_family()`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md), [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md) |
| Joint tests (slices) and slopes | [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md), [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md) |
| Resampling tests | [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md), [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md), [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/permanova.lmer.md), [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md) |
| Semiparametric models | [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md), [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md), [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md) |
| R², variance components, diagnostics | [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md), [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md), [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md), [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md), [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md) |
| Reference grid / K matrix | [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md) |

## Documentation

The [package website](https://dongwenluo.github.io/predictmeans/) has
the help pages and vignettes: getting started; contrasts, joint tests
and slopes; link and response scale; permutation and bootstrap tests;
smooth curves in mixed models; missing cells and non-estimable means; R²
for mixed models; and migrating from version 1.x.

For AI assistants and scripted use, a compact description of the result
tables, the defaults and common recipes is installed with the package:
`system.file("llms-context.md", package = "predictmeans")` ([on
GitHub](https://github.com/dongwenluo/predictmeans/blob/main/inst/llms-context.md)).
