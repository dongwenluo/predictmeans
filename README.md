# predictmeans

Predicted (least-squares) means, standard errors, LSDs, pairwise comparisons and
contrasts for linear, generalized linear and mixed models — `aov`, `lm`, `glm`, `gls`,
`lme`, `lmer`, `glmer`, `glmmTMB`, `gam` and semiparametric `semireg` models — with
Kenward–Roger, Satterthwaite or containment degrees of freedom, permutation and bootstrap
tests, R² for mixed models and diagnostics. Results match SAS `LSMEANS` and emmeans.

## Installation

```r
install.packages("predictmeans")                       # CRAN
# remotes::install_github("dongwenluo/predictmeans")    # development version
```

## Example

```r
library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)

pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
pm$means     # one row per mean: estimate, std.error, df, conf.low, conf.high, group, n, ...
pm$pairs     # one row per comparison: estimate, std.error, df, statistic, p.value, p.adj, ...
ggplot2::autoplot(pm)
```

Every `predictmeans()` result has three stable data frames, `means`, `pairs` and
`settings`, with fixed column names; `as.data.frame()`, `tidy()`, `coef()` and `vcov()`
work on it.

## Main functions

| Task | Function |
|---|---|
| Means and pairwise comparisons | `predictmeans()`, `covariatemeans()` |
| Contrasts | `specify_contrast()`, `contrast_family()`, `contrastmeans()` |
| Joint tests (slices) and slopes | `pm_joint()`, `pm_slopes()` |
| Resampling tests | `permmodels()`, `permlmer()`, `permanova.lmer()`, `bootmodels()` |
| Semiparametric models | `semireg()`, `semireg_tmb()`, `semipred()` |
| R², variance components, diagnostics | `R2_model()`, `R2_glmm()`, `varcomp()`, `pm_check()`, `residplot()`, `CookD()` |
| Reference grid / K matrix | `Kmatrix()` |

## Documentation

The [package website](https://dongwenluo.github.io/predictmeans/) has the help pages
and vignettes: getting started; contrasts, joint tests and slopes; link and response
scale; permutation and bootstrap tests; smooth curves in mixed models; missing cells and
non-estimable means; R² for mixed models; and migrating from version 1.x.

For AI assistants and scripted use, a compact description of the result tables, the
defaults and common recipes is installed with the package:
`system.file("llms-context.md", package = "predictmeans")`
([on GitHub](https://github.com/dongwenluo/predictmeans/blob/main/inst/llms-context.md)).
