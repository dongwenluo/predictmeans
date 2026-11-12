# Migrating from predictmeans 1.x

Most scripts written for predictmeans 1.x run unchanged. This page lists
what gives different numbers, what is new, and how to get the old
results back.

``` r

library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
```

## Results that change

### Kenward–Roger adjusted standard errors are the default

For `lmer` models fitted by REML, standard errors, confidence limits,
LSDs and p-values now use the Kenward–Roger adjusted covariance, as SAS
`DDFM=KR` and emmeans do. In 1.x the unadjusted covariance was used with
Kenward–Roger df. For balanced designs the two agree; for unbalanced
ones the adjusted SEs are slightly larger.

To reproduce 1.x numbers, use `vcov_adjust = FALSE`:

``` r

Om <- Oats[-c(1, 2, 30), ]   # unbalanced
fm_u <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Om)
new <- predictmeans(fm_u, "Variety", plot = FALSE, prtnum = FALSE)
old <- predictmeans(fm_u, "Variety", vcov_adjust = FALSE, plot = FALSE, prtnum = FALSE)
cbind(SE_2.0 = new$means$std.error, SE_1.x = old$means$std.error)
#>        SE_2.0   SE_1.x
#> [1,] 8.086277 8.085432
#> [2,] 8.057580 8.057580
#> [3,] 8.137792 8.132226
```

### Means that cannot be estimated are NA

A mean that depends on aliased (non-estimable) coefficients used to be
printed as a number that depended on which coefficient R happened to
drop. It is now `NA`, with a message; pairwise differences that are
estimable are still shown. Means over empty cells keep the 1.x rule: the
average of the observed cells, with a warning.

### Permutation tests

[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
now uses Freedman–Lane permutations by default (`method = "auto"`): each
term, coefficient or pair is tested with the model without it, so the
other terms keep their effects. It also works for `glm`, `glmer` and
`glmmTMB` models, and needs no refits. Mixed-model residuals are
decorrelated with the symmetric square root of their covariance instead
of the Cholesky factor, whose results depended on the order of the data.
For the 1.x scheme use `method = "raw"`; its p-values for mixed models
still differ slightly from 1.x because of the new square root.

### Bugs fixed

Some 1.x results were wrong and are now correct, for example `permlist`
together with `atvar`, a user-supplied `Df` for lmer models, factor
levels containing `":"`, and several
[`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
cases with transformed covariates. See `NEWS.md` for the full list.

## New in 2.0

- **Stable tables** `means`, `pairs` and `settings` in every result,
  with [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  `tidy()`, [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).
- **Degrees of freedom**: `df_method = "Kenward-Roger"`,
  `"Satterthwaite"` or `"containment"` (SAS `DDFM=CONTAIN`).
- **Averaging weights**: `weights = "proportional"` (SAS `OM`),
  `"outer"`, `"cells"`, `"flat"`.
- **[`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md)**
  to write contrasts by the names of the means.
- **Covariates**:
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  handles [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`,
  [`log()`](https://rdrr.io/r/base/Log.html), smooth terms of `gam` and
  covariates that enter only through interactions; `cov_reduce` and
  `exclude` pass through
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  and
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md).
- **More comparisons**: `adj = "sidak"`, `"scheffe"` and `"mvt"` (exact
  single-step, as SAS `ADJUST=SIMULATE`);
  [`contrast_family()`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md)
  (Dunnett, polynomial trends, …);
  [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
  for joint F tests and slices;
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md)
  for slopes by group; [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) of the means. See
  [`vignette("contrasts")`](https://dongwenluo.github.io/predictmeans/articles/contrasts.md).
- **Response scale**: delta-method standard errors of back-transformed
  means, and ratios or odds ratios for log and logit models.
- **[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)**
  (R² as numbers, Piepho or Nakagawa),
  **[`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)**
  (parametric bootstrap), and
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  on
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fits.
- **Permutation tests for GLMs and GLMMs**:
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  on the working scale,
  [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  (also
  [`permglmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md),
  [`permglmmTMB()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md))
  for nested GLMMs, and `perms` for permutations restricted by the
  design
  ([`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)).
  See
  [`vignette("resampling")`](https://dongwenluo.github.io/predictmeans/articles/resampling.md).
- **Other model classes** can be added with a few methods (see
  [`?pm_coef`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)).
- **Fewer dependencies** and faster Kenward–Roger and Satterthwaite
  computations.

## Old and new output side by side

The 1.x list elements are still there, so code that reads
`pm$"Predicted Means"` or `pm$"Pairwise p-value"` keeps working:

``` r

pm <- predictmeans(fm, "nitro", adj = "BH", plot = FALSE, prtnum = FALSE)
names(pm)
#>  [1] "Predicted Means"               "Standard Error of Means"      
#>  [3] "Standard Error of Differences" "LSD"                          
#>  [5] "Pairwise p-value"              "predictmeansPlot"             
#>  [7] "predictmeansBarPlot"           "mean_table"                   
#>  [9] "p_valueMatrix"                 "means"                        
#> [11] "pairs"                         "settings"
pm$`Predicted Means`
#> nitro
#>        0      0.2      0.4      0.6 
#>  79.3889  98.8889 114.2222 123.3889
pm$means
#>   nitro  estimate std.error      df  conf.low conf.high group  n estimable
#> 1     0  79.38889  7.174754 6.79202  62.31741  96.46037  A    18      TRUE
#> 2   0.2  98.88889  7.174754 6.79202  81.81741 115.96037   B   18      TRUE
#> 3   0.4 114.22222  7.174754 6.79202  97.15074 131.29370    C  18      TRUE
#> 4   0.6 123.38889  7.174754 6.79202 106.31741 140.46037     D 18      TRUE
```
