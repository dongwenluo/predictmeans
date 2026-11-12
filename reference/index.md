# Package index

## Predicted means and comparisons

- [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  : Predicted Means of a Linear Model
- [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  : Predicted Means of a Linear Model with Covariate Variable(s)
- [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
  : Linear Contrast Tests for a Linear Model
- [`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md)
  : Contrasts Written with Level Names
- [`contrast_family()`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md)
  : Families of Contrasts among Means
- [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
  : Joint F Tests of Contrasts among Means (Slices)
- [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md)
  : Slopes of a Covariate by Group
- [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  : Matrix of Coefficients in a Linear Model

## Results

- [`tidy(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`as.data.frame(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`coef(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`vcov(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  : Predicted Means and Pairwise Comparisons as Data Frames
- [`autoplot(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  [`summary(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  [`print(`*`<summary.predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  : Plots and Summary of Predicted Means
- [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md)
  : Level Plot of a Matrix of p-values.
- [`ci_mcp()`](https://dongwenluo.github.io/predictmeans/reference/ci_mcp.md)
  : Multiple Comparisons Based on the Confidence Intervals

## Degrees of freedom, tests and diagnostics

- [`df_term()`](https://dongwenluo.github.io/predictmeans/reference/df_term.md)
  : Calculate degree of freedom of a modelterm (contrast) for a lmer
  model

- [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md)
  : An adjusted coefficient of determination (R2) for generalized linear
  mixed models

- [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
  : R-squared for Mixed Models as Numbers

- [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  :

  Calculate SE and CI of variance components for `lmer`, `glmer`, `lme`,
  `glmmTMB` model

- [`se_ranef()`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  : Extract Standard Errors of Model Random Effects

- [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  : Calculates and plots Cook's distances for a Linear (Mixed) Model

- [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)
  : Diagnostic Plots for a Linear (Mixed) Model

## Permutation tests

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  : Permutation Test of Linear Model

- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  [`permglmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  [`permglmmTMB()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  :

  Permutation Test of random or fixed effects for `lmer`, `glmer` and
  `glmmTMB` models.

- [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/permanova.lmer.md)
  :

  Permutation ANOVA for `lmer` Model

- [`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)
  : Permutation Index

- [`permstats()`](https://dongwenluo.github.io/predictmeans/reference/permstats.md)
  : Effect Sizes and Model Checks from Permutations

- [`permcombine()`](https://dongwenluo.github.io/predictmeans/reference/permcombine.md)
  : Combine Permutation Tests of Several Responses

- [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  : Parametric Bootstrap Tests for Mixed Models

## Semiparametric regression

- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  : Fitting Semi Parametric Models Using lme4 Ecosystem
- [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  : Fitting Semi Parametric Models Using glmmTMB
- [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  : Predicted Means of a Semi Paramatric Model with Covariate
  Variable(s)
- [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/smZ.md)
  : Generate Sparse Matrix Z for penalized spline smoothing

## Extending to new model classes

- [`pm_coef()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_vcov()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_basis()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_df()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_linkinv()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  : Supporting a New Model Class

## Data

- [`ATP`](https://dongwenluo.github.io/predictmeans/reference/ATP.md) :

  `ATP` containing data

- [`Body_length`](https://dongwenluo.github.io/predictmeans/reference/Body_length.md)
  :

  `Body_length` data

- [`Clinical`](https://dongwenluo.github.io/predictmeans/reference/Clinical.md)
  :

  `Clinical` data

- [`copper`](https://dongwenluo.github.io/predictmeans/reference/copper.md)
  :

  `copper` data

- [`Drug`](https://dongwenluo.github.io/predictmeans/reference/Drug.md)
  :

  `Drug` data

- [`Poisson_ind`](https://dongwenluo.github.io/predictmeans/reference/Poisson_ind.md)
  :

  `Poisson_ind` data

## Package

- [`predictmeans-package`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-package.md)
  : Predicted Means for Linear and Semiparametric Models
