# Package index

## Predicted means and comparisons

Means of factor terms with their standard errors, intervals and pairwise
comparisons; contrasts, joint tests and slopes.

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

## Results and plots

Stable tables (tidy, as.data.frame), plots drawn from results, and the
advice kept with every result.

- [`tidy(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`glance(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`as.data.frame(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`coef(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  [`vcov(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md)
  : Predicted Means and Pairwise Comparisons as Data Frames
- [`tidy(`*`<pm_curves>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md)
  [`as.data.frame(`*`<pm_curves>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md)
  [`print(`*`<pm_curves>`*`)`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md)
  : Curves and Their Differences as Data Frames
- [`autoplot(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  [`plot(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  [`summary(`*`<predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  [`print(`*`<summary.predictmeans>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.predictmeans.md)
  : Plots and Summary of Predicted Means
- [`autoplot(`*`<pm_curves>`*`)`](https://dongwenluo.github.io/predictmeans/reference/autoplot.pm_curves.md)
  : Plots of Curves, Surfaces and Their Differences
- [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md)
  : Level Plot of a Matrix of p-values
- [`advice()`](https://dongwenluo.github.io/predictmeans/reference/advice.md)
  : Advice on a Result

## Degrees of freedom, tests and diagnostics

F tests with Kenward-Roger df, R-squared, variance components and model
checks.

- [`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md)
  : ANOVA Tables for lmer Fits with Kenward-Roger or Satterthwaite df
- [`df_term()`](https://dongwenluo.github.io/predictmeans/reference/df_term.md)
  : Denominator Degrees of Freedom of the Means or Contrasts of an lmer
  Model
- [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
  : R-squared of a Model as Numbers
- [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  : Variance Components with Standard Errors and Confidence Intervals
- [`se_ranef`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  [`print.pm_ranef`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  [`tidy.pm_ranef`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  [`as.data.frame.pm_ranef`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  [`as.list.pm_ranef`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)
  : Standard Errors of the Predicted Random Effects
- [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)
  : Diagnostic Plots for a Fitted Model
- [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  : Cook's Distances of a Linear (Mixed) Model
- [`pm_check()`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md)
  : Checks of a Fitted Model
- [`pm_boxcox()`](https://dongwenluo.github.io/predictmeans/reference/pm_boxcox.md)
  [`pm_yeojohnson()`](https://dongwenluo.github.io/predictmeans/reference/pm_boxcox.md)
  : Box-Cox and Yeo-Johnson Transformations of a Response

## Permutation and bootstrap tests

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  : Permutation Tests of the Terms and Coefficients of a Model
- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  : Permutation Test of Random or Fixed Effects of Mixed Models
- [`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)
  : Permutation Index
- [`permstats()`](https://dongwenluo.github.io/predictmeans/reference/permstats.md)
  : Effect Sizes and Model Checks from Permutations
- [`permcombine()`](https://dongwenluo.github.io/predictmeans/reference/permcombine.md)
  : Combine Permutation Tests of Several Responses
- [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  : Parametric Bootstrap Tests

## Semiparametric regression

Penalized splines and surfaces fitted as mixed models, with predictions,
differences and derivatives.

- [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  : Semiparametric Mixed Models Fitted with lme4
- [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  : Semiparametric Mixed Models Fitted with glmmTMB
- [`anova.semireg`](https://dongwenluo.github.io/predictmeans/reference/anova.semireg.md)
  : Tests for a Semi-Parametric Model
- [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  : Predicted Curves and Surfaces of a Semiparametric Model
- [`semiboundary()`](https://dongwenluo.github.io/predictmeans/reference/semiboundary.md)
  : Boundary of the region covered by two-dimensional data

## Extending to new model classes

- [`pm_coef()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_vcov()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_basis()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_df()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_linkinv()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  [`pm_contrasts()`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)
  : Supporting a New Model Class
- [`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md)
  : Data of a Fitted Model, with the Response on Both Scales

## Data

- [`ATP`](https://dongwenluo.github.io/predictmeans/reference/ATP.md) :
  ATP Content of Preserved Dog Hearts
- [`Body_length`](https://dongwenluo.github.io/predictmeans/reference/Body_length.md)
  : Body Length of Beetles (Simulated)
- [`Clinical`](https://dongwenluo.github.io/predictmeans/reference/Clinical.md)
  : Multicentre Clinical Trial of a Skin Treatment
- [`Caschool`](https://dongwenluo.github.io/predictmeans/reference/Caschool.md)
  : California Test Score Data
- [`copper`](https://dongwenluo.github.io/predictmeans/reference/copper.md)
  [`copperKnots`](https://dongwenluo.github.io/predictmeans/reference/copper.md)
  : Copper Grade at Sampling Points
- [`DataExam5.2`](https://dongwenluo.github.io/predictmeans/reference/DataExam5.2.md)
  : Heights of Trees: 37 Seedlots at 6 Sites
- [`Drug`](https://dongwenluo.github.io/predictmeans/reference/Drug.md)
  : Potency of Analgesic Drugs in Mice
- [`growthIndiana`](https://dongwenluo.github.io/predictmeans/reference/growthIndiana.md)
  : Growth of Indiana Children
- [`ozoneCalif`](https://dongwenluo.github.io/predictmeans/reference/ozoneCalif.md)
  : Ozone Concentration and Meteorology in California
- [`Poisson_ind`](https://dongwenluo.github.io/predictmeans/reference/Poisson_ind.md)
  : Egg Counts of Beetles (Simulated)
- [`SBMDfem`](https://dongwenluo.github.io/predictmeans/reference/SBMDfem.md)
  : Spinal Bone Mineral Density of Girls
- [`scallop`](https://dongwenluo.github.io/predictmeans/reference/scallop.md)
  [`scallopKnots`](https://dongwenluo.github.io/predictmeans/reference/scallop.md)
  [`scallopBdry`](https://dongwenluo.github.io/predictmeans/reference/scallop.md)
  : Scallop Abundance off Long Island
- [`wageTUnion`](https://dongwenluo.github.io/predictmeans/reference/wageTUnion.md)
  : Trade Union Membership and Wages

## Package

- [`predictmeans-package`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-package.md)
  : Predicted Means for Linear and Semiparametric Models
- [`predictmeans-defunct`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  [`smZ`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  [`permglmer`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  [`permglmmTMB`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  [`permanova.lmer`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  [`ci_mcp`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  [`R2_glmm`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  : Functions removed in predictmeans 2.0
