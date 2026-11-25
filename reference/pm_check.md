# Checks of a Fitted Model

A table of checks of a model's assumptions and structure, grouped in
sections, with a plain note on each: the fit itself (singular or not
converged); the residuals (outliers, influence, normality, a
transformation); the variance (equal across factor levels, constant with
the fitted values); the distribution of generalized models (over- and
underdispersion, excess zeros, separation); the design (balance and
empty cells, aliasing, collinearity, few levels of a random effect); and
the terms (interactions, straight lines). For any model that
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
handles.

## Usage

``` r
pm_check(model, plot = FALSE, row = NULL, col = NULL, coords = NULL, by = NULL, data = NULL)
```

## Arguments

- model:

  A fitted model, as for
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md).

- plot:

  `TRUE`: also draw the diagnostic plots of
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md).

- row, col, coords, by, data:

  The layout of the plots, as in
  [`residvariogram`](https://dongwenluo.github.io/predictmeans/reference/residvariogram.md):
  when given, a row checks for spatial pattern in the residuals (the
  global test of the residual variogram, 99 simulations).

## Value

A data frame of class `"pm_check"`: `section`, `check`, `result`,
`p.value`, `look`, `status` (`"!"` worth a look, `"i"` information,
`"OK"`) and `note`; printed by section with the notes below.
[`summary()`](https://rdrr.io/r/base/summary.html) keeps only the rows
worth a look. The rows worth a look are also kept as advice
([`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).

## Details

The checks inform; they do not decide. Rows marked `!` are worth a look
at the plots of
[`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md);
rows marked `i` are information. Keep the terms the design calls for
(blocks, whole plots, the treatment structure): dropping a term because
its test is not significant makes the later tests and intervals too
optimistic. Tests of normality and of equal variances reject small
departures in large data sets and miss large ones in small data sets;
the plots say how large a departure is.

- Fit:

  singular fits (a variance component estimated as zero) and convergence
  warnings.

- Outliers:

  residuals beyond 3 (of the residuals
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)
  uses: standardized, normalized or randomized quantile), against the
  number expected for normal residuals; worth a look when that many, or
  one as large as the largest (Bonferroni over the observations), has
  probability below 0.01.

- Influence:

  (`lm`, `aov`, `glm`) the largest Cook's distance; above 0.5 is worth a
  look
  ([`CookD`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  for mixed models).

- Normality:

  Shapiro-Wilk tests of the residuals and of each term's random effects
  (intercepts).

- Transformation of the response:

  (`lm` and `aov` with a positive response that is not transformed in
  the formula) the Box-Cox power the residuals favour: the profile
  likelihood estimate of \\\lambda\\ and its 95% interval (Box and Cox,
  1964; as
  [`MASS::boxcox()`](https://rdrr.io/pkg/MASS/man/boxcox.html)). Worth a
  look when the interval leaves out 1 (no transformation); the note
  names the usual transformation within the interval.

- Equal variances:

  for each factor: the ratio of the largest to the smallest standard
  deviation of the residuals within its levels, and the Brown-Forsythe
  test (ANOVA of the absolute deviations from the medians).

- Spread vs fitted values:

  (gaussian models) Spearman correlation of the absolute residuals with
  the fitted values; worth a look when \\\|\rho\| \ge 0.2\\ and p \<
  0.01 (the variance grows with the mean).

- Dispersion:

  (Poisson, binomial with more than one trial, negative binomial; `glm`,
  `glmer`, `glmmTMB`) Pearson chi-square / residual df, about 1 when the
  variance function is right; overdispersion (ratio \> 1.5, p \< 0.01)
  makes SEs too small, underdispersion (ratio \< 0.7, p \< 0.01) too
  large. For mixed models the residual df counts the variance
  parameters, an approximation. A 0/1 response carries no information on
  dispersion; for families with an estimated dispersion (Gamma, quasi)
  its estimate is shown.

- Zero counts:

  (Poisson, negative binomial) the observed zeros against those
  simulated from the fitted model (200 simulations, fixed seed; the
  user's random numbers are not changed); worth a look when there are
  more (p \< 0.01).

- Separation:

  (binomial) fitted probabilities of 0 or 1, or coefficients beyond 10
  on the logit scale or with standard errors above 50: the estimates and
  Wald tests of these terms are unreliable.

- Balance of cells:

  of the crossed fixed factors: balanced, unbalanced or with empty
  cells; unbalanced data make Type I, II and III tests differ (see
  [`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md),
  section 'Which type of test').

- Estimability:

  coefficients dropped as aliased (empty cells or confounded factors).

- Collinearity:

  (models of main effects only) the generalized variance inflation
  factor of each term, \\GVIF^{1/(2 df)}\\ (Fox and Monette, 1992);
  above \\\sqrt{10}\\ is worth a look.

- Random-effect levels:

  a grouping factor with fewer than 5 levels: its variance is poorly
  estimated.

- Interaction:

  the Type III test of each interaction term (as
  `permmodels(type = "III")`).

- Straight line:

  for a factor with numeric levels and at least three of them: the joint
  test of the polynomial contrasts above the linear one (with the actual
  spacing of the levels), from
  [`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md).

## References

Fox, J. and Monette, G. (1992). Generalized collinearity diagnostics.
*Journal of the American Statistical Association*, 87, 178-183.

## See also

[`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md),
[`CookD`](https://dongwenluo.github.io/predictmeans/reference/CookD.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md),
[`pm_boxcox`](https://dongwenluo.github.io/predictmeans/reference/pm_boxcox.md),
[`residvariogram`](https://dongwenluo.github.io/predictmeans/reference/residvariogram.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
pm_check(fm)
#> Checks of a lmerMod model (72 observations)
#> 
#> Fit
#>   OK Fit                                         no warnings                   
#> 
#> Residuals
#>   OK Outliers (|residual| > 3)                   0 of 72 (expected 0.2)        
#>   OK Normality of residuals                      Shapiro-Wilk             0.270
#>   OK Normality of random effects: Variety:Block  Shapiro-Wilk, 18 le...   0.564
#>   OK Normality of random effects: Block          Shapiro-Wilk, 6 levels   0.031
#> 
#> Variance
#>   OK Equal variances: nitro                      SD ratio 1.32 (0.4 ...   0.340
#>   OK Equal variances: Variety                    SD ratio 1.37 (Gold...   0.060
#>   OK Spread vs fitted values                     Spearman rho 0.04        0.751
#> 
#> Design
#>   OK Balance of cells                            balanced (6 per cell)         
#> 
#> Terms
#>   OK Interaction: nitro:Variety                  Type III test            0.932
#>   OK Straight line in nitro                      curvature, 2 df          0.265
#> 
#> Notes
#> - Normality of random effects: Variety:Block: Shapiro-Wilk, 18 levels.
#> - Normality of random effects: Block: Few levels: little power.
#> - Equal variances: nitro: SD ratio 1.32 (0.4 / 0).
#> - Equal variances: Variety: SD ratio 1.37 (Golden Rain / Marvellous).
#> - Interaction: nitro:Variety: No evidence of interaction; keep the term if
#>     the design calls for it.
#> - Straight line in nitro: The means of nitro are near a straight line:
#>     nitro could also enter as a number (a response line); as a factor it
#>     compares the levels.
#> 
#> Nothing flagged.
#> (! worth a look, i information, OK no concern.) Checks inform, they do not decide; see residplot().

# counts: overdispersion
gm <- glm(breaks ~ wool * tension, family = poisson, data = warpbreaks)
pm_check(gm)
#> Checks of a glm model (54 observations)
#> 
#> Fit
#>   OK Fit                                   no warnings                         
#> 
#> Residuals
#>   !  Outliers (|residual| > 3)             6 of 54 (expected 0.1): 5...        
#>   OK Influence (Cook's distance)           largest 0.341 (5)                   
#>   OK Normality of residuals (quantile)     Shapiro-Wilk                   0.296
#> 
#> Variance
#>   OK Equal variances: wool                 SD ratio 1.41 (A / B)          0.023
#>   OK Equal variances: tension              SD ratio 1.42 (L / H)          0.171
#> 
#> Distribution
#>   !  Dispersion (Pearson chi-square / df)  3.76 on 48 df                 <1e-04
#> 
#> Design
#>   OK Balance of cells                      balanced (9 per cell)               
#> 
#> Terms
#>   OK Interaction: wool:tension             Type III test                 <1e-04
#> 
#> Notes
#> ! Outliers (|residual| > 3): 6 of 54 (expected 0.1): 5 (3.6), 24 (3.4), 4
#>     (-3.2), .... Check these observations (data entry, unusual units)
#>     before anything else.
#> ! Dispersion (Pearson chi-square / df): Overdispersion: more variation than
#>     the poisson distribution allows, so SEs are too small and tests too
#>     liberal. Consider a negative binomial model, an observation-level
#>     random effect, or a quasi family (glm).
#> - Interaction: wool:tension: The effect of one factor depends on the other:
#>     compare the levels of one within each level of the other (atvar or
#>     pm_joint(by = ...)).
#> 
#> 2 worth a look: Outliers (|residual| > 3); Dispersion (Pearson chi-square / df).
#> (! worth a look, i information, OK no concern.) Checks inform, they do not decide; see residplot().
summary(pm_check(gm))
#> Checks of a glm model (54 observations)
#> 
#> Residuals
#>   !  Outliers (|residual| > 3)             6 of 54 (expected 0.1): 5...        
#> 
#> Distribution
#>   !  Dispersion (Pearson chi-square / df)  3.76 on 48 df                 <1e-04
#> 
#> Notes
#> ! Outliers (|residual| > 3): 6 of 54 (expected 0.1): 5 (3.6), 24 (3.4), 4
#>     (-3.2), .... Check these observations (data entry, unusual units)
#>     before anything else.
#> ! Dispersion (Pearson chi-square / df): Overdispersion: more variation than
#>     the poisson distribution allows, so SEs are too small and tests too
#>     liberal. Consider a negative binomial model, an observation-level
#>     random effect, or a quasi family (glm).
#> 
#> 2 worth a look: Outliers (|residual| > 3); Dispersion (Pearson chi-square / df).
#> (! worth a look, i information, OK no concern.) Checks inform, they do not decide; see residplot().
```
