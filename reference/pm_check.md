# Checks of a Fitted Model

A table of checks of a model's assumptions and structure, with a plain
note on each: the fit itself (singular or not converged), outliers and
influential observations, normality of the residuals and of the random
effects, equal variances across the levels of each factor,
overdispersion (Poisson and binomial models), interactions, whether a
factor with numeric levels follows a straight line, and coefficients
that cannot be estimated. For any model that
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
handles.

## Usage

``` r
pm_check(model, plot = FALSE)
```

## Arguments

- model:

  A fitted model, as for
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md).

- plot:

  `TRUE`: also draw the diagnostic plots of
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md).

## Value

A data frame of class `"pm_check"`: `check`, `result`, `p.value`, `look`
and `note`; printed as a table with the notes below it.

## Details

The checks inform; they do not decide. Rows marked `look = TRUE` are
worth a look at the plots of
[`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md).
Keep the terms the design calls for (blocks, whole plots, the treatment
structure): dropping a term because its test is not significant makes
the later tests and intervals too optimistic. Tests of normality and of
equal variances reject small departures in large data sets and miss
large ones in small data sets; the plots say how large a departure is.

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

- Equal variances:

  for each factor: the ratio of the largest to the smallest standard
  deviation of the residuals within its levels, and the Brown-Forsythe
  test (ANOVA of the absolute deviations from the medians).

- Overdispersion:

  Poisson and binomial (with more than one trial) models: Pearson
  chi-square / residual df, about 1 when the variance function is right.

- Interaction:

  the Type III test of each interaction term (as
  `permmodels(type = "III")`).

- Straight line:

  for a factor with numeric levels (e.g. nitrogen rates 0, 0.2, 0.4,
  0.6) and at least three of them: the joint test of the polynomial
  contrasts above the linear one (with the actual spacing of the
  levels), from
  [`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md).

- Estimability:

  coefficients dropped as aliased (empty cells or confounded factors).

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
#>    check                                     
#>    Fit                                       
#>    Outliers (|residual| > 3)                 
#>    Normality of residuals                    
#>    Normality of random effects: Variety:Block
#>    Normality of random effects: Block        
#>    Equal variances: nitro                    
#>    Equal variances: Variety                  
#>    Interaction: nitro:Variety                
#>    Straight line in nitro                    
#>  result                                   p.value
#>  no warnings                                     
#>  0 of 72 (expected 0.2)                          
#>  Shapiro-Wilk                             0.270  
#>  Shapiro-Wilk, 18 levels                  0.564  
#>  Shapiro-Wilk, 6 levels                   0.031  
#>  SD ratio 1.32 (0.4 / 0)                  0.340  
#>  SD ratio 1.37 (Golden Rain / Marvellous) 0.060  
#>  Type III test                            0.932  
#>  curvature, 2 df                          0.265  
#> 
#> - Normality of random effects: Block: Few levels: little power.
#> - Interaction: nitro:Variety: No evidence of interaction; keep the term
#>     if the design calls for it.
#> - Straight line in nitro: The means of nitro are near a straight line:
#>     nitro could also enter as a number (a response line); as a factor
#>     it compares the levels.
#> 
#> Checks inform, they do not decide; see residplot() for the plots.

# counts: overdispersion
gm <- glm(breaks ~ wool * tension, family = poisson, data = warpbreaks)
pm_check(gm)
#> Checks of a glm model (54 observations)
#> 
#>    check                                   
#>    Fit                                     
#>  ! Outliers (|residual| > 3)               
#>    Influence (Cook's distance)             
#>    Normality of residuals (quantile)       
#>    Equal variances: wool                   
#>    Equal variances: tension                
#>  ! Overdispersion (Pearson chi-square / df)
#>    Interaction: wool:tension               
#>  result                                                   p.value
#>  no warnings                                                     
#>  6 of 54 (expected 0.1): 5 (3.6), 24 (3.4), 4 (-3.2), ...        
#>  largest 0.341 (5)                                               
#>  Shapiro-Wilk                                             0.296  
#>  SD ratio 1.41 (A / B)                                    0.023  
#>  SD ratio 1.42 (L / H)                                    0.171  
#>  3.76                                                     <1e-04 
#>  Type III test                                            <1e-04 
#> 
#> ! Outliers (|residual| > 3): Check these observations (data entry,
#>     unusual units) before anything else.
#> ! Overdispersion (Pearson chi-square / df): More variation than the
#>     distribution allows: SEs are too small. Consider a negative
#>     binomial model or an observation-level random effect.
#> - Interaction: wool:tension: The effect of one factor depends on the
#>     other: compare the levels of one within each level of the other
#>     (atvar or pm_joint(by = ...)).
#> 
#> Checks inform, they do not decide; see residplot() for the plots.
```
