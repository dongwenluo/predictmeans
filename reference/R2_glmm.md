# An adjusted coefficient of determination (R2) for generalized linear mixed models

This function produces adjusted R2 for generalized linear mixed models
which was crafted following the guidance provided by Professor
Hans-Peter Piepho.

## Usage

``` r
R2_glmm(model, over_disp=FALSE)
```

## Arguments

- model:

  An object returned by `lmer`, `glmer` or `glmmTMB`.

- over_disp:

  A logical scalar to indicate whether `model` with over-dispersion or
  not. The default value is FALSE.

## Value

Adjusted R2 in percentage for Total (fixed + random), Fiexd, Random and
individual random term.

## References

Piepho HP. An adjusted coefficient of determination (R2 ) for
generalized linear mixed models in one go. Biom J. 2023
Oct;65(7):e2200290. doi: 10.1002/bimj.202200290. Epub 2023 May 1. PMID:
37127864.

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  (fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats))
#> Linear mixed model fit by REML ['lmerModLmerTest']
#> Formula: yield ~ nitro * Variety + (1 | Block/Variety)
#>    Data: Oats
#> REML criterion at convergence: 529.0285
#> Random effects:
#>  Groups        Name        Std.Dev.
#>  Variety:Block (Intercept) 10.30   
#>  Block         (Intercept) 14.65   
#>  Residual                  13.31   
#> Number of obs: 72, groups:  Variety:Block, 18; Block, 6
#> Fixed Effects:
#>                (Intercept)                    nitro0.2  
#>                    80.0000                     18.5000  
#>                   nitro0.4                    nitro0.6  
#>                    34.6667                     44.8333  
#>          VarietyMarvellous              VarietyVictory  
#>                     6.6667                     -8.5000  
#> nitro0.2:VarietyMarvellous  nitro0.4:VarietyMarvellous  
#>                     3.3333                     -4.1667  
#> nitro0.6:VarietyMarvellous     nitro0.2:VarietyVictory  
#>                    -4.6667                     -0.3333  
#>    nitro0.4:VarietyVictory     nitro0.6:VarietyVictory  
#>                     4.6667                      2.1667  
  R2_glmm(fm)
#> # Adjusted R2 for Mixed Models
#>         Total         Fixed        Random Variety:Block         Block 
#>      "75.81%"      "37.19%"      "38.63%"      "13.87%"      "24.75%" 
  (gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
#> Generalized linear mixed model fit by maximum likelihood (Laplace
#>   Approximation) [glmerMod]
#>  Family: binomial  ( logit )
#> Formula: cbind(incidence, size - incidence) ~ period + (1 | herd)
#>    Data: cbpp
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>  194.0526  204.1793  -92.0263  184.0526        51 
#> Random effects:
#>  Groups Name        Std.Dev.
#>  herd   (Intercept) 0.6423  
#> Number of obs: 56, groups:  herd, 15
#> Fixed Effects:
#> (Intercept)      period2      period3      period4  
#>     -1.3985      -0.9923      -1.1287      -1.5803  
  R2_glmm(gm)
#> # Adjusted R2 for Mixed Models
#>    Total    Fixed   Random     herd 
#>  "32.1%" "13.52%" "18.58%" "18.58%" 
```
