# Calculate SE and CI of variance components for `lmer`, `glmer`, `lme`, `glmmTMB` model

This function calculates SE and CI of variance components for `lmer`,
`glmer`, `lme`, `glmmTMB` model. For `lme` models it needs package
lmeInfo.

## Usage

``` r
varcomp(model, ci=TRUE, level=0.95)
```

## Arguments

- model:

  Model object returned by `lmer`, `glmer`, `lme`, `glmmTMB`.

- ci:

  a logical value to indicates wheather or not to simulate a confidence
  interval for `lmer` and `glmmTMB` models, the default value is TRUE.

- level:

  level of confidence of CI, the default value is 0.95.

## Value

Variance components table.

## Details

For `glmmTMB` models the table holds the variances and covariances of
each random-effect term (conditional model, then the zero-inflation
model with prefix `zi.`) and, for gaussian models, the residual
variance. SEs are from the delta method. CIs are Wald intervals on the
log scale for variances and on the Fisher z (atanh) scale for
correlations, times the two standard deviations for covariances.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
if (FALSE) varcomp(fm) # \dontrun{}
fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
varcomp(fm1)
#>                                      vcov        SE     lower     upper
#> Block.Block.var((Intercept))     214.4745 168.83242  45.84595 1003.3449
#> Variety.Variety.var((Intercept)) 106.0618  67.87554  30.25561  371.8022
#> sigma_sq                         177.0835  37.33248 117.14643  267.6870

data(Orthodont, package="nlme")
mod <- lmer(distance ~ age + (age|Subject), data=Orthodont)
if (FALSE) varcomp(mod) # \dontrun{}
mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
varcomp(mod1)
#>                                             vcov        SE       lower
#> Subject.Subject.var((Intercept))      5.41508791 4.9583411  0.89986361
#> Subject.Subject.cov(age,(Intercept)) -0.32106066 0.4243794 -0.49433350
#> Subject.Subject.var(age)              0.05126955 0.0414508  0.01051371
#> sigma_sq                              1.71620400 0.3302836  1.17675496
#>                                           upper
#> Subject.Subject.var((Intercept))     32.5862463
#> Subject.Subject.cov(age,(Intercept))  0.1568360
#> Subject.Subject.var(age)              0.2500132
#> sigma_sq                              2.5029477
```
