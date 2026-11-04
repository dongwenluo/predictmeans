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
  interval for `lmer` model, the default value is TRUE.

- level:

  level of confidence of CI, the default value is 0.95.

## Value

Variance components table.

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
#> Block.Block.var((Intercept))     214.4748 168.83266  45.84586 1003.3500
#> Variety.Variety.var((Intercept)) 106.0620  67.87561  30.25561  371.8034
#> sigma_sq                         177.0835  37.33247 117.14638  267.6869

data(Orthodont, package="nlme")
mod <- lmer(distance ~ age + (age|Subject), data=Orthodont)
if (FALSE) varcomp(mod) # \dontrun{}
mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
varcomp(mod1)
#>                                             vcov        SE       lower
#> Subject.Subject.var((Intercept))      5.41508652 4.9583408  0.89980808
#> Subject.Subject.cov(age,(Intercept)) -0.32106059 0.4243794 -0.49435952
#> Subject.Subject.var(age)              0.05126954 0.0414508  0.01050857
#> sigma_sq                              1.71620402 0.3302836  1.17696609
#>                                           upper
#> Subject.Subject.var((Intercept))     32.5882404
#> Subject.Subject.cov(age,(Intercept))  0.1570346
#> Subject.Subject.var(age)              0.2501354
#> sigma_sq                              2.5024988
```
