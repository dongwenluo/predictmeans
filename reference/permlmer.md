# Permutation Test of random or fixed effects for `lmer` model.

This function provides permutation tests for the terms in a linear mixed
model of `lmer`.

## Usage

``` r
permlmer(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed)
```

## Arguments

- lmer0:

  `lmer` model under H0, note that `lmer0` model must nest within
  `lmer1` model.

- lmer1:

  `lmer` model under H1, note that `lmer0` model must nest within
  `lmer1` model.

- nperm:

  Number of permutation, the default value is 999.

- ncore:

  Number of core for parallel computing, the default value is 3.

- plot:

  Plot permutation distribution or not, the default value is FALSE.

- seed:

  Specify a random number generator seed, for reproducible results.

## Value

Permutation p-value.

## References

Oliver E. Lee and Thomas M. Braun (2012), *Permutation Tests for Random
Effects in Linear Mixed Models. Biometrics*, Journal 68(2).

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
## Test random effects (99 permutations to keep the example quick;
## use 999 or more in practice)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
anova(fm1, fm2, fm3)
#> refitting model(s) with ML (instead of REML)
#> Data: sleepstudy
#> Models:
#> fm3: Reaction ~ Days + (1 | Subject)
#> fm2: Reaction ~ Days + (Days || Subject)
#> fm1: Reaction ~ Days + (Days | Subject)
#>     npar    AIC    BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
#> fm3    4 1802.1 1814.8 -897.04    1794.1                          
#> fm2    5 1762.0 1778.0 -876.00    1752.0 42.0754  1  8.782e-11 ***
#> fm1    6 1763.9 1783.1 -875.97    1751.9  0.0639  1     0.8004    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# \donttest{
permlmer(fm3, fm2, nperm=99, ncore=1)
#> Data: sleepstudy
#> Models:
#> lmer0: Reaction ~ Days + (1 | Subject)
#> lmer1: Reaction ~ Days + (Days || Subject)
#>       npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> lmer0    4 1794.5 1807.2 -893.23    1786.5                            
#> lmer1    5 1753.7 1769.6 -871.83    1743.7 42.796  1 6.0763e-11   0.01
permlmer(fm2, fm1, nperm=99, ncore=1)
#> Warning: Model failed to converge with max|grad| = 0.00228659 (tol = 0.002, component 1)
#>   See ?lme4::convergence and ?lme4::troubleshooting.
#> Data: sleepstudy
#> Models:
#> lmer0: Reaction ~ Days + (Days || Subject)
#> lmer1: Reaction ~ Days + (Days | Subject)
#>       npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq) Perm-p
#> lmer0    5 1753.7 1769.6 -871.83    1743.7                           
#> lmer1    6 1755.6 1774.8 -871.81    1743.6 0.041  1     0.8395   0.81
# }

## Test fixed effects
Oats$nitro <- factor(Oats$nitro)
fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# \donttest{
permlmer(fm0, fm, nperm=99, ncore=1)
#> refitting model(s) with ML (instead of REML)
#> Data: Oats
#> Models:
#> lmer0: yield ~ nitro + Variety + (1 | Block/Variety)
#> lmer1: yield ~ nitro * Variety + (1 | Block/Variety)
#>       npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> lmer0    9 616.04 636.53 -299.02    598.04                            
#> lmer1   15 625.91 660.06 -297.95    595.91 2.1375  6    0.90663   0.91
# }
```
