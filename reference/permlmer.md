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
# library(predictmeans)
## Test random effects
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
# fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
# fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
# anova(fm1, fm2, fm3)
# permlmer(fm3, fm2)
# permlmer(fm2, fm1)

## Test fixed effects
# Oats$nitro <- factor(Oats$nitro)
# fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# permlmer(fm0, fm)
```
