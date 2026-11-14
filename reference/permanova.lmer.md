# Permutation ANOVA for `lmer` Model

This function provides permutation ANOVA for `lmer` model (and, since
predictmeans 2.0, `glmer` model): each term is tested by
[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
comparing the model with and without it.

## Usage

``` r
permanova.lmer(model, nperm = 999, ncore=1, type = c("I", "II", "III",  "1", "2", "3"),
 ...)
```

## Arguments

- model:

  Model object returned by `lmer` or `glmer`.

- nperm:

  Number of permutation, the default value is 999.

- ncore:

  Number of cores for parallel computing, the default value is 1 (before
  predictmeans 2.0: 3).

- type:

  The type of ANOVA table requested (using SAS terminology) with Type I
  being the familiar sequential ANOVA table.

- ...:

  Use to setup option: seed – Specify a random number generator seed,
  for reproducible results.

## Value

Permutation ANOVA table.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## Permutation test for model terms (99 permutations to keep the example
## quick; use 999 or more in practice)
# \donttest{
permanova.lmer(fm, nperm=99, ncore=1)
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> Analysis of Variance Table
#>               npar  Sum Sq Mean Sq F value Perm.p
#> nitro            3 20020.5  6673.5 37.6857   0.01
#> Variety          2   526.1   263.0  1.4853   0.29
#> nitro:Variety    6   321.7    53.6  0.3028   0.96
permanova.lmer(fm, nperm=99, ncore=1, type=2)
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> refitting model(s) with ML (instead of REML)
#> Analysis of Variance Table
#>               npar  Sum Sq Mean Sq F value Perm.p
#> nitro            3 20020.5  6673.5 37.6857   0.01
#> Variety          2   526.1   263.0  1.4853   0.23
#> nitro:Variety    6   321.7    53.6  0.3028   0.94
# }
## Compare to F test
fm0 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
anova(fm0)
#>               numDF denDF   F-value p-value
#> (Intercept)       1    45 245.14333  <.0001
#> nitro             3    45  37.68561  <.0001
#> Variety           2    10   1.48534  0.2724
#> nitro:Variety     6    45   0.30282  0.9322
```
