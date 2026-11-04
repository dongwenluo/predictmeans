# Permutation ANOVA for `lmer` Model

This function provides permutation ANOVA for `lmer` model.

## Usage

``` r
permanova.lmer(model, nperm = 999, ncore=3, type = c("I", "II", "III",  "1", "2", "3"),
 ...)
```

## Arguments

- model:

  Model object returned by `lmer`.

- nperm:

  Number of permutation, the default value is 999.

- ncore:

  Number of core for parallel computing, the default value is 3.

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
# library(predictmeans)
# Oats$nitro <- factor(Oats$nitro) 
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## Permutation Test for model terms
# permanova.lmer(fm)
# permanova.lmer(fm, type=2)
## Compare to F test
# fm0 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# anova(fm0)
```
