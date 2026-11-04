# Permutation Test of Linear Model

This function provides permutation t-tests for coefficients of (fixed)
effects and permutation F-tests for the terms in a linear model such as
`aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

## Usage

``` r
permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
  test.statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, 
  fo=NULL, prt=TRUE, ncore=3, seed)
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

- nperm:

  The number of permutations. The default is 4999.

- type:

  type of ANOVA test, "I", "II", "III", 1, 2, or 3. Roman numerals are
  equivalent to the corresponding Arabic numerals.

- test.statistic:

  For type I ANOVA, F test is applied to all models, while for type II
  and III ANOVA, F test is performed for `lm`, Chisq test for `lme` and
  `gls` model, Chisq or F tests for `lmer` model and Likelihood ratio,
  Wald or F tests for `glm` model.

- exact:

  A logical variable to indicate whether or not exact no. of
  permutations will be used (applicable only to free the permutation
  case). The default is FALSE.

- data:

  In some cases, you need to provide the data set used in model fitting,
  especially when you have applied some variable trnasformation in the
  model.

- fo:

  A model formula used in the `model`; `fo!=NULL` when the formula is
  specified by function `formula`.

- prt:

  A logical variable to indicate whether or not to print output on the
  screen. The default is TRUE.

- ncore:

  Number of core for parallel computing, the default value is 3.

- seed:

  Specify a random number generator seed, for reproducible results.

## Value

The function produces permutation t-test table for coefficients of
(fixed) effects, permutation ANOVA table for model terms and a model
parameter list `permlist`, a list containing `nsim=4999` times
permutation refitted `model` parameters which are used in functions
`predictmeans` and `contrastmeans`.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
## Not run for simplifying process of submiting pkg to CRAN
#library(predictmeans)
#Oats$nitro <- factor(Oats$nitro) 
#fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
## library(lme4)
## fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#
## Permutation Test for model terms
#system.time(
#  permlme <- permmodels(model=fm, nperm=999)
#)  
#
## Permutation Test for multiple comparisons
#predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH", 
#  permlist=permlme, plot=FALSE)
#
## Permutation Test for specified contrasts
#cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
```
