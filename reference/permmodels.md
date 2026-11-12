# Permutation Test of Linear Model

This function provides permutation t-tests for coefficients of (fixed)
effects and permutation F-tests for the terms in a linear model such as
`aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

## Usage

``` r
permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
  test.statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, 
  fo=NULL, prt=TRUE, ncore=3, seed, method=c("raw", "freedman-lane"))
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

- method:

  How the data are permuted: `"raw"` (the default) or `"freedman-lane"`
  (`lm` and `aov` models only). See Details.

## Value

The function produces permutation t-test table for coefficients of
(fixed) effects, permutation ANOVA table for model terms and a model
parameter list `permlist`, a list containing `nsim=4999` times
permutation refitted `model` parameters which are used in functions
`predictmeans` and `contrastmeans`. With `method = "freedman-lane"`,
`permlist` is `NULL`.

## Details

With `method = "raw"`, the response is permuted for `lm`, `aov` and
`glm` models; for `gls`, `lme` and `lmer` models the fitted values and
the decorrelated residuals of the model are permuted together. Each
permutation refits the model, and the refits are kept in `permlist`. The
residuals are decorrelated with the symmetric square root of their
covariance, so the results do not depend on the order of the data
(predictmeans 1.x used the Cholesky factor, which made tests of terms
that vary within blocks reject too often).

With `method = "freedman-lane"` (Freedman and Lane, 1983), each test has
its own permutations: the residuals of the model without the tested term
(for type I and II tests also without the terms that contain it; for a t
test, without that coefficient) are permuted and added back to that
model's fitted values, and the statistic is computed again. The
statistics come from QR projections, so no model is refitted and it is
fast; weights and offsets are allowed. Use it when the model has a
covariate with extreme values (skewed, or with outliers) that is related
to the tested term: in a simulation study, permuting the response then
rejected 6 to 7% of true null hypotheses at the 5% level, and
Freedman-Lane 4 to 5% (Anderson and Legendre (1999) found the same). In
other designs the two methods agreed.

For mixed models with a covariate with extreme values related to the
tested term, the same study found the `"raw"` scheme slightly liberal
(6.3% at the 5% level), while the Satterthwaite F test kept the level
(4.6%); prefer the Kenward-Roger or Satterthwaite F test there (e.g.
[`anova()`](https://rdrr.io/r/stats/anova.html) of an lmerTest fit, or
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)).

## References

Anderson, M. J. and Legendre, P. (1999). An empirical comparison of
permutation methods for tests of partial regression coefficients in a
linear model. *Journal of Statistical Computation and Simulation*, 62,
271-303.

Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of
reported significance levels. *Journal of Business & Economic
Statistics*, 1, 292-298.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

# \donttest{
## Permutation test for model terms (99 permutations to keep the example
## quick; use the default or more in practice)
permlme <- permmodels(model=fm, nperm=99, ncore=1)
#> 
#> Coefficients of (fixed) effects:
#>                              Value Std.Error DF t-value p-value Perm_p_value
#> (Intercept)                80.0000    9.1070 45  8.7845  0.0000           NA
#> nitro0.2                   18.5000    7.6830 45  2.4079  0.0202         0.01
#> nitro0.4                   34.6667    7.6830 45  4.5122  0.0000         0.01
#> nitro0.6                   44.8333    7.6830 45  5.8354  0.0000         0.01
#> VarietyMarvellous           6.6667    9.7150 10  0.6862  0.5082         0.44
#> VarietyVictory             -8.5000    9.7150 10 -0.8749  0.4021         0.37
#> nitro0.2:VarietyMarvellous  3.3333   10.8653 45  0.3068  0.7604         0.71
#> nitro0.4:VarietyMarvellous -4.1667   10.8653 45 -0.3835  0.7032         0.69
#> nitro0.6:VarietyMarvellous -4.6667   10.8653 45 -0.4295  0.6696         0.68
#> nitro0.2:VarietyVictory    -0.3333   10.8653 45 -0.0307  0.9757         0.96
#> nitro0.4:VarietyVictory     4.6667   10.8653 45  0.4295  0.6696         0.66
#> nitro0.6:VarietyVictory     2.1667   10.8653 45  0.1994  0.8428         0.84
#> 
#> Note: Perm_p_value of t test is obtained using ‘99’ permutations.
#> 
#> ANOVA:
#>               numDF denDF  F-value p-value Perm_p_value
#> (Intercept)       1    45 245.1433  0.0000             
#> nitro             3    45  37.6856  0.0000         0.01
#> Variety           2    10   1.4853  0.2724         0.24
#> nitro:Variety     6    45   0.3028  0.9322         0.88
#> 
#> Note: Perm_p_value of F test is obtained using ‘99’ permutations.
#> 

## Permutation test for multiple comparisons
predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH",
  permlist=permlme, plot=FALSE, ncore=1)
#> $`Predicted Means`
#>       Variety Golden Rain Marvellous  Victory
#> nitro                                        
#> 0                 80.0000    86.6667  71.5000
#> 0.2               98.5000   108.5000  89.6667
#> 0.4              114.6667   117.1667 110.8333
#> 0.6              124.8333   126.8333 118.5000
#> 
#> $`Standard Error of Means`
#> All means have the same SE 
#>                    9.10696 
#> 
#> $`Standard Error of Differences`
#>  Max.SED  Min.SED Aveg.SED 
#> 9.715028 7.682958 9.160827 
#> attr(,"For the Same Level of Factor")
#>             nitro  Variety
#> Aveg.SED 9.715028 7.682958
#> Min.SED  9.715028 7.682958
#> Max.SED  9.715028 7.682958
#> 
#> $`Approximated LSD`
#>  Max.LSD  Min.LSD Aveg.LSD 
#> 19.43006 15.36592 18.32165 
#> attr(,"For the Same Level of Factor")
#>             nitro  Variety
#> Aveg.LSD 19.43006 15.36592
#> Min.LSD  19.43006 15.36592
#> Max.LSD  19.43006 15.36592
#> attr(,"Note")
#> [1] "This is a approximate LSD (i.e. 2*SED) at 0.05 level."
#> 
#> $`Pairwise '99' times permuted p-value (adjusted by 'BH' method)`
#> [1] "For variable 'nitro' at each level of 'Variety'"
#> 
#> $`Golden Rain`
#>     0     0.2   0.4  0.6 Group
#> 0   1                    A    
#> 0.2 0.015 1               B   
#> 0.4 0.015 0.06  1         BC  
#> 0.6 0.015 0.015 0.17 1     C  
#> 
#> $Marvellous
#>     0    0.2  0.4  0.6 Group
#> 0   1                  A    
#> 0.2 0.02 1              B   
#> 0.4 0.02 0.25 1         B   
#> 0.6 0.02 0.06 0.24 1    B   
#> 
#> $Victory
#>     0     0.2   0.4  0.6 Group
#> 0   1                    A    
#> 0.2 0.036 1               B   
#> 0.4 0.015 0.015 1          C  
#> 0.6 0.015 0.015 0.28 1     C  
#> 
#> $mean_table
#> $mean_table$Table
#>        Variety nitro     Mean      SE Df  LL(95%)  UL(95%) LetterGrp
#> 1  Golden Rain     0  80.0000 9.10696 NA  61.7861  98.2139       A  
#> 4  Golden Rain   0.2  98.5000 9.10696 NA  80.2861 116.7139        B 
#> 7  Golden Rain   0.4 114.6667 9.10696 NA  96.4528 132.8806        BC
#> 10 Golden Rain   0.6 124.8333 9.10696 NA 106.6194 143.0472         C
#> 2   Marvellous     0  86.6667 9.10696 NA  68.4528 104.8806        A 
#> 5   Marvellous   0.2 108.5000 9.10696 NA  90.2861 126.7139         B
#> 8   Marvellous   0.4 117.1667 9.10696 NA  98.9528 135.3806         B
#> 11  Marvellous   0.6 126.8333 9.10696 NA 108.6194 145.0472         B
#> 3      Victory     0  71.5000 9.10696 NA  53.2861  89.7139       A  
#> 6      Victory   0.2  89.6667 9.10696 NA  71.4528 107.8806        B 
#> 9      Victory   0.4 110.8333 9.10696 NA  92.6194 129.0472         C
#> 12     Victory   0.6 118.5000 9.10696 NA 100.2861 136.7139         C
#> 
#> $mean_table$Note
#> [1] "Letter-based representation of pairwise comparisons at significant level '0.05' at each level of Variety"
#> 
#> 

## Permutation test for specified contrasts
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
#> $`The t tests of the specified contrasts`
#>      Estimate Std. Error t value Permuted Pr(>|t|)
#> [1,]  -8.5000      9.715 -0.8749              0.37
#> [2,] -45.6667      9.715 -4.7006              0.01
#> attr(,"Note")
#> [1] "The permuted p-value is obtained using ‘99’ permutations."
#> 
#> $K
#>      0:Golden Rain 0:Marvellous 0:Victory 0.2:Golden Rain 0.2:Marvellous
#> [1,]            -1            0         1               0              0
#> [2,]             0            0         1               0              0
#>      0.2:Victory 0.4:Golden Rain 0.4:Marvellous 0.4:Victory 0.6:Golden Rain
#> [1,]           0               0              0           0               0
#> [2,]           0               0             -1           0               0
#>      0.6:Marvellous 0.6:Victory
#> [1,]              0           0
#> [2,]              0           0
#> 
# }

## Freedman-Lane permutations for a linear model with a covariate
fm2 <- lm(mpg ~ wt + factor(cyl), data = mtcars)
permmodels(fm2, nperm = 999, type = "II", seed = 1, method = "freedman-lane")
#> 
#> Coefficients of (fixed) effects:
#>              Estimate Std. Error t value Pr(>|t|) Perm_p_value
#> (Intercept)   33.9908     1.8878 18.0056   0.0000           NA
#> wt            -3.2056     0.7539 -4.2521   0.0002        0.001
#> factor(cyl)6  -4.2556     1.3861 -3.0702   0.0047        0.009
#> factor(cyl)8  -6.0709     1.6523 -3.6742   0.0010        0.003
#> 
#> Note: Perm_p_value of t test is obtained using ‘999’ Freedman-Lane permutations.
#> 
#> ANOVA:
#> Anova Table (Type II tests)
#> 
#>  Response: mpg
#>               Sum Sq Df F value Pr(>F) Perm_p_value
#> wt          118.2039  1 18.0801 0.0002        0.001
#> factor(cyl)  95.2633  2  7.2856 0.0028        0.005
#> Residuals   183.0586 28                            
#> 
#> Note: Perm_p_value of F test is obtained using ‘999’ Freedman-Lane permutations.
#> 
```
