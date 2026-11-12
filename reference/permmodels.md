# Permutation Test of Linear Model

This function provides permutation t-tests for coefficients of (fixed)
effects and permutation F-tests for the terms in a linear model such as
`aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer` and `glmmTMB`.

## Usage

``` r
permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
  test.statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, 
  fo=NULL, prt=TRUE, ncore=1, seed, method=c("auto", "raw", "freedman-lane"), perms=NULL,
  refit=TRUE, robust=FALSE)
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`,
  `glmer` or `glmmTMB`.

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

  Number of cores for parallel computing, the default value is 1 (before
  predictmeans 2.0: 3).

- seed:

  Specify a random number generator seed, for reproducible results.

- method:

  How the data are permuted: `"freedman-lane"` or `"raw"`; `"auto"` (the
  default) chooses the more robust one for the model class, which is
  `"freedman-lane"` for every class at present. `"raw"` is the scheme of
  predictmeans 1.x. See Details.

- perms:

  Optional matrix of permutations to use instead of free ones, one row
  per observation used in the model and one column per permutation (e.g.
  from
  [`permindex`](https://dongwenluo.github.io/predictmeans/reference/permindex.md),
  for permutations within blocks or of whole plots); `nperm` is then its
  number of columns. Use `exact = TRUE` when its columns are all the
  permutations the design allows.

- refit:

  For `method = "freedman-lane"` with mixed models (`lmer`, `lme`,
  `gls`, and `glmer`, `glmmTMB` through their working linear mixed
  model): refit the model to each permuted data set (the default), so
  the variance parameters are estimated again each time. `FALSE` holds
  the covariance of the data at the fit (no refits, much faster), which
  rejects too often when the variance parameters rest on few units. In a
  simulation (test of the whole-plot factor of a split plot, true null,
  5% level) `refit = FALSE` rejected 8.6%, 7.2%, 6.0% and 5.6% with 6,
  12, 24 and 48 blocks (whole-plot error df 15, 33, 69 and 141), and
  refits 5.0%, 6.2%, 5.4% and 4.6% (Kenward-Roger 4.6% to 6.4%); at the
  1% level `refit = FALSE` rejected up to 2.8%. So `refit = FALSE` suits
  a quick look, or a final test when the error df of the tested terms
  exceed about 100; otherwise keep the refits (for 199 permutations of
  the split plot's three terms on one core: 8 s with 6 blocks, 25 s with
  48).

- robust:

  For `lm`, `aov` and `glm` models with `method = "freedman-lane"`:
  tests for error variances that are unequal while the model assumes
  them equal. The statistics are heteroscedasticity-robust (HC3
  sandwich) t and Wald F statistics (DiCiccio and Romano, 2017), and the
  signs of the reduced model's residuals are flipped at random instead
  of permuted (a wild bootstrap; Davidson and Flachaire, 2008), so each
  residual stays with its own unit and keeps its variance. In a
  simulation with group SDs 3, 1, 1 (2000 data sets per case) the usual
  and permutation F tests rejected 19-24% of true null hypotheses at the
  5% level, these tests 5-7% (normal errors) or about 9% (strongly
  skewed errors, against 11-14% for the HC3 F test), and 5.4% with equal
  variances. The tables show the robust statistics; pairs and contrasts
  from the result use them too (robust t, and confidence limits from
  robust standard errors). Type II tests match
  `car::Anova(model, type = "II", white.adjust = "hc3")`. Not with
  `perms`. The default is `FALSE`.

## Value

The function produces permutation t-test table for coefficients of
(fixed) effects, permutation ANOVA table for model terms and (with
`method = "raw"`) a model parameter list `permlist`, a list containing
`nsim=4999` times permutation refitted `model` parameters which are used
in functions `predictmeans` and `contrastmeans`. With
`method = "freedman-lane"`, `permlist` is `NULL` and the result has a
fourth component, `perm`: what is needed to test any comparison with the
same permutations. Give the whole result to
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
or
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
(argument `permlist`) for Freedman-Lane p-values of pairwise comparisons
and contrasts. If refits fail, the number left out is in attribute
`"failures"`.

## Details

With `method = "auto"` (the default since predictmeans 2.0) the method
is chosen by model class from a simulation study
(`dev/permutation-plan.md` in the package sources): at present
Freedman-Lane for every class. Use `method = "raw"` for the results of
predictmeans 1.x.

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
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)),
or `method = "freedman-lane"` (4.6%).

Since predictmeans 2.0, `method = "freedman-lane"` also handles `lmer`,
`lme` and `gls` models, and `glm`, `glmer` and `glmmTMB` models. For
`lmer`, `lme` and `gls` models (with `refit = TRUE`, the default) the
model without the tested term is fitted, its residuals are decorrelated
with the symmetric square root of its covariance, permuted, correlated
again and added to its fitted values, and the model is refitted to each
permuted data set; the statistics are the model's Wald t and F
statistics (unadjusted covariance). One set of refits serves a term and
its coefficients; pairs and contrasts in
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
get one set of refits for their family (the model with all the compared
means equal). With `refit = FALSE` the covariance of the data, V, is
fixed at the fitted model: the data and the design are decorrelated with
the symmetric square root of V, so least squares on them gives the same
statistics, and Freedman-Lane permutations are done as for `lm`, without
refits. That is much faster, but the permuted statistics then behave as
if V were known, so comparisons whose error rests on few units rejected
too many true null hypotheses in a simulation study (whole-plot factor
of a split plot with 6 blocks: 8.6% at the 5% level, against 5.0% with
refits; pairs of whole-plot means 6.2%). GLMs and GLMMs are first
linearised at the fit (working response and working variances on the
link scale, as in
[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md));
GLMMs are then treated as the working linear mixed model (weights from
the working variances), refitted to each permuted working response when
`refit = TRUE`. Their t and F statistics are Wald statistics on the
working scale with the dispersion estimated from the working residuals,
so overdispersion affects them much less than likelihood-ratio tests (in
a pilot study with overdispersed counts: 8-10% against 35% at the 5%
level). For these classes the tables show these statistics and the
permutation p-values; type I tests are sequential in the GLS metric (as
[`anova()`](https://rdrr.io/r/stats/anova.html) of lme4 and nlme). The
result can be given to
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
and
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
as for `lm`. For `glmer` and `glmmTMB` models only
`method = "freedman-lane"` is available.

## References

Davidson, R. and Flachaire, E. (2008). The wild bootstrap, tamed at
last. *Journal of Econometrics*, 146, 162-169.

DiCiccio, C. J. and Romano, J. P. (2017). Robust permutation tests for
correlation and regression coefficients. *Journal of the American
Statistical Association*, 112, 1211-1220.

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
#>                            Estimate Std. Error t value Perm_p_value
#> (Intercept)                 80.0000     9.1070  8.7845           NA
#> nitro0.2                    18.5000     7.6830  2.4079         0.01
#> nitro0.4                    34.6667     7.6830  4.5122         0.01
#> nitro0.6                    44.8333     7.6830  5.8354         0.01
#> VarietyMarvellous            6.6667     9.7150  0.6862         0.48
#> VarietyVictory              -8.5000     9.7150 -0.8749         0.39
#> nitro0.2:VarietyMarvellous   3.3333    10.8653  0.3068         0.73
#> nitro0.4:VarietyMarvellous  -4.1667    10.8653 -0.3835         0.67
#> nitro0.6:VarietyMarvellous  -4.6667    10.8653 -0.4295         0.69
#> nitro0.2:VarietyVictory     -0.3333    10.8653 -0.0307         0.96
#> nitro0.4:VarietyVictory      4.6667    10.8653  0.4295         0.70
#> nitro0.6:VarietyVictory      2.1667    10.8653  0.1994         0.84
#> 
#> Note: Perm_p_value of t test is obtained using ‘99’ Freedman-Lane permutations.
#> 
#> ANOVA:
#> Type I F tests (Freedman-Lane permutations; the model refitted to each permuted data set)
#>               Df F value Perm_p_value
#> nitro          3 37.6856         0.01
#> Variety        2  1.4853         0.30
#> nitro:Variety  6  0.3028         0.96
#> 
#> Note: Perm_p_value of F test is obtained using ‘99’ Freedman-Lane permutations.
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
#> 9.715030 7.682956 9.160828 
#> attr(,"For the Same Level of Factor")
#>            nitro  Variety
#> Aveg.SED 9.71503 7.682956
#> Min.SED  9.71503 7.682956
#> Max.SED  9.71503 7.682956
#> 
#> $`Approximated LSD`
#>  Max.LSD  Min.LSD Aveg.LSD 
#> 19.43006 15.36591 18.32166 
#> attr(,"For the Same Level of Factor")
#>             nitro  Variety
#> Aveg.LSD 19.43006 15.36591
#> Min.LSD  19.43006 15.36591
#> Max.LSD  19.43006 15.36591
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
#> 0.4 0.015 0.048 1          C  
#> 0.6 0.015 0.015 0.15 1     C  
#> 
#> $Marvellous
#>     0    0.2   0.4  0.6 Group
#> 0   1                   A    
#> 0.2 0.02 1               B   
#> 0.4 0.02 0.23  1         BC  
#> 0.6 0.02 0.045 0.23 1     C  
#> 
#> $Victory
#>     0     0.2   0.4  0.6 Group
#> 0   1                    A    
#> 0.2 0.036 1               B   
#> 0.4 0.015 0.015 1          C  
#> 0.6 0.015 0.015 0.27 1     C  
#> 
#> $mean_table
#> $mean_table$Table
#>        Variety nitro     Mean      SE Df  LL(95%)  UL(95%) LetterGrp
#> 1  Golden Rain     0  80.0000 9.10696 NA  61.7861  98.2139       A  
#> 4  Golden Rain   0.2  98.5000 9.10696 NA  80.2861 116.7139        B 
#> 7  Golden Rain   0.4 114.6667 9.10696 NA  96.4527 132.8806         C
#> 10 Golden Rain   0.6 124.8333 9.10696 NA 106.6194 143.0473         C
#> 2   Marvellous     0  86.6667 9.10696 NA  68.4527 104.8806       A  
#> 5   Marvellous   0.2 108.5000 9.10696 NA  90.2861 126.7139        B 
#> 8   Marvellous   0.4 117.1667 9.10696 NA  98.9527 135.3806        BC
#> 11  Marvellous   0.6 126.8333 9.10696 NA 108.6194 145.0473         C
#> 3      Victory     0  71.5000 9.10696 NA  53.2861  89.7139       A  
#> 6      Victory   0.2  89.6667 9.10696 NA  71.4527 107.8806        B 
#> 9      Victory   0.4 110.8333 9.10696 NA  92.6194 129.0473         C
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
#>      Estimate Std. Error t value Permuted Pr(>|t|) lower.CL upper.CL
#> [1,]  -8.5000      9.715 -0.8749              0.37 -31.6879  14.6879
#> [2,] -45.6667      9.715 -4.7006              0.01 -68.1443 -23.1890
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

## The same by lmer, type II tests; refit = FALSE is much faster (no refits)
fm3 <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
permmodels(fm3, nperm = 99, type = "II", seed = 1, ncore = 1)
#> 
#> Coefficients of (fixed) effects:
#>                            Estimate Std. Error t value Perm_p_value
#> (Intercept)                 80.0000     9.1070  8.7844           NA
#> nitro0.2                    18.5000     7.6829  2.4079         0.02
#> nitro0.4                    34.6667     7.6829  4.5122         0.01
#> nitro0.6                    44.8333     7.6829  5.8354         0.01
#> VarietyMarvellous            6.6667     9.7150  0.6862         0.57
#> VarietyVictory              -8.5000     9.7150 -0.8749         0.43
#> nitro0.2:VarietyMarvellous   3.3333    10.8653  0.3068         0.78
#> nitro0.4:VarietyMarvellous  -4.1667    10.8653 -0.3835         0.71
#> nitro0.6:VarietyMarvellous  -4.6667    10.8653 -0.4295         0.71
#> nitro0.2:VarietyVictory     -0.3333    10.8653 -0.0307         0.99
#> nitro0.4:VarietyVictory      4.6667    10.8653  0.4295         0.66
#> nitro0.6:VarietyVictory      2.1667    10.8653  0.1994         0.92
#> 
#> Note: Perm_p_value of t test is obtained using ‘99’ Freedman-Lane permutations.
#> 
#> ANOVA:
#> Type II F tests (Freedman-Lane permutations; the model refitted to each permuted data set)
#>               Df F value Perm_p_value
#> nitro          3 37.6857         0.01
#> Variety        2  1.4853         0.24
#> nitro:Variety  6  0.3028         0.93
#> 
#> Note: Perm_p_value of F test is obtained using ‘99’ Freedman-Lane permutations.
#> 
permmodels(fm3, nperm = 999, type = "II", seed = 1, refit = FALSE)
#> 
#> Coefficients of (fixed) effects:
#>                            Estimate Std. Error t value Perm_p_value
#> (Intercept)                 80.0000     9.1070  8.7844           NA
#> nitro0.2                    18.5000     7.6829  2.4079        0.019
#> nitro0.4                    34.6667     7.6829  4.5122        0.001
#> nitro0.6                    44.8333     7.6829  5.8354        0.001
#> VarietyMarvellous            6.6667     9.7150  0.6862        0.494
#> VarietyVictory              -8.5000     9.7150 -0.8749        0.390
#> nitro0.2:VarietyMarvellous   3.3333    10.8653  0.3068        0.753
#> nitro0.4:VarietyMarvellous  -4.1667    10.8653 -0.3835        0.712
#> nitro0.6:VarietyMarvellous  -4.6667    10.8653 -0.4295        0.654
#> nitro0.2:VarietyVictory     -0.3333    10.8653 -0.0307        0.975
#> nitro0.4:VarietyVictory      4.6667    10.8653  0.4295        0.688
#> nitro0.6:VarietyVictory      2.1667    10.8653  0.1994        0.849
#> 
#> Note: Perm_p_value of t test is obtained using ‘999’ Freedman-Lane permutations.
#> 
#> ANOVA:
#> Type II F tests with the covariance of the data fixed at the fitted model
#>               Df F value Perm_p_value
#> nitro          3 37.6857        0.001
#> Variety        2  1.4853        0.230
#> nitro:Variety  6  0.3028        0.939
#> 
#> Note: Perm_p_value of F test is obtained using ‘999’ Freedman-Lane permutations.
#> 
# }

## Freedman-Lane permutations for a linear model with a covariate
mtcars$cyl <- factor(mtcars$cyl)
fm2 <- lm(mpg ~ wt + cyl, data = mtcars)
pfl <- permmodels(fm2, nperm = 999, type = "II", seed = 1, method = "freedman-lane")
#> 
#> Coefficients of (fixed) effects:
#>             Estimate Std. Error t value Pr(>|t|) Perm_p_value
#> (Intercept)  33.9908     1.8878 18.0056   0.0000           NA
#> wt           -3.2056     0.7539 -4.2521   0.0002        0.001
#> cyl6         -4.2556     1.3861 -3.0702   0.0047        0.009
#> cyl8         -6.0709     1.6523 -3.6742   0.0010        0.003
#> 
#> Note: Perm_p_value of t test is obtained using ‘999’ Freedman-Lane permutations.
#> 
#> ANOVA:
#> Anova Table (Type II tests)
#> 
#>  Response: mpg
#>             Sum Sq Df F value Pr(>F) Perm_p_value
#> wt        118.2039  1 18.0801 0.0002        0.001
#> cyl        95.2633  2  7.2856 0.0028        0.005
#> Residuals 183.0586 28                            
#> 
#> Note: Perm_p_value of F test is obtained using ‘999’ Freedman-Lane permutations.
#> 
## ... and for the pairwise comparisons of the cyl means
predictmeans(fm2, "cyl", permlist = pfl, plot = FALSE)$pairs
#> 
#>  The predicted means are estimated at 
#>     wt 
#> 3.2172 
#>   level1 level2 estimate std.error df statistic p.value p.adj  conf.low
#> 1      4      6 4.255582  1.386073 NA  3.070244   0.009 0.009  1.281825
#> 2      4      8 6.070860  1.652288 NA  3.674214   0.003 0.003  2.625644
#> 3      6      8 1.815277  1.357634 NA  1.337089   0.185 0.185 -1.036960
#>   conf.high      lsd    p.method
#> 1  7.229339 2.973757 permutation
#> 2  9.516075 3.445215 permutation
#> 3  4.667515 2.852238 permutation
## Robust (HC3) statistics, for error variances that differ between the cyl groups
permmodels(fm2, nperm = 999, type = "II", seed = 1, robust = TRUE)
#> 
#> Coefficients of (fixed) effects (robust standard errors):
#>             Estimate Robust SE t value Perm_p_value
#> (Intercept)  33.9908    2.1414 15.8734           NA
#> wt           -3.2056    0.7548 -4.2467        0.001
#> cyl6         -4.2556    1.2893 -3.3007        0.003
#> cyl8         -6.0709    1.6077 -3.7761        0.004
#> 
#> Note: Perm_p_value of robust t test is obtained using ‘999’ sign flips of the Freedman-Lane residuals (wild bootstrap).
#> 
#> ANOVA:
#> Type II robust (HC3) Wald F tests (Freedman-Lane residuals, sign flips)
#>     Df F value Perm_p_value
#> wt   1 18.0348        0.001
#> cyl  2  7.1814        0.007
#> 
#> Note: Perm_p_value of robust F test is obtained using ‘999’ sign flips of the Freedman-Lane residuals (wild bootstrap).
#> 
```
