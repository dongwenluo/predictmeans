# Permutation Tests of the Terms and Coefficients of a Model

Permutation t tests of the coefficients of the fixed effects and
permutation F (or chi-square) tests of the terms of a model fitted by
`aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer` or `glmmTMB`. p-values
are (b + 1)/(B + 1): b of the B permuted statistics at least as large as
the observed one. Results are the same for any `ncore`.

## Usage

``` r
permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
  test_statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, ncore=1, seed,
  method=c("auto", "raw", "freedman-lane"), perms=NULL, refit=TRUE, robust=FALSE, ...)
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`,
  `glmer` or `glmmTMB`.

- nperm:

  The number of permutations. The default is 4999.

- type:

  The type of the tests: "I", "II" or "III" (or 1, 2, 3). Type III tests
  are those of sum-to-zero contrasts whatever contrasts the model was
  fitted with (as SAS Type 3);
  [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) gives them
  only for models fitted with such contrasts (with R's default treatment
  contrasts it tests each main effect at the first level of the other
  factors).

- test_statistic:

  The statistic of Type II and III tests with `method = "raw"`: F for
  `lm`; chi-square for `lme` and `gls`; `"Chisq"` or `"F"` for `lmer`;
  `"LR"`, `"Wald"` or `"F"` for `glm`. Type I tests use F. A statistic
  the model does not have is replaced by the one it has, with a message;
  with `method = "freedman-lane"` the F (or t) statistics are permuted
  and a different `test_statistic` is ignored with a message. The F
  tests of an `lmer` model are Kenward-Roger F tests with the observed
  information of the variance parameters, as in the rest of predictmeans
  ([`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) uses the
  expected information; the two agree for balanced data).

- exact:

  `TRUE`: `perms` holds every permutation the design allows (a complete
  enumeration), so p-values are the proportion of them with a statistic
  at least as large as the observed one. Otherwise (the default `FALSE`)
  the observed data count as one of the permutations: p = (b + 1)/(B +
  1), never 0. `exact = TRUE` without `perms`, or with repeated columns,
  is ignored with a warning.

- data:

  The data the model was fitted to, when they cannot be found from the
  model's call (e.g. after transforming variables outside the data).

- ...:

  Not used: `prt` and `fo` (before 2.0) stop with their 2.0 spelling.
  Nothing is printed during the call:
  [`print()`](https://rdrr.io/r/base/print.html) of the result shows the
  tables. A model whose formula is held in a variable (`lm(f, data)`;
  1.x's `fo`) is detected and refitted in this session only.

- ncore:

  Number of cores (default 1; before predictmeans 2.0: 3). Each
  permutation has its own seed, so the results do not depend on `ncore`.

- seed:

  A random number seed, for reproducible results; the user's random
  number stream is left as it was.

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

A list of class `"pm_permmodels"`: `COEFFICENTS` (permutation t tests of
the coefficients), `ANOVA` (permutation tests of the terms, of the
`type` asked for) and `permlist` (with `method = "raw"`: the
coefficients and covariance of each refit, for
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
and
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)).
Values are stored rounded to 4 decimals; the notes give the smallest
attainable p-value, 1/(B + 1). Advice on the number of permutations is
kept in `attr(, "advice")`
([`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).
With `method = "freedman-lane"`, `permlist` is `NULL` and the result has
a fourth component, `perm`: what is needed to test any comparison with
the same permutations. Give the whole result to
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
(4.6%); prefer the Kenward-Roger or Satterthwaite F test there
([`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md)
or
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

## See also

[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md),
[`permindex`](https://dongwenluo.github.io/predictmeans/reference/permindex.md),
[`permstats`](https://dongwenluo.github.io/predictmeans/reference/permstats.md),
[`permcombine`](https://dongwenluo.github.io/predictmeans/reference/permcombine.md),
[`bootmodels`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

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

## Permutation test for multiple comparisons
predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH",
  permlist=permlme)
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
#> 9.715029 7.682957 9.160827 
#> attr(,"For the Same Level of Factor")
#>             nitro  Variety
#> Aveg.SED 9.715029 7.682957
#> Min.SED  9.715029 7.682957
#> Max.SED  9.715029 7.682957
#> 
#> $`Approximated LSD`
#>  Max.LSD  Min.LSD Aveg.LSD 
#> 19.43006 15.36591 18.32165 
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
#> 0.6 0.015 0.015 0.13 1     C  
#> 
#> $Marvellous
#>     0    0.2  0.4  0.6 Group
#> 0   1                  A    
#> 0.2 0.06 1             AB   
#> 0.4 0.03 0.23 1         B   
#> 0.6 0.03 0.06 0.23 1    B   
#> 
#> $Victory
#>     0     0.2   0.4  0.6 Group
#> 0   1                    A    
#> 0.2 0.015 1               B   
#> 0.4 0.015 0.024 1          C  
#> 0.6 0.015 0.015 0.33 1     C  
#> 
#> $mean_table
#> $mean_table$Table
#>        Variety nitro     Mean      SE Df  LL(95%)  UL(95%) LetterGrp
#> 1  Golden Rain     0  80.0000 9.10696 NA  60.7021  99.2979       A  
#> 4  Golden Rain   0.2  98.5000 9.10696 NA  79.2021 117.7979        B 
#> 7  Golden Rain   0.4 114.6667 9.10696 NA  95.3688 133.9645         C
#> 10 Golden Rain   0.6 124.8333 9.10696 NA 105.5355 144.1312         C
#> 2   Marvellous     0  86.6667 9.10696 NA  67.3688 105.9645        A 
#> 5   Marvellous   0.2 108.5000 9.10696 NA  89.2021 127.7979        AB
#> 8   Marvellous   0.4 117.1667 9.10696 NA  97.8688 136.4645         B
#> 11  Marvellous   0.6 126.8333 9.10696 NA 107.5355 146.1312         B
#> 3      Victory     0  71.5000 9.10696 NA  52.2021  90.7979       A  
#> 6      Victory   0.2  89.6667 9.10696 NA  70.3688 108.9645        B 
#> 9      Victory   0.4 110.8333 9.10696 NA  91.5355 130.1312         C
#> 12     Victory   0.6 118.5000 9.10696 NA  99.2021 137.7979         C
#> 
#> $mean_table$Note
#> [1] "Letter-based representation of pairwise comparisons at significant level '0.05' at each level of Variety"
#> 
#> 

## Permutation test for specified contrasts
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
#> Contrasts of 'nitro:Variety' (95% limits; p.adj: none; 99 permutations)
#>    contrast estimate std.error df statistic p.value p.adj conf.low conf.high
#>  contrast 1    -8.50     9.715 NA   -0.8749    0.38  0.38   -28.04     11.04
#>  contrast 2   -45.67     9.715 NA   -4.7006    0.01  0.01   -65.22    -26.12

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
#> Note: Perm_p_value of t test is obtained using ‘99’ Freedman-Lane permutations. The smallest attainable p-value is 1/(99 + 1) = 0.01.
#> 
#> ANOVA:
#> Type II F tests (Freedman-Lane permutations; the model refitted to each permuted data set)
#>               Df F value Perm_p_value
#> nitro          3 37.6857         0.01
#> Variety        2  1.4853         0.24
#> nitro:Variety  6  0.3028         0.93
#> 
#> Note: Perm_p_value of F test is obtained using ‘99’ Freedman-Lane permutations. The smallest attainable p-value is 1/(99 + 1) = 0.01.
#> 
#> Advice:
#> ! Only 99 resamples: p-values at the floor (0.01) only say 'as small as can
#>   be shown' (Monte Carlo SE about 0.022 at p = 0.05). Use nperm = 999 or
#>   more (4999 for reporting); 'seed' makes it reproducible.
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
#> Note: Perm_p_value of t test is obtained using ‘999’ Freedman-Lane permutations. The smallest attainable p-value is 1/(999 + 1) = 0.001.
#> 
#> ANOVA:
#> Type II F tests with the covariance of the data fixed at the fitted model
#>               Df F value Perm_p_value
#> nitro          3 37.6857        0.001
#> Variety        2  1.4853        0.230
#> nitro:Variety  6  0.3028        0.939
#> 
#> Note: Perm_p_value of F test is obtained using ‘999’ Freedman-Lane permutations. The smallest attainable p-value is 1/(999 + 1) = 0.001.
#> 
# }

## Freedman-Lane permutations for a linear model with a covariate
mtcars$cyl <- factor(mtcars$cyl)
fm2 <- lm(mpg ~ wt + cyl, data = mtcars)
pfl <- permmodels(fm2, nperm = 999, type = "II", seed = 1, method = "freedman-lane")
## ... and for the pairwise comparisons of the cyl means
predictmeans(fm2, "cyl", permlist = pfl)$pairs
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
#> Note: Perm_p_value of robust t test is obtained using ‘999’ sign flips of the Freedman-Lane residuals (wild bootstrap). The smallest attainable p-value is 1/(999 + 1) = 0.001.
#> 
#> ANOVA:
#> Type II robust (HC3) Wald F tests (Freedman-Lane residuals, sign flips)
#>     Df F value Perm_p_value
#> wt   1 18.0348        0.001
#> cyl  2  7.1814        0.007
#> 
#> Note: Perm_p_value of robust F test is obtained using ‘999’ sign flips of the Freedman-Lane residuals (wild bootstrap). The smallest attainable p-value is 1/(999 + 1) = 0.001.
#> 
```
