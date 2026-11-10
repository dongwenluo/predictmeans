# Permutation and bootstrap tests

F and t tests rest on the model’s distributional assumptions and, for
mixed models, on approximate degrees of freedom. With few experimental
units, skewed residuals or a GLMM, resampling tests are a useful check.
predictmeans has:

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md):
  permutation tests of the terms and coefficients of `lm`, `glm`, `aov`,
  `gls`, `lme` and `lmer` models;
- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  and
  [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/permanova.lmer.md):
  permutation likelihood-ratio tests for `lmer` models, including tests
  of random effects;
- [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md):
  parametric bootstrap tests for `lmer` and `glmer` models.

The results of
[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
and
[`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
can be passed to
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
and
[`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
(argument `permlist`), which then give resampling p-values for pairwise
comparisons and contrasts.

The examples below use 99 or 199 resamples so that they run quickly; use
999 or more in practice. The smallest p-value possible is 1/(number of
resamples + 1). `ncore` sets the number of cores used for the refits
(`ncore = 1` here, as for CRAN checks); `seed` makes the results
reproducible.

``` r

library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
```

## Permutation tests of terms and coefficients

For linear models the response is permuted. For mixed models the
residuals are first decorrelated with the Cholesky factor of their
estimated covariance, permuted, and correlated again (Lee and Braun
2012), so each permuted data set keeps the model’s correlation
structure; the model is then refitted to each.

``` r

perm <- permmodels(fm, nperm = 199, ncore = 1, seed = 1)
#> 
#> Coefficients of (fixed) effects:
#>                            Estimate Std. Error      df t value Pr(>|t|)
#> (Intercept)                 80.0000     9.1070 16.0816  8.7844   0.0000
#> nitro0.2                    18.5000     7.6829 45.0000  2.4079   0.0202
#> nitro0.4                    34.6667     7.6829 45.0000  4.5122   0.0000
#> nitro0.6                    44.8333     7.6829 45.0000  5.8354   0.0000
#> VarietyMarvellous            6.6667     9.7150 30.2308  0.6862   0.4978
#> VarietyVictory              -8.5000     9.7150 30.2308 -0.8749   0.3885
#> nitro0.2:VarietyMarvellous   3.3333    10.8653 45.0000  0.3068   0.7604
#> nitro0.4:VarietyMarvellous  -4.1667    10.8653 45.0000 -0.3835   0.7032
#> nitro0.6:VarietyMarvellous  -4.6667    10.8653 45.0000 -0.4295   0.6696
#> nitro0.2:VarietyVictory     -0.3333    10.8653 45.0000 -0.0307   0.9757
#> nitro0.4:VarietyVictory      4.6667    10.8653 45.0000  0.4295   0.6696
#> nitro0.6:VarietyVictory      2.1667    10.8653 45.0000  0.1994   0.8428
#>                            Perm_p_value
#> (Intercept)                          NA
#> nitro0.2                          0.020
#> nitro0.4                          0.005
#> nitro0.6                          0.005
#> VarietyMarvellous                 0.495
#> VarietyVictory                    0.350
#> nitro0.2:VarietyMarvellous        0.775
#> nitro0.4:VarietyMarvellous        0.640
#> nitro0.6:VarietyMarvellous        0.680
#> nitro0.2:VarietyVictory           0.965
#> nitro0.4:VarietyVictory           0.660
#> nitro0.6:VarietyVictory           0.885
#> 
#> Note: Perm_p_value of t test is obtained using '199' permutations.
#> 
#> ANOVA:
#>                   Sum Sq   Mean Sq NumDF DenDF F value Pr(>F) Perm_p_value
#> nitro         20020.5001 6673.5000     3    45 37.6857 0.0000        0.005
#> Variety         526.0575  263.0288     2    10  1.4853 0.2724        0.200
#> nitro:Variety   321.7500   53.6250     6    45  0.3028 0.9322        0.930
#> 
#> Note: Perm_p_value of F test is obtained using '199' permutations.
```

The tables repeat the model’s tests with a `Perm_p_value` column: the
proportion of permuted statistics at least as large as the observed one.

## Pairwise comparisons and contrasts

The same refits give permutation p-values for pairwise comparisons, in
the `pairs` table and the legacy p-value matrix:

``` r

pm <- predictmeans(fm, "nitro", permlist = perm, plot = FALSE, prtnum = FALSE, ncore = 1)
pm$pairs[, c("level1", "level2", "estimate", "statistic", "p.value", "p.method")]
#>   level1 level2   estimate statistic p.value    p.method
#> 1      0    0.2 -19.500000 -4.396098   0.005 permutation
#> 2      0    0.4 -34.833333 -7.852859   0.005 permutation
#> 3    0.2    0.4 -15.333333 -3.456761   0.005 permutation
#> 4      0    0.6 -44.000000 -9.919400   0.005 permutation
#> 5    0.2    0.6 -24.500000 -5.523302   0.005 permutation
#> 6    0.4    0.6  -9.166667 -2.066542   0.050 permutation
```

and for contrasts:

``` r

cm <- contrast_family(Kmatrix(fm, "nitro"), "trt.vs.ctrl")
contrastmeans(fm, "nitro", ctrmatrix = cm, permlist = perm)[[1]]
#>         Estimate Std. Error t value Permuted Pr(>|t|)
#> 0.2 - 0  19.5000     4.4358  4.3961             0.005
#> 0.4 - 0  34.8333     4.4358  7.8529             0.005
#> 0.6 - 0  44.0000     4.4358  9.9194             0.005
#> attr(,"Note")
#> [1] "The permuted p-value is obtained using '199' permutations."
```

`adj` adjusts the permutation p-values with any
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method or
`"sidak"`. `"tukey"`, `"scheffe"` and `"mvt"` need the t distribution
and are not available with `permlist`.

## Testing random effects: `permlmer()`

[`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
compares two nested `lmer` models with a permutation likelihood-ratio
test (Lee and Braun 2012). This is useful for random effects, where the
usual chi-square reference distribution is wrong because the null value
of a variance lies on the boundary. Does the slope for `Days` vary
between subjects?

``` r

fm1 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
fm0 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
permlmer(fm0, fm1, nperm = 99, ncore = 1, seed = 1)
#> Data: sleepstudy
#> Models:
#> lmer0: Reaction ~ Days + (1 | Subject)
#> lmer1: Reaction ~ Days + (Days | Subject)
#>       npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> lmer0    4 1794.5 1807.2 -893.23    1786.5                            
#> lmer1    6 1755.6 1774.8 -871.81    1743.6 42.837  2   4.99e-10   0.01
```

[`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/permanova.lmer.md)
applies the same test to each fixed term of an `lmer` model in turn
(Type I, II or III), for an ANOVA table with permutation p-values.

## Parametric bootstrap: `bootmodels()`

[`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
simulates new responses from the model fitted without its fixed effects
(an intercept and the same random effects: the global null hypothesis),
refits the full model to each and compares the statistics. Unlike
permutation, it works for GLMMs. For the clinical trial data (a binomial
GLMM):

``` r

data(Clinical)
fm_clin <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic),
                 family = binomial, data = Clinical)
boot <- bootmodels(fm_clin, nsim = 99, seed = 1)
#> 
#> Coefficients of (fixed) effects:
#>               Estimate Std. Error z value Pr(>|z|) Boot_p_value
#> (Intercept)    -1.1965     0.5529 -2.1641   0.0305           NA
#> Treatmentdrug   0.7382     0.3004  2.4575   0.0140         0.03
#> 
#> ANOVA:
#>           npar Sum Sq Mean Sq F value Boot_p_value
#> Treatment    1 6.2054  6.2054  6.2054         0.03
#> 
#> Note: Boot_p_value is from '99' parametric bootstrap samples under the global null model.
```

Its result is used as a `permlist`; the `pairs` table then says
`"bootstrap"`:

``` r

pmb <- predictmeans(fm_clin, "Treatment", permlist = boot, plot = FALSE, prtnum = FALSE,
                    ncore = 1)
pmb$pairs[, c("level1", "level2", "estimate", "p.value", "p.method")]
#>    level1 level2   estimate p.value  p.method
#> 1 control   drug -0.7382041    0.03 bootstrap
```

## Which to use

- **Terms and pairwise comparisons in `lm`, `lme` and `lmer` models**:
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md).
- **A variance component, or any comparison of two nested `lmer`
  models**:
  [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md).
- **GLMMs**:
  [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md).
- With large samples and well-behaved residuals the F and t tests (with
  Kenward–Roger df for `lmer`) agree closely with these, at a fraction
  of the time.

## References

Lee, O. E. and Braun, T. M. (2012). Permutation tests for random effects
in linear mixed models. *Biometrics* 68, 486–493.
