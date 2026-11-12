# Permutation and bootstrap tests

F and t tests rest on the model’s distributional assumptions and, for
mixed models, on approximate degrees of freedom. With few experimental
units, skewed residuals or a GLMM, resampling tests are a useful check.
predictmeans has:

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md):
  permutation tests of the terms and coefficients of `lm`, `aov`, `glm`,
  `gls`, `lme`, `lmer`, `glmer` and `glmmTMB` models;
- [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  (also as
  [`permglmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
  and
  [`permglmmTMB()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md))
  and
  [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/permanova.lmer.md):
  permutation likelihood-ratio tests comparing two nested mixed models,
  including tests of random effects;
- [`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md):
  permutations restricted by the design (within blocks, of whole plots),
  for randomisation tests;
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

The examples below use 99 to 999 resamples so that they run quickly; use
999 or more in practice. The smallest p-value possible is 1/(number of
resamples + 1). `ncore` sets the number of cores used for refits
(`ncore = 1` here, as for CRAN checks); `seed` makes the results
reproducible without changing your own random number stream.

``` r

library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
```

## How the data are permuted

[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
has two schemes.

- **Freedman–Lane** (`method = "freedman-lane"`; the default, through
  `method = "auto"`). For each test, the model without the tested term
  (for a pair of means: with the two means equal) gives fitted values
  and residuals; the residuals are permuted and added back to those
  fitted values, and the statistic is computed again (Freedman and Lane
  1983). The other terms keep their effects, so each test is of its own
  term. For linear mixed models the residuals are decorrelated with the
  symmetric square root of their covariance before they are permuted and
  correlated again after, and the model is refitted to each permuted
  data set, so its variance parameters are estimated again each time
  (one set of refits per term). `refit = FALSE` instead fixes the
  covariance of the data at the fit, which turns the model into a linear
  model with no refits: much faster, but slightly liberal when the
  variance parameters rest on few units (whole plots with few blocks:
  8.6% instead of 5% in a simulation study, against 5.0% with refits).
  GLMs and GLMMs are first linearised at the fit (working response and
  working variances on the link scale); a GLMM is then treated as its
  working linear mixed model and refitted in the same way.
- **Raw** (`method = "raw"`, the scheme of predictmeans 1.x). The
  response (`lm`, `glm`) or the fitted values together with the
  decorrelated residuals (mixed models) are permuted, and the model is
  refitted to each permuted data set. The refits are kept in `permlist`.
  This scrambles the effects of all the terms at once, so it tests the
  global null hypothesis; with a covariate that has extreme values and
  is related to the tested term it rejected true null hypotheses too
  often (6–7% at the 5% level, against 4–5% for Freedman–Lane, in a
  simulation study; Anderson and Legendre 1999 found the same for linear
  models).

``` r

perm <- permmodels(fm, nperm = 199, seed = 1, ncore = 1)
#> 
#> Coefficients of (fixed) effects:
#>                            Estimate Std. Error t value Perm_p_value
#> (Intercept)                 80.0000     9.1070  8.7844           NA
#> nitro0.2                    18.5000     7.6829  2.4079        0.025
#> nitro0.4                    34.6667     7.6829  4.5122        0.005
#> nitro0.6                    44.8333     7.6829  5.8354        0.005
#> VarietyMarvellous            6.6667     9.7150  0.6862        0.545
#> VarietyVictory              -8.5000     9.7150 -0.8749        0.410
#> nitro0.2:VarietyMarvellous   3.3333    10.8653  0.3068        0.765
#> nitro0.4:VarietyMarvellous  -4.1667    10.8653 -0.3835        0.715
#> nitro0.6:VarietyMarvellous  -4.6667    10.8653 -0.4295        0.645
#> nitro0.2:VarietyVictory     -0.3333    10.8653 -0.0307        0.990
#> nitro0.4:VarietyVictory      4.6667    10.8653  0.4295        0.670
#> nitro0.6:VarietyVictory      2.1667    10.8653  0.1994        0.885
#> 
#> Note: Perm_p_value of t test is obtained using '199' Freedman-Lane permutations.
#> 
#> ANOVA:
#> Type I F tests (Freedman-Lane permutations; the model refitted to each permuted data set)
#>               Df F value Perm_p_value
#> nitro          3 37.6857        0.005
#> Variety        2  1.4853        0.260
#> nitro:Variety  6  0.3028        0.925
#> 
#> Note: Perm_p_value of F test is obtained using '199' Freedman-Lane permutations.
```

For `lm`, `aov`, `lme` and `gls` models the tables are the model’s own
tests with a `Perm_p_value` column; for the other classes they show the
statistics that were permuted (Wald t and F with the covariance of the
data fixed at the fit).

``` r

fm_cars <- lm(mpg ~ wt + factor(cyl), data = mtcars)
permmodels(fm_cars, nperm = 999, type = "II", seed = 1)
#> 
#> Coefficients of (fixed) effects:
#>              Estimate Std. Error t value Pr(>|t|) Perm_p_value
#> (Intercept)   33.9908     1.8878 18.0056   0.0000           NA
#> wt            -3.2056     0.7539 -4.2521   0.0002        0.001
#> factor(cyl)6  -4.2556     1.3861 -3.0702   0.0047        0.009
#> factor(cyl)8  -6.0709     1.6523 -3.6742   0.0010        0.003
#> 
#> Note: Perm_p_value of t test is obtained using '999' Freedman-Lane permutations.
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
#> Note: Perm_p_value of F test is obtained using '999' Freedman-Lane permutations.
```

## Pairwise comparisons and contrasts

The same permutations give p-values for pairwise comparisons, in the
`pairs` table and the legacy p-value matrix. With Freedman–Lane each
pair is tested with its own reduced model, in which that pair’s means
are equal, so its p-value does not depend on the other means:

``` r

pm <- predictmeans(fm, "nitro", permlist = perm, plot = FALSE, prtnum = FALSE, ncore = 1)
pm$pairs[, c("level1", "level2", "estimate", "statistic", "p.value", "p.method")]
#>   level1 level2   estimate statistic p.value    p.method
#> 1      0    0.2 -19.500000 -4.396098   0.005 permutation
#> 2      0    0.4 -34.833333 -7.852859   0.005 permutation
#> 3    0.2    0.4 -15.333333 -3.456761   0.005 permutation
#> 4      0    0.6 -44.000000 -9.919400   0.005 permutation
#> 5    0.2    0.6 -24.500000 -5.523302   0.005 permutation
#> 6    0.4    0.6  -9.166667 -2.066542   0.055 permutation
```

and for contrasts:

``` r

cm <- contrast_family(Kmatrix(fm, "nitro"), "trt.vs.ctrl")
contrastmeans(fm, "nitro", ctrmatrix = cm, permlist = perm)[[1]]
#>         Estimate Std. Error t value Permuted Pr(>|t|) lower.CL upper.CL
#> 0.2 - 0  19.5000     4.4358  4.3961             0.005  11.2629  27.7371
#> 0.4 - 0  34.8333     4.4358  7.8529             0.005  26.8245  42.8422
#> 0.6 - 0  44.0000     4.4358  9.9194             0.005  35.5097  52.4903
#> attr(,"Note")
#> [1] "The permuted p-value is obtained using '199' permutations."
```

With `method = "raw"` all pairs use the same refits. `adj` adjusts the
permutation p-values with any
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method or
`"sidak"`, and with the max-\|t\| method of Westfall and Young (1993),
which uses the correlation between the comparisons: `adj = "tukey"` (or
`"mvt"`) compares each pair’s \|t\| with the largest \|t\| of the family
in each permuted data set (a permutation Tukey test), and
`adj = "westfall"` does the same step by step, which is more powerful.
In a simulation study this kept the familywise error rate as Tukey’s
test does.

``` r

pmw <- predictmeans(fm, "nitro", permlist = perm, adj = "westfall", plot = FALSE,
                    prtnum = FALSE, ncore = 1)
pmw$pairs[, c("level1", "level2", "p.value", "p.adj")]
#>   level1 level2 p.value p.adj
#> 1      0    0.2   0.005 0.005
#> 2      0    0.4   0.005 0.005
#> 3    0.2    0.4   0.005 0.010
#> 4      0    0.6   0.005 0.005
#> 5    0.2    0.6   0.005 0.005
#> 6    0.4    0.6   0.055 0.055
```

With a `permlist`, the observed t statistics use the unadjusted
covariance (no Kenward–Roger adjustment), as the permuted ones do.

## GLMs and GLMMs

For GLMs and GLMMs
[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
works on the working (link) scale, with the dispersion estimated from
the working residuals, which makes it much less sensitive to
overdispersion than the likelihood-ratio test. For the clinical trial
data (a binomial GLMM):

``` r

data(Clinical)
fm_clin <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic),
                 family = binomial, data = Clinical)
perm_clin <- permmodels(fm_clin, nperm = 199, seed = 1, ncore = 1)
#> 
#> Coefficients of (fixed) effects:
#>               Estimate Std. Error t value Perm_p_value
#> (Intercept)    -1.1399     0.5591 -2.0389           NA
#> Treatmentdrug   0.7207     0.3070  2.3474        0.075
#> 
#> Note: Perm_p_value of t test is obtained using '199' Freedman-Lane permutations.
#> 
#> ANOVA:
#> Type I F tests on the working (link) scale (Freedman-Lane permutations; the working linear mixed model refitted to each permuted data set)
#>           Df F value Perm_p_value
#> Treatment  1  5.5102        0.075
#> 
#> Note: Perm_p_value of F test is obtained using '199' Freedman-Lane permutations.
```

## Testing random effects: `permlmer()`

[`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
compares two nested mixed models with a permutation likelihood-ratio
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

For GLMMs (`glmer` or `glmmTMB`, also as
[`permglmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
and
[`permglmmTMB()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md))
the same is done on the working scale. For 0/1 responses, where the
working scale fails with rare events or tiny clusters, randomised
quantile residuals are permuted instead and the GLMMs refitted
(`method = "quantile"`, chosen automatically). Do tick counts vary
between locations, beyond the variation between broods?

``` r

gm1 <- glmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD), family = poisson,
             data = lme4::grouseticks)
gm2 <- update(gm1, . ~ . + (1 | LOCATION))
permglmer(gm1, gm2, nperm = 99, ncore = 1, seed = 1)
#> Data: lme4::grouseticks
#> Models:
#> m0: TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD)
#> m1: TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION)
#>    npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> m0    5 1988.1 2008.1 -989.04    1978.1                            
#> m1    6 1987.9 2011.9 -987.94    1975.9 2.1992  1    0.13809   0.08
```

[`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/permanova.lmer.md)
applies the same test to each fixed term of an `lmer` or `glmer` model
in turn (Type I, II or III), for an ANOVA table with permutation
p-values.

## Permutations restricted by the design: `permindex()`

In a designed experiment the randomisation defines which permutations
are possible: treatments are randomised within blocks, whole-plot
treatments to whole plots.
[`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)
gives such permutations, and
[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
and
[`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
use them through `perms` (a randomisation test):

``` r

fm_rcb <- lm(yield ~ Block + nitro, data = Oats)
idx <- permindex(Oats, block = "Block", nsim = 999, seed = 1)
permmodels(fm_rcb, perms = idx)$ANOVA
#> 
#> Coefficients of (fixed) effects:
#>             Estimate Std. Error t value Pr(>|t|) Perm_p_value
#> (Intercept)  79.3889     3.7668 21.0759   0.0000           NA
#> Block.L      29.4823     4.6134  6.3906   0.0000        0.069
#> Block.Q      20.0124     4.6134  4.3379   0.0001        1.000
#> Block.C       5.3728     4.6134  1.1646   0.2486        0.001
#> Block^4       4.7718     4.6134  1.0343   0.3049        1.000
#> Block^5      -1.2651     4.6134 -0.2742   0.7848        0.001
#> nitro0.2     19.5000     5.3271  3.6605   0.0005        0.002
#> nitro0.4     34.8333     5.3271  6.5389   0.0000        0.001
#> nitro0.6     44.0000     5.3271  8.2597   0.0000        0.001
#> 
#> Note: Perm_p_value of t test is obtained using '999' Freedman-Lane permutations.
#> 
#> ANOVA:
#> Analysis of Variance Table
#> 
#>  Response: yield
#>           Df   Sum Sq   Mean Sq F value Pr(>F) Perm_p_value
#> Block      5 15875.28 3175.0556 12.4317      0            1
#> nitro      3 20020.50 6673.5000 26.1297      0        0.001
#> Residuals 63 16090.17  255.3995                            
#> 
#> Note: Perm_p_value of F test is obtained using '999' Freedman-Lane permutations.
#>           Df   Sum Sq   Mean Sq F value Pr(>F) Perm_p_value
#> Block      5 15875.28 3175.0556 12.4317      0            1
#> nitro      3 20020.50 6673.5000 26.1297      0        0.001
#> Residuals 63 16090.17  255.3995
```

## Parametric bootstrap: `bootmodels()`

[`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
simulates new responses from the model fitted without its fixed effects
(an intercept and the same random effects: the global null hypothesis),
refits the full model to each and compares the statistics. For the
clinical trial data:

``` r

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

## Unequal variances: `robust = TRUE`

No permutation scheme repairs a variance model that is wrong: when
groups have very different error variances and the model assumes them
equal, the F test and the permutation tests above all reject true null
hypotheses far too often (about 20% at the 5% level in our simulations
with group SDs 3, 1, 1). For `lm`, `aov` and `glm` models,
`robust = TRUE` uses heteroscedasticity-robust (HC3) t and F statistics
and flips the signs of the reduced model’s residuals at random instead
of permuting them (a wild bootstrap), so each residual keeps its own
variance:

``` r

mt <- transform(mtcars, cyl = factor(cyl))
fm_mt <- lm(mpg ~ wt + cyl, data = mt)
permmodels(fm_mt, nperm = 999, type = "II", seed = 1, robust = TRUE)
#> 
#> Coefficients of (fixed) effects (robust standard errors):
#>             Estimate Robust SE t value Perm_p_value
#> (Intercept)  33.9908    2.1414 15.8734           NA
#> wt           -3.2056    0.7548 -4.2467        0.001
#> cyl6         -4.2556    1.2893 -3.3007        0.003
#> cyl8         -6.0709    1.6077 -3.7761        0.004
#> 
#> Note: Perm_p_value of robust t test is obtained using '999' sign flips of the Freedman-Lane residuals (wild bootstrap).
#> 
#> ANOVA:
#> Type II robust (HC3) Wald F tests (Freedman-Lane residuals, sign flips)
#>     Df F value Perm_p_value
#> wt   1 18.0348        0.001
#> cyl  2  7.1814        0.007
#> 
#> Note: Perm_p_value of robust F test is obtained using '999' sign flips of the Freedman-Lane residuals (wild bootstrap).
```

Pairs and contrasts from this result use the robust statistics and
robust standard errors. For mixed models, model the variances instead
(e.g. `lme(..., weights = varIdent())`) and use
[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
on that model.

## More from the permutations: `permstats()` and `permcombine()`

[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
keeps the permuted statistics.
[`permstats()`](https://dongwenluo.github.io/predictmeans/reference/permstats.md)
turns them into a standardised effect size for each term and
coefficient, and checks the parametric test against them: `param.size`
is the proportion of permuted data sets (in which the term has no
effect) where the parametric test rejects at `level`. Near `level`, the
parametric test suits these data; well above it, its p-values are too
small here.

``` r

p_mpg <- permmodels(fm_mt, nperm = 999, type = "II", seed = 1, prt = FALSE)
permstats(p_mpg)
#>           statistic df1 df2      p.param p.perm       ses param.size
#> wt        18.080056   1  28 0.0002130435  0.001 11.695665 0.04904905
#> cyl        7.285567   2  28 0.0028353022  0.005  4.622323 0.05705706
#> coef:wt    4.252065   1  28 0.0002130435  0.001  5.599425 0.04904905
#> coef:cyl6  3.070244   1  28 0.0047178339  0.009  3.291703 0.05905906
#> coef:cyl8  3.674214   1  28 0.0009991893  0.003  4.186066 0.05405405
```

When the same model is fitted to several responses measured on the same
units (traits, time points), results made with the same permutations
(same rows, `seed` and `nperm`) can be combined.
[`permcombine()`](https://dongwenluo.github.io/predictmeans/reference/permcombine.md)
gives, for each term, one p-value for “no effect on any response”
(Fisher’s or Tippett’s combination; Pesarin and Salmaso 2010) and
p-values adjusted for testing every response (min-p; Westfall and Young
1993). The responses’ correlation is kept by the shared permutations,
without a multivariate model.

``` r

p_qsec <- permmodels(lm(qsec ~ wt + cyl, data = mt), nperm = 999, type = "II", seed = 1,
                     prt = FALSE)
permcombine(mpg = p_mpg, qsec = p_qsec)
#>             mpg  qsec combined adj.mpg adj.qsec
#> wt        0.001 0.002    0.001   0.002    0.003
#> cyl       0.005 0.002    0.002   0.008    0.004
#> coef:wt   0.001 0.002    0.001   0.002    0.003
#> coef:cyl6 0.009 0.002    0.002   0.016    0.003
#> coef:cyl8 0.003 0.002    0.003   0.004    0.003
```

For many comparisons within one response, such as a treatment difference
at each of several time points,
`contrastmeans(..., permlist = , adj = "westfall")` gives the
max-statistic (step-down max-\|t\|) adjustment over the time points.

## Which to use

| Model | Terms, coefficients, pairs, contrasts | Random effects, nested models |
|----|----|----|
| `lm`, `aov` | [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) (Freedman–Lane) | — |
| `glm` | [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) (Freedman–Lane on the working scale) | — |
| `lmer`, `lme`, `gls` | [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) (Freedman–Lane with refits); or the Kenward–Roger F test | [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md) |
| `glmer`, `glmmTMB` | [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) (Freedman–Lane on the working scale, with refits) | [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md): working scale; quantile residuals for 0/1 responses |
| designed experiments | `permmodels(perms = permindex(...))` for a randomisation test |  |
| `lm`, `glm` with unequal variances | `permmodels(robust = TRUE)` (HC3 statistics, sign flips) |  |
| several responses, same units | [`permcombine()`](https://dongwenluo.github.io/predictmeans/reference/permcombine.md) of [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) results with the same `seed` |  |

- With large samples and well-behaved residuals the F and t tests (with
  Kenward–Roger df for `lmer`) agree closely with these, at a fraction
  of the time.
- [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  simulates from the global null model; it is useful as a check for
  GLMMs and very small samples.
- The choices come from simulation studies of the type I error and power
  of each scheme (in the package sources: `dev/permutation-plan.md` and
  `dev/sim-study/`).

## References

Davidson, R. and Flachaire, E. (2008). The wild bootstrap, tamed at
last. *Journal of Econometrics* 146, 162–169.

DiCiccio, C. J. and Romano, J. P. (2017). Robust permutation tests for
correlation and regression coefficients. *Journal of the American
Statistical Association* 112, 1211–1220.

Pesarin, F. and Salmaso, L. (2010). *Permutation Tests for Complex Data:
Theory, Applications and Software*. Wiley.

Anderson, M. J. and Legendre, P. (1999). An empirical comparison of
permutation methods for tests of partial regression coefficients in a
linear model. *Journal of Statistical Computation and Simulation* 62,
271–303.

Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of
reported significance levels. *Journal of Business & Economic
Statistics* 1, 292–298.

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics* 5, 236–244.

Lee, O. E. and Braun, T. M. (2012). Permutation tests for random effects
in linear mixed models. *Biometrics* 68, 486–493.

Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.
