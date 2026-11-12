# Linear Contrast Tests for a Linear Model

Performs t-tests (or permuted t-tests) of specified contrasts for linear
models obtained from functions `aov`, `lm`, `glm`, `gls`, `lme`, or
`lmer`.

## Usage

``` r
contrastmeans(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", Df, permlist,
                vcov_adjust=TRUE, weights=c("equal", "proportional", "outer", "cells", "flat"),
                df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment"), joint=FALSE,
                level=0.05)
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

- modelterm:

  Name (in "quotes") for indicating which factor term's contrast to be
  calculated. The `modelterm` must be given exactly as it appears in the
  printed model, e.g. "A" or "A:B".

- ctrmatrix:

  A specified contrast matrix. If `ctrmatrix` is missing, the programe
  will ask user to enter it. The columns can be named by level (e.g.
  from
  [`specify_contrast`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md));
  levels that are not given get zero weight.

- ctrnames:

  Names of the specified contrasts, e.g. c("A vs D", "C vs B", ...)

- adj:

  Name (in "quote") for indicating a method for adjusting p-values of
  pairwise comparisons. The choices are "none", "holm", "hochberg",
  "hommel", "bonferroni", "BH", "BY", "fdr", "sidak", "scheffe", "mvt"
  and "westfall", in any letter case. "sidak", "scheffe" (rank of
  `ctrmatrix`) and "mvt" (single-step multivariate t, the exact
  counterpart of SAS `ADJUST=SIMULATE`; needs package mvtnorm) are as in
  emmeans. With `permlist`, "mvt" gives the single-step and "westfall"
  the step-down max-\|t\| adjustment of the permutation p-values
  (Westfall and Young, 1993); "westfall" needs `permlist`, and "scheffe"
  is not available with it. The default method is "none".

- Df:

  A denominator degree of freedom for `modelterm`. (For `glmer` models
  the `Df` needs to be specified, while for the other models, `Df` is
  obtained from the fitted model automatically).

- permlist:

  A model parameter list containing `nsim` parameters produced by the
  function `permmodels`. When `permlist != NULL`, the option `Df` will
  be non-functional. This is a key option for the permutation test. With
  the result of `permmodels(method = "freedman-lane")` (`lm` and `aov`
  models), each contrast is tested by Freedman-Lane permutations of the
  model in which that contrast is 0. With a `permlist`, `vcov_adjust` is
  not used (the permuted refits use the unadjusted covariance).

- vcov_adjust:

  For `lmer` models fitted by REML: if `TRUE` (the default), standard
  errors use the Kenward-Roger adjusted covariance matrix of the fixed
  effects, as SAS `DDFM=KR` and emmeans do, to go with the Kenward-Roger
  degrees of freedom. `FALSE` uses the unadjusted covariance and
  reproduces results of predictmeans 1.x. It has no effect for other
  models.

- weights:

  Passed to
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md):
  how the cells of factors not in `modelterm` are weighted when
  averaging over them (`"equal"`, the default, `"proportional"`,
  `"outer"`, `"cells"` or `"flat"`).

- df_method:

  Denominator df for `lmer` models: `"auto"` (Kenward-Roger for REML
  fits, Satterthwaite for ML fits), `"Kenward-Roger"`, `"Satterthwaite"`
  or `"containment"` (SAS `DDFM=CONTAIN`; see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).

- joint:

  Logical: also test all the contrasts together (all of them zero), as
  [`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
  does: an F test with Kenward-Roger (REML `lmer` fits, the scaled F of
  Kenward and Roger, 1997) or Satterthwaite (ML fits) denominator df,
  nlme's df for `lme`, the residual df for other models (normal-based,
  i.e. chi-square / df1, when there is none). Contrasts that are linear
  combinations of others are dropped, so `df1` is the rank of
  `ctrmatrix`. With `permlist`, also a permutation p-value of the joint
  Wald F statistic, from the same permutations (for Freedman-Lane, the
  model with all the contrasts zero). The default is `FALSE`.

- level:

  Significance level for the confidence limits of the contrasts
  (`lower.CL`, `upper.CL`); the default 0.05 gives 95% limits. They are
  t limits, or with `permlist` permutation-t limits: estimate +- c SE,
  with c the critical value of the permutation test (the k-th largest
  permuted \|t\|, k = floor(level (nperm + 1))); with `adj` "mvt" or
  "westfall" c comes from the largest \|t\| of the contrasts, so the
  limits are simultaneous.

## Value

There are two components in the output which are

- Table:

  A table showing t-test results for the specified linear contrasts.

- K:

  A contrast matrix.

and, with `joint = TRUE`, a third,

- Joint F test of the contrasts:

  `df1`, `df2`, `F.ratio` and `p.value` of the joint test (and
  `Perm_p_value` with `permlist`).

## References

Torsten Hothorn, Frank Bretz and Peter Westfall (2008), *Simultaneous
Inference in General Parametric Models. Biometrical*, Journal 50(3),
346–363.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
# ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## The contrast has a contrast matrix as follows:
#     0:Golden Rain 0:Marvellous 0:Victory 
#[1,]            -1            0         1 
#[2,]             0            0         1 
#     0.2:Golden Rain 0.2:Marvellous 0.2:Victory 
#[1,]               0              0           0 
#[2,]               0              0           0 
#     0.4:Golden Rain  0.4:Marvellous 0.4:Victory
#[1,]               0               0           0
#[2,]               0              -1           0
#      0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#[1,]                0              0           0
#[2,]                0              0           0

# 1. Enter above contrast matrix into a pop up window, then close the window
if (interactive()) contrastmeans(fm, "nitro:Variety")
 
# 2. Construct the contrast matrix directly
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)
#> 
#>  Denominator degree of freedom for ‘nitro:Variety’ and its marginal terms vary between ‘10’ and ‘45’ .
#>  Probabilities will be calculated using ‘10’ Df. 
#> $`The t tests of the specified contrasts`
#>      Estimate Std. Error t value df Pr(>|t|) lower.CL upper.CL
#> [1,]  -8.5000      9.715 -0.8749 10   0.4021 -30.1464  13.1464
#> [2,] -45.6667      9.715 -4.7006 10   0.0008 -67.3131 -24.0202
#> attr(,"Note")
#> [1] "The p-value is adjusted by ‘none’ method, if p-value = 0 means p-value < 0.0001."
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

# 3. ... and both contrasts tested together (an F test with 2 numerator df)
contrastmeans(fm, "nitro:Variety", ctrmatrix=cm, joint=TRUE)
#> 
#>  Denominator degree of freedom for ‘nitro:Variety’ and its marginal terms vary between ‘10’ and ‘45’ .
#>  Probabilities will be calculated using ‘10’ Df. 
#> $`The t tests of the specified contrasts`
#>      Estimate Std. Error t value df Pr(>|t|) lower.CL upper.CL
#> [1,]  -8.5000      9.715 -0.8749 10   0.4021 -30.1464  13.1464
#> [2,] -45.6667      9.715 -4.7006 10   0.0008 -67.3131 -24.0202
#> attr(,"Note")
#> [1] "The p-value is adjusted by ‘none’ method, if p-value = 0 means p-value < 0.0001."
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
#> $`Joint F test of the contrasts`
#>               df1 df2 F.ratio p.value
#> all contrasts   2  10 12.4991  0.0019
#> 
```
