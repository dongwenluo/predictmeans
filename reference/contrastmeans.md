# Linear Contrast Tests for a Linear Model

t tests (or permutation t tests) of given contrasts among the predicted
means of a model term, for the models of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
(`aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer`, `glmmTMB`, `gam`).

## Usage

``` r
contrastmeans(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", permlist,
              vcov_adjust=TRUE,
              weights=c("equal", "proportional", "outer", "cells", "flat"),
              df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment",
      "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
              joint=FALSE, conf_level=0.95, vcov=NULL, equivalence=NULL, ...)
```

## Arguments

- model:

  A fitted model, as for
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- modelterm:

  Name (in "quotes") for indicating which factor term's contrast to be
  calculated. The `modelterm` must be given exactly as it appears in the
  printed model, e.g. "A" or "A:B".

- ctrmatrix:

  A specified contrast matrix. If it is missing, an edit window asks for
  it (interactive sessions only). The columns can be named by level
  (e.g. from
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

- permlist:

  A list from
  [`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md):
  each contrast is then tested by permutation, with p-value \\(b +
  1)/(B + 1)\\ (\\b\\ of the \\B\\ permuted \|t\| at least the observed
  one), and `df_method` is not used. With the result of
  `permmodels(method = "freedman-lane")` (`lm` and `aov` models), each
  contrast is tested by Freedman-Lane permutations of the model in which
  that contrast is 0. With a `permlist`, `vcov_adjust` is not used (the
  permuted refits use the unadjusted covariance).

- vcov_adjust:

  For models fitted by REML (`lmer`, `lme`, `gls`) with Kenward-Roger
  df: if `TRUE` (the default), standard errors use the Kenward-Roger
  adjusted covariance matrix of the fixed effects, as SAS `DDFM=KR` and
  emmeans do, to go with the Kenward-Roger degrees of freedom. `FALSE`
  uses the unadjusted covariance and reproduces results of predictmeans
  1.x. It has no effect for other models.

- weights:

  Passed to
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md):
  how the cells of factors not in `modelterm` are weighted when
  averaging over them (`"equal"`, the default, `"proportional"`,
  `"outer"`, `"cells"` or `"flat"`).

- df_method:

  Denominator df, by method or as one number (a fixed df, `Df` before
  2.0): `"auto"` (Kenward-Roger for REML fits of `lmer`, `lme` and
  `gls`, Satterthwaite for ML fits), `"Kenward-Roger"`,
  `"Satterthwaite"` or `"containment"` (SAS `DDFM=CONTAIN`), and for
  `lmer`, `lme` and `gls` models `"between-within"` (SAS
  `DDFM=BETWITHIN`), and `"KR-expected"` (Kenward-Roger with the
  expected information, as pbkrtest; `"Kenward-Roger"` uses the observed
  information, as SAS `DDFM=KR`), `"KR-average"` (the average
  information, as Genstat); for `lme` and `gls` also `"Satterthwaite"`,
  `"Kenward-Roger"`, `"KR-expected"`, `"KR-average"` and `"containment"`
  (nlme's or the residual df, as predictmeans 1.x); see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- joint:

  Logical: also test all the contrasts together (all of them zero), as
  [`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
  does: an F test with Kenward-Roger (REML fits, the scaled F of Kenward
  and Roger, 1997) or Satterthwaite (ML fits) denominator df, the
  residual df for other models (normal-based, i.e. chi-square / df1,
  when there is none). Contrasts that are linear combinations of others
  are dropped, so `df1` is the rank of `ctrmatrix`. With `permlist`,
  also a permutation p-value of the joint Wald F statistic, from the
  same permutations (for Freedman-Lane, the model with all the contrasts
  zero). The default is `FALSE`.

- conf_level:

  Confidence level of the limits of the contrasts (`level = 0.05` before
  2.0); the tests are at the significance level `1 - conf_level`. The
  limits are t limits of each contrast, or with `permlist` permutation-t
  limits: estimate +- c SE, with c the critical value of the permutation
  test (the k-th largest permuted \|t\|, k = floor((1 - conf_level)
  (nperm + 1))). With `adj` "bonferroni", "sidak", "scheffe" or "mvt"
  (with `permlist`: "bonferroni", "sidak", "mvt" or "westfall", from the
  largest permuted \|t\| of the contrasts) the table also has
  `adj.lower.CL` and `adj.upper.CL`: simultaneous limits that exclude 0
  exactly when the adjusted p-value is below `1 - conf_level`, as
  emmeans' [`confint()`](https://rdrr.io/r/stats/confint.html) and SAS
  `AdjLower`/`AdjUpper`. The step-down and false discovery rate methods
  have none.

- vcov:

  Optional covariance of the coefficients to use instead of the model's.
  For `lm`, `aov` and `glm` models: `"HC0"`, `"HC1"`, `"HC2"` or `"HC3"`
  (heteroscedasticity-consistent, as
  [`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)),
  or a one-sided formula naming a cluster variable of the data, e.g.
  `~ Block` (cluster-robust CR1, as
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)).
  For any model: a covariance matrix of the coefficients, or a function
  that returns one from the model (for mixed models e.g. from the
  clubSandwich package). The degrees of freedom stay those of the model
  (the residual df; for `lmer` Kenward-Roger or Satterthwaite df
  computed from the model), and `vcov_adjust` is not used. With a
  `permlist`, `vcov` is not used (for robust permutation tests see
  `permmodels(robust = TRUE)`). The default `NULL` uses the model's
  covariance.

- equivalence:

  Optional equivalence margins for equivalence tests (two one-sided t
  tests, TOST) of the contrasts: a positive number \\d\\ (equivalence
  region \\-d\\ to \\d\\) or `c(lower, upper)`, on the scale of the
  comparisons (the link scale for GLMs: e.g. log ratios). The TOST
  p-value is the larger of the two one-sided p-values; below
  `1 - conf_level`, the difference lies within the margins (as emmeans'
  `test(..., side = "equivalence")`). Not adjusted for multiplicity;
  with a `permlist` it is t-based (normal-based in `contrastmeans`).

- ...:

  Not used: arguments removed in 2.0 (`Df`, `level`) stop with their 2.0
  spelling.

## Value

An object of class `"contrastmeans"`, a list.
[`print()`](https://rdrr.io/r/base/print.html) shows the contrasts,
`tidy()` and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
the table `contrasts` (`what = "joint"`: the joint test). Its first
components are those of 1.x:

- The t tests of the specified contrasts:

  A table (rounded to 4 decimals) of the t tests of the contrasts, with
  the adjusted p-value in `Pr(>|t|)`.

- K:

  A contrast matrix.

- Joint F test of the contrasts:

  With `joint = TRUE`: `df1`, `df2`, `F.ratio` and `p.value` of the
  joint test (and `Perm_p_value` with `permlist`).

and then the tables in the layout of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
unrounded:

- contrasts:

  One row per contrast: `contrast` (its name), `over` (contrasts of
  means over empty cells: the cells used), `estimate`, `std.error`,
  `df`, `statistic` (t), `p.value` (unadjusted; with `permlist` the
  permutation p-value), `p.adj` (adjusted by `adj`), `p.equiv` (with
  `equivalence`), `conf.low`, `conf.high`, `adj.conf.low`,
  `adj.conf.high` (simultaneous, with the adjustments that have them),
  `p.method`, and for a log or logit scale (link, or a log response)
  `ratio` or `odds.ratio` with its limits when the weights of every
  contrast sum to zero.

- joint:

  With `joint = TRUE`: the joint test, unrounded.

- settings:

  `modelterm`, `level` (significance level, `1 - conf_level`), `adjust`,
  `df_method`, `vcov_adjust`, `vcov`, `equivalence`, `nperm`,
  `model_class`, `version`.

## References

Torsten Hothorn, Frank Bretz and Peter Westfall (2008), *Simultaneous
Inference in General Parametric Models. Biometrical*, Journal 50(3),
346–363.

## See also

[`specify_contrast`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
[`contrast_family`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md),
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md),
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`tidy.predictmeans`](https://dongwenluo.github.io/predictmeans/reference/tidy.predictmeans.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)

# 1. Enter a contrast matrix in an edit window (interactive sessions only)
if (interactive()) contrastmeans(fm, "nitro:Variety")
 
# 2. Give the contrast matrix: Victory - Golden Rain at nitro 0, Victory at 0 - Marvellous at 0.4
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)
#> Contrasts of 'nitro:Variety' (95% limits; p.adj: none; Kenward-Roger df)
#>    contrast estimate std.error    df statistic   p.value     p.adj conf.low
#>  contrast 1    -8.50     9.715 30.23   -0.8749 3.885e-01 3.885e-01   -28.33
#>  contrast 2   -45.67     9.715 30.23   -4.7006 5.332e-05 5.332e-05   -65.50
#>  conf.high
#>      11.33
#>     -25.83

# 3. ... and both contrasts tested together (an F test with 2 numerator df)
ct <- contrastmeans(fm, "nitro:Variety", ctrmatrix=cm, joint=TRUE)
ct
#> Contrasts of 'nitro:Variety' (95% limits; p.adj: none; Kenward-Roger df)
#>    contrast estimate std.error    df statistic   p.value     p.adj conf.low
#>  contrast 1    -8.50     9.715 30.23   -0.8749 3.885e-01 3.885e-01   -28.33
#>  contrast 2   -45.67     9.715 30.23   -4.7006 5.332e-05 5.332e-05   -65.50
#>  conf.high
#>      11.33
#>     -25.83
#> Joint F test: F(2, 30.2) = 12.499, p = 0.000111
as.data.frame(ct)
#>     contrast  estimate std.error       df  statistic      p.value        p.adj
#> 1 contrast 1  -8.50000  9.715029 30.23074 -0.8749331 3.885088e-01 3.885088e-01
#> 2 contrast 2 -45.66667  9.715029 30.23074 -4.7006208 5.331539e-05 5.331539e-05
#>    conf.low conf.high p.method
#> 1 -28.33439  11.33439        t
#> 2 -65.50105 -25.83228        t
```
