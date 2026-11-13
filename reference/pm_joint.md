# Joint F Tests of Contrasts among Means (Slices)

F tests of sets of contrasts among the predicted means of a model term.
With `by`, the means of the term are compared within each level of the
`by` factor(s): the test of `nitro` within each `Variety`, as SAS
`LSMEANS nitro*Variety / SLICE=Variety` or
`emmeans::joint_tests(model, by = "Variety")`. With `ctrmatrix`, all its
rows are tested together.

## Usage

``` r
pm_joint(model, modelterm, by = NULL, ctrmatrix = NULL,
         df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "inner-outer", "residual"),
         weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL)
```

## Arguments

- model:

  A fitted model, as for
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- modelterm:

  Name (in "quotes") of the model term whose means are compared, e.g.
  `"nitro:Variety"`.

- by:

  Optional name(s) of factor(s) of `modelterm`: test the differences
  among the means within each of their levels. Without `by` and
  `ctrmatrix`, all differences among the means of `modelterm` are
  tested.

- ctrmatrix:

  Optional contrast matrix over the means of `modelterm` (as for
  [`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
  e.g. from
  [`contrast_family`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md)
  or
  [`specify_contrast`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md));
  its rows are tested together.

- df_method:

  For `lmer` models (`"between-within"` also for `lme` and `gls`): the
  denominator df (see Details).

- weights:

  How the means average over other factors, as in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

- vcov:

  Optional covariance of the coefficients, as in
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md):
  the F statistic is then the Wald F with this covariance, and the
  denominator df those of the model's own test.

## Value

A data frame with one row per test: the `by` levels (if any), `df1`,
`df2`, `F.ratio` and `p.value`. Attribute `"df_method"` says how `df2`
was found.

## Details

The F statistic and its denominator degrees of freedom follow the model:

- `lmer`:

  Kenward-Roger (the default for REML fits): the scaled F test of
  Kenward and Roger (1997) with the observed information of the variance
  parameters, as SAS `DDFM=KR`; Satterthwaite (the default for ML fits):
  as `lmerTest::contest(joint = TRUE)`; containment: the Wald F with the
  containment df of `modelterm`.

- KR-expected:

  (`lmer`, `lme`, `gls`) Kenward-Roger with the expected information, as
  [`pbkrtest::KRmodcomp()`](https://rdrr.io/pkg/pbkrtest/man/kr__modcomp.html)
  (see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).

- between-within:

  (`lmer`, `lme`, `gls`) the Wald F with the between-within df of
  `modelterm` (SAS `DDFM=BETWITHIN`; see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).

- `lme`, `gls` with `df_method = "Satterthwaite"`, `"Kenward-Roger"` or
  `"KR-expected"`:

  as for `lmer` (Satterthwaite df of the F test, lmerTest's rule).

- `lme`, `gls`:

  Kenward-Roger for REML fits (the default), as for `lmer`; with
  `df_method = "containment"` (and for ML fits) the Wald F with nlme's
  df of `modelterm` (as [`anova()`](https://rdrr.io/r/stats/anova.html);
  `"containment (nlme)"`) or the residual df (`gls`).

- other models:

  the Wald F with the residual df (normal-based when there is none).

Contrasts that are linear combinations of others are dropped, so `df1`
is the rank of the set.
[`emmeans::joint_tests()`](https://rvlenth.github.io/emmeans/reference/joint_tests.html)
gives a slightly different Kenward-Roger test (the unscaled F with its
own df).

## References

Kenward, M. G. and Roger, J. H. (1997), Small sample inference for fixed
effects from restricted maximum likelihood, *Biometrics*, 53, 983-997.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  # nitrogen within each variety (SAS SLICE)
  pm_joint(fm, "nitro:Variety", by = "Variety")
#>       Variety df1      df2   F.ratio      p.value
#> 1 Golden Rain   3 45.00014 13.022765 2.969017e-06
#> 2  Marvellous   3 45.00014  9.952878 3.769239e-05
#> 3     Victory   3 45.00014 15.315709 5.256861e-07
  # varieties within each nitrogen level
  pm_joint(fm, "nitro:Variety", by = "nitro")
#>   nitro df1     df2   F.ratio   p.value
#> 1     0   2 30.2308 1.2245400 0.3080589
#> 2   0.2   2 30.2308 1.8814476 0.1697812
#> 3   0.4   2 30.2308 0.2156336 0.8072596
#> 4   0.6   2 30.2308 0.4010511 0.6731259
  # a set of contrasts: linear and quadratic nitrogen trends together
  pm_joint(fm, "nitro", ctrmatrix = contrast_family(Kmatrix(fm, "nitro"), "poly", max.degree = 2))
#>   df1      df2  F.ratio     p.value
#> 1   2 45.00014 56.51839 5.31137e-13
```
