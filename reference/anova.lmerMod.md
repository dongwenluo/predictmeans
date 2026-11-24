# ANOVA Tables for lmer Fits with Kenward-Roger or Satterthwaite df

[`anova()`](https://rdrr.io/r/stats/anova.html) of a single `lmer` fit
(and of a gaussian
[`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
fit, `fit$semer`) gives F tests with denominator degrees of freedom,
computed by the package's own code: Type III tests with Kenward-Roger df
by default, as SAS PROC MIXED `DDFM=KR`. This replaces the tables of
lmerTest, which predictmeans attached before 2.0. With several models,
`anova(m1, m2)` is lme4's likelihood ratio comparison, unchanged.

## Usage

``` r
# S3 method for class 'lmerMod'
anova(
  object,
  ...,
  type = c("III", "II", "I"),
  ddf = c("Kenward-Roger", "Satterthwaite", "lme4"),
  advice = TRUE
)
```

## Arguments

- object:

  An `lmer` fit (class `lmerMod`).

- ...:

  Other fitted models: then lme4's comparison of the models.

- type:

  `"III"` (default), `"II"` or `"I"`; also `3`, `2`, `1`.

- ddf:

  `"Kenward-Roger"` (default), `"Satterthwaite"` or `"lme4"`.

- advice:

  Logical: add a note on the choice of type when the data are unbalanced
  or have empty cells (default `TRUE`).

## Value

An `anova` table: for each term the numerator df, the denominator df
(column heading names the method), the F value and its p-value. With
`advice = TRUE` a note on the type is kept in `attr(, "advice")` (read
it with
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).

## Details

Type III hypotheses do not depend on the contrasts the model was fitted
with (sum-to-zero coding, as SAS Type 3). With empty cells (aliased
coefficients) each term is tested on the cells with data: an interaction
by its contrasts that use no empty cell (one df fewer per empty cell), a
main effect whose means average over empty cells on the cells its means
share, as
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md);
the table says so, and the intercept is not tested; Type II are those of
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html); Type I
(sequential) are lmerTest's, always with Satterthwaite df. Kenward-Roger
F tests use the adjusted covariance of the fixed effects (`pbkrtest`'s,
observed information); a variance component estimated as zero is left
out of the adjustment, as in SAS. `ddf = "lme4"` gives lme4's own table
(sums of squares and F values without df or p-values).

With lmerTest attached,
[`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) is lmerTest's and its
fits are of class `lmerModLmerTest`: their
[`anova()`](https://rdrr.io/r/stats/anova.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) are lmerTest's
(Satterthwaite by default), and this method applies to fits of
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). Everything
else in predictmeans works the same with either class.

## Which type of test

With balanced data (equal numbers in all cells) Types I, II and III give
the same tests. With unbalanced data:

- Type I:

  sequential: each term adjusted for the terms before it, so the order
  in the formula matters. Use it when the order is meaningful, e.g. a
  block or a covariate to be removed first.

- Type II:

  each term adjusted for all terms that do not contain it (the principle
  of marginality). The most powerful tests of main effects when the
  interactions are negligible; with an interaction present the main
  effects are not well defined.

- Type III:

  each term adjusted for all others, testing equally weighted cell
  means; meaningful with an interaction present, and here independent of
  the model's contrasts. A covariate in an interaction makes the other
  terms be tested at covariate = 0: centre it.

With empty cells no type answers the question of a main effect over all
cells: Type III uses the cells with data (and says so). State the
comparisons of interest on the cell means instead, with
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
or
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).

## See also

[`anova.semireg`](https://dongwenluo.github.io/predictmeans/reference/anova.semireg.md),
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md),
[`df_term`](https://dongwenluo.github.io/predictmeans/reference/df_term.md),
[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
anova(fm)
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: yield
#>                      F Df Df.res    Pr(>F)    
#> (Intercept)   245.1374  1      5 1.932e-05 ***
#> nitro          37.6857  3     45 2.458e-12 ***
#> Variety         1.4853  2     10    0.2724    
#> nitro:Variety   0.3028  6     45    0.9322    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(fm, ddf = "Satterthwaite")
#> Analysis of Deviance Table (Type III Wald F tests with Satterthwaite df)
#> 
#> Response: yield
#>                      F Df Df.res    Pr(>F)    
#> (Intercept)   245.1374  1      5 1.932e-05 ***
#> nitro          37.6857  3     45 2.458e-12 ***
#> Variety         1.4853  2     10    0.2724    
#> nitro:Variety   0.3028  6     45    0.9322    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
anova(fm, type = "I")
#> Type I Analysis of Variance Table with Satterthwaite's method
#>                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#> nitro         20020.5  6673.5     3    45 37.6857 2.458e-12 ***
#> Variety         526.1   263.0     2    10  1.4853    0.2724    
#> nitro:Variety   321.7    53.6     6    45  0.3028    0.9322    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
