# ANOVA Tables for lmer Fits with Kenward-Roger or Satterthwaite df

[`anova()`](https://rdrr.io/r/stats/anova.html) of a single `lmer` fit
(and of a gaussian
[`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
fit, `fit$semer`) gives F tests with denominator degrees of freedom,
computed by the package's own code: Type III tests with Kenward-Roger df
by default, as SAS PROC MIXED `DDFM=KR`. This replaces the tables of
lmerTest, which predictmeans attached before 2.0. With several models,
`anova(m1, m2)` is lme4's likelihood ratio comparison, unchanged.

## Arguments

- object:

  An `lmer` fit (class `lmerMod`).

- ...:

  Other fitted models: then lme4's comparison of the models.

- type:

  `"III"` (default), `"II"` or `"I"`.

- ddf:

  `"Kenward-Roger"` (default), `"Satterthwaite"` or `"lme4"`.

## Value

An `anova` table.

## Details

Type III hypotheses do not depend on the contrasts the model was fitted
with (sum-to-zero coding, as SAS Type 3); Type II are those of
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
