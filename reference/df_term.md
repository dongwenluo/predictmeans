# Denominator Degrees of Freedom of the Means or Contrasts of an lmer Model

The denominator df of each mean of a model term, or of each row of a
contrast matrix, for an `lmer` fit: Kenward-Roger (default),
Satterthwaite, containment or between-within.

## Usage

``` r
df_term(model, modelterm, covariate=NULL, ctrmatrix=NULL, ctrnames=NULL,
        df_method=c("Kenward-Roger", "Satterthwaite", "containment",
          "between-within", "KR-expected", "KR-average"), ...)
```

## Arguments

- model:

  Model object returned by `lmer`.

- modelterm:

  Name (in "quotes") of the model term whose means get df, as written in
  the model formula, e.g. "A" or "A:B".

- covariate:

  Name (in "quotes") of a covariate in the `model`, passed to
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).

- ctrmatrix:

  A contrast matrix on the fixed-effect coefficients (one column per
  coefficient of `fixef(model)`). When given, `modelterm` is ignored
  (except for "containment" and "between-within", which need it).

- ctrnames:

  Names of the rows of `ctrmatrix`, e.g. c("A vs D", "C vs B").

- df_method:

  The method for the degrees of freedom (called `type` before 2.0):
  "Kenward-Roger" (the default), "Satterthwaite", "containment",
  "between-within", "KR-expected" or "KR-average". Kenward-Roger needs a
  REML fit; for models fitted by ML, Satterthwaite is used with a
  warning. "containment" gives the SAS `DDFM=CONTAIN` df of `modelterm`
  (needed then) for every row, and "between-within" the SAS
  `DDFM=BETWITHIN` df (see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).
  "Kenward-Roger" uses the observed information of the variance
  parameters, as SAS `DDFM=KR`; "KR-expected" the expected information,
  as pbkrtest and emmeans; "KR-average" the average information, as
  Genstat.

- ...:

  Arguments removed in 2.0 (`type`, `vcov_adj`) stop with their 2.0
  spelling.

## Value

A numeric vector of df, one per mean of `modelterm` (rows of
`Kmatrix(model, modelterm)$K`) or per row of `ctrmatrix`.

## Details

Kenward-Roger df that cannot be computed fall back to Satterthwaite with
a warning. The Satterthwaite df equal those of lmerTest and the
Kenward-Roger df those of pbkrtest (`"KR-expected"`) and SAS `DDFM=KR`
(`tests/testthat/test-satterthwaite.R`, `test-kr-sas.R`).

## See also

[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
(argument `df_method`),
[`anova.lmerMod`](https://dongwenluo.github.io/predictmeans/reference/anova.lmerMod.md),
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
df_term(fm, "nitro:Variety")
#>   0:Golden Rain    0:Marvellous       0:Victory 0.2:Golden Rain  0.2:Marvellous 
#>        16.08144        16.08144        16.08144        16.08144        16.08144 
#>     0.2:Victory 0.4:Golden Rain  0.4:Marvellous     0.4:Victory 0.6:Golden Rain 
#>        16.08144        16.08144        16.08144        16.08144        16.08144 
#>  0.6:Marvellous     0.6:Victory 
#>        16.08144        16.08144 
df_term(fm, "nitro:Variety", df_method = "containment")
#>  [1] 45 45 45 45 45 45 45 45 45 45 45 45
# two contrasts of the 12 coefficients
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
df_term(fm, ctrmatrix=cm, df_method="Satterthwaite")
#> [1] 47.32750 45.00004
```
