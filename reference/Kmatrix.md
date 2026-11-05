# Matrix of Coefficients in a Linear Model

This function obtains a matrix of coefficients (the linear functions
\\K\\ such that the predicted means are \\K\hat\beta\\) for models such
as `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer`, `glmmTMB`, `gam`
(package mgcv) and `semireg`.

The predicted means are built on a reference grid: every combination of
the factor levels in the model, with covariates fixed at their means (or
other chosen values). Each row of `K` is the equally weighted average of
the grid rows for one level combination of `modelterm`. The design
matrix of the grid is evaluated with the model's own stored basis
information, so covariates may enter the model through transformations
such as `poly(x, 2)`, `ns(x, 3)`, `bs(x)`, `log(x)`, `scale(x)` or
`I(x^2)`, only through interactions (e.g. `A:x` without `x`, or `x:z`),
or through smooth terms (`s(x)`, `te(x, z)`, `s(x, by = A)`) of a `gam`.

## Usage

``` r
Kmatrix(model, modelterm, covariate=NULL, covariateV=NULL, data=NULL, prtnum=FALSE,
        cov_reduce=mean, exclude=NULL,
        weights=c("equal", "proportional", "outer", "cells", "flat"))
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`,
  `glmer`, `glmmTMB`, `gam`/`bam` (package mgcv; for `gamm` use its
  `$gam` component) or `semireg`.

- modelterm:

  Name (in "quotes") for indicating which model term's predicted mean to
  be calculated. The `modelterm` must be given exactly as it appears in
  the printed model, e.g. "A" or "A:B". A factor created in the model
  formula, e.g. `factor(dose)`, can be given either as "factor(dose)" or
  as "dose".

- covariate:

  A numerical vector to specify values of covariates for calculating
  predicted means, default values are the means of the associated
  covariates. It also can be the name of one or more covariates in the
  model, in which case predictions are made over a range of values of
  those covariates (see `covariateV`). For a `semireg` model this may
  also be a covariate that only appears in the smooth terms.

- covariateV:

  Values for the covariate(s) named in `covariate`. For one covariate, a
  numeric vector. For several covariates, either a list of numeric
  vectors (all combinations are used) or a matrix with one named column
  per covariate, where each row is one prediction point (e.g. a set of
  locations for a two-dimensional smooth). When `covariateV` is `NULL`,
  50 (one covariate) or 30 (several covariates) equally spaced values
  over the observed range are used.

- data:

  In some cases, you need to provide the data set used in model fitting,
  for example when the data used to fit the model can no longer be
  found. Usually the data are recovered from the model.

- prtnum:

  An option for printing covariate info on the screen or not. The
  default is FALSE.

- cov_reduce:

  How covariates that are not varied are handled. A function returning a
  single value, applied to each covariate (the default `mean`; `median`
  is another choice), or the string `"average"`, which averages the
  predictions over the observed values of those covariates instead of
  fixing them at one value. With `"average"`, nonlinear terms such as
  `poly(x, 2)` or `s(x)` are averaged over the data, which for a `gam`
  gives predicted means adjusted for (rather than evaluated at the
  centre of) a spatial or temporal trend.

- exclude:

  For `gam` models only: labels of smooth terms to leave out of the
  predictions (their columns of `K` are set to zero), e.g. `"s(Block)"`.
  The default `NULL` excludes random-effect smooths (`bs = "re"` and
  `bs = "fs"`), so the predicted means are population-level means, as
  for `lmer`. Use `exclude = character(0)` to keep all smooth terms.

- weights:

  How the cells of the factors that are not in `modelterm` are weighted
  when averaging over them: `"equal"` (the default; equal weights),
  `"proportional"` (weights proportional to the frequencies of those
  cells in the data, the same for every mean, as SAS `LSMEANS / OM`),
  `"outer"` (products of each factor's own frequencies), `"cells"` (each
  mean weights the cells by their own frequencies) or `"flat"` (equal
  weights over the non-empty cells). The definitions are those of
  emmeans.

## Value

- K:

  Coefficients matrix. It has an attribute `"estimable"`, a logical
  vector with one value per row that is `FALSE` for rows that are not
  estimable functions of the model parameters (e.g. means for a missing
  cell in a rank-deficient model).

- fctnames:

  A model frame contains factor(s) info in the model.

- response:

  The name of response variable in the model.

- preddf:

  The data of the predictor variables used in model fitting.

- refgrid:

  The reference grid used to build `K`: a list with the grid data
  (`grid`), the sparse weight matrix `W` that averages grid rows into
  rows of `K`, and the rows of `fctnames` (`target`). Any other basis
  evaluated on `grid`, e.g. the smooth terms of a `semireg` model, can
  be averaged in the same way as `W %*% basis`.

## References

Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
*Prediction in linear mixed models*, Australian and New Zealand Journal
of Statistics, 46(3), 325-347.

## Author

This function heavily depends on the codes from package "lsmeans".

## Examples

``` r
  library(predictmeans)
  data(Oats, package="nlme")
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  Kmatrix(fm, "Variety", prtnum=TRUE)$K
#> 
#>  The predicted means are estimated at 
#> nitro 
#>   0.3 
#>             (Intercept) nitro VarietyMarvellous VarietyVictory
#> Golden Rain           1   0.3                 0              0
#> Marvellous            1   0.3                 1              0
#> Victory               1   0.3                 0              1
#>             nitro:VarietyMarvellous nitro:VarietyVictory
#> Golden Rain                     0.0                  0.0
#> Marvellous                      0.3                  0.0
#> Victory                         0.0                  0.3
#> attr(,"estimable")
#> Golden Rain  Marvellous     Victory 
#>        TRUE        TRUE        TRUE 
  Kmatrix(fm, "Variety", 0.5, prtnum=TRUE)$K
#> 
#>  The predicted means are estimated at 
#> nitro 
#>   0.5 
#>             (Intercept) nitro VarietyMarvellous VarietyVictory
#> Golden Rain           1   0.5                 0              0
#> Marvellous            1   0.5                 1              0
#> Victory               1   0.5                 0              1
#>             nitro:VarietyMarvellous nitro:VarietyVictory
#> Golden Rain                     0.0                  0.0
#> Marvellous                      0.5                  0.0
#> Victory                         0.0                  0.5
#> attr(,"estimable")
#> Golden Rain  Marvellous     Victory 
#>        TRUE        TRUE        TRUE 
 # Kmatrix(fm, "Variety", "nitro")$K
  Kmatrix(fm, "Variety", "nitro", covariateV=seq(0, 0.6, 0.1))$K
#>                 (Intercept) nitro VarietyMarvellous VarietyVictory
#> Golden Rain:0             1   0.0                 0              0
#> Golden Rain:0.1           1   0.1                 0              0
#> Golden Rain:0.2           1   0.2                 0              0
#> Golden Rain:0.3           1   0.3                 0              0
#> Golden Rain:0.4           1   0.4                 0              0
#> Golden Rain:0.5           1   0.5                 0              0
#> Golden Rain:0.6           1   0.6                 0              0
#> Marvellous:0              1   0.0                 1              0
#> Marvellous:0.1            1   0.1                 1              0
#> Marvellous:0.2            1   0.2                 1              0
#> Marvellous:0.3            1   0.3                 1              0
#> Marvellous:0.4            1   0.4                 1              0
#> Marvellous:0.5            1   0.5                 1              0
#> Marvellous:0.6            1   0.6                 1              0
#> Victory:0                 1   0.0                 0              1
#> Victory:0.1               1   0.1                 0              1
#> Victory:0.2               1   0.2                 0              1
#> Victory:0.3               1   0.3                 0              1
#> Victory:0.4               1   0.4                 0              1
#> Victory:0.5               1   0.5                 0              1
#> Victory:0.6               1   0.6                 0              1
#>                 nitro:VarietyMarvellous nitro:VarietyVictory
#> Golden Rain:0                       0.0                  0.0
#> Golden Rain:0.1                     0.0                  0.0
#> Golden Rain:0.2                     0.0                  0.0
#> Golden Rain:0.3                     0.0                  0.0
#> Golden Rain:0.4                     0.0                  0.0
#> Golden Rain:0.5                     0.0                  0.0
#> Golden Rain:0.6                     0.0                  0.0
#> Marvellous:0                        0.0                  0.0
#> Marvellous:0.1                      0.1                  0.0
#> Marvellous:0.2                      0.2                  0.0
#> Marvellous:0.3                      0.3                  0.0
#> Marvellous:0.4                      0.4                  0.0
#> Marvellous:0.5                      0.5                  0.0
#> Marvellous:0.6                      0.6                  0.0
#> Victory:0                           0.0                  0.0
#> Victory:0.1                         0.0                  0.1
#> Victory:0.2                         0.0                  0.2
#> Victory:0.3                         0.0                  0.3
#> Victory:0.4                         0.0                  0.4
#> Victory:0.5                         0.0                  0.5
#> Victory:0.6                         0.0                  0.6
#> attr(,"estimable")
#>   Golden Rain:0 Golden Rain:0.1 Golden Rain:0.2 Golden Rain:0.3 Golden Rain:0.4 
#>            TRUE            TRUE            TRUE            TRUE            TRUE 
#> Golden Rain:0.5 Golden Rain:0.6    Marvellous:0  Marvellous:0.1  Marvellous:0.2 
#>            TRUE            TRUE            TRUE            TRUE            TRUE 
#>  Marvellous:0.3  Marvellous:0.4  Marvellous:0.5  Marvellous:0.6       Victory:0 
#>            TRUE            TRUE            TRUE            TRUE            TRUE 
#>     Victory:0.1     Victory:0.2     Victory:0.3     Victory:0.4     Victory:0.5 
#>            TRUE            TRUE            TRUE            TRUE            TRUE 
#>     Victory:0.6 
#>            TRUE 

  # Covariates entering through transformations or only through interactions
  fm2 <- lm(yield ~ Block + Variety + poly(nitro, 2), data=Oats)
  Kmatrix(fm2, "Variety")$K
#>             (Intercept)       Block.L       Block.Q       Block.C      Block^4
#> Golden Rain           1 -6.938894e-17 -2.775558e-17 -4.163336e-17 2.775558e-17
#> Marvellous            1 -6.938894e-17 -2.775558e-17 -4.163336e-17 2.775558e-17
#> Victory               1 -6.938894e-17 -2.775558e-17 -4.163336e-17 2.775558e-17
#>                   Block^5 VarietyMarvellous VarietyVictory poly(nitro, 2)1
#> Golden Rain -1.387779e-17                 0              0               0
#> Marvellous  -1.387779e-17                 1              0               0
#> Victory     -1.387779e-17                 0              1               0
#>             poly(nitro, 2)2
#> Golden Rain      -0.1473139
#> Marvellous       -0.1473139
#> Victory          -0.1473139
#> attr(,"estimable")
#> Golden Rain  Marvellous     Victory 
#>        TRUE        TRUE        TRUE 
  fm3 <- lm(yield ~ Variety + Variety:nitro, data=Oats)
  Kmatrix(fm3, "Variety", "nitro", covariateV=c(0, 0.3, 0.6))$K
#>                 (Intercept) VarietyMarvellous VarietyVictory
#> Golden Rain:0             1                 0              0
#> Golden Rain:0.3           1                 0              0
#> Golden Rain:0.6           1                 0              0
#> Marvellous:0              1                 1              0
#> Marvellous:0.3            1                 1              0
#> Marvellous:0.6            1                 1              0
#> Victory:0                 1                 0              1
#> Victory:0.3               1                 0              1
#> Victory:0.6               1                 0              1
#>                 VarietyGolden Rain:nitro VarietyMarvellous:nitro
#> Golden Rain:0                        0.0                     0.0
#> Golden Rain:0.3                      0.3                     0.0
#> Golden Rain:0.6                      0.6                     0.0
#> Marvellous:0                         0.0                     0.0
#> Marvellous:0.3                       0.0                     0.3
#> Marvellous:0.6                       0.0                     0.6
#> Victory:0                            0.0                     0.0
#> Victory:0.3                          0.0                     0.0
#> Victory:0.6                          0.0                     0.0
#>                 VarietyVictory:nitro
#> Golden Rain:0                    0.0
#> Golden Rain:0.3                  0.0
#> Golden Rain:0.6                  0.0
#> Marvellous:0                     0.0
#> Marvellous:0.3                   0.0
#> Marvellous:0.6                   0.0
#> Victory:0                        0.0
#> Victory:0.3                      0.3
#> Victory:0.6                      0.6
#> attr(,"estimable")
#>   Golden Rain:0 Golden Rain:0.3 Golden Rain:0.6    Marvellous:0  Marvellous:0.3 
#>            TRUE            TRUE            TRUE            TRUE            TRUE 
#>  Marvellous:0.6       Victory:0     Victory:0.3     Victory:0.6 
#>            TRUE            TRUE            TRUE            TRUE 

  # A smooth term in an mgcv model, averaged over the observed nitro values
  if (requireNamespace("mgcv", quietly = TRUE)) {
    fm4 <- mgcv::gam(yield ~ Variety + s(nitro, k = 4) + s(Block, bs = "re"), data = Oats)
    Kmatrix(fm4, "Variety", cov_reduce = "average")$K
  }
#>             (Intercept) VarietyMarvellous VarietyVictory    s(nitro).1
#> Golden Rain           1                 0              0 -8.326673e-17
#> Marvellous            1                 1              0 -8.326673e-17
#> Victory               1                 0              1 -8.326673e-17
#>               s(nitro).2   s(nitro).3 s(Block).1 s(Block).2 s(Block).3
#> Golden Rain 2.428613e-17 5.551115e-17          0          0          0
#> Marvellous  2.428613e-17 5.551115e-17          0          0          0
#> Victory     2.428613e-17 5.551115e-17          0          0          0
#>             s(Block).4 s(Block).5 s(Block).6
#> Golden Rain          0          0          0
#> Marvellous           0          0          0
#> Victory              0          0          0
#> attr(,"estimable")
#> Golden Rain  Marvellous     Victory 
#>        TRUE        TRUE        TRUE 
```
