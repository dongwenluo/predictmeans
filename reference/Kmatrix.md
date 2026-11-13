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
  rows of `K`, and the rows of `fctnames` (`target`), and `basis`, the
  function that gives the model's basis for grid rows. Any other basis
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
  Kmatrix(fm, "Variety", "nitro")$K
#>                                (Intercept)      nitro VarietyMarvellous
#> Golden Rain:0                            1 0.00000000                 0
#> Golden Rain:0.0122448979591837           1 0.01224490                 0
#> Golden Rain:0.0244897959183673           1 0.02448980                 0
#> Golden Rain:0.036734693877551            1 0.03673469                 0
#> Golden Rain:0.0489795918367347           1 0.04897959                 0
#> Golden Rain:0.0612244897959184           1 0.06122449                 0
#> Golden Rain:0.073469387755102            1 0.07346939                 0
#> Golden Rain:0.0857142857142857           1 0.08571429                 0
#> Golden Rain:0.0979591836734694           1 0.09795918                 0
#> Golden Rain:0.110204081632653            1 0.11020408                 0
#> Golden Rain:0.122448979591837            1 0.12244898                 0
#> Golden Rain:0.13469387755102             1 0.13469388                 0
#> Golden Rain:0.146938775510204            1 0.14693878                 0
#> Golden Rain:0.159183673469388            1 0.15918367                 0
#> Golden Rain:0.171428571428571            1 0.17142857                 0
#> Golden Rain:0.183673469387755            1 0.18367347                 0
#> Golden Rain:0.195918367346939            1 0.19591837                 0
#> Golden Rain:0.2                          1 0.20000000                 0
#> Golden Rain:0.208163265306122            1 0.20816327                 0
#> Golden Rain:0.220408163265306            1 0.22040816                 0
#> Golden Rain:0.23265306122449             1 0.23265306                 0
#> Golden Rain:0.244897959183673            1 0.24489796                 0
#> Golden Rain:0.257142857142857            1 0.25714286                 0
#> Golden Rain:0.269387755102041            1 0.26938776                 0
#> Golden Rain:0.281632653061224            1 0.28163265                 0
#> Golden Rain:0.293877551020408            1 0.29387755                 0
#> Golden Rain:0.306122448979592            1 0.30612245                 0
#> Golden Rain:0.318367346938775            1 0.31836735                 0
#> Golden Rain:0.330612244897959            1 0.33061224                 0
#> Golden Rain:0.342857142857143            1 0.34285714                 0
#> Golden Rain:0.355102040816327            1 0.35510204                 0
#> Golden Rain:0.36734693877551             1 0.36734694                 0
#> Golden Rain:0.379591836734694            1 0.37959184                 0
#> Golden Rain:0.391836734693878            1 0.39183673                 0
#> Golden Rain:0.4                          1 0.40000000                 0
#> Golden Rain:0.404081632653061            1 0.40408163                 0
#> Golden Rain:0.416326530612245            1 0.41632653                 0
#> Golden Rain:0.428571428571429            1 0.42857143                 0
#> Golden Rain:0.440816326530612            1 0.44081633                 0
#> Golden Rain:0.453061224489796            1 0.45306122                 0
#> Golden Rain:0.46530612244898             1 0.46530612                 0
#> Golden Rain:0.477551020408163            1 0.47755102                 0
#> Golden Rain:0.489795918367347            1 0.48979592                 0
#> Golden Rain:0.502040816326531            1 0.50204082                 0
#> Golden Rain:0.514285714285714            1 0.51428571                 0
#> Golden Rain:0.526530612244898            1 0.52653061                 0
#> Golden Rain:0.538775510204082            1 0.53877551                 0
#> Golden Rain:0.551020408163265            1 0.55102041                 0
#> Golden Rain:0.563265306122449            1 0.56326531                 0
#> Golden Rain:0.575510204081633            1 0.57551020                 0
#> Golden Rain:0.587755102040816            1 0.58775510                 0
#> Golden Rain:0.6                          1 0.60000000                 0
#> Marvellous:0                             1 0.00000000                 1
#> Marvellous:0.0122448979591837            1 0.01224490                 1
#> Marvellous:0.0244897959183673            1 0.02448980                 1
#> Marvellous:0.036734693877551             1 0.03673469                 1
#> Marvellous:0.0489795918367347            1 0.04897959                 1
#> Marvellous:0.0612244897959184            1 0.06122449                 1
#> Marvellous:0.073469387755102             1 0.07346939                 1
#> Marvellous:0.0857142857142857            1 0.08571429                 1
#> Marvellous:0.0979591836734694            1 0.09795918                 1
#> Marvellous:0.110204081632653             1 0.11020408                 1
#> Marvellous:0.122448979591837             1 0.12244898                 1
#> Marvellous:0.13469387755102              1 0.13469388                 1
#> Marvellous:0.146938775510204             1 0.14693878                 1
#> Marvellous:0.159183673469388             1 0.15918367                 1
#> Marvellous:0.171428571428571             1 0.17142857                 1
#> Marvellous:0.183673469387755             1 0.18367347                 1
#> Marvellous:0.195918367346939             1 0.19591837                 1
#> Marvellous:0.2                           1 0.20000000                 1
#> Marvellous:0.208163265306122             1 0.20816327                 1
#> Marvellous:0.220408163265306             1 0.22040816                 1
#> Marvellous:0.23265306122449              1 0.23265306                 1
#> Marvellous:0.244897959183673             1 0.24489796                 1
#> Marvellous:0.257142857142857             1 0.25714286                 1
#> Marvellous:0.269387755102041             1 0.26938776                 1
#> Marvellous:0.281632653061224             1 0.28163265                 1
#> Marvellous:0.293877551020408             1 0.29387755                 1
#> Marvellous:0.306122448979592             1 0.30612245                 1
#> Marvellous:0.318367346938775             1 0.31836735                 1
#> Marvellous:0.330612244897959             1 0.33061224                 1
#> Marvellous:0.342857142857143             1 0.34285714                 1
#> Marvellous:0.355102040816327             1 0.35510204                 1
#> Marvellous:0.36734693877551              1 0.36734694                 1
#> Marvellous:0.379591836734694             1 0.37959184                 1
#> Marvellous:0.391836734693878             1 0.39183673                 1
#> Marvellous:0.4                           1 0.40000000                 1
#> Marvellous:0.404081632653061             1 0.40408163                 1
#> Marvellous:0.416326530612245             1 0.41632653                 1
#> Marvellous:0.428571428571429             1 0.42857143                 1
#> Marvellous:0.440816326530612             1 0.44081633                 1
#> Marvellous:0.453061224489796             1 0.45306122                 1
#> Marvellous:0.46530612244898              1 0.46530612                 1
#> Marvellous:0.477551020408163             1 0.47755102                 1
#> Marvellous:0.489795918367347             1 0.48979592                 1
#> Marvellous:0.502040816326531             1 0.50204082                 1
#> Marvellous:0.514285714285714             1 0.51428571                 1
#> Marvellous:0.526530612244898             1 0.52653061                 1
#> Marvellous:0.538775510204082             1 0.53877551                 1
#> Marvellous:0.551020408163265             1 0.55102041                 1
#> Marvellous:0.563265306122449             1 0.56326531                 1
#> Marvellous:0.575510204081633             1 0.57551020                 1
#> Marvellous:0.587755102040816             1 0.58775510                 1
#> Marvellous:0.6                           1 0.60000000                 1
#> Victory:0                                1 0.00000000                 0
#> Victory:0.0122448979591837               1 0.01224490                 0
#> Victory:0.0244897959183673               1 0.02448980                 0
#> Victory:0.036734693877551                1 0.03673469                 0
#> Victory:0.0489795918367347               1 0.04897959                 0
#> Victory:0.0612244897959184               1 0.06122449                 0
#> Victory:0.073469387755102                1 0.07346939                 0
#> Victory:0.0857142857142857               1 0.08571429                 0
#> Victory:0.0979591836734694               1 0.09795918                 0
#> Victory:0.110204081632653                1 0.11020408                 0
#> Victory:0.122448979591837                1 0.12244898                 0
#> Victory:0.13469387755102                 1 0.13469388                 0
#> Victory:0.146938775510204                1 0.14693878                 0
#> Victory:0.159183673469388                1 0.15918367                 0
#> Victory:0.171428571428571                1 0.17142857                 0
#> Victory:0.183673469387755                1 0.18367347                 0
#> Victory:0.195918367346939                1 0.19591837                 0
#> Victory:0.2                              1 0.20000000                 0
#> Victory:0.208163265306122                1 0.20816327                 0
#> Victory:0.220408163265306                1 0.22040816                 0
#> Victory:0.23265306122449                 1 0.23265306                 0
#> Victory:0.244897959183673                1 0.24489796                 0
#> Victory:0.257142857142857                1 0.25714286                 0
#> Victory:0.269387755102041                1 0.26938776                 0
#> Victory:0.281632653061224                1 0.28163265                 0
#> Victory:0.293877551020408                1 0.29387755                 0
#> Victory:0.306122448979592                1 0.30612245                 0
#> Victory:0.318367346938775                1 0.31836735                 0
#> Victory:0.330612244897959                1 0.33061224                 0
#> Victory:0.342857142857143                1 0.34285714                 0
#> Victory:0.355102040816327                1 0.35510204                 0
#> Victory:0.36734693877551                 1 0.36734694                 0
#> Victory:0.379591836734694                1 0.37959184                 0
#> Victory:0.391836734693878                1 0.39183673                 0
#> Victory:0.4                              1 0.40000000                 0
#> Victory:0.404081632653061                1 0.40408163                 0
#> Victory:0.416326530612245                1 0.41632653                 0
#> Victory:0.428571428571429                1 0.42857143                 0
#> Victory:0.440816326530612                1 0.44081633                 0
#> Victory:0.453061224489796                1 0.45306122                 0
#> Victory:0.46530612244898                 1 0.46530612                 0
#> Victory:0.477551020408163                1 0.47755102                 0
#> Victory:0.489795918367347                1 0.48979592                 0
#> Victory:0.502040816326531                1 0.50204082                 0
#> Victory:0.514285714285714                1 0.51428571                 0
#> Victory:0.526530612244898                1 0.52653061                 0
#> Victory:0.538775510204082                1 0.53877551                 0
#> Victory:0.551020408163265                1 0.55102041                 0
#> Victory:0.563265306122449                1 0.56326531                 0
#> Victory:0.575510204081633                1 0.57551020                 0
#> Victory:0.587755102040816                1 0.58775510                 0
#> Victory:0.6                              1 0.60000000                 0
#>                                VarietyVictory nitro:VarietyMarvellous
#> Golden Rain:0                               0              0.00000000
#> Golden Rain:0.0122448979591837              0              0.00000000
#> Golden Rain:0.0244897959183673              0              0.00000000
#> Golden Rain:0.036734693877551               0              0.00000000
#> Golden Rain:0.0489795918367347              0              0.00000000
#> Golden Rain:0.0612244897959184              0              0.00000000
#> Golden Rain:0.073469387755102               0              0.00000000
#> Golden Rain:0.0857142857142857              0              0.00000000
#> Golden Rain:0.0979591836734694              0              0.00000000
#> Golden Rain:0.110204081632653               0              0.00000000
#> Golden Rain:0.122448979591837               0              0.00000000
#> Golden Rain:0.13469387755102                0              0.00000000
#> Golden Rain:0.146938775510204               0              0.00000000
#> Golden Rain:0.159183673469388               0              0.00000000
#> Golden Rain:0.171428571428571               0              0.00000000
#> Golden Rain:0.183673469387755               0              0.00000000
#> Golden Rain:0.195918367346939               0              0.00000000
#> Golden Rain:0.2                             0              0.00000000
#> Golden Rain:0.208163265306122               0              0.00000000
#> Golden Rain:0.220408163265306               0              0.00000000
#> Golden Rain:0.23265306122449                0              0.00000000
#> Golden Rain:0.244897959183673               0              0.00000000
#> Golden Rain:0.257142857142857               0              0.00000000
#> Golden Rain:0.269387755102041               0              0.00000000
#> Golden Rain:0.281632653061224               0              0.00000000
#> Golden Rain:0.293877551020408               0              0.00000000
#> Golden Rain:0.306122448979592               0              0.00000000
#> Golden Rain:0.318367346938775               0              0.00000000
#> Golden Rain:0.330612244897959               0              0.00000000
#> Golden Rain:0.342857142857143               0              0.00000000
#> Golden Rain:0.355102040816327               0              0.00000000
#> Golden Rain:0.36734693877551                0              0.00000000
#> Golden Rain:0.379591836734694               0              0.00000000
#> Golden Rain:0.391836734693878               0              0.00000000
#> Golden Rain:0.4                             0              0.00000000
#> Golden Rain:0.404081632653061               0              0.00000000
#> Golden Rain:0.416326530612245               0              0.00000000
#> Golden Rain:0.428571428571429               0              0.00000000
#> Golden Rain:0.440816326530612               0              0.00000000
#> Golden Rain:0.453061224489796               0              0.00000000
#> Golden Rain:0.46530612244898                0              0.00000000
#> Golden Rain:0.477551020408163               0              0.00000000
#> Golden Rain:0.489795918367347               0              0.00000000
#> Golden Rain:0.502040816326531               0              0.00000000
#> Golden Rain:0.514285714285714               0              0.00000000
#> Golden Rain:0.526530612244898               0              0.00000000
#> Golden Rain:0.538775510204082               0              0.00000000
#> Golden Rain:0.551020408163265               0              0.00000000
#> Golden Rain:0.563265306122449               0              0.00000000
#> Golden Rain:0.575510204081633               0              0.00000000
#> Golden Rain:0.587755102040816               0              0.00000000
#> Golden Rain:0.6                             0              0.00000000
#> Marvellous:0                                0              0.00000000
#> Marvellous:0.0122448979591837               0              0.01224490
#> Marvellous:0.0244897959183673               0              0.02448980
#> Marvellous:0.036734693877551                0              0.03673469
#> Marvellous:0.0489795918367347               0              0.04897959
#> Marvellous:0.0612244897959184               0              0.06122449
#> Marvellous:0.073469387755102                0              0.07346939
#> Marvellous:0.0857142857142857               0              0.08571429
#> Marvellous:0.0979591836734694               0              0.09795918
#> Marvellous:0.110204081632653                0              0.11020408
#> Marvellous:0.122448979591837                0              0.12244898
#> Marvellous:0.13469387755102                 0              0.13469388
#> Marvellous:0.146938775510204                0              0.14693878
#> Marvellous:0.159183673469388                0              0.15918367
#> Marvellous:0.171428571428571                0              0.17142857
#> Marvellous:0.183673469387755                0              0.18367347
#> Marvellous:0.195918367346939                0              0.19591837
#> Marvellous:0.2                              0              0.20000000
#> Marvellous:0.208163265306122                0              0.20816327
#> Marvellous:0.220408163265306                0              0.22040816
#> Marvellous:0.23265306122449                 0              0.23265306
#> Marvellous:0.244897959183673                0              0.24489796
#> Marvellous:0.257142857142857                0              0.25714286
#> Marvellous:0.269387755102041                0              0.26938776
#> Marvellous:0.281632653061224                0              0.28163265
#> Marvellous:0.293877551020408                0              0.29387755
#> Marvellous:0.306122448979592                0              0.30612245
#> Marvellous:0.318367346938775                0              0.31836735
#> Marvellous:0.330612244897959                0              0.33061224
#> Marvellous:0.342857142857143                0              0.34285714
#> Marvellous:0.355102040816327                0              0.35510204
#> Marvellous:0.36734693877551                 0              0.36734694
#> Marvellous:0.379591836734694                0              0.37959184
#> Marvellous:0.391836734693878                0              0.39183673
#> Marvellous:0.4                              0              0.40000000
#> Marvellous:0.404081632653061                0              0.40408163
#> Marvellous:0.416326530612245                0              0.41632653
#> Marvellous:0.428571428571429                0              0.42857143
#> Marvellous:0.440816326530612                0              0.44081633
#> Marvellous:0.453061224489796                0              0.45306122
#> Marvellous:0.46530612244898                 0              0.46530612
#> Marvellous:0.477551020408163                0              0.47755102
#> Marvellous:0.489795918367347                0              0.48979592
#> Marvellous:0.502040816326531                0              0.50204082
#> Marvellous:0.514285714285714                0              0.51428571
#> Marvellous:0.526530612244898                0              0.52653061
#> Marvellous:0.538775510204082                0              0.53877551
#> Marvellous:0.551020408163265                0              0.55102041
#> Marvellous:0.563265306122449                0              0.56326531
#> Marvellous:0.575510204081633                0              0.57551020
#> Marvellous:0.587755102040816                0              0.58775510
#> Marvellous:0.6                              0              0.60000000
#> Victory:0                                   1              0.00000000
#> Victory:0.0122448979591837                  1              0.00000000
#> Victory:0.0244897959183673                  1              0.00000000
#> Victory:0.036734693877551                   1              0.00000000
#> Victory:0.0489795918367347                  1              0.00000000
#> Victory:0.0612244897959184                  1              0.00000000
#> Victory:0.073469387755102                   1              0.00000000
#> Victory:0.0857142857142857                  1              0.00000000
#> Victory:0.0979591836734694                  1              0.00000000
#> Victory:0.110204081632653                   1              0.00000000
#> Victory:0.122448979591837                   1              0.00000000
#> Victory:0.13469387755102                    1              0.00000000
#> Victory:0.146938775510204                   1              0.00000000
#> Victory:0.159183673469388                   1              0.00000000
#> Victory:0.171428571428571                   1              0.00000000
#> Victory:0.183673469387755                   1              0.00000000
#> Victory:0.195918367346939                   1              0.00000000
#> Victory:0.2                                 1              0.00000000
#> Victory:0.208163265306122                   1              0.00000000
#> Victory:0.220408163265306                   1              0.00000000
#> Victory:0.23265306122449                    1              0.00000000
#> Victory:0.244897959183673                   1              0.00000000
#> Victory:0.257142857142857                   1              0.00000000
#> Victory:0.269387755102041                   1              0.00000000
#> Victory:0.281632653061224                   1              0.00000000
#> Victory:0.293877551020408                   1              0.00000000
#> Victory:0.306122448979592                   1              0.00000000
#> Victory:0.318367346938775                   1              0.00000000
#> Victory:0.330612244897959                   1              0.00000000
#> Victory:0.342857142857143                   1              0.00000000
#> Victory:0.355102040816327                   1              0.00000000
#> Victory:0.36734693877551                    1              0.00000000
#> Victory:0.379591836734694                   1              0.00000000
#> Victory:0.391836734693878                   1              0.00000000
#> Victory:0.4                                 1              0.00000000
#> Victory:0.404081632653061                   1              0.00000000
#> Victory:0.416326530612245                   1              0.00000000
#> Victory:0.428571428571429                   1              0.00000000
#> Victory:0.440816326530612                   1              0.00000000
#> Victory:0.453061224489796                   1              0.00000000
#> Victory:0.46530612244898                    1              0.00000000
#> Victory:0.477551020408163                   1              0.00000000
#> Victory:0.489795918367347                   1              0.00000000
#> Victory:0.502040816326531                   1              0.00000000
#> Victory:0.514285714285714                   1              0.00000000
#> Victory:0.526530612244898                   1              0.00000000
#> Victory:0.538775510204082                   1              0.00000000
#> Victory:0.551020408163265                   1              0.00000000
#> Victory:0.563265306122449                   1              0.00000000
#> Victory:0.575510204081633                   1              0.00000000
#> Victory:0.587755102040816                   1              0.00000000
#> Victory:0.6                                 1              0.00000000
#>                                nitro:VarietyVictory
#> Golden Rain:0                            0.00000000
#> Golden Rain:0.0122448979591837           0.00000000
#> Golden Rain:0.0244897959183673           0.00000000
#> Golden Rain:0.036734693877551            0.00000000
#> Golden Rain:0.0489795918367347           0.00000000
#> Golden Rain:0.0612244897959184           0.00000000
#> Golden Rain:0.073469387755102            0.00000000
#> Golden Rain:0.0857142857142857           0.00000000
#> Golden Rain:0.0979591836734694           0.00000000
#> Golden Rain:0.110204081632653            0.00000000
#> Golden Rain:0.122448979591837            0.00000000
#> Golden Rain:0.13469387755102             0.00000000
#> Golden Rain:0.146938775510204            0.00000000
#> Golden Rain:0.159183673469388            0.00000000
#> Golden Rain:0.171428571428571            0.00000000
#> Golden Rain:0.183673469387755            0.00000000
#> Golden Rain:0.195918367346939            0.00000000
#> Golden Rain:0.2                          0.00000000
#> Golden Rain:0.208163265306122            0.00000000
#> Golden Rain:0.220408163265306            0.00000000
#> Golden Rain:0.23265306122449             0.00000000
#> Golden Rain:0.244897959183673            0.00000000
#> Golden Rain:0.257142857142857            0.00000000
#> Golden Rain:0.269387755102041            0.00000000
#> Golden Rain:0.281632653061224            0.00000000
#> Golden Rain:0.293877551020408            0.00000000
#> Golden Rain:0.306122448979592            0.00000000
#> Golden Rain:0.318367346938775            0.00000000
#> Golden Rain:0.330612244897959            0.00000000
#> Golden Rain:0.342857142857143            0.00000000
#> Golden Rain:0.355102040816327            0.00000000
#> Golden Rain:0.36734693877551             0.00000000
#> Golden Rain:0.379591836734694            0.00000000
#> Golden Rain:0.391836734693878            0.00000000
#> Golden Rain:0.4                          0.00000000
#> Golden Rain:0.404081632653061            0.00000000
#> Golden Rain:0.416326530612245            0.00000000
#> Golden Rain:0.428571428571429            0.00000000
#> Golden Rain:0.440816326530612            0.00000000
#> Golden Rain:0.453061224489796            0.00000000
#> Golden Rain:0.46530612244898             0.00000000
#> Golden Rain:0.477551020408163            0.00000000
#> Golden Rain:0.489795918367347            0.00000000
#> Golden Rain:0.502040816326531            0.00000000
#> Golden Rain:0.514285714285714            0.00000000
#> Golden Rain:0.526530612244898            0.00000000
#> Golden Rain:0.538775510204082            0.00000000
#> Golden Rain:0.551020408163265            0.00000000
#> Golden Rain:0.563265306122449            0.00000000
#> Golden Rain:0.575510204081633            0.00000000
#> Golden Rain:0.587755102040816            0.00000000
#> Golden Rain:0.6                          0.00000000
#> Marvellous:0                             0.00000000
#> Marvellous:0.0122448979591837            0.00000000
#> Marvellous:0.0244897959183673            0.00000000
#> Marvellous:0.036734693877551             0.00000000
#> Marvellous:0.0489795918367347            0.00000000
#> Marvellous:0.0612244897959184            0.00000000
#> Marvellous:0.073469387755102             0.00000000
#> Marvellous:0.0857142857142857            0.00000000
#> Marvellous:0.0979591836734694            0.00000000
#> Marvellous:0.110204081632653             0.00000000
#> Marvellous:0.122448979591837             0.00000000
#> Marvellous:0.13469387755102              0.00000000
#> Marvellous:0.146938775510204             0.00000000
#> Marvellous:0.159183673469388             0.00000000
#> Marvellous:0.171428571428571             0.00000000
#> Marvellous:0.183673469387755             0.00000000
#> Marvellous:0.195918367346939             0.00000000
#> Marvellous:0.2                           0.00000000
#> Marvellous:0.208163265306122             0.00000000
#> Marvellous:0.220408163265306             0.00000000
#> Marvellous:0.23265306122449              0.00000000
#> Marvellous:0.244897959183673             0.00000000
#> Marvellous:0.257142857142857             0.00000000
#> Marvellous:0.269387755102041             0.00000000
#> Marvellous:0.281632653061224             0.00000000
#> Marvellous:0.293877551020408             0.00000000
#> Marvellous:0.306122448979592             0.00000000
#> Marvellous:0.318367346938775             0.00000000
#> Marvellous:0.330612244897959             0.00000000
#> Marvellous:0.342857142857143             0.00000000
#> Marvellous:0.355102040816327             0.00000000
#> Marvellous:0.36734693877551              0.00000000
#> Marvellous:0.379591836734694             0.00000000
#> Marvellous:0.391836734693878             0.00000000
#> Marvellous:0.4                           0.00000000
#> Marvellous:0.404081632653061             0.00000000
#> Marvellous:0.416326530612245             0.00000000
#> Marvellous:0.428571428571429             0.00000000
#> Marvellous:0.440816326530612             0.00000000
#> Marvellous:0.453061224489796             0.00000000
#> Marvellous:0.46530612244898              0.00000000
#> Marvellous:0.477551020408163             0.00000000
#> Marvellous:0.489795918367347             0.00000000
#> Marvellous:0.502040816326531             0.00000000
#> Marvellous:0.514285714285714             0.00000000
#> Marvellous:0.526530612244898             0.00000000
#> Marvellous:0.538775510204082             0.00000000
#> Marvellous:0.551020408163265             0.00000000
#> Marvellous:0.563265306122449             0.00000000
#> Marvellous:0.575510204081633             0.00000000
#> Marvellous:0.587755102040816             0.00000000
#> Marvellous:0.6                           0.00000000
#> Victory:0                                0.00000000
#> Victory:0.0122448979591837               0.01224490
#> Victory:0.0244897959183673               0.02448980
#> Victory:0.036734693877551                0.03673469
#> Victory:0.0489795918367347               0.04897959
#> Victory:0.0612244897959184               0.06122449
#> Victory:0.073469387755102                0.07346939
#> Victory:0.0857142857142857               0.08571429
#> Victory:0.0979591836734694               0.09795918
#> Victory:0.110204081632653                0.11020408
#> Victory:0.122448979591837                0.12244898
#> Victory:0.13469387755102                 0.13469388
#> Victory:0.146938775510204                0.14693878
#> Victory:0.159183673469388                0.15918367
#> Victory:0.171428571428571                0.17142857
#> Victory:0.183673469387755                0.18367347
#> Victory:0.195918367346939                0.19591837
#> Victory:0.2                              0.20000000
#> Victory:0.208163265306122                0.20816327
#> Victory:0.220408163265306                0.22040816
#> Victory:0.23265306122449                 0.23265306
#> Victory:0.244897959183673                0.24489796
#> Victory:0.257142857142857                0.25714286
#> Victory:0.269387755102041                0.26938776
#> Victory:0.281632653061224                0.28163265
#> Victory:0.293877551020408                0.29387755
#> Victory:0.306122448979592                0.30612245
#> Victory:0.318367346938775                0.31836735
#> Victory:0.330612244897959                0.33061224
#> Victory:0.342857142857143                0.34285714
#> Victory:0.355102040816327                0.35510204
#> Victory:0.36734693877551                 0.36734694
#> Victory:0.379591836734694                0.37959184
#> Victory:0.391836734693878                0.39183673
#> Victory:0.4                              0.40000000
#> Victory:0.404081632653061                0.40408163
#> Victory:0.416326530612245                0.41632653
#> Victory:0.428571428571429                0.42857143
#> Victory:0.440816326530612                0.44081633
#> Victory:0.453061224489796                0.45306122
#> Victory:0.46530612244898                 0.46530612
#> Victory:0.477551020408163                0.47755102
#> Victory:0.489795918367347                0.48979592
#> Victory:0.502040816326531                0.50204082
#> Victory:0.514285714285714                0.51428571
#> Victory:0.526530612244898                0.52653061
#> Victory:0.538775510204082                0.53877551
#> Victory:0.551020408163265                0.55102041
#> Victory:0.563265306122449                0.56326531
#> Victory:0.575510204081633                0.57551020
#> Victory:0.587755102040816                0.58775510
#> Victory:0.6                              0.60000000
#> attr(,"estimable")
#>                  Golden Rain:0 Golden Rain:0.0122448979591837 
#>                           TRUE                           TRUE 
#> Golden Rain:0.0244897959183673  Golden Rain:0.036734693877551 
#>                           TRUE                           TRUE 
#> Golden Rain:0.0489795918367347 Golden Rain:0.0612244897959184 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.073469387755102 Golden Rain:0.0857142857142857 
#>                           TRUE                           TRUE 
#> Golden Rain:0.0979591836734694  Golden Rain:0.110204081632653 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.122448979591837   Golden Rain:0.13469387755102 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.146938775510204  Golden Rain:0.159183673469388 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.171428571428571  Golden Rain:0.183673469387755 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.195918367346939                Golden Rain:0.2 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.208163265306122  Golden Rain:0.220408163265306 
#>                           TRUE                           TRUE 
#>   Golden Rain:0.23265306122449  Golden Rain:0.244897959183673 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.257142857142857  Golden Rain:0.269387755102041 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.281632653061224  Golden Rain:0.293877551020408 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.306122448979592  Golden Rain:0.318367346938775 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.330612244897959  Golden Rain:0.342857142857143 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.355102040816327   Golden Rain:0.36734693877551 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.379591836734694  Golden Rain:0.391836734693878 
#>                           TRUE                           TRUE 
#>                Golden Rain:0.4  Golden Rain:0.404081632653061 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.416326530612245  Golden Rain:0.428571428571429 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.440816326530612  Golden Rain:0.453061224489796 
#>                           TRUE                           TRUE 
#>   Golden Rain:0.46530612244898  Golden Rain:0.477551020408163 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.489795918367347  Golden Rain:0.502040816326531 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.514285714285714  Golden Rain:0.526530612244898 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.538775510204082  Golden Rain:0.551020408163265 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.563265306122449  Golden Rain:0.575510204081633 
#>                           TRUE                           TRUE 
#>  Golden Rain:0.587755102040816                Golden Rain:0.6 
#>                           TRUE                           TRUE 
#>                   Marvellous:0  Marvellous:0.0122448979591837 
#>                           TRUE                           TRUE 
#>  Marvellous:0.0244897959183673   Marvellous:0.036734693877551 
#>                           TRUE                           TRUE 
#>  Marvellous:0.0489795918367347  Marvellous:0.0612244897959184 
#>                           TRUE                           TRUE 
#>   Marvellous:0.073469387755102  Marvellous:0.0857142857142857 
#>                           TRUE                           TRUE 
#>  Marvellous:0.0979591836734694   Marvellous:0.110204081632653 
#>                           TRUE                           TRUE 
#>   Marvellous:0.122448979591837    Marvellous:0.13469387755102 
#>                           TRUE                           TRUE 
#>   Marvellous:0.146938775510204   Marvellous:0.159183673469388 
#>                           TRUE                           TRUE 
#>   Marvellous:0.171428571428571   Marvellous:0.183673469387755 
#>                           TRUE                           TRUE 
#>   Marvellous:0.195918367346939                 Marvellous:0.2 
#>                           TRUE                           TRUE 
#>   Marvellous:0.208163265306122   Marvellous:0.220408163265306 
#>                           TRUE                           TRUE 
#>    Marvellous:0.23265306122449   Marvellous:0.244897959183673 
#>                           TRUE                           TRUE 
#>   Marvellous:0.257142857142857   Marvellous:0.269387755102041 
#>                           TRUE                           TRUE 
#>   Marvellous:0.281632653061224   Marvellous:0.293877551020408 
#>                           TRUE                           TRUE 
#>   Marvellous:0.306122448979592   Marvellous:0.318367346938775 
#>                           TRUE                           TRUE 
#>   Marvellous:0.330612244897959   Marvellous:0.342857142857143 
#>                           TRUE                           TRUE 
#>   Marvellous:0.355102040816327    Marvellous:0.36734693877551 
#>                           TRUE                           TRUE 
#>   Marvellous:0.379591836734694   Marvellous:0.391836734693878 
#>                           TRUE                           TRUE 
#>                 Marvellous:0.4   Marvellous:0.404081632653061 
#>                           TRUE                           TRUE 
#>   Marvellous:0.416326530612245   Marvellous:0.428571428571429 
#>                           TRUE                           TRUE 
#>   Marvellous:0.440816326530612   Marvellous:0.453061224489796 
#>                           TRUE                           TRUE 
#>    Marvellous:0.46530612244898   Marvellous:0.477551020408163 
#>                           TRUE                           TRUE 
#>   Marvellous:0.489795918367347   Marvellous:0.502040816326531 
#>                           TRUE                           TRUE 
#>   Marvellous:0.514285714285714   Marvellous:0.526530612244898 
#>                           TRUE                           TRUE 
#>   Marvellous:0.538775510204082   Marvellous:0.551020408163265 
#>                           TRUE                           TRUE 
#>   Marvellous:0.563265306122449   Marvellous:0.575510204081633 
#>                           TRUE                           TRUE 
#>   Marvellous:0.587755102040816                 Marvellous:0.6 
#>                           TRUE                           TRUE 
#>                      Victory:0     Victory:0.0122448979591837 
#>                           TRUE                           TRUE 
#>     Victory:0.0244897959183673      Victory:0.036734693877551 
#>                           TRUE                           TRUE 
#>     Victory:0.0489795918367347     Victory:0.0612244897959184 
#>                           TRUE                           TRUE 
#>      Victory:0.073469387755102     Victory:0.0857142857142857 
#>                           TRUE                           TRUE 
#>     Victory:0.0979591836734694      Victory:0.110204081632653 
#>                           TRUE                           TRUE 
#>      Victory:0.122448979591837       Victory:0.13469387755102 
#>                           TRUE                           TRUE 
#>      Victory:0.146938775510204      Victory:0.159183673469388 
#>                           TRUE                           TRUE 
#>      Victory:0.171428571428571      Victory:0.183673469387755 
#>                           TRUE                           TRUE 
#>      Victory:0.195918367346939                    Victory:0.2 
#>                           TRUE                           TRUE 
#>      Victory:0.208163265306122      Victory:0.220408163265306 
#>                           TRUE                           TRUE 
#>       Victory:0.23265306122449      Victory:0.244897959183673 
#>                           TRUE                           TRUE 
#>      Victory:0.257142857142857      Victory:0.269387755102041 
#>                           TRUE                           TRUE 
#>      Victory:0.281632653061224      Victory:0.293877551020408 
#>                           TRUE                           TRUE 
#>      Victory:0.306122448979592      Victory:0.318367346938775 
#>                           TRUE                           TRUE 
#>      Victory:0.330612244897959      Victory:0.342857142857143 
#>                           TRUE                           TRUE 
#>      Victory:0.355102040816327       Victory:0.36734693877551 
#>                           TRUE                           TRUE 
#>      Victory:0.379591836734694      Victory:0.391836734693878 
#>                           TRUE                           TRUE 
#>                    Victory:0.4      Victory:0.404081632653061 
#>                           TRUE                           TRUE 
#>      Victory:0.416326530612245      Victory:0.428571428571429 
#>                           TRUE                           TRUE 
#>      Victory:0.440816326530612      Victory:0.453061224489796 
#>                           TRUE                           TRUE 
#>       Victory:0.46530612244898      Victory:0.477551020408163 
#>                           TRUE                           TRUE 
#>      Victory:0.489795918367347      Victory:0.502040816326531 
#>                           TRUE                           TRUE 
#>      Victory:0.514285714285714      Victory:0.526530612244898 
#>                           TRUE                           TRUE 
#>      Victory:0.538775510204082      Victory:0.551020408163265 
#>                           TRUE                           TRUE 
#>      Victory:0.563265306122449      Victory:0.575510204081633 
#>                           TRUE                           TRUE 
#>      Victory:0.587755102040816                    Victory:0.6 
#>                           TRUE                           TRUE 
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
#>             (Intercept)       Block.L       Block.Q       Block.C       Block^4
#> Golden Rain           1 -6.938894e-17 -2.775558e-17 -6.938894e-18 -2.081668e-17
#> Marvellous            1 -6.938894e-17 -2.775558e-17 -6.938894e-18 -2.081668e-17
#> Victory               1 -6.938894e-17 -2.775558e-17 -6.938894e-18 -2.081668e-17
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
#> Golden Rain           1                 0              0 -6.938894e-17
#> Marvellous            1                 1              0 -6.938894e-17
#> Victory               1                 0              1 -6.938894e-17
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
