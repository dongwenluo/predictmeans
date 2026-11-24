# Standard Errors of the Predicted Random Effects

The predicted random effects (BLUPs, conditional modes) of a mixed model
with their standard errors and limits, and optionally the predicted
value of each level (fixed part plus random effect).

Two standard errors are available:

- `type = "prediction"`:

  (the default) the prediction error SE, \\\sqrt{Var(\hat u - u)}\\,
  from the inverse of the mixed model equations (Henderson's
  \\C\_{22}\\; SAS PROC MIXED `RANDOM / SOLUTION`, "Std Err Pred"). It
  includes the uncertainty of the estimated fixed effects, which matters
  most for a term that contains a fixed factor (e.g. `Block:Variety`
  with `Variety` fixed).

- `type = "conditional"`:

  the conditional SE, \\\sqrt{Var(u \mid y)}\\ with the fixed effects
  treated as known: lme4's `ranef(condVar = TRUE)`, the numbers of
  `se_ranef()` in predictmeans 1.x. It is smaller.

Both treat the variance parameters as known (the usual plug-in), except
for `glmmTMB`: there `"conditional"` is `ranef(condVar = TRUE)` of
glmmTMB, which is TMB's `sdreport()` and includes the uncertainty of the
fixed effects *and* of the variance parameters (so it can exceed
`"prediction"`, which inverts the \\(\beta, b)\\ block of TMB's joint
precision, variance parameters known). Correlated terms, nested and
crossed terms, glmmTMB structures (`us`, `cs`, `ar1`, ...) and `lme`
correlation and variance structures are all covered, as the SEs come
from the joint covariance of all random effects. For GLMMs, all values
are on the link scale.

With `predict = TRUE`, each level gets its predicted value: for a random
term with no fixed factor in it (`Block`), the overall mean (the fixed
factors averaged with equal weights, as
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md);
covariates at their means) plus \\\hat u\\; for a term that contains
fixed factors (`Block:Variety`, `Variety` fixed), the predicted mean of
that fixed level plus \\\hat u\\. Its SE is from the joint covariance of
\\(\hat\beta, \hat u - u)\\, \\L C L'\\ with \\L = \[k', e_i\]\\
(`"prediction"` always). Only intercept rows get a prediction (random
slopes: `NA`). For GLMMs the response scale is added by the inverse link
(SE by the delta method).

## Usage

``` r
se_ranef(
  object,
  rand_term = NULL,
  type = c("prediction", "conditional"),
  predict = FALSE,
  conf_level = 0.95
)

# S3 method for class 'pm_ranef'
print(x, n = 10, digits = 4, ...)
```

## Arguments

- object:

  An `lmer`, `glmer`, `glmmTMB` (conditional model), `lme` or
  [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fit.

- rand_term:

  Optional name of one grouping factor, e.g. `"Block"` or
  `"Variety:Block"` (`lme`: the level name, e.g. `"Variety"` for
  `~ 1 | Block/Variety`).

- type:

  `"prediction"` (default) or `"conditional"` (1.x); see Description.

- predict:

  Add the predicted value of each level (see Description)?

- conf_level:

  Confidence level of the limits.

- x:

  A `pm_ranef` result.

- n:

  Rows printed per group.

- digits:

  Digits printed.

- ...:

  Not used.

## Value

A data frame of class `pm_ranef`: `group`, `level`, `term`, `estimate`
(the BLUP), `std.error`, `conf.low`, `conf.high` (normal limits \\\hat u
\pm z\\SE\\: for `"prediction"` a prediction interval of \\u\\); with
`predict = TRUE` also `prediction`, `prediction.std.error`,
`prediction.conf.low`, `prediction.conf.high`, and for GLMMs `response`,
`response.std.error`, `response.conf.low`, `response.conf.high`.
`attr(, "settings")` records `type`, `conf_level`, `limits`,
`model_class`, `scale` and `conditional` (what the conditional SE is for
this class). [`as.list()`](https://rdrr.io/r/base/list.html) gives the
1.x result: one matrix of SEs (levels x terms) per group, or that matrix
only when `rand_term` was given.

## References

Henderson, C. R. (1975). Best linear unbiased estimation and prediction
under a selection model. *Biometrics*, 31, 423-447.

## See also

[`varcomp`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
[`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
se_ranef(fm, "Block")
#> Predicted random effects (BLUPs), prediction error SE; 95% normal prediction intervals of u
#> 
#> Block (6 levels)
#>  level        term estimate std.error conf.low conf.high
#>     VI (Intercept)   -6.260     8.342  -22.611    10.091
#>      V (Intercept)  -10.583     8.342  -26.934     5.768
#>    III (Intercept)   -6.530     8.342  -22.881     9.821
#>     IV (Intercept)   -4.706     8.342  -21.057    11.645
#>     II (Intercept)    2.657     8.342  -13.694    19.008
#>      I (Intercept)   25.422     8.342    9.071    41.773
se_ranef(fm, "Block", type = "conditional")     # 1.x numbers
#> Predicted random effects (BLUPs), conditional SE; 95% normal intervals from the conditional SE
#> 
#> Block (6 levels)
#>  level        term estimate std.error conf.low conf.high
#>     VI (Intercept)   -6.260     6.373  -18.751     6.232
#>      V (Intercept)  -10.583     6.373  -23.075     1.909
#>    III (Intercept)   -6.530     6.373  -19.022     5.962
#>     IV (Intercept)   -4.706     6.373  -17.198     7.786
#>     II (Intercept)    2.657     6.373   -9.835    15.149
#>      I (Intercept)   25.422     6.373   12.930    37.913
se_ranef(fm, "Variety:Block", predict = TRUE)   # cell values
#> Predicted random effects (BLUPs), prediction error SE; 95% normal prediction intervals of u
#> 
#> Variety:Block (18 levels)
#>            level        term estimate std.error conf.low conf.high prediction
#>   Golden Rain:VI (Intercept)   -5.637     7.782  -20.889     9.614      98.86
#>    Golden Rain:V (Intercept)    1.117     7.782  -14.135    16.369     105.62
#>  Golden Rain:III (Intercept)   -7.916     7.782  -23.168     7.336      96.58
#>   Golden Rain:IV (Intercept)    5.789     7.782   -9.462    21.041     110.29
#>   Golden Rain:II (Intercept)    4.299     7.782  -10.953    19.550     108.80
#>    Golden Rain:I (Intercept)    2.348     7.782  -12.904    17.600     106.85
#>    Marvellous:VI (Intercept)    3.858     7.782  -11.394    19.110     113.65
#>     Marvellous:V (Intercept)   -9.848     7.782  -25.100     5.404      99.94
#>   Marvellous:III (Intercept)   10.751     7.782   -4.501    26.003     120.54
#>    Marvellous:IV (Intercept)   -7.116     7.782  -22.367     8.136     102.68
#>  prediction.std.error prediction.conf.low prediction.conf.high
#>                 9.274               80.69                117.0
#>                 9.274               87.44                123.8
#>                 9.274               78.41                114.8
#>                 9.274               92.11                128.5
#>                 9.274               90.62                127.0
#>                 9.274               88.67                125.0
#>                 9.274               95.47                131.8
#>                 9.274               81.77                118.1
#>                 9.274              102.37                138.7
#>                 9.274               84.50                120.9
#> ... 8 more rows: as.data.frame(x) gives the table.
as.list(se_ranef(fm, type = "conditional"))     # 1.x list of matrices
#> $`Variety:Block`
#>                 (Intercept)
#> Golden Rain:VI     7.173039
#> Golden Rain:V      7.173039
#> Golden Rain:III    7.173039
#> Golden Rain:IV     7.173039
#> Golden Rain:II     7.173039
#> Golden Rain:I      7.173039
#> Marvellous:VI      7.173039
#> Marvellous:V       7.173039
#> Marvellous:III     7.173039
#> Marvellous:IV      7.173039
#> Marvellous:II      7.173039
#> Marvellous:I       7.173039
#> Victory:VI         7.173039
#> Victory:V          7.173039
#> Victory:III        7.173039
#> Victory:IV         7.173039
#> Victory:II         7.173039
#> Victory:I          7.173039
#> 
#> $Block
#>     (Intercept)
#> VI     6.373413
#> V      6.373413
#> III    6.373413
#> IV     6.373413
#> II     6.373413
#> I      6.373413
#> 
```
