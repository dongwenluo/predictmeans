# Predicted Means of a Linear Model with Covariate Variable(s)

This function obtains predicted means with graph for a new set of
covariate values.

## Usage

``` r
covariatemeans(model, modelterm=NULL, covariate, as.is=FALSE, covariateV=NULL, 
  data=NULL, level=0.05, Df=NULL, trans=NULL, transOff=0, responsen=NULL, trellis=TRUE, 
  plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=FALSE,
  cov_reduce=mean, exclude=NULL,
  weights=c("equal", "proportional", "outer", "cells", "flat"))
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

- modelterm:

  Name (in "quotes") for indicating which factor term's predicted mean
  to be calculated. The `modelterm` must be given exactly as it appears
  in the printed model, e.g. "A" or "A:B".

- covariate:

  Name (in "quotes") of one the covariate variables in the `model`.

- as.is:

  A logic value to specify wheather or not using original covariate
  values' rage for graph, the default is FALSE.

- covariateV:

  A numeric vector when as.is is FALSE, then covariatemeans will produce
  the result for `covariate` at value of `covariateV`.

- data:

  In some cases, you need to provide the data set used in model fitting,
  especially when you have applied some variable trnasformation in the
  model.

- level:

  A significant level for calculating confident interval. The default
  value is 0.05.

- Df:

  A degree of freedom for calculating CI of predicted means (you can
  manually specified ddf here). For the above models, ddf is obtained
  from the function automatically.

- trans:

  A function object for calculating the back transformed means, e.g.
  `trans=exp`.

- transOff:

  When you use `trans=exp(x+1)`, then `transOff=1`, the default is 0.

- responsen:

  Name (in "quotes") of the back transformed response variable in the
  `model`.

- trellis:

  A logical variable. If set to TRUE (default), a trellis plots of
  predicted means with CI will be drawn.

- plotord:

  A numeric (or character) vector specifying the order of plotting for
  two or three way interaction (e.g. `plotord = c(2, 1, 3)` or
  `plotord = c("B", "A", "C")` will put the second variable 'B' in
  `modelterm` on the `X` axis, the first variable 'A' as the grouping
  variable, and the third one 'C' as the panel variable). The defaults
  are `c(1, 2)` and `c(1, 2, 3)` for two and three way interactions.

- mtitle:

  The main title in the graph.

- ci:

  A logical variable to indicate whether to print confidence interval.
  The default value is TRUE.

- point:

  A logical variable to indicate whether to print raw data points. The
  default value is TRUE.

- jitterv:

  A degree of jitter in x and y direction in the graph. The default is
  zero.

- newwd:

  A logical variable to indicate whether to print graph in a new window.
  The default value is FALSE.

- cov_reduce:

  How covariates that are not varied are handled when building the
  predicted means: a function such as `mean` (the default) or `median`,
  or `"average"` to average the predictions over the observed covariate
  values (e.g. means adjusted for a spatial smooth in a `gam`). See
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).

- exclude:

  For `gam` models: smooth terms to leave out of the predictions. The
  default excludes random-effect smooths. See
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md).

- weights:

  Passed to
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md):
  how the cells of the factors that are not in `modelterm` are weighted
  when averaging over them: `"equal"` (the default; equal weights),
  `"proportional"` (weights proportional to the frequencies of those
  cells in the data, the same for every mean, as SAS `LSMEANS / OM`),
  `"outer"` (products of each factor's own frequencies), `"cells"` (each
  mean weights the cells by their own frequencies) or `"flat"` (equal
  weights over the non-empty cells). The definitions are those of
  emmeans.

## Value

- Predicted Means:

  A table of predicted means.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  data(Oats, package="nlme")
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  covariatemeans(fm, "Variety", covariate="nitro")

  covariatemeans(fm, "Variety", covariate="nitro", covariateV=seq(0, 0.6, 0.1))$pltdf

#>        Variety xvar      ses      Mean        LL        UL     factors
#> 1  Golden Rain  0.0 8.570716  81.90000  62.80325 100.99675 Golden Rain
#> 2  Golden Rain  0.1 8.150210  89.43333  71.27353 107.59313 Golden Rain
#> 3  Golden Rain  0.2 7.887152  96.96667  79.39300 114.54034 Golden Rain
#> 4  Golden Rain  0.3 7.797494 104.50000  87.12610 121.87390 Golden Rain
#> 5  Golden Rain  0.4 7.887152 112.03333  94.45966 129.60700 Golden Rain
#> 6  Golden Rain  0.5 8.150210 119.56667 101.40687 137.72647 Golden Rain
#> 7  Golden Rain  0.6 8.570716 127.10000 108.00325 146.19675 Golden Rain
#> 8   Marvellous  0.0 8.570716  90.41667  71.31992 109.51341  Marvellous
#> 9   Marvellous  0.1 8.150210  96.87500  78.71520 115.03480  Marvellous
#> 10  Marvellous  0.2 7.887152 103.33333  85.75966 120.90700  Marvellous
#> 11  Marvellous  0.3 7.797494 109.79167  92.41777 127.16557  Marvellous
#> 12  Marvellous  0.4 7.887152 116.25000  98.67633 133.82367  Marvellous
#> 13  Marvellous  0.5 8.150210 122.70833 104.54853 140.86813  Marvellous
#> 14  Marvellous  0.6 8.570716 129.16667 110.06992 148.26341  Marvellous
#> 15     Victory  0.0 8.570716  73.30000  54.20325  92.39675     Victory
#> 16     Victory  0.1 8.150210  81.40833  63.24853  99.56813     Victory
#> 17     Victory  0.2 7.887152  89.51667  71.94300 107.09034     Victory
#> 18     Victory  0.3 7.797494  97.62500  80.25110 114.99890     Victory
#> 19     Victory  0.4 7.887152 105.73333  88.15966 123.30700     Victory
#> 20     Victory  0.5 8.150210 113.84167  95.68187 132.00147     Victory
#> 21     Victory  0.6 8.570716 121.95000 102.85325 141.04675     Victory
```
