# Supporting a New Model Class

[`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
and
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
get everything they need from a model through these seven generics.
Their default methods work for any model class with the usual methods
(`coef`, `vcov`, `terms`, `df.residual`, `family`), so many classes need
no extra code. For a class that lacks one of these, or needs something
special, write a method for the corresponding generic only; nothing
inside predictmeans needs to change.

## Usage

``` r
pm_coef(model, ...)
pm_vcov(model, ...)
pm_basis(model, newdata, ...)
pm_df(model, L, ...)
pm_linkinv(model, ...)
pm_contrasts(model, ...)
```

## Arguments

- model:

  A fitted model.

- newdata:

  A data frame of predictor values (the reference grid).

- L:

  A matrix of linear functions of the coefficients, one row per contrast
  or mean.

- ...:

  Further arguments. Kmatrix() passes `terms`, `contrasts` and `xlev` to
  `pm_basis()`, and `exclude` for `gam` models.

## Value

As described in Details.

## Details

- `pm_coef(model)`:

  Named fixed-effect estimates, with `NA` for aliased coefficients.
  Default: `coef(model)`; methods for `lme`, `merMod` and `glmmTMB`
  return the fixed effects.

- `pm_vcov(model)`:

  Covariance matrix of the non-aliased estimates. Default:
  `vcov(model)`.

- `pm_basis(model, newdata, ...)`:

  Design matrix of the fixed effects for `newdata`, with columns named
  as in `pm_coef(model)`: the columns are matched to the coefficients by
  name, and a coefficient without a column stops with a message (a
  column for an aliased coefficient may be there or not). Default:
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) of the
  model's terms (with their stored `predvars`, so
  [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`,
  [`scale()`](https://rdrr.io/r/base/scale.html) etc. use the fitted
  basis), contrasts and factor levels. The `gam` method uses the linear
  predictor matrix of package mgcv.

- `pm_data(model)`:

  The data used to fit the model. The default method (see
  [`pm_data`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md))
  finds them from the `data` argument of the call,
  [`getData()`](https://rdrr.io/pkg/nlme/man/getData.html) or the model
  frame;
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  uses a
  [`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md)
  method only when the class has its own.

- `pm_df(model, L)`:

  Denominator degrees of freedom for each row of `L`, or `NULL` (the
  default) to use the model's single residual df: `df.residual(model)`,
  or `Inf` (normal-based) when there is none. Kenward-Roger df for
  `lmer` and the rules of nlme for `lme` are built in.

- `pm_linkinv(model)`:

  Inverse link function used to back-transform the means, or `NULL` for
  none. Default: the inverse link of `family(model)` unless it is the
  identity.

- `pm_contrasts(model)`:

  The contrasts of the factors used in the fit, as for the
  `contrasts.arg` argument of
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  Default: those of `model.matrix(model)`, else `model$contrasts`; a
  method is only needed for a class that has neither.

The model also needs a [`terms()`](https://rdrr.io/r/stats/terms.html)
method for its fixed effects (used to find the variables in the model),
and its call ([`getCall()`](https://rdrr.io/r/stats/update.html)) or a
[`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md)
method so that the data can be found.

## What a new class gets

- Means, pairs and contrasts:

  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
  [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md),
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md),
  [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  need [`coef()`](https://rdrr.io/r/stats/coef.html) (or `pm_coef()`),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) (or `pm_vcov()`) and
  [`terms()`](https://rdrr.io/r/stats/terms.html), plus the data (found
  from the call, or a
  [`pm_data()`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md)
  method). The columns of `pm_basis()` are matched to the coefficients
  by name.

- Residual plots:

  [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)
  and
  [`pm_check()`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md)
  also need [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) methods
  (residuals divided by their SD); without them they stop and say so.

- The rest:

  [`anova()`](https://rdrr.io/r/stats/anova.html),
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md),
  [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md),
  [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
  [`se_ranef()`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md),
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  and
  [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
  have their own code for the classes they list, and stop with that list
  for any other class.

## See also

[`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
[`pm_data`](https://dongwenluo.github.io/predictmeans/reference/pm_data.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  # A made-up model class that wraps lm() and has only four methods. (Methods are
  # found where predictmeans looks for them: in a package's namespace or, as
  # here, in the global environment; registerS3method() also works.)
  fit_mymodel <- function(formula, data) {
    structure(list(fit = lm(formula, data), call = match.call()), class = "mymodel")
  }
  registerS3method("coef", "mymodel", function(object, ...) coef(object$fit))
  registerS3method("vcov", "mymodel", function(object, ...) vcov(object$fit))
  registerS3method("terms", "mymodel", function(x, ...) terms(x$fit))
  registerS3method("df.residual", "mymodel", function(object, ...) object$fit$df.residual)

  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  m <- fit_mymodel(yield ~ Block + nitro*Variety, data = Oats)
  predictmeans(m, "nitro")$means
#>   nitro  estimate std.error df conf.low conf.high comparison.low
#> 1     0  79.38889  3.758096 55  71.8575  86.92028       74.06339
#> 2   0.2  98.88889  3.758096 55  91.3575 106.42028       93.56339
#> 3   0.4 114.22222  3.758096 55 106.6908 121.75362      108.89672
#> 4   0.6 123.38889  3.758096 55 115.8575 130.92028      118.06339
#>   comparison.high  n estimable
#> 1        84.71439 18      TRUE
#> 2       104.21439 18      TRUE
#> 3       119.54772 18      TRUE
#> 4       128.71439 18      TRUE

  # Degrees of freedom of its own, e.g. normal-based:
  registerS3method("pm_df", "mymodel", function(model, L, ...) rep(Inf, nrow(L)))
  predictmeans(m, "nitro")$means
#>   nitro  estimate std.error  df  conf.low conf.high comparison.low
#> 1     0  79.38889  3.758096 Inf  72.02316  86.75462       74.18053
#> 2   0.2  98.88889  3.758096 Inf  91.52316 106.25462       93.68053
#> 3   0.4 114.22222  3.758096 Inf 106.85649 121.58796      109.01386
#> 4   0.6 123.38889  3.758096 Inf 116.02316 130.75462      118.18053
#>   comparison.high  n estimable
#> 1        84.59725 18      TRUE
#> 2       104.09725 18      TRUE
#> 3       119.43058 18      TRUE
#> 4       128.59725 18      TRUE
```
