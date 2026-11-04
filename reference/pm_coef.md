# Supporting a New Model Class

[`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md),
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
and
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
get everything they need from a model through these six generics. Their
default methods work for any model class with the usual methods (`coef`,
`vcov`, `terms`, `df.residual`, `family`), so many classes need no extra
code. For a class that lacks one of these, or needs something special,
write a method for the corresponding generic only; nothing inside
predictmeans needs to change.

## Usage

``` r
pm_coef(model, ...)
pm_vcov(model, ...)
pm_basis(model, newdata, ...)
pm_data(model, ...)
pm_df(model, L, ...)
pm_linkinv(model, ...)
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
  as in `pm_coef(model)`. Default:
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) of the
  model's terms (with their stored `predvars`, so
  [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`,
  [`scale()`](https://rdrr.io/r/base/scale.html) etc. use the fitted
  basis), contrasts and factor levels. The `gam` method uses the linear
  predictor matrix of package mgcv.

- `pm_data(model)`:

  The data used to fit the model, or `NULL` (the default), in which case
  they are found from the `data` argument of the call,
  [`getData()`](https://rdrr.io/pkg/nlme/man/getData.html), the model
  frame or the formula environment.

- `pm_df(model, L)`:

  Denominator degrees of freedom for each row of `L`, or `NULL` (the
  default) to use the model's single residual df: `df.residual(model)`,
  or `Inf` (normal-based) when there is none. Kenward-Roger df for
  `lmer` and the rules of nlme for `lme` are built in.

- `pm_linkinv(model)`:

  Inverse link function used to back-transform the means, or `NULL` for
  none. Default: the inverse link of `family(model)` unless it is the
  identity.

The model also needs a [`terms()`](https://rdrr.io/r/stats/terms.html)
method for its fixed effects (used to find the variables in the model),
and its call ([`getCall()`](https://rdrr.io/r/stats/update.html)) or a
`pm_data()` method so that the data can be found.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  # A made-up model class that wraps lm() and has only four methods.
  fit_mymodel <- function(formula, data) {
    structure(list(fit = lm(formula, data), call = match.call()), class = "mymodel")
  }
  coef.mymodel <- function(object, ...) coef(object$fit)
  vcov.mymodel <- function(object, ...) vcov(object$fit)
  terms.mymodel <- function(x, ...) terms(x$fit)
  df.residual.mymodel <- function(object, ...) object$fit$df.residual

  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  m <- fit_mymodel(yield ~ Block + nitro*Variety, data = Oats)
  predictmeans(m, "nitro", plot = FALSE)$means
#> Error in x$terms %||% attr(x, "terms") %||% stop("no terms component nor attribute"): no terms component nor attribute

  # Degrees of freedom of its own, e.g. normal-based:
  pm_df.mymodel <- function(model, L, ...) rep(Inf, nrow(L))
  predictmeans(m, "nitro", plot = FALSE)$means
#> Error in x$terms %||% attr(x, "terms") %||% stop("no terms component nor attribute"): no terms component nor attribute
```
