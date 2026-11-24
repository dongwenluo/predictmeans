# Data of a Fitted Model, with the Response on Both Scales

`pm_data()` returns the rows used to fit a model (after `subset` and
missing values), with the response on the scale of the model and on the
original scale: `.y` is the response as modelled (e.g. `log(y + 1)`; for
a GLM the observed response, as a proportion for
`cbind(successes, failures)`), and `.y_orig` the response on the
original scale (`y`), taken from the data or, when only the transformed
response is available, computed with the inverse transformation. A
transformation given to
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
as `trans` (for a response transformed before the fit) is used the same
way. It is also the generic through which a new model class can say
where its data are (see
[`pm_coef`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md)).

## Usage

``` r
pm_data(model, ...)

# Default S3 method
pm_data(model, data = NULL, trans = NULL, ...)
```

## Arguments

- model:

  A fitted model, or a
  [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fit.

- ...:

  Not used.

- data:

  Optional: the data the model was fitted to, when they cannot be found
  from the model.

- trans:

  Optional: a function that back-transforms the response (as `trans` of
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)),
  used when the model formula has an untransformed response.

## Value

A data frame with the variables of the data, `.y`, `.y_orig` and
`.weights` (prior weights, when there are any), and attribute
`"response"`: a list with `expr` (the response as written, e.g.
`"log(y + 1)"`), `name` (the response variable, e.g. `"y"`), `scale` of
the model's means (`"identity"`, `"transformed"` or `"link"`), `link`,
`family` and `back` (what back-transforms the means: `"link"`,
`"formula"`, `"trans"` or `NA`).

## Examples

``` r
  oats <- as.data.frame(nlme::Oats)
  fm <- lm(log(yield + 1) ~ Variety, data = oats)
  d <- pm_data(fm)
  head(d[c("yield", ".y", ".y_orig")])
#>   yield       .y .y_orig
#> 1   111 4.718499     111
#> 2   130 4.875197     130
#> 3   157 5.062595     157
#> 4   174 5.164786     174
#> 5   117 4.770685     117
#> 6   114 4.744932     114
  attr(d, "response")
#> $expr
#> [1] "log(yield + 1)"
#> 
#> $name
#> [1] "yield"
#> 
#> $scale
#> [1] "transformed"
#> 
#> $link
#> NULL
#> 
#> $family
#> [1] "gaussian"
#> 
#> $back
#> [1] "formula"
#> 
```
