# Plot a Residual Variogram

`autoplot()` of
[`residvariogram`](https://dongwenluo.github.io/predictmeans/reference/residvariogram.md):
the semivariance against the lag with the global envelope shaded and the
lags outside it marked in red. Grid layouts: `type = "edges"` (default)
the lags along a column and along a row; `"surface"` every row x column
lag as tiles.

## Usage

``` r
# S3 method for class 'pm_variogram'
autoplot(object, type = c("edges", "surface"), ...)
```

## Arguments

- object:

  A result of
  [`residvariogram()`](https://dongwenluo.github.io/predictmeans/reference/residvariogram.md).

- type:

  `"edges"` or `"surface"` (grid layouts only).

- ...:

  Not used.

## Value

A `ggplot` object.

## See also

[`residvariogram`](https://dongwenluo.github.io/predictmeans/reference/residvariogram.md)
