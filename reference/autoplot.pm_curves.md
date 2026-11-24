# Plots of Curves, Surfaces and Their Differences

`autoplot()` draws a result of
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
or
[`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
as a ggplot2 object (a plotly widget for `type = "3d"`): the curves with
their bands and the data, surfaces as filled contours, slopes
(`semipred(deriv = )`) with the line at zero, or the difference of two
curves (`contr`). It reads the tables of the result only, so the layout
can be changed without computing again. `plot = TRUE` in
[`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
and
[`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
prints it, and `$plt` of the result holds it; their `plot_args`
(superseded) give the same arguments under their 1.x names. Those
settings are kept in the result (`attr(x, "plot_args")`) and are the
defaults of `autoplot(x)`; arguments given to `autoplot()` win
(`autoplot(x, interval = "ci")` draws the bands of a result made with
`ci = FALSE`).

The curves and their bands share one legend. When each panel is one
group (one factor, the default layout) the strips name the groups and
there is no legend. With more than 8 groups the colours are evenly
spaced hues, and with more than 12 the legend is set in columns with
small keys; `panel = character(0)` puts all curves in one panel.

## Usage

``` r
# S3 method for class 'pm_curves'
autoplot(object, type = NULL, group = NULL, panel = NULL,
         scales = c("fixed", "free", "free_x", "free_y"), scale = NULL,
         data = TRUE, jitter = 0, interval = c("ci", "none"), title = NULL,
         style = list(), ...)
```

## Arguments

- object:

  A result of
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  or
  [`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md).

- type:

  `"lines"` (one covariate: curves or slopes), `"surface"` (two
  covariates: filled contours, one panel per level), `"3d"` (two
  covariates: a plotly 3-D surface, one widget with a scene per level;
  needs plotly) or `"difference"` (a result with `contr`). The default
  is `"difference"` with `contr`, else `"lines"` or `"surface"` by the
  number of covariates.

- group, panel:

  Names of factors of `modelterm`: `group` colours the curves, `panel`
  makes the panels (one factor: `facet_wrap()`; two: `facet_grid()` rows
  by columns; more: the first gives the rows, the others together the
  columns). Factors in neither are joined to `group`, so each
  combination has its own curve. The defaults: with one factor, colour
  and one panel per level (as 1.x); with more, the first factor colours
  and the others make the panels. `panel = character(0)` draws one
  panel. These replace `plotord` and `trellis` of 1.x.

- scales:

  The axes of the panels: `"fixed"` (default), `"free"`, `"free_x"`,
  `"free_y"`.

- scale:

  `"model"` (the link or the transformed response, column `estimate`) or
  `"response"` (the back-transformed curves, columns `response...`).
  Default: `"response"` when the curves are back-transformed. For
  `contr`: `"model"` shows the difference (as a ratio or odds ratio for
  a log or logit), `"response"` the difference of the back-transformed
  curves; the default is `contr_scale` of the call. Slopes are on the
  model scale.

- data:

  Draw the data (`TRUE`, jittered by `jitter`): the original response on
  the response scale, the modelled one on the model scale (none for a
  link on the model scale). Not drawn for slopes and differences.

- jitter:

  Jitter of the data points, in x and y (units of the axes).

- interval:

  `"ci"`: the band the result has (simultaneous with
  `band = "simultaneous"`); `"none"`.

- title:

  The title; `NULL` gives the default
  ([`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md):
  that of 1.x), `FALSE` none.

- style:

  Colours and sizes, any of: `line` and `band` (one curve), `diff_line`
  and `diff_band` (a difference), `ref` (the reference line and the
  marks where a band leaves it out), `alpha`, `linewidth`, `points`,
  `point_alpha`, `point_size`, `palette` (curves by group; ggplot2's own
  beyond its length), `surface_points`.

- ...:

  Not used.

## Value

A `ggplot` object, or for `type = "3d"` a plotly widget. Add to it with
ggplot2 (`+ labs(y = "Rate")`).

## See also

[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
[`semipred`](https://dongwenluo.github.io/predictmeans/reference/semipred.md),
[`tidy.pm_curves`](https://dongwenluo.github.io/predictmeans/reference/tidy.pm_curves.md)

## Examples

``` r
  library(predictmeans)
  data(Oats, package = "nlme")
  Oats$B <- factor(Oats$Block %in% c("I", "II", "III"), labels = c("late", "early"))
  fm <- lm(yield ~ nitro * Variety * B, data = Oats)
  cm <- covariatemeans(fm, "Variety:B", covariate = "nitro")
  ggplot2::autoplot(cm)                                   # colour Variety, panels B

  ggplot2::autoplot(cm, group = "B", panel = "Variety")

  ggplot2::autoplot(cm, panel = character(0), data = FALSE)
```
