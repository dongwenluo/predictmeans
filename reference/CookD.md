# Cook's Distances of a Linear (Mixed) Model

Cook's distances of the observations (or groups of observations) of a
model fitted by `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer`,
`glmmTMB` or [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), with
a plot.

## Usage

``` r
CookD(model, group=NULL, plot=TRUE, id=3, ...)
```

## Arguments

- model:

  Model object returned by `aov` (without `Error()`), `lm`, `glm`,
  `gls`, `lme`, `lmer`, `glmer`, `glmmTMB` or `gam`.

- group:

  Optional name (in "quotes") of a variable whose levels are deleted one
  at a time (e.g. the clusters of a mixed model); the default `NULL`
  deletes single observations.

- plot:

  `TRUE` (the default): draw the distances (ggplot2, as the Cook's
  distance panel of
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)).

- id:

  The number of the largest Cook's distances labelled in the plot (`idn`
  before 2.0, now as in
  [`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)).
  The default value is 3.

- ...:

  Not used: `idn`, `newwd` (before 2.0) stop with their 2.0 spelling.

## Value

Invisibly, the Cook's distances (named by observation or group).

## Details

Cook's distance of an observation (or of a group of observations,
`group`) is \\D = (\hat\beta - \hat\beta\_{(-i)})' V^{-1} (\hat\beta -
\hat\beta\_{(-i)}) / p\\, where \\\hat\beta\_{(-i)}\\ are the fixed
effects of the model refitted without it, \\V\\ the covariance of
\\\hat\beta\\ and \\p\\ their number (for a `gam`, the effective degrees
of freedom, as `mgcv`'s
[`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html);
the smoothing parameters are re-estimated for each deletion). For `lm`
and `aov` models and single observations this is the exact closed form
of
[`stats::cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html).
Other models are refitted for each deletion: exact also for `glm`
([`stats::cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
for a `glm` is a one-step approximation), and for mixed models and `gls`
fitted by maximum likelihood (the REML criterion changes with the fixed
effects). A deletion that leaves a coefficient without data gives `NA`.
Refitting a generalized linear mixed model for every observation takes
time; `group` deletes whole clusters, usually the more useful check for
mixed models.

## See also

[`residplot`](https://dongwenluo.github.io/predictmeans/reference/residplot.md),
[`pm_check`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
CookD(fm)

CookD(fm, group = "Block")   # delete whole blocks
```
