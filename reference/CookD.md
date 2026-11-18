# Calculates and plots Cook's distances for a Linear (Mixed) Model

This function produces Cook's distance plots for a model obtained from
functions `aov`, `lm`, `glm`, `gls`, `lme`, `lmer`, `glmer`, `glmmTMB`
or [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

## Usage

``` r
CookD(model, group=NULL, plot=TRUE, id=3, ...)
```

## Arguments

- model:

  Model object returned by `aov` (without `Error()`), `lm`, `glm`,
  `gls`, `lme`, `lmer`, `glmer`, `glmmTMB` or `gam`.

- group:

  Name (in "quotes") for indicating how observations are deleted for
  Cook's distance calculation. If `group!=NULL` then deletions will be
  along levels of `group` variable, otherwise, will be along individual
  observations.

- plot:

  A logical variable; if it is true, a plot of Cook's distance will be
  presented. The default is TRUE.

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

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
CookD(fm)
```
