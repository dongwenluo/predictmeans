# Calculates and plots Cook's distances for a Linear (Mixed) Model

This function produces Cook's distance plots for a linear model obtained
from functions `aov`, `lm`, `glm`, `gls`, `lme`, or `lmer`.

## Usage

``` r
CookD(model, group=NULL, plot=TRUE, idn=3, newwd=FALSE)
```

## Arguments

- model:

  Model object returned by `aov`, `lm`, `glm`, `gls`, `lme`, and `lmer`.

- group:

  Name (in "quotes") for indicating how observations are deleted for
  Cook's distance calculation. If `group!=NULL` then deletions will be
  along levels of `group` variable, otherwise, will be along individual
  observations.

- plot:

  A logical variable; if it is true, a plot of Cook's distance will be
  presented. The default is TRUE.

- idn:

  An integer indicating the number of top Cook's distances to be
  labelled in the plot. The default value is 3.

- newwd:

  A logical variable to indicate whether to print graph in a new window.
  The default value is FALSE.

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
