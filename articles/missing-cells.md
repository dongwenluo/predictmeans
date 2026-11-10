# Missing cells and non-estimable means

A predicted mean is a linear combination of the model’s coefficients.
When some treatment combinations have no data, or two factors are
confounded, some of these combinations cannot be estimated from the
data: emmeans prints them as `nonEst` and SAS as `Non-est`. This
vignette shows what predictmeans does in these cases and how to get
means and comparisons that are estimable.

In short:

- **Cell means** of combinations with no data are not shown.
- **A mean that averages over an empty cell** (e.g. the mean of a
  nitrogen level over the varieties when one variety lacks that level)
  keeps the rule of predictmeans 1.x: it is the average of the observed
  cells, with a warning. Such a mean is not comparable with the means
  that average over all cells.
- **Means that the data cannot determine** (because factors are
  confounded) are `NA`, with a message. Their pairwise differences are
  still given when they can be estimated, and letters are given among
  those.
- The `means` table has an `estimable` column and the number of
  observations, `n`.

``` r

library(predictmeans)
```

## An empty cell

Suppose the plots of Victory at the highest nitrogen level were lost in
the oats split plot:

``` r

Oats <- as.data.frame(nlme::Oats)
Oats$N <- factor(Oats$nitro)
Om <- subset(Oats, !(Variety == "Victory" & N == "0.6"))
with(Om, table(Variety, N))
#>              N
#> Variety       0 0.2 0.4 0.6
#>   Golden Rain 6   6   6   6
#>   Marvellous  6   6   6   6
#>   Victory     6   6   6   0
fm <- lmer(yield ~ N * Variety + (1 | Block/Variety), data = Om)
#> fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
```

The model with the interaction has one coefficient fewer (lme4 drops
it). The cell means are those of the 11 combinations with data:

``` r

pm_cells <- predictmeans(fm, "N:Variety", plot = FALSE, prtnum = FALSE)
nrow(pm_cells$means)
#> [1] 11
```

The nitrogen means average over the varieties. For N = 0.6 only two
varieties have data, so its mean is their average, and predictmeans
warns:

``` r

pm_n <- predictmeans(fm, "N", plot = FALSE, prtnum = FALSE)
#> Warning in Kmatrix(model, modelterm, covariate, data = data, prtnum = prtnum, :
#> Missing treatments' combination appeared, predicted means maybe misleading!
pm_n$means[, c("N", "estimate", "std.error", "n")]
#>     N  estimate std.error  n
#> 1   0  79.38889  6.972348 18
#> 2 0.2  98.88889  6.972348 18
#> 3 0.4 114.22222  6.972348 18
#> 4 0.6 125.83333  7.471667 12
```

That mean leaves out Victory, the lowest-yielding variety at the other
nitrogen levels, so it is too high compared with the others; emmeans
calls it non-estimable. Comparisons that stay within the observed cells
are sound, e.g. nitrogen levels within each variety:

``` r

pm_v <- predictmeans(fm, "N:Variety", atvar = "Variety", adj = "tukey", plot = FALSE,
                     prtnum = FALSE)
within <- subset(pm_v$pairs, !is.na(by))
head(within[, c("level1", "level2", "by", "estimate", "p.adj")])
#>             level1          level2          by  estimate        p.adj
#> 4    0:Golden Rain 0.2:Golden Rain Golden Rain -18.50000 0.0954637136
#> 8     0:Marvellous  0.2:Marvellous  Marvellous -21.83333 0.0358850523
#> 13       0:Victory     0.2:Victory     Victory -18.16667 0.0607840836
#> 16   0:Golden Rain 0.4:Golden Rain Golden Rain -34.66667 0.0003460119
#> 19 0.2:Golden Rain 0.4:Golden Rain Golden Rain -16.16667 0.1743878815
#> 23    0:Marvellous  0.4:Marvellous  Marvellous -30.50000 0.0017447497
```

or contrasts written with
[`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
which refuses means that cannot be estimated:

``` r

cm <- specify_contrast(Kmatrix(fm, "N:Variety"),
  N_effect_GR_and_M = (`0.6:Golden Rain` + `0.6:Marvellous`)/2 - (`0:Golden Rain` + `0:Marvellous`)/2)
contrastmeans(fm, "N:Variety", ctrmatrix = cm)[[1]]
#>                   Estimate Std. Error t value      df Pr(>|t|)
#> N_effect_GR_and_M     42.5     5.4742  7.7637 40.0194        0
#> attr(,"Note")
#> [1] "The p-value is adjusted by 'none' method, if p-value = 0 means p-value < 0.0001."
```

## Confounded factors

In this spray trial the highest spray volume was used only with fine
droplets, and the other two volumes only with medium droplets:

``` r

set.seed(3)
spray <- data.frame(
  droplet = factor(rep(c("medium", "fine"), c(40, 10)), levels = c("medium", "fine")),
  vol = factor(c(rep(c("130", "200"), 20), rep("400", 10))),
  conc = factor(rep(c("none", "half"), 25)))
spray$y <- 5 + 0.3 * (spray$conc == "half") - 0.4 * (spray$droplet == "fine") +
  rnorm(50, sd = 0.2)
with(spray, table(droplet, vol))
#>         vol
#> droplet  130 200 400
#>   medium  20  20   0
#>   fine     0   0  10
```

In a model with `droplet` and `vol` as separate factors, the effect of
volume 400 cannot be told apart from that of fine droplets. A volume
mean averaged over both droplet sizes needs both, so none of the volume
means can be estimated:

``` r

fs <- lm(y ~ droplet + conc + vol, data = spray)
pv <- predictmeans(fs, "vol", pairwise = TRUE, plot = FALSE, prtnum = FALSE)
#> Means that cannot be estimated are shown as NA: 130, 200, 400. Their pairwise differences are shown when estimable. If factors are confounded or cells are empty, weights = "cells" averages over the observed cells only.
pv$means[, c("vol", "estimate", "n", "estimable", "group")]
#>   vol estimate  n estimable group
#> 1 130       NA 20     FALSE     A
#> 2 200       NA 20     FALSE     A
#> 3 400       NA 10     FALSE  <NA>
pv$pairs[, c("level1", "level2", "estimate", "std.error", "p.value")]
#>   level1 level2   estimate std.error    p.value
#> 1    130    200 -0.2166482 0.1249636 0.08967159
#> 2    130    400         NA        NA         NA
#> 3    200    400         NA        NA         NA
```

The difference between volumes 130 and 200 (both with medium droplets)
can still be estimated, and the letters compare those two only.

There are two ways to get means. `weights = "cells"` averages over the
observed combinations only, weighting them by their numbers of
observations; these means depend on how the data happen to be spread
over the other factors:

``` r

predictmeans(fs, "vol", weights = "cells", plot = FALSE, prtnum = FALSE)$means[, c("vol", "estimate", "std.error")]
#>   vol estimate  std.error
#> 1 130 4.987406 0.03951695
#> 2 200 5.289860 0.03951695
#> 3 400 4.731546 0.05588541
```

Better, state the structure in the model: volume is nested in droplet
size. The means of the observed droplet–volume combinations, averaged
equally over `conc`, are then estimable (and equal to what emmeans gives
after detecting the nesting):

``` r

fs2 <- lm(y ~ conc + droplet/vol, data = spray)
pn <- predictmeans(fs2, "droplet:vol", pairwise = TRUE, plot = FALSE, prtnum = FALSE)
#> Warning in Kmatrix(model, modelterm, covariate, data = data, prtnum = prtnum, :
#> Missing treatments' combination appeared, predicted means maybe misleading!
pn$means[, c("droplet", "vol", "estimate", "std.error", "n", "group")]
#>   droplet vol estimate  std.error  n group
#> 1  medium 130 5.030309 0.06844537 20    A 
#> 2  medium 200 5.246957 0.06844537 20    A 
#> 3    fine 400 4.731546 0.05588541 10     B
```

The warning about missing combinations refers to the three combinations
that cannot occur in this design.

## Checklist

- Look at the `n` and `estimable` columns of the `means` table, and at
  warnings and messages.
- Compare means that average over the same cells; use `atvar` or
  contrasts within the observed cells when some are empty.
- When factors are confounded, write the structure into the model
  (nesting, or one factor for the observed combinations) rather than
  averaging over combinations that do not exist.
