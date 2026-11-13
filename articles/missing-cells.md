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
  cells, with a warning naming it.
- **Its comparisons** with other means are made on the cells both means
  have (column `over` of `pairs`), not as the difference of two averages
  over different cells.
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
#> Warning: No data for N:Variety 0.6:Victory. The mean of N = 0.6 averages only
#> its observed cells; comparisons with it use only the cells that the compared
#> means all have (column 'over' of the pairs).
pm_n$means[, c("N", "estimate", "std.error", "n")]
#>     N  estimate std.error  n
#> 1   0  79.38889  6.972348 18
#> 2 0.2  98.88889  6.972348 18
#> 3 0.4 114.22222  6.972348 18
#> 4 0.6 125.83333  7.471667 12
```

That mean leaves out Victory, the lowest-yielding variety at the other
nitrogen levels, so it is too high compared with the others; emmeans
calls it non-estimable. The difference 125.8 − 79.4 = 46.4 would mix the
nitrogen effect with the variety effect (Victory is in the mean of N = 0
but not in that of 0.6). So predictmeans compares N = 0.6 with each
other level on the varieties both have, Golden Rain and Marvellous
(column `over`): 0.6 − 0 is the average of their two differences, 42.5.
Pairs among the complete levels are the differences of the means as
usual:

``` r

pm_n2 <- suppressWarnings(predictmeans(fm, "N", pairwise = TRUE, plot = FALSE, prtnum = FALSE))
pm_n2$pairs[, c("level1", "level2", "over", "estimate", "std.error", "p.value")]
#>   level1 level2                            over   estimate std.error
#> 1      0    0.2                            <NA> -19.500000  4.469636
#> 2      0    0.4                            <NA> -34.833333  4.469636
#> 3    0.2    0.4                            <NA> -15.333333  4.469636
#> 4      0    0.6 Variety Golden Rain, Marvellous -42.500000  5.474164
#> 5    0.2    0.6 Variety Golden Rain, Marvellous -22.333333  5.474164
#> 6    0.4    0.6 Variety Golden Rain, Marvellous  -9.916667  5.474164
#>        p.value
#> 1 8.760703e-05
#> 2 1.494073e-09
#> 3 1.409747e-03
#> 4 1.639811e-09
#> 5 2.085466e-04
#> 6 7.755869e-02
```

Comparing the two averages instead would wrongly call 0.6 and 0.4
different here (p = 0.03 against 0.08 on the common varieties). In
simulations (`dev/review-defaults.md` in the package sources) comparing
averages over different cells rejected 8–78% of true null hypotheses
when the other factor had an effect; comparing on the common cells 4–6%.
When two means have no cell in common their difference is NA
(`over = "none"`).
[`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
tests a term on the cells common to all its means, and
[`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
compares on the cells the means of each contrast have in common.

Other comparisons that stay within the observed cells are sound too,
e.g. nitrogen levels within each variety:

``` r

pm_v <- predictmeans(fm, "N:Variety", atvar = "Variety", adj = "tukey", plot = FALSE,
                     prtnum = FALSE)
within <- subset(pm_v$pairs, !is.na(by))
head(within[, c("level1", "level2", "by", "estimate", "p.adj")])
#>             level1          level2          by  estimate        p.adj
#> 4    0:Golden Rain 0.2:Golden Rain Golden Rain -18.50000 0.0954457441
#> 8     0:Marvellous  0.2:Marvellous  Marvellous -21.83333 0.0358722237
#> 13       0:Victory     0.2:Victory     Victory -18.16667 0.0607716581
#> 16   0:Golden Rain 0.4:Golden Rain Golden Rain -34.66667 0.0003454196
#> 19 0.2:Golden Rain 0.4:Golden Rain Golden Rain -16.16667 0.1743690902
#> 23    0:Marvellous  0.4:Marvellous  Marvellous -30.50000 0.0017427519
```

or contrasts written with
[`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
which refuses means that cannot be estimated:

``` r

cm <- specify_contrast(Kmatrix(fm, "N:Variety"),
  N_effect_GR_and_M = (`0.6:Golden Rain` + `0.6:Marvellous`)/2 - (`0:Golden Rain` + `0:Marvellous`)/2)
contrastmeans(fm, "N:Variety", ctrmatrix = cm)[[1]]
#>                   Estimate Std. Error t value      df Pr(>|t|) lower.CL
#> N_effect_GR_and_M     42.5     5.4742  7.7637 40.0636        0  31.4368
#>                   upper.CL
#> N_effect_GR_and_M  53.5632
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
