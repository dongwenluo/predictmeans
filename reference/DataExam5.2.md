# Heights of Trees: 37 Seedlots at 6 Sites

Mean heights of 37 seedlots (26 Acacia, Allocasuarina and Casuarina
species) grown at six sites in Thailand, the data of Table 5.2 of
Williams, Harwood and Matheson (2024), used there for a regression of
each seedlot on the site mean (Finlay-Wilkinson regression). 30
site-seedlot combinations have no height.

## Usage

``` r
data(DataExam5.2)
```

## Format

A data frame with 222 rows and 4 columns:

- site:

  factor, 6 sites.

- seedlot:

  factor, 37 seedlot numbers.

- ht:

  height (cm); NA when the seedlot was not measured at the site.

- sitemean:

  the mean of `ht` at the site (one value per site).

## Source

Williams, E. R., Harwood, C. E. and Matheson, A. C. (2024).
*Experimental Design and Analysis for Tree Improvement* (3rd ed.). CSIRO
Publishing, Table 5.2.
[doi:10.1071/9781486317110](https://doi.org/10.1071/9781486317110)

## Details

`sitemean` takes one value per site, so in
`ht ~ site + sitemean * seedlot` it is aliased with `site` (its
coefficient is NA); see the section "Covariates aliased with a factor"
of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).
On the four sites and 27 seedlots with complete data (and `sitemean`
recomputed for them)
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
reproduces the means and slopes of the book's Table 5.7
(`tests/testthat/test-level-covariate.R`).

## Examples

``` r
data(DataExam5.2)
fit <- lm(ht ~ site + sitemean * seedlot, data = DataExam5.2)
head(predictmeans(fit, "seedlot", plot = FALSE)$mean_table)
#> Covariate sitemean is constant within the levels of site: each mean uses the value of each level averaged over.
#>   seedlot     Mean       SE  Df  LL(95%)  UL(95%)
#> 1   13877 291.0000 30.20294 114 231.1682 350.8318
#> 2   13866 302.0000 30.20294 114 242.1682 361.8318
#> 3   13689 462.8333 30.20294 114 403.0015 522.6651
#> 4   13688 458.1667 30.20294 114 398.3349 517.9985
#> 5   13861 538.5000 30.20294 114 478.6682 598.3318
#> 6   13854 536.5000 30.20294 114 476.6682 596.3318
# each seedlot averaged over the sites where it was measured
head(predictmeans(fit, "seedlot", plot = FALSE, weights = "cells")$mean_table)
#> Covariate sitemean is constant within the levels of site: each mean uses the value of each level averaged over.
#>   seedlot     Mean       SE  Df  LL(95%)  UL(95%)
#> 1   13877 291.0000 30.20294 114 231.1682 350.8318
#> 2   13866 302.0000 30.20294 114 242.1682 361.8318
#> 3   13689 462.8333 30.20294 114 403.0015 522.6651
#> 4   13688 458.1667 30.20294 114 398.3349 517.9985
#> 5   13861 538.5000 30.20294 114 478.6682 598.3318
#> 6   13854 536.5000 30.20294 114 476.6682 596.3318
```
