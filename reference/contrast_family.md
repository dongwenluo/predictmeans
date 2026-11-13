# Families of Contrasts among Means

Build a standard family of contrasts among the means of a model term, as
a matrix for the `ctrmatrix` argument of
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md).
The families are those of emmeans, with the same coefficients and names:

- `"pairwise"`, `"revpairwise"`:

  all pairs, earlier minus later level (or later minus earlier);

- `"trt.vs.ctrl"`:

  each level minus the reference level `ref` (Dunnett's comparisons; use
  `adj = "mvt"` in `contrastmeans` for the exact Dunnett adjustment);

- `"consec"`:

  each level minus the one before;

- `"poly"`:

  orthogonal polynomials (linear, quadratic, ...) with integer
  coefficients, for equally spaced levels;

- `"eff"`:

  each level minus the average of all levels;

- `"del.eff"`:

  each level minus the average of the other levels.

With `by`, the family is built within each level (or combination of
levels) of the `by` factor(s), as `atvar` does for pairwise comparisons.
Note that `contrastmeans` adjusts the p-values of all rows as one
family.

## Usage

``` r
contrast_family(.x, family = c("pairwise", "revpairwise", "trt.vs.ctrl", "consec", "poly",
                "eff", "del.eff"), ref = 1, by = NULL, max.degree = NULL)
```

## Arguments

- .x:

  A
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  or
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  result for the model term, or a character vector of the level names
  (then `by` can't be used).

- family:

  The family of contrasts (see Description).

- ref:

  For `"trt.vs.ctrl"`: the reference level, by position or name (within
  each `by` group), the default is the first level.

- by:

  Optional name(s) of factor(s) of the term: build the family within
  each of their levels.

- max.degree:

  For `"poly"`: the highest degree, the default is the number of
  levels - 1 (at most 6).

## Value

A matrix with one row per contrast (named as in emmeans, e.g.
`"0.2 - 0"`, `"linear"`; with `by` the group comes first, e.g.
`"Victory: 0.2 - 0"`) and one column per level of the term. Levels whose
means can't be estimated are left out of the families.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  # Dunnett: each nitrogen level against no nitrogen, exact adjustment
  K <- Kmatrix(fm, "nitro")
  contrastmeans(fm, "nitro", ctrmatrix = contrast_family(K, "trt.vs.ctrl"), adj = "mvt")
#> $`The t tests of the specified contrasts`
#>         Estimate Std. Error t value      df Pr(>|t|) lower.CL upper.CL
#> 0.2 - 0  19.5000     4.4358  4.3961 45.0001    2e-04  10.5659  28.4341
#> 0.4 - 0  34.8333     4.4358  7.8529 45.0001    0e+00  25.8993  43.7674
#> 0.6 - 0  44.0000     4.4358  9.9194 45.0001    0e+00  35.0659  52.9341
#>         adj.lower.CL adj.upper.CL
#> 0.2 - 0       8.7182      30.2818
#> 0.4 - 0      24.0516      45.6151
#> 0.6 - 0      33.2182      54.7818
#> attr(,"Note")
#> [1] "The p-value is adjusted by ‘mvt’ method, if p-value = 0 means p-value < 0.0001."
#> 
#> $K
#>          0 0.2 0.4 0.6
#> 0.2 - 0 -1   1   0   0
#> 0.4 - 0 -1   0   1   0
#> 0.6 - 0 -1   0   0   1
#> 
  # linear and quadratic trends in nitrogen within each variety
  K2 <- Kmatrix(fm, "nitro:Variety")
  cf <- contrast_family(K2, "poly", by = "Variety", max.degree = 2)
  contrastmeans(fm, "nitro:Variety", ctrmatrix = cf)
#> $`The t tests of the specified contrasts`
#>                        Estimate Std. Error t value      df Pr(>|t|) lower.CL
#> Golden Rain: linear    150.6667    24.2956  6.2014 45.0001   0.0000 101.7328
#> Golden Rain: quadratic  -8.3333    10.8653 -0.7670 45.0001   0.4471 -30.2172
#> Marvellous: linear     129.1667    24.2956  5.3165 45.0001   0.0000  80.2328
#> Marvellous: quadratic  -12.1667    10.8653 -1.1198 45.0001   0.2688 -34.0506
#> Victory: linear        162.1667    24.2956  6.6747 45.0001   0.0000 113.2328
#> Victory: quadratic     -10.5000    10.8653 -0.9664 45.0001   0.3390 -32.3839
#>                        upper.CL
#> Golden Rain: linear    199.6005
#> Golden Rain: quadratic  13.5506
#> Marvellous: linear     178.1005
#> Marvellous: quadratic    9.7172
#> Victory: linear        211.1005
#> Victory: quadratic      11.3839
#> attr(,"Note")
#> [1] "The p-value is adjusted by ‘none’ method, if p-value = 0 means p-value < 0.0001."
#> 
#> $K
#>                        0:Golden Rain 0:Marvellous 0:Victory 0.2:Golden Rain
#> Golden Rain: linear               -3            0         0              -1
#> Golden Rain: quadratic             1            0         0              -1
#> Marvellous: linear                 0           -3         0               0
#> Marvellous: quadratic              0            1         0               0
#> Victory: linear                    0            0        -3               0
#> Victory: quadratic                 0            0         1               0
#>                        0.2:Marvellous 0.2:Victory 0.4:Golden Rain
#> Golden Rain: linear                 0           0               1
#> Golden Rain: quadratic              0           0              -1
#> Marvellous: linear                 -1           0               0
#> Marvellous: quadratic              -1           0               0
#> Victory: linear                     0          -1               0
#> Victory: quadratic                  0          -1               0
#>                        0.4:Marvellous 0.4:Victory 0.6:Golden Rain
#> Golden Rain: linear                 0           0               3
#> Golden Rain: quadratic              0           0               1
#> Marvellous: linear                  1           0               0
#> Marvellous: quadratic              -1           0               0
#> Victory: linear                     0           1               0
#> Victory: quadratic                  0          -1               0
#>                        0.6:Marvellous 0.6:Victory
#> Golden Rain: linear                 0           0
#> Golden Rain: quadratic              0           0
#> Marvellous: linear                  3           0
#> Marvellous: quadratic               1           0
#> Victory: linear                     0           3
#> Victory: quadratic                  0           1
#> 
```
