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
                "eff", "del.eff"), ref = 1, by = NULL, max_degree = NULL, ...)
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

- max_degree:

  For `"poly"`: the highest degree, the default is the number of
  levels - 1 (at most 6).

- ...:

  Not used: `max.degree` (before 2.0) stops with its 2.0 spelling.

## Value

A matrix with one row per contrast (named as in emmeans, e.g.
`"0.2 - 0"`, `"linear"`; with `by` the group comes first, e.g.
`"Victory: 0.2 - 0"`) and one column per level of the term. Levels whose
means can't be estimated are left out of the families.

## See also

[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
[`specify_contrast`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  # Dunnett: each nitrogen level against no nitrogen, exact adjustment
  K <- Kmatrix(fm, "nitro")
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
    contrastmeans(fm, "nitro", ctrmatrix = contrast_family(K, "trt.vs.ctrl"),
                  adj = "mvt")
  }
#> Contrasts of 'nitro' (95% limits; p.adj: mvt; Kenward-Roger df)
#>  contrast estimate std.error df statistic   p.value     p.adj conf.low
#>   0.2 - 0    19.50     4.436 45     4.396 6.657e-05 2.200e-04    10.57
#>   0.4 - 0    34.83     4.436 45     7.853 5.649e-10 3.368e-09    25.90
#>   0.6 - 0    44.00     4.436 45     9.919 6.697e-13 8.626e-14    35.07
#>  conf.high adj.conf.low adj.conf.high
#>      28.43        8.718         30.28
#>      43.77       24.052         45.62
#>      52.93       33.218         54.78
  # linear and quadratic trends in nitrogen within each variety
  K2 <- Kmatrix(fm, "nitro:Variety")
  cf <- contrast_family(K2, "poly", by = "Variety", max_degree = 2)
  contrastmeans(fm, "nitro:Variety", ctrmatrix = cf)
#> Contrasts of 'nitro:Variety' (95% limits; p.adj: none; Kenward-Roger df)
#>                contrast estimate std.error df statistic   p.value     p.adj
#>     Golden Rain: linear  150.667     24.30 45    6.2014 1.568e-07 1.568e-07
#>  Golden Rain: quadratic   -8.333     10.87 45   -0.7670 4.471e-01 4.471e-01
#>      Marvellous: linear  129.167     24.30 45    5.3165 3.186e-06 3.186e-06
#>   Marvellous: quadratic  -12.167     10.87 45   -1.1198 2.688e-01 2.688e-01
#>         Victory: linear  162.167     24.30 45    6.6747 3.099e-08 3.099e-08
#>      Victory: quadratic  -10.500     10.87 45   -0.9664 3.390e-01 3.390e-01
#>  conf.low conf.high
#>    101.73   199.601
#>    -30.22    13.551
#>     80.23   178.101
#>    -34.05     9.717
#>    113.23   211.101
#>    -32.38    11.384
#> 
#> Advice:
#> i 6 contrasts tested without adjustment: some will be significant by
#>   chance. If they are a family, use adj = "holm" (any design) or
#>   "tukey"/"mvt" (all pairs).
```
