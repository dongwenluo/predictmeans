# Permutation Index

Permutations allowed by the design: free, within blocks, of whole groups
(e.g. whole plots), or of whole groups within blocks. Give the result to
[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
or
[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
(argument `perms`) for restricted (design-based) permutation tests; when
its columns are every permutation the design allows, use `exact = TRUE`
in
[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md).

## Usage

``` r
permindex(data, block=NULL, group=NULL, nsim=4999, seed)
```

## Arguments

- data:

  The data the model was fitted to.

- block:

  Optional name (in "quotes") of the blocking factor in `data`.

- group:

  Optional name (in "quotes") of the factor whose levels are permuted as
  whole units.

- nsim:

  The number of permutations. The default is 4999.

- seed:

  Specify a random number generator seed, for reproducible results. The
  user's random number stream is left as it was.

## Value

A matrix with one row per row of `data` and up to `nsim` columns, each a
permutation of the row numbers of `data` (row i of the permuted data is
row `index[i]` of `data`). Duplicated permutations are dropped, so there
are fewer columns when the design allows fewer distinct permutations.

## Details

Groups are permuted as whole units, so the groups that change places
(within a block, if `block` is given) must have the same number of rows.
The data need not be sorted: before predictmeans 2.0 the indices
referred to the data sorted by `block` and `group`, so they were right
only for data already in that order (the permutations are the same as
before for sorted data and the same `seed`).

## See also

[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md),
[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  block <- rep(1:3, each=12)
  group <- rep(rep(1:3, each=4), 3)
  data <- data.frame(block, group)
  cbind(data, permindex(data, block="block", group="group", nsim=5))  
#>    block group  1  2  3  4  5
#> 1      1     1 11  3 10  1  2
#> 2      1     1 12  2 12  3  4
#> 3      1     1  9  4  9  4  1
#> 4      1     1 10  1 11  2  3
#> 5      1     2  6  9  4  5 10
#> 6      1     2  8 11  3  8 11
#> 7      1     2  5 12  1  7 12
#> 8      1     2  7 10  2  6  9
#> 9      1     3  3  7  5 12  8
#> 10     1     3  2  6  8  9  7
#> 11     1     3  4  5  7 11  6
#> 12     1     3  1  8  6 10  5
#> 13     2     1 16 19 19 22 20
#> 14     2     1 15 17 17 23 19
#> 15     2     1 13 18 18 24 18
#> 16     2     1 14 20 20 21 17
#> 17     2     2 23 21 16 17 16
#> 18     2     2 22 22 14 19 13
#> 19     2     2 21 24 13 20 15
#> 20     2     2 24 23 15 18 14
#> 21     2     3 18 13 23 13 22
#> 22     2     3 17 15 21 14 21
#> 23     2     3 19 14 22 16 24
#> 24     2     3 20 16 24 15 23
#> 25     3     1 31 26 30 29 33
#> 26     3     1 29 28 29 30 34
#> 27     3     1 32 27 32 31 35
#> 28     3     1 30 25 31 32 36
#> 29     3     2 36 35 28 33 27
#> 30     3     2 35 34 26 34 25
#> 31     3     2 34 33 25 35 26
#> 32     3     2 33 36 27 36 28
#> 33     3     3 26 32 34 26 30
#> 34     3     3 25 30 33 28 29
#> 35     3     3 28 31 35 27 31
#> 36     3     3 27 29 36 25 32
                        # Permute group as a whole within each block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, block="block",  nsim=5)) 
#>    block group  1  2  3  4  5
#> 1      1     1  1  2  6  5  3
#> 2      1     1  9  5  5  9  5
#> 3      1     1  8  1  9  7  6
#> 4      1     1  5 12  1  2 12
#> 5      1     2  2  7  7 11  4
#> 6      1     2 11 11  2 10 11
#> 7      1     2 10  6 11  8  7
#> 8      1     2  4 10 12  3  1
#> 9      1     3  6  4  4  1  8
#> 10     1     3 12  3 10 12  2
#> 11     1     3  7  8  8  4  9
#> 12     1     3  3  9  3  6 10
#> 13     2     1 24 21 16 18 21
#> 14     2     1 23 15 21 21 13
#> 15     2     1 14 23 14 19 23
#> 16     2     1 20 16 17 23 17
#> 17     2     2 21 24 23 17 18
#> 18     2     2 16 22 13 20 24
#> 19     2     2 19 20 19 13 16
#> 20     2     2 22 17 22 14 20
#> 21     2     3 13 18 15 16 22
#> 22     2     3 18 19 20 24 19
#> 23     2     3 15 13 24 22 15
#> 24     2     3 17 14 18 15 14
#> 25     3     1 30 36 35 34 31
#> 26     3     1 36 27 31 33 35
#> 27     3     1 32 28 27 26 27
#> 28     3     1 26 32 29 31 33
#> 29     3     2 33 29 30 36 28
#> 30     3     2 31 35 26 27 36
#> 31     3     2 25 30 36 35 34
#> 32     3     2 28 34 33 25 25
#> 33     3     3 27 25 28 29 30
#> 34     3     3 34 31 32 28 29
#> 35     3     3 29 26 34 32 32
#> 36     3     3 35 33 25 30 26
                        # Permute obs within each block only.
  cbind(data, permindex(data, group="group", nsim=5)) 
#>    block group  1  2  3  4  5
#> 1      1     1 34 15 35 27 33
#> 2      1     1 35 28 11  4 12
#> 3      1     1 23  1 23 26 35
#> 4      1     1 21  2 12  3 11
#> 5      1     2  6 10 32 20  8
#> 6      1     2 31  9 20  8 17
#> 7      1     2 20 36  5  7 20
#> 8      1     2  7 11  6 31 29
#> 9      1     3 13  8  1 23 13
#> 10     1     3  4 29 27  9 26
#> 11     1     3  1 20 26 11  4
#> 12     1     3 15 32 15 36 27
#> 13     2     1 22 27 34 14 24
#> 14     2     1 24  4  9 13  9
#> 15     2     1 11 25 10  1 36
#> 16     2     1 10 13 22 15 22
#> 17     2     2  8 21  8  6 18
#> 18     2     2 30 34 19 32  7
#> 19     2     2 29 23 18 19 30
#> 20     2     2  5 24 17 30  5
#> 21     2     3 25  7 28 33 15
#> 22     2     3  2 18 25 22  1
#> 23     2     3  3 17  2 10 14
#> 24     2     3 16  5  4 24  2
#> 25     3     1 12 26 24 28 34
#> 26     3     1  9  3 21  2 23
#> 27     3     1 33 14 33 16 10
#> 28     3     1 36 16 36 25 21
#> 29     3     2 17 12  7  5  6
#> 30     3     2 18 22 30 29 31
#> 31     3     2 19 33 31 17 32
#> 32     3     2 32 35 29 18 19
#> 33     3     3 27 30 13 21 28
#> 34     3     3 28 31  3 34  3
#> 35     3     3 14  6 16 12 25
#> 36     3     3 26 19 14 35 16
                        # Permute groups as a whole block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, nsim=5))  # Free permutation.
#>    block group  1  2  3  4  5
#> 1      1     1 28 28  8  7 13
#> 2      1     1 17 14 33 29 27
#> 3      1     1  7 15 14 26  6
#> 4      1     1 13 30  3 33 11
#> 5      1     2 15 26 29 17 32
#> 6      1     2 29  9 32 14 33
#> 7      1     2  4 11 27 28 25
#> 8      1     2 35 22 11 22  1
#> 9      1     3 22 10 22 21  3
#> 10     1     3 21 31  4  8 31
#> 11     1     3 10 20 12 24 23
#> 12     1     3  6  8 19 10 28
#> 13     2     1 23 13  9 18 20
#> 14     2     1 20 19 28 31  4
#> 15     2     1 12 29  6  1 36
#> 16     2     1  8 24 21 32 21
#> 17     2     2 27  1 34 19  8
#> 18     2     2 32 33 26 35 35
#> 19     2     2 33 32 17 15 16
#> 20     2     2 34 36 20  9 26
#> 21     2     3  9 21 31  4 15
#> 22     2     3 31 18  2 36 12
#> 23     2     3 36 23 23 30 34
#> 24     2     3  5  7 30 34 18
#> 25     3     1  1  2 13 27 17
#> 26     3     1 11 27 25 12 29
#> 27     3     1 18  6 10 13 24
#> 28     3     1 16 34 15  2  2
#> 29     3     2 24 17 18 16 14
#> 30     3     2  3 16 36 11 10
#> 31     3     2  2 35  5 23 30
#> 32     3     2 25  3 24 25 19
#> 33     3     3 26  4  1  3  5
#> 34     3     3 14 25 16  5  7
#> 35     3     3 30  5  7 20  9
#> 36     3     3 19 12 35  6 22

  ## A randomisation test of nitrogen in a randomised block design:
  ## permutations within blocks only
  Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro, data = Oats)
  idx <- permindex(Oats, block = "Block", nsim = 199, seed = 1)
  permmodels(fm, perms = idx, ncore = 1)
#> 
#> Coefficients of (fixed) effects:
#>             Estimate Std. Error t value Pr(>|t|) Perm_p_value
#> (Intercept)  79.3889     3.7668 21.0759   <1e-04           NA
#> Block.L      29.4823     4.6134  6.3906   <1e-04        0.070
#> Block.Q      20.0124     4.6134  4.3379   0.0001        1.000
#> Block.C       5.3728     4.6134  1.1646   0.2486        0.005
#> Block^4       4.7718     4.6134  1.0343   0.3049        1.000
#> Block^5      -1.2651     4.6134 -0.2742   0.7848        0.005
#> nitro0.2     19.5000     5.3271  3.6605   0.0005        0.005
#> nitro0.4     34.8333     5.3271  6.5389   <1e-04        0.005
#> nitro0.6     44.0000     5.3271  8.2597   <1e-04        0.005
#> 
#> Note: Perm_p_value of t test is obtained using ‘199’ Freedman-Lane permutations. The smallest attainable p-value is 1/(199 + 1) = 0.005.
#> 
#> ANOVA:
#> Analysis of Variance Table
#> Response: yield
#>           Df   Sum Sq   Mean Sq F value Pr(>F) Perm_p_value
#> Block      5 15875.28 3175.0556 12.4317 <1e-04            1
#> nitro      3 20020.50 6673.5000 26.1297 <1e-04        0.005
#> Residuals 63 16090.17  255.3995                            
#> 
#> Note: Perm_p_value of F test is obtained using ‘199’ Freedman-Lane permutations. The smallest attainable p-value is 1/(199 + 1) = 0.005.
#> 
#> Advice:
#> ! Only 199 resamples: p-values at the floor (0.005) only say 'as small as
#>   can be shown' (Monte Carlo SE about 0.015 at p = 0.05). Use nperm = 999
#>   or more (4999 for reporting); 'seed' makes it reproducible.
#> 
```
