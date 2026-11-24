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
#> 1      1     1  1  6  3  3  1
#> 2      1     1  2  7  2  2  4
#> 3      1     1  3  8  4  4  2
#> 4      1     1  4  5  1  1  3
#> 5      1     2  7 11  8  8 11
#> 6      1     2  8  9  7  5 10
#> 7      1     2  5 10  5  7  9
#> 8      1     2  6 12  6  6 12
#> 9      1     3  9  2  9 10  7
#> 10     1     3 10  3 10 12  5
#> 11     1     3 11  1 12  9  8
#> 12     1     3 12  4 11 11  6
#> 13     2     1 17 15 16 17 22
#> 14     2     1 19 13 15 18 21
#> 15     2     1 18 14 13 20 24
#> 16     2     1 20 16 14 19 23
#> 17     2     2 23 24 20 16 15
#> 18     2     2 24 23 19 13 13
#> 19     2     2 22 22 17 14 14
#> 20     2     2 21 21 18 15 16
#> 21     2     3 13 20 21 22 18
#> 22     2     3 16 19 24 24 19
#> 23     2     3 14 18 23 21 20
#> 24     2     3 15 17 22 23 17
#> 25     3     1 34 35 27 30 29
#> 26     3     1 35 33 26 31 30
#> 27     3     1 36 34 25 29 31
#> 28     3     1 33 36 28 32 32
#> 29     3     2 27 28 36 27 33
#> 30     3     2 25 27 33 25 35
#> 31     3     2 28 26 34 26 34
#> 32     3     2 26 25 35 28 36
#> 33     3     3 32 30 32 36 28
#> 34     3     3 31 31 31 34 25
#> 35     3     3 30 32 29 35 27
#> 36     3     3 29 29 30 33 26
                        # Permute group as a whole within each block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, block="block",  nsim=5)) 
#>    block group  1  2  3  4  5
#> 1      1     1  7  8  5 12  6
#> 2      1     1  6  4  2  8 10
#> 3      1     1  4  2  3  6  3
#> 4      1     1  8  9  7  5 12
#> 5      1     2  2  5 10  3  4
#> 6      1     2  3 11  4 11  5
#> 7      1     2  5  1  6  7 11
#> 8      1     2  9  3 12 10  2
#> 9      1     3  1  6  8  1  8
#> 10     1     3 10 12  1  9  7
#> 11     1     3 11  7  9  2  9
#> 12     1     3 12 10 11  4  1
#> 13     2     1 20 24 18 23 13
#> 14     2     1 22 14 22 16 23
#> 15     2     1 18 13 24 18 22
#> 16     2     1 17 17 19 22 16
#> 17     2     2 21 15 21 15 20
#> 18     2     2 16 23 13 17 15
#> 19     2     2 15 22 16 14 18
#> 20     2     2 13 18 15 19 21
#> 21     2     3 14 21 20 13 24
#> 22     2     3 24 20 23 21 19
#> 23     2     3 19 19 17 20 14
#> 24     2     3 23 16 14 24 17
#> 25     3     1 28 33 31 26 26
#> 26     3     1 33 34 28 36 29
#> 27     3     1 29 26 36 28 35
#> 28     3     1 32 25 29 27 25
#> 29     3     2 30 30 25 32 31
#> 30     3     2 36 35 27 34 36
#> 31     3     2 31 29 33 29 33
#> 32     3     2 25 28 34 33 32
#> 33     3     3 27 32 35 31 30
#> 34     3     3 26 36 32 25 28
#> 35     3     3 34 31 30 35 27
#> 36     3     3 35 27 26 30 34
                        # Permute obs within each block only.
  cbind(data, permindex(data, group="group", nsim=5)) 
#>    block group  1  2  3  4  5
#> 1      1     1 30 18 16  6 35
#> 2      1     1 20  7 26  7 36
#> 3      1     1 18 20 13 29  9
#> 4      1     1 17 17 27  5 21
#> 5      1     2 28 11 35 25  8
#> 6      1     2 13 36 34 14  7
#> 7      1     2 26 33 22 15 31
#> 8      1     2 25  9 12  2 32
#> 9      1     3 12 13  6 35 13
#> 10     1     3 35 26  7  9  2
#> 11     1     3 10 27 32 33 15
#> 12     1     3 33  4 30 34 25
#> 13     2     1 32 29  3  8 33
#> 14     2     1 19 30 14 17 12
#> 15     2     1  5 32  1 30 22
#> 16     2     1  6  8  4 18 11
#> 17     2     2  3 10 10  1 29
#> 18     2     2 15 23 23 27 17
#> 19     2     2  4 24 33  3 18
#> 20     2     2  1 34 11 13  5
#> 21     2     3 36 25 18 23 28
#> 22     2     3 24 15 20 12 16
#> 23     2     3 34 28 31 10 14
#> 24     2     3 23  2 29 24 26
#> 25     3     1 29  6 15 19 10
#> 26     3     1  8  5  2 20 34
#> 27     3     1 31 31 28 31 23
#> 28     3     1  7 19 25 32 24
#> 29     3     2 27 35 24  4 20
#> 30     3     2 16 21  9 28 19
#> 31     3     2  2 12 21 26 30
#> 32     3     2 14 22 36 16  6
#> 33     3     3 22  3 17 21  1
#> 34     3     3 11 14  5 36  4
#> 35     3     3  9  1  8 11  3
#> 36     3     3 21 16 19 22 27
                        # Permute groups as a whole block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, nsim=5))  # Free permutation.
#>    block group  1  2  3  4  5
#> 1      1     1 17 23 24  5 14
#> 2      1     1  9 28  5 19 16
#> 3      1     1 10 22 21  1  4
#> 4      1     1 11  3 20 28 19
#> 5      1     2 27 13 22 35 21
#> 6      1     2 13 17 23 23 34
#> 7      1     2 32 34 27 30 23
#> 8      1     2 28 29 31 31  1
#> 9      1     3 29 32 15 27 36
#> 10     1     3  2  4 17 13 20
#> 11     1     3 31 16 33 24  9
#> 12     1     3 19 33 30 25 28
#> 13     2     1  8  8 28 36 13
#> 14     2     1  6 36 10 14 24
#> 15     2     1 21 10 13  9 22
#> 16     2     1 24  2 26 29 31
#> 17     2     2 34  9 12  4 25
#> 18     2     2  4  6 14 20 27
#> 19     2     2  1 21 29  3 29
#> 20     2     2 16 35 34 12 18
#> 21     2     3  7 30  2  6  3
#> 22     2     3 30  7  9 33  6
#> 23     2     3 25 27 35 11 33
#> 24     2     3  5 24  3 32  8
#> 25     3     1 23 31 25 26 30
#> 26     3     1  3  1  8 15  7
#> 27     3     1 20 25 16 16 26
#> 28     3     1 35 26  7 18 35
#> 29     3     2 36 14 19  8  2
#> 30     3     2 15  5 11 22  5
#> 31     3     2 12 11  6 10 11
#> 32     3     2 14 15 36  2 17
#> 33     3     3 26 19 18 17 12
#> 34     3     3 22 20  4 21 32
#> 35     3     3 18 18  1 34 15
#> 36     3     3 33 12 32  7 10

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
