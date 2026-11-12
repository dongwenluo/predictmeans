# Permutation Index

This function obtains permutation index for a dataset: free
permutations, permutations within blocks, permutations of whole groups
(e.g. whole plots), or of whole groups within blocks. The result can be
given to
[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
or
[`permlmer`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md)
(argument `perms`) for tests that permute only as the design allows
(restricted, or design-based, permutation tests).

## Usage

``` r
permindex(data, block=NULL, group=NULL, nsim=4999, seed)
```

## Arguments

- data:

  Data object used in the `model` fitting.

- block:

  Name (in "quotes") for the blocking factor in the `data`.

- group:

  Name (in "quotes") for the group factor in the `data`.

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
#> 1      1     1  9  8 12  9  3
#> 2      1     1 11  5  9 12  1
#> 3      1     1 12  7 10 11  2
#> 4      1     1 10  6 11 10  4
#> 5      1     2  7  4  5  7 10
#> 6      1     2  6  1  8  5 11
#> 7      1     2  5  3  6  8  9
#> 8      1     2  8  2  7  6 12
#> 9      1     3  2  9  3  3  7
#> 10     1     3  1 11  2  2  6
#> 11     1     3  4 12  4  1  8
#> 12     1     3  3 10  1  4  5
#> 13     2     1 23 18 17 20 20
#> 14     2     1 22 17 18 17 18
#> 15     2     1 21 20 19 19 19
#> 16     2     1 24 19 20 18 17
#> 17     2     2 19 22 16 15 15
#> 18     2     2 18 23 15 13 13
#> 19     2     2 20 24 13 16 14
#> 20     2     2 17 21 14 14 16
#> 21     2     3 13 13 24 22 22
#> 22     2     3 16 14 23 23 24
#> 23     2     3 15 16 21 24 21
#> 24     2     3 14 15 22 21 23
#> 25     3     1 31 31 32 26 27
#> 26     3     1 32 30 31 25 26
#> 27     3     1 29 32 30 27 25
#> 28     3     1 30 29 29 28 28
#> 29     3     2 28 25 25 36 35
#> 30     3     2 26 28 28 34 34
#> 31     3     2 25 26 26 33 33
#> 32     3     2 27 27 27 35 36
#> 33     3     3 34 35 33 29 32
#> 34     3     3 35 33 36 31 30
#> 35     3     3 33 36 35 32 31
#> 36     3     3 36 34 34 30 29
                        # Permute group as a whole within each block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, block="block",  nsim=5)) 
#>    block group  1  2  3  4  5
#> 1      1     1  7 11  6 11  8
#> 2      1     1 12  4  7  3  7
#> 3      1     1  1  5  4 10  5
#> 4      1     1  4  1 12 12  6
#> 5      1     2 11  8  8  1  4
#> 6      1     2  8 12  1  2 12
#> 7      1     2  5  9  5  8  2
#> 8      1     2  9 10  2  5 11
#> 9      1     3  2  3 11  7  9
#> 10     1     3  3  6  3  6 10
#> 11     1     3  6  2  9  4  1
#> 12     1     3 10  7 10  9  3
#> 13     2     1 19 13 13 16 22
#> 14     2     1 23 24 21 18 13
#> 15     2     1 18 20 20 15 21
#> 16     2     1 21 17 17 19 19
#> 17     2     2 20 14 18 20 24
#> 18     2     2 16 22 23 14 23
#> 19     2     2 13 18 24 17 18
#> 20     2     2 15 16 19 23 16
#> 21     2     3 24 15 15 22 17
#> 22     2     3 22 19 22 13 14
#> 23     2     3 14 21 14 21 20
#> 24     2     3 17 23 16 24 15
#> 25     3     1 33 27 27 36 26
#> 26     3     1 27 31 35 30 35
#> 27     3     1 35 30 36 35 32
#> 28     3     1 26 29 29 34 36
#> 29     3     2 36 33 34 28 29
#> 30     3     2 28 32 28 31 31
#> 31     3     2 29 34 33 26 28
#> 32     3     2 30 28 30 32 30
#> 33     3     3 25 36 31 27 27
#> 34     3     3 34 25 32 33 25
#> 35     3     3 32 26 26 29 33
#> 36     3     3 31 35 25 25 34
                        # Permute obs within each block only.
  cbind(data, permindex(data, group="group", nsim=5)) 
#>    block group  1  2  3  4  5
#> 1      1     1 26  8  8 27 16
#> 2      1     1 28  5  5 14 13
#> 3      1     1  4 17  7 13  4
#> 4      1     1 13 29 18 16  3
#> 5      1     2 23 25 35  7 30
#> 6      1     2 36 28 10  8 18
#> 7      1     2 22  2 12 20  6
#> 8      1     2 24  3 11  6 31
#> 9      1     3 20  9 28 22 35
#> 10     1     3 30 21 27 10 12
#> 11     1     3 17 12 15 11 23
#> 12     1     3  7 35 25 34 21
#> 13     2     1 15 32 29 26  2
#> 14     2     1 14 18 20 28 28
#> 15     2     1  2  6 19 15 27
#> 16     2     1  3 19 17  4  1
#> 17     2     2  9 14 33 17 32
#> 18     2     2 34 27 34 30 20
#> 19     2     2 35  1  9 18 19
#> 20     2     2 10  4 24 32 29
#> 21     2     3 29 36 16 35 36
#> 22     2     3 19 22 26 36 22
#> 23     2     3 18 23 13 21 10
#> 24     2     3 32 34  2 12 24
#> 25     3     1 25 31  6 25 25
#> 26     3     1 27 30 31  2 26
#> 27     3     1 16 20 32  3 14
#> 28     3     1  1  7 30  1 15
#> 29     3     2 21 15 36 31  8
#> 30     3     2 33 13 23 19  5
#> 31     3     2 11 16 21  5 17
#> 32     3     2 12 26 22 29  7
#> 33     3     3 31 24  1 24 11
#> 34     3     3  5 10  3  9  9
#> 35     3     3  8 11 14 33 33
#> 36     3     3  6 33  4 23 34
                        # Permute groups as a whole block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, nsim=5))  # Free permutation.
#>    block group  1  2  3  4  5
#> 1      1     1  9 10 34  3 26
#> 2      1     1  6 29 32 27 34
#> 3      1     1 10 18 16 26 13
#> 4      1     1 27 19  2 21 36
#> 5      1     2 20 32  6 20 22
#> 6      1     2 30  8 27 25 31
#> 7      1     2 33  3  3 15  3
#> 8      1     2 17 28 21 10 14
#> 9      1     3 14 30  4 16  6
#> 10     1     3 31  2  7 29 15
#> 11     1     3  8  1 24 24  1
#> 12     1     3 34 16 29  4 27
#> 13     2     1 32 24  8  8 33
#> 14     2     1  7 26 12 19 32
#> 15     2     1  5  6 30 28  5
#> 16     2     1  2  7 14 23 12
#> 17     2     2 36 36 19 12 10
#> 18     2     2 24  4 11 11  8
#> 19     2     2 18 12 13  5  7
#> 20     2     2 13 33 20 13 21
#> 21     2     3 35  5 15  7  2
#> 22     2     3 26 35 36 22 17
#> 23     2     3 28 22 10 34 24
#> 24     2     3  4  9 18 33 30
#> 25     3     1 29 31 31 35 11
#> 26     3     1 16 17 22 30 23
#> 27     3     1 12 14 35  6 19
#> 28     3     1 23 23 17 17  9
#> 29     3     2 25 11 23  1 20
#> 30     3     2  1 21  1 14 18
#> 31     3     2 22 15  9 36  4
#> 32     3     2 11 13 26 31 28
#> 33     3     3  3 25 28  9 25
#> 34     3     3 19 20 25 32 29
#> 35     3     3 15 34  5 18 16
#> 36     3     3 21 27 33  2 35

  ## A randomisation test of nitrogen in a randomised block design:
  ## permutations within blocks only
  Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro, data = Oats)
  idx <- permindex(Oats, block = "Block", nsim = 199, seed = 1)
  permmodels(fm, perms = idx, ncore = 1)
#> 
#> Coefficients of (fixed) effects:
#>             Estimate Std. Error t value Pr(>|t|) Perm_p_value
#> (Intercept)  79.3889     3.7668 21.0759   0.0000           NA
#> Block.L      29.4823     4.6134  6.3906   0.0000        0.070
#> Block.Q      20.0124     4.6134  4.3379   0.0001        1.000
#> Block.C       5.3728     4.6134  1.1646   0.2486        0.005
#> Block^4       4.7718     4.6134  1.0343   0.3049        1.000
#> Block^5      -1.2651     4.6134 -0.2742   0.7848        0.005
#> nitro0.2     19.5000     5.3271  3.6605   0.0005        0.005
#> nitro0.4     34.8333     5.3271  6.5389   0.0000        0.005
#> nitro0.6     44.0000     5.3271  8.2597   0.0000        0.005
#> 
#> Note: Perm_p_value of t test is obtained using ‘199’ Freedman-Lane permutations.
#> 
#> ANOVA:
#> Analysis of Variance Table
#> 
#>  Response: yield
#>           Df   Sum Sq   Mean Sq F value Pr(>F) Perm_p_value
#> Block      5 15875.28 3175.0556 12.4317      0            1
#> nitro      3 20020.50 6673.5000 26.1297      0        0.005
#> Residuals 63 16090.17  255.3995                            
#> 
#> Note: Perm_p_value of F test is obtained using ‘199’ Freedman-Lane permutations.
#> 
```
