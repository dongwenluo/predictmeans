# Contrasts Written with Level Names

Build a contrast matrix for
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
by writing each contrast as an expression in the level names of a model
term, so there is no need to count columns, even when some treatment
combinations are missing. Each level name (in backticks when it is not a
syntactic name, e.g. `` `Wet:Bio` ``) stands for that mean.

## Usage

``` r
specify_contrast(.x, ..., .zero_sum = TRUE)
```

## Arguments

- .x:

  A result of
  [`Kmatrix`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
  or
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  for the term, or a character vector of its level names (as in the row
  names of `Kmatrix(model, modelterm)$K`).

- ...:

  Named expressions, one per contrast, e.g.
  `` Bio_vs_None = (`Wet:Bio` + `Dry:Bio`)/2 - (`Wet:None` + `Dry:None`)/2 ``.

- .zero_sum:

  If `TRUE` (the default), every contrast's weights must sum to zero.
  Use `FALSE` for other linear combinations of the means, such as the
  average of two means.

## Value

A matrix with one row per contrast (named as given) and one column per
level. Columns of levels that are not mentioned are zero. Means that
cannot be estimated (e.g. missing cells) cannot be used.

## Details

The arguments start with a dot so that contrasts can be given any name,
including `x`.

## See also

[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  cm <- specify_contrast(Kmatrix(fm, "nitro:Variety"),
    N_effect_Victory = `0.6:Victory` - `0:Victory`,
    N_effect_vs_others = (`0.6:Victory` - `0:Victory`) -
                         ((`0.6:Golden Rain` - `0:Golden Rain`) + (`0.6:Marvellous` - `0:Marvellous`))/2)
  contrastmeans(fm, "nitro:Variety", ctrmatrix = cm)
#> $`The t tests of the specified contrasts`
#>                    Estimate Std. Error t value df Pr(>|t|) lower.CL upper.CL
#> N_effect_Victory       47.0     7.6829  6.1174 45   0.0000  31.5257  62.4743
#> N_effect_vs_others      4.5     9.4097  0.4782 45   0.6348 -14.4520  23.4520
#> attr(,"Note")
#> [1] "The p-value is adjusted by ‘none’ method, if p-value = 0 means p-value < 0.0001."
#> 
#> $K
#>                    0:Golden Rain 0:Marvellous 0:Victory 0.2:Golden Rain
#> N_effect_Victory             0.0          0.0        -1               0
#> N_effect_vs_others           0.5          0.5        -1               0
#>                    0.2:Marvellous 0.2:Victory 0.4:Golden Rain 0.4:Marvellous
#> N_effect_Victory                0           0               0              0
#> N_effect_vs_others              0           0               0              0
#>                    0.4:Victory 0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#> N_effect_Victory             0             0.0            0.0           1
#> N_effect_vs_others           0            -0.5           -0.5           1
#> 
```
