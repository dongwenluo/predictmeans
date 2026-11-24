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

[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
[`contrast_family`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md),
[`pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)

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
                         ((`0.6:Golden Rain` - `0:Golden Rain`) +
                          (`0.6:Marvellous` - `0:Marvellous`))/2)
  contrastmeans(fm, "nitro:Variety", ctrmatrix = cm)
#> Contrasts of 'nitro:Variety' (95% limits; p.adj: none; Kenward-Roger df)
#>            contrast estimate std.error df statistic   p.value     p.adj
#>    N_effect_Victory     47.0     7.683 45    6.1174 2.090e-07 2.090e-07
#>  N_effect_vs_others      4.5     9.410 45    0.4782 6.348e-01 6.348e-01
#>  conf.low conf.high
#>     31.53     62.47
#>    -14.45     23.45
```
