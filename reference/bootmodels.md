# Parametric Bootstrap Tests for Mixed Models

Parametric bootstrap p-values for the terms and coefficients of an
`lmer` or `glmer` model. Responses are simulated from the model fitted
without its fixed effects (only an intercept and the same random
effects: the global null hypothesis), the full model is refitted to
each, and a p-value is the proportion of simulated statistics at least
as large as the observed one. Useful for GLMMs and small samples, where
Wald or F tests can be poor.

## Usage

``` r
bootmodels(model, nsim = 999, type = c("I", "II", "III", "1", "2", "3"), seed = NULL,
           ncore = 1, prt = TRUE)
```

## Arguments

- model:

  A model fitted by `lmer` or `glmer`.

- nsim:

  Number of simulated data sets, the default is 999.

- type:

  The type of ANOVA table (as in
  [`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)):
  Type I gives sequential F statistics, Types II and III Wald chi-square
  statistics.

- seed:

  Optional random number seed, for reproducible results; the user's
  random number stream is left as it was.

- ncore:

  Number of cores for refitting, the default is 1.

- prt:

  Logical: print the tables? The default is TRUE.

## Value

Invisibly, a list with `permlist` (the coefficients and covariance
matrices of the refits), `COEFFICENTS` and `ANOVA` (the model's tables
with a `Boot_p_value` column). Attribute `"resampling"` is
`"parametric bootstrap"`.

## Details

The result has the same form as that of
[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md):
its `permlist` can be given to
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
or
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
(argument `permlist`) for bootstrap p-values of pairwise comparisons or
contrasts (`p.method` `"bootstrap"` in the `pairs` table).

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  # \donttest{
  bt <- bootmodels(fm, nsim = 99, seed = 1)   # use 999 or more in practice
#> 
#> Coefficients of (fixed) effects:
#>                            Estimate Std. Error t value Boot_p_value
#> (Intercept)                 80.0000     9.1070  8.7844           NA
#> nitro0.2                    18.5000     7.6829  2.4079         0.03
#> nitro0.4                    34.6667     7.6829  4.5122         0.01
#> nitro0.6                    44.8333     7.6829  5.8354         0.01
#> VarietyMarvellous            6.6667     9.7150  0.6862         0.44
#> VarietyVictory              -8.5000     9.7150 -0.8749         0.34
#> nitro0.2:VarietyMarvellous   3.3333    10.8653  0.3068         0.76
#> nitro0.4:VarietyMarvellous  -4.1667    10.8653 -0.3835         0.67
#> nitro0.6:VarietyMarvellous  -4.6667    10.8653 -0.4295         0.63
#> nitro0.2:VarietyVictory     -0.3333    10.8653 -0.0307         0.97
#> nitro0.4:VarietyVictory      4.6667    10.8653  0.4295         0.72
#> nitro0.6:VarietyVictory      2.1667    10.8653  0.1994         0.85
#> 
#> ANOVA:
#>               npar     Sum Sq   Mean Sq F value Boot_p_value
#> nitro            3 20020.5000 6673.5000 37.6857         0.01
#> Variety          2   526.0575  263.0288  1.4853         0.30
#> nitro:Variety    6   321.7500   53.6250  0.3028         0.89
#> 
#> Note: Boot_p_value is from ‘99’ parametric bootstrap samples under the global null model.
#> 
  predictmeans(fm, "nitro", permlist = bt, plot = FALSE, ncore = 1)$pairs
#>   level1 level2   estimate std.error df statistic p.value p.adj  conf.low
#> 1      0    0.2 -19.500000  4.435752 NA -4.396098    0.01  0.01 -29.32361
#> 2      0    0.4 -34.833333  4.435752 NA -7.852859    0.01  0.01 -43.98725
#> 3    0.2    0.4 -15.333333  4.435752 NA -3.456761    0.01  0.01 -23.69685
#> 4      0    0.6 -44.000000  4.435752 NA -9.919400    0.01  0.01 -53.65647
#> 5    0.2    0.6 -24.500000  4.435752 NA -5.523302    0.01  0.01 -32.96682
#> 6    0.4    0.6  -9.166667  4.435752 NA -2.066542    0.06  0.06 -18.51942
#>    conf.high      lsd  p.method
#> 1  -9.676390 9.823610 bootstrap
#> 2 -25.679414 9.153919 bootstrap
#> 3  -6.969817 8.363516 bootstrap
#> 4 -34.343534 9.656466 bootstrap
#> 5 -16.033178 8.466822 bootstrap
#> 6   0.186089 9.352756 bootstrap
  # }
```
