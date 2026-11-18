# Parametric Bootstrap Tests

Parametric bootstrap p-values for the terms and coefficients of a model.
Responses are simulated from the model fitted without its fixed effects
(only an intercept, any offset, and the same random effects and
covariance structure: the global null hypothesis), the full model is
refitted to each, and a p-value is the proportion of simulated
statistics at least as large as the observed one. Useful for GLMs and
GLMMs and small samples, where Wald or F tests can be poor.

## Usage

``` r
bootmodels(model, nsim = 999, type = c("I", "II", "III", "1", "2", "3"), seed = NULL,
           ncore = 1, ...)
```

## Arguments

- model:

  A model fitted by `lm`, `aov` (without `Error()`), `glm` (a family
  with a distribution: gaussian, binomial, Poisson, Gamma; also
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html)), `gls`,
  `lme`, `lmer`, `glmer` or `glmmTMB`. Responses of `lme` and `gls`
  models are simulated from the null model's covariance of the data
  (random effects, correlation and variance structures).

- nsim:

  Number of simulated data sets, the default is 999.

- type:

  The type of ANOVA table (as in
  [`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)):
  Type I gives sequential F statistics (for a `glm` with a known
  dispersion, the deviance; `glmmTMB` has none and uses Type II), Types
  II and III Wald chi-square statistics (`lm`: F; `glm`: likelihood
  ratio).

- seed:

  Optional random number seed, for reproducible results; the user's
  random number stream is left as it was.

- ncore:

  Number of cores for refitting, the default is 1.

- ...:

  Not used: `prt` (before 2.0) stops with its 2.0 spelling;
  [`print()`](https://rdrr.io/r/base/print.html) of the result shows the
  tables.

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

In a simulation study (500 data sets per case, 5% level;
`dev/sim-boot.R` in the package sources) the bootstrap corrected the
liberal Wald chi-square test of a whole-plot factor in a split plot with
6 blocks (10% to 6-7%, as Kenward-Roger) and kept the level in Poisson
and binomial GLMMs, also when other terms had effects. Like the
likelihood-ratio test it relies on the model: with overdispersed counts
fitted as Poisson it rejected 37% of true null hypotheses, where
[`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
(Freedman-Lane on the working scale) gave 7%.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  # \donttest{
  bt <- bootmodels(fm, nsim = 99, seed = 1)   # use 999 or more in practice
  predictmeans(fm, "nitro", permlist = bt, plot = FALSE, ncore = 1)$pairs
#>   level1 level2   estimate std.error df statistic p.value p.adj  conf.low
#> 1      0    0.2 -19.500000  4.435752 NA -4.396098    0.01  0.01 -29.32367
#> 2      0    0.4 -34.833333  4.435752 NA -7.852859    0.01  0.01 -43.98725
#> 3    0.2    0.4 -15.333333  4.435752 NA -3.456761    0.01  0.01 -23.69685
#> 4      0    0.6 -44.000000  4.435752 NA -9.919400    0.01  0.01 -53.65647
#> 5    0.2    0.6 -24.500000  4.435752 NA -5.523302    0.01  0.01 -32.96682
#> 6    0.4    0.6  -9.166667  4.435752 NA -2.066542    0.06  0.06 -18.51942
#>    conf.high      lsd  p.method
#> 1  -9.676325 9.823675 bootstrap
#> 2 -25.679414 9.153919 bootstrap
#> 3  -6.969817 8.363516 bootstrap
#> 4 -34.343534 9.656466 bootstrap
#> 5 -16.033178 8.466822 bootstrap
#> 6   0.186089 9.352756 bootstrap
  # }
```
