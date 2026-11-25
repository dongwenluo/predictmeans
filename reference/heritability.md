# Heritability on a Line-Mean (Entry-Mean) Basis

The share of the variance of genotype means that is genetic, from a
mixed model with random genotype effects: the standard form for balanced
trials and the mixed-model forms of Cullis et al. (2006), Piepho and
Möhring (2007) and Oakey et al. (2006) for unbalanced trials, incomplete
blocks or correlated errors.

## Usage

``` r
heritability(
  model,
  genotype,
  method = c("cullis", "piepho", "standard", "oakey", "reliability", "simulated",
    "family", "repeatability", "all"),
  relationship = NULL,
  within_variance = NULL,
  within_n = NULL,
  n_env = NULL,
  n_rep = NULL,
  piepho = c("fixed", "refit"),
  nsim = 10000,
  nboot = 0,
  ncore = 1,
  seed = NULL,
  conf_level = 0.95
)
```

## Arguments

- model:

  An `lmer` fit, or a gaussian `semireg` fit, with a random genotype
  term `(1 | genotype)`. A genotypic variance estimated as 0 gives 0 for
  every form, with a note.

- genotype:

  Name of the genotype (line, entry, family) factor.

- method:

  `"cullis"` (default), `"piepho"`, `"standard"`, `"oakey"`,
  `"reliability"`, `"simulated"`, `"family"`, `"repeatability"`, or
  `"all"` (the line-mean forms: standard to simulated).

- relationship:

  For `"family"`: the coefficient of relationship among members of a
  family.

- within_variance, within_n:

  For `"family"` and `"repeatability"` from plot means: the variance
  among plants within plots and the number of plants per plot (numbers,
  or names of columns of the data: the mean of the variances, the
  harmonic mean of the numbers). The single-plant variance is then
  \\\sigma^2_f + \sum_t \sigma^2_t + (\sigma^2_e - \sigma^2_w/n) +
  \sigma^2_w\\, as Williams et al. (2023, Example 6.2), with
  \\\sigma^2_e - \sigma^2_w/n\\ the plot variance; \\\sigma^2_w\\ is
  taken as known in the standard error.

- n_env, n_rep:

  For `"standard"`, a model with one genotype x environment term: the
  number of environments and of replicates per environment (`n_rep`
  alone: replicates, with the environments from the data); by default
  both from the data (harmonic means), the only way with several such
  terms.

- piepho:

  For `"piepho"`: `"fixed"` (default) or `"refit"`, see Details.

- nsim:

  For `"simulated"`: the number of simulated draws.

- nboot:

  For the mixed-model forms (Cullis to simulated): parametric bootstrap
  draws for their standard error and percentile limits (each draw refits
  the model; 0, the default, gives none). The other forms get theirs by
  the delta method from the covariance of the variance parameters, as
  [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md).

- ncore:

  Cores for the bootstrap; results are the same for any number.

- seed:

  Seed for the simulated and bootstrap draws.

- conf_level:

  Confidence level of the limits (logit-normal for the delta method, so
  they stay in (0, 1)).

## Value

A data frame of class `"pm_heritability"`: `method`, `estimate`,
`std.error`, `conf.low`, `conf.high`; `attr(, "components")` holds the
variance components and the mean variances of differences used.

## Details

With \\\sigma^2_g\\ the genotypic variance:

- `"standard"`:

  \\\sigma^2_g / (\sigma^2_g + \sum_t \sigma^2_t / n_t + \sigma^2_e /
  n)\\, over the random terms \\t\\ that contain the genotype (genotype
  x environment, ...), with \\n_t\\ the number of their levels per
  genotype and \\n\\ the number of plots per genotype (harmonic means
  over genotypes; give `n_env`, `n_rep` to set them). With one
  environment and \\r\\ replicates, \\\sigma^2_g/(\sigma^2_g +
  \sigma^2_e/r)\\; with \\e\\ environments, \\\sigma^2_g/(\sigma^2_g +
  \sigma^2\_{ge}/e + \sigma^2_e/(re))\\. Exact for balanced data only.

- `"cullis"`:

  \\1 - \bar v\_{BLUP}/(2\sigma^2_g)\\, \\\bar v\_{BLUP}\\ the mean
  variance of a difference of two genotype BLUPs (prediction error
  variance, from the inverse of the mixed model equations, as
  `se_ranef(type = "prediction")`).

- `"piepho"`:

  \\\sigma^2_g / (\sigma^2_g + \bar v\_{BLUE}/2)\\, \\\bar v\_{BLUE}\\
  the mean variance of a difference of two genotype BLUEs from the same
  model with genotype as a fixed effect (other random terms kept). By
  default (`piepho = "fixed"`) its variance components are those of the
  random-genotype fit, as Schmidt et al. (2019) do, so both terms of the
  ratio rest on the same estimates (no refit); `piepho = "refit"`
  estimates them again, as Piepho and Möhring (2007). On the oats alpha
  design of John and Williams (1995): 0.808 (Schmidt et al. 2019,
  Example 2) and 0.803 (Piepho and Möhring 2007, Example 1).
  (`inti::H2cal(emmeans = TRUE)` and heritable take emmeans's
  Kenward-Roger adjusted covariance: 0.797 there.)

- `"oakey"`:

  the mean of the non-zero eigenvalues of \\I - C\_{gg}/\sigma^2_g\\,
  \\C\_{gg}\\ the prediction error covariance of the genotype BLUPs (the
  zero ones are directions confounded with fixed effects, as the mean of
  all genotypes with the intercept). With independent genotypes it
  equals `"cullis"` when the intercept is the only fixed effect the
  genotypes are confounded with; with genotypes nested in a fixed factor
  (families in provenances) the two differ (a note says so), and with a
  relationship matrix.

- `"reliability"`:

  the mean reliability of the genotype BLUPs, \\\bar r^2 = m^{-1}\sum_i
  (1 - C\_{ii}/\sigma^2_g)\\ (Schmidt et al. 2019, Eq. 44); it ignores
  the covariances of the prediction errors, so it is a little below the
  others for small \\m\\.

- `"simulated"`:

  the expected squared correlation between true genotypic values and
  their BLUPs, simulated from their joint covariance (Piepho and Möhring
  2007; Schmidt et al. 2019, Eq. 41), with `nsim` draws.

- `"family"`:

  narrow-sense heritability of a single tree (plant, animal) from
  families: \\(\sigma^2_f/\rho) / (\sigma^2_f + \sum_t \sigma^2_t +
  \sigma^2_e)\\, \\\sigma^2_f\\ the family variance (`genotype` is the
  family factor), \\\rho\\ the `relationship` among members of a family
  (1/4 for half-sibs, 1/2 for full-sibs; often taken as 1/3 for
  open-pollinated progeny, which are partly full-sibs or selfs), and
  \\\sigma^2_t\\ the random terms with the family in them (family x
  block, i.e. plot; family x site). The other design terms (blocks,
  rows, columns) are left out of the phenotypic variance. Can exceed 1
  when \\\rho\\ is too small for the material. Family-mean heritability
  is the line-mean forms above with the family factor as `genotype`.

- `"repeatability"`:

  \\\sigma^2_s / (\sigma^2_s + \sum_t \sigma^2_t + \sigma^2_e)\\ with
  `genotype` the subject that has repeated records: the share of the
  variance of a single record due to the subject.

The four line-mean forms agree for balanced data with independent
genotypes; with imbalance, incomplete blocks or spatial errors they
differ, and Schmidt et al. (2019) recommend the mixed-model forms.
Narrow-sense heritability from a pedigree or genomic relationship matrix
needs a model with that matrix, which lme4 cannot fit: heritable
computes these forms (and the per-genotype \\H^2\_\Delta\\ of Schmidt et
al. 2019) from asreml fits with a relationship matrix. Its standard,
Cullis and Oakey values for lmer fits are those of this function
(tested); this function adds standard errors, the family and
repeatability forms, and advice.

## References

Cullis, B. R., Smith, A. B. and Coombes, N. E. (2006). On the design of
early generation variety trials with correlated data. *Journal of
Agricultural, Biological, and Environmental Statistics*, 11, 381-393.

Oakey, H., Verbyla, A., Pitchford, W., Cullis, B. and Kuchel, H. (2006).
Joint modeling of additive and non-additive genetic line effects in
single field trials. *Theoretical and Applied Genetics*, 113, 809-819.

Piepho, H.-P. and Möhring, J. (2007). Computing heritability and
selection response from unbalanced plant breeding trials. *Genetics*,
177, 1881-1888.

Schmidt, P., Hartung, J., Bennewitz, J. and Piepho, H.-P. (2019).
Heritability in plant breeding on a genotype-difference basis.
*Genetics*, 212(4), 991-1008.
[doi:10.1534/genetics.119.302134](https://doi.org/10.1534/genetics.119.302134)

## See also

[`varcomp`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
[`se_ranef`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md)

## Examples

``` r
library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro + (1 | Block) + (1 | Variety), data = Oats)
heritability(fm, "Variety", method = "all")
#> Heritability on a line-mean basis (95% limits)
#> Variance of Variety 27.45; residual 234.5
#> 
#>       method estimate std.error conf.low conf.high
#>     standard   0.7375    0.2668   0.1587    0.9767
#>       cullis   0.7375        NA       NA        NA
#>       piepho   0.7375        NA       NA        NA
#>        oakey   0.7375        NA       NA        NA
#>  reliability   0.4916        NA       NA        NA
#>    simulated   0.7553        NA       NA        NA
#> 
#> Simulated form: Monte Carlo SE 0.0029 (10000 draws; nsim for more).
# \donttest{
heritability(fm, "Variety", nboot = 199, seed = 1)   # SE of the Cullis form by bootstrap
#> Heritability on a line-mean basis (95% limits)
#> Variance of Variety 27.45; residual 234.5
#> 
#>  method estimate std.error conf.low conf.high
#>  cullis   0.7375    0.3436        0    0.9395
#> 
#> SE and limits of the mixed-model forms from 199 parametric bootstrap draws (percentile).
# }
```
