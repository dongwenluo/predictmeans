# R-squared for Mixed Models as Numbers

The R-squared of a mixed model as a numeric result. `method = "piepho"`
(the default) gives the same values as
[`R2_glmm`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md)
(Piepho, 2023), as proportions: total, fixed, random and each
random-effect component. `method = "nakagawa"` gives the marginal (fixed
effects) and conditional (fixed and random effects) R-squared of
Nakagawa and Schielzeth (2013), with random slopes as Johnson (2014):
the variance of each random-effect term is the mean of \\z_i' \Sigma
z_i\\ over the observations.

## Usage

``` r
R2_model(model, method = c("piepho", "nakagawa"), over_disp = FALSE)
```

## Arguments

- model:

  A model fitted by `lmer`, `glmer` or `glmmTMB`.

- method:

  `"piepho"` or `"nakagawa"`.

- over_disp:

  For `"piepho"`: as in
  [`R2_glmm`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md).

## Value

A named numeric vector of proportions (class `"pm_r2"`, printed as
percentages). For `"nakagawa"`, `Marginal` and `Conditional`, with
attribute `"variances"`: the fixed-effect, random-effect and
observation-level variances.

## Details

For `"nakagawa"` the observation-level variance is the residual variance
for gaussian models and the variance of the latent distribution for
binomial models: \\\pi^2/3\\ (logit), 1 (probit) or \\\pi^2/6\\
(cloglog). Other families have several competing approximations; use
`method = "piepho"` for them.

## References

Piepho, H.-P. (2023), An adjusted coefficient of determination (R2) for
generalized linear mixed models in one go, *Biometrical Journal*, 65(7),
2200290.

Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for
obtaining R2 from generalized linear mixed-effects models, *Methods in
Ecology and Evolution*, 4, 133-142.

Johnson, P. C. D. (2014), Extension of Nakagawa & Schielzeth's R2GLMM to
random slopes models, *Methods in Ecology and Evolution*, 5, 944-946.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  R2_model(fm)
#> # Adjusted R2 for mixed models (Piepho)
#>        Total        Fixed       Random      Subject Subject.Days 
#>       79.36%       28.02%       51.34%       18.32%       30.43% 
  R2_model(fm, "nakagawa")
#> # R2 for mixed models (Nakagawa and Schielzeth)
#>    Marginal Conditional 
#>      27.87%      79.92% 
```
