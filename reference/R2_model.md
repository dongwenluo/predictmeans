# R-squared of a Model as Numbers

The R-squared of a model as a numeric result, for every model class of
[`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md).
`method = "piepho"` (the default) gives the same values as
[`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
of predictmeans 1.x (removed in 2.0; Piepho, 2023), as proportions:
total, fixed, random and each random-effect component.
`method = "nakagawa"` gives the marginal (fixed effects) and conditional
(fixed and random effects) R-squared of Nakagawa and Schielzeth (2013),
with random slopes as Johnson (2014): the variance of each random-effect
term is the mean of \\z_i' \Sigma z_i\\ over the observations.

## Usage

``` r
R2_model(model, method = c("piepho", "nakagawa"), over_disp = FALSE)
```

## Arguments

- model:

  A model of any class
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  handles: `lm`, `aov` (also with `Error()`), `glm`, `gls`, `lme`,
  `lmer`, `glmer`, `glmmTMB`, `gam` or a
  [`semireg`](https://dongwenluo.github.io/predictmeans/reference/semireg.md)
  fit.

- method:

  `"piepho"` (every class) or `"nakagawa"` (`lmer`, `glmer`, `glmmTMB`
  and `lme`; for `lme` with a variance or correlation structure the
  observation-level variance is the mean residual variance).

- over_disp:

  For `"piepho"` with `lmer`, `glmer` and `glmmTMB` models: `TRUE`
  leaves an observation-level random effect (one level per observation,
  modelling overdispersion) out of the random share; the model must have
  one. Default `FALSE`.

## Value

A named numeric vector of proportions (class `"pm_r2"`, printed as
percentages). For `"nakagawa"`, `Marginal` and `Conditional`, with
attribute `"variances"`: the fixed-effect, random-effect and
observation-level variances. The Nakagawa values agree with
`performance::r2_nakagawa()` (`tests/testthat/test-r2-model.R`). Notes
on the choice of method are kept in `attr(, "advice")`
([`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)).

## Details

Piepho's R-squared is one definition for every class: the variance of
the fixed part of the model (less the bias from estimating the
coefficients) and of the random effects, as shares of the average
variance of the data. For models without random effects (`lm`, `aov`,
`glm`, `gls`) `Total` equals `Fixed` and `Random` is 0; for `lm` and
`aov` it is the usual adjusted R-squared. `glm` models are measured on
the link scale (working variances; the latent distribution for a binary
response); `gls` and `lme` models with any correlation and variance
structure; `aov` models with `Error()` from the stratum variance
components of
[`varcomp`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md).
In `gam` models the smooths are part of the fixed share, except
random-effect smooths (`bs = "re"`).

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

## See also

[`varcomp`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
[`advice`](https://dongwenluo.github.io/predictmeans/reference/advice.md)

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
  R2_model(lm(Reaction ~ Days, data = sleepstudy))   # the adjusted R-squared
#> # Adjusted R2 (Piepho)
#>  Total  Fixed Random 
#> 28.25% 28.25%  0.00% 
  R2_model(glm(breaks ~ wool*tension, family = poisson, data = warpbreaks))
#> # Adjusted R2 (Piepho)
#>  Total  Fixed Random 
#> 63.54% 63.54%  0.00% 
#> 
#> Advice:
#> i For a non-gaussian model the R2 is on the link (latent) scale: not
#>   comparable with an lm R2 or across families and links; use it to compare
#>   models of the same family.
```
