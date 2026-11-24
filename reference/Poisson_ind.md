# Egg Counts of Beetles (Simulated)

Simulated egg counts of individual female beetles in containers within
populations, under two treatments in two habitats, of the kind used by
Nakagawa and Schielzeth (2013). Used in the `R2_glmm` vignette for a
Poisson GLMM with an individual-level (overdispersion) random effect.

## Usage

``` r
data(Poisson_ind)
```

## Format

A data frame with 480 rows and 7 columns:

- Population:

  factor, 12 populations.

- Container:

  factor, 60 containers (nested in populations).

- Sex:

  factor, `Female` only.

- Habitat:

  factor, `B` or `A`.

- Treatment:

  factor, `Exp` or `Cont`.

- Egg:

  integer, number of eggs.

- individual:

  factor, individual within container.

## References

Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for
obtaining \\R^2\\ from generalized linear mixed-effects models. *Methods
in Ecology and Evolution*, 4, 133–142.

## Examples

``` r
data(Poisson_ind)
gm <- glmer(Egg ~ Treatment + Habitat + (1|Population/Container) + (1|individual:Container),
            family = poisson, data = Poisson_ind)
R2_model(gm)
#> # Adjusted R2 for mixed models (Piepho)
#>                Total                Fixed               Random 
#>               64.78%                9.13%               55.65% 
#> individual:Container Container:Population           Population 
#>               13.86%                3.09%               38.71% 
#> 
#> Advice:
#> i For a non-gaussian model the R2 is on the link (latent) scale: not
#>   comparable with an lm R2 or across families and links; use it to compare
#>   models of the same family.
R2_model(gm, over_disp = TRUE)   # without the overdispersion term
#> # Adjusted R2 for mixed models (Piepho)
#>  Total  Fixed Random 
#> 50.93%  9.13% 41.79% 
#> 
#> Advice:
#> i For a non-gaussian model the R2 is on the link (latent) scale: not
#>   comparable with an lm R2 or across families and links; use it to compare
#>   models of the same family.
```
