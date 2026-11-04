# `Poisson_ind` data

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
