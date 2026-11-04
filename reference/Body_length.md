# `Body_length` data

Simulated body lengths of beetles reared in containers within
populations, under two treatments in two habitats, of the kind used by
Nakagawa and Schielzeth (2013) to illustrate \\R^2\\ for mixed models.
Used in the `R2_glmm` vignette.

## Usage

``` r
data(Body_length)
```

## Format

A data frame with 960 rows and 6 columns:

- Population:

  factor, 12 populations.

- Container:

  factor, 120 containers (nested in populations).

- Sex:

  factor, `Female` or `Male`.

- Habitat:

  factor, `A` or `B`.

- Treatment:

  factor, `Cont` or `Exp`.

- BodyL:

  numeric, body length.

## References

Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for
obtaining \\R^2\\ from generalized linear mixed-effects models. *Methods
in Ecology and Evolution*, 4, 133–142.

## Examples

``` r
data(Body_length)
fm <- lmer(BodyL ~ Treatment + Habitat + Sex + (1|Population/Container), data = Body_length)
R2_glmm(fm)
#> # Adjusted R2 for Mixed Models
#>                Total                Fixed               Random 
#>             "73.39%"             "40.09%"              "33.3%" 
#> Container:Population           Population 
#>              "5.18%"             "28.13%" 
```
