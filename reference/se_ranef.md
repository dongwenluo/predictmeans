# Extract Standard Errors of Model Random Effects

These functions extract standard errors of model random effects from
objects returned by modeling functions.

## Usage

``` r
se_ranef(object, rand_term=NULL)
```

## Arguments

- object:

  object of `merMod` and `glmmTMB` fit

- rand_term:

  a name of random term in the model.

## Value

`se_ranef` gives a list of standard errors for `ranef`.

## Details

`se_ranef` extracts standard errors of the random effects from objects
returned by lmer, glmer and glmmTMB functions.

## References

This function is modified from function 'se.ranef' at package 'arm'.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard
