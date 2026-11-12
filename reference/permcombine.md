# Combine Permutation Tests of Several Responses

Nonparametric combination (NPC) of permutation tests made with the same
permutations, for example the same model fitted to several responses
(traits, time points) measured on the same units: for each term or
coefficient, one p-value for the hypothesis that it has no effect on any
of the responses, and p-values adjusted for testing it on every
response.

## Usage

``` r
permcombine(..., method = c("fisher", "tippett"), tests = NULL)
```

## Arguments

- ...:

  Results of
  [`permmodels`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  with Freedman-Lane permutations (the default), all made with the same
  permutations: the same observations (rows of the data) and the same
  `seed` (or the same `perms`) and `nperm`. Names given to the arguments
  name the responses. A single list of such results can be given
  instead.

- method:

  How the p-values of the responses are combined: `"fisher"` (the
  default), \\-2 \sum \log p\\, good when several responses have
  effects; `"tippett"`, the smallest p-value, good when one response has
  a large effect.

- tests:

  Names of the terms and coefficients (`"coef:<name>"`) to combine; the
  default is every one that all the results have.

## Value

A data frame with one row per test: the permutation p-value of each
response, the combined p-value (`combined`), and the p-values adjusted
for the number of responses (columns `adj.<response>`): single-step
min-p adjusted p-values (Westfall and Young 1993), which control the
familywise error rate over the responses and allow for their
correlation.

## Details

Each response's statistic is turned into a p-value within each permuted
data set (its rank among all the permutations, the observed data counted
as one of them); the p-values of the responses are combined in every
permuted data set, and the combined p-value is the proportion of
permuted data sets whose combined statistic is at least the observed one
(Pesarin and Salmaso 2010). Because the same permutation of the units is
applied to every response, the correlation between the responses is kept
and needs no model. Permuted data sets in which a refit failed for any
response are left out.

Several tests within one family over time or over levels of a factor (a
max-statistic test over time points) are available from
[`contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
with `adj = "westfall"` and a `permlist`.

## References

Pesarin, F. and Salmaso, L. (2010). *Permutation Tests for Complex Data:
Theory, Applications and Software*. Wiley.

Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley.

Winkler, A. M., Webster, M. A., Brooks, J. C., Tracey, I., Smith, S. M.
and Nichols, T. E. (2016). Non-parametric combination and related
permutation tests for neuroimaging. *Human Brain Mapping*, 37,
1486-1511.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
mtcars$cyl <- factor(mtcars$cyl)
p1 <- permmodels(lm(mpg ~ wt + cyl, data = mtcars), nperm = 999, type = "II", seed = 1, prt = FALSE)
p2 <- permmodels(lm(qsec ~ wt + cyl, data = mtcars), nperm = 999, type = "II", seed = 1, prt = FALSE)
permcombine(mpg = p1, qsec = p2)
#>             mpg  qsec combined adj.mpg adj.qsec
#> wt        0.001 0.002    0.001   0.002    0.003
#> cyl       0.005 0.002    0.002   0.008    0.004
#> coef:wt   0.001 0.002    0.001   0.002    0.003
#> coef:cyl6 0.009 0.002    0.002   0.016    0.003
#> coef:cyl8 0.003 0.002    0.003   0.004    0.003
```
