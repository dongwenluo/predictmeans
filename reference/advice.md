# Advice on a Result

Results of predictmeans functions carry short advice when this fit and
these data call for a choice or a caution the user may not know about:
an unbalanced design, empty cells, a covariate aliased with a factor, a
singular fit, an extrapolated mean, few degrees of freedom, many tests,
too few permutations, and so on. Advice is printed at the end of the
result as an "Advice:" block (`!` a caution: the result may mislead as
it stands; `i` a note: a choice worth knowing). It is never generic, and
no test is run that was not asked for. `advice()` returns it as a data
frame, for code and AI assistants.

`options(predictmeans.advice = FALSE)` stops the printing; the advice is
still stored.

## Usage

``` r
advice(x)
```

## Arguments

- x:

  A result of
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
  [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
  [`anova()`](https://rdrr.io/r/stats/anova.html) of an `lmer` or
  `semireg` fit,
  [`pm_check()`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md),
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md),
  [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md),
  [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md),
  [`residvariogram()`](https://dongwenluo.github.io/predictmeans/reference/residvariogram.md),
  [`heritability()`](https://dongwenluo.github.io/predictmeans/reference/heritability.md),
  [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md),
  [`se_ranef()`](https://dongwenluo.github.io/predictmeans/reference/se_ranef.md),
  [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md),
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md),
  [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  or
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md).

## Value

A data frame with columns `id` (stable, listed below), `level`
(`"caution"` or `"note"`) and `text`; zero rows when there is no advice.

## Advice ids

|  |  |  |
|----|----|----|
| **id** | **from** | **when** |
| `nonestimable_means` | predictmeans, covariatemeans | some means cannot be estimated (NA) |
| `empty_cells` | predictmeans | means average over the observed cells only |
| `aliased_covariate` | predictmeans | a covariate is constant within a factor's levels |
| `large_se` | predictmeans | a mean's SE is over 3 times the median SE |
| `lsd_varies` | predictmeans | the average LSD judges some pair differently from its own test |
| `bias_adjust` | predictmeans, covariatemeans | back-transformed means are a typical cluster's or a median, not the mean |
| `contrast_common_cells` | contrastmeans | contrasts are estimated on the cells their means share |
| `small_df` | predictmeans, contrastmeans | degrees of freedom below 5 |
| `many_tests_unadjusted` | predictmeans, contrastmeans | 6 or more requested tests with `adj = "none"` |
| `singular_fit` | predictmeans, contrastmeans, R2_model, semireg | a variance at zero or a correlation at +-1 |
| `anova_unbalanced_typeI`, `_typeII`, `_typeIII` | anova | unbalanced data: what the type tests |
| `anova_nested_in_data` | anova | a factor nested in another in the data, without the nesting in the model |
| `h2_standard_unbalanced` | heritability | the standard form with unequal plots per genotype or incomplete blocks |
| `h2_forms_differ` | heritability | the forms differ by more than 0.1 |
| `h2_above_one` | heritability | a family heritability above 1 (beyond or within its sampling error) |
| `h2_genotype_zero` | heritability | the genotypic variance estimated as 0 |
| `h2_genotype_nested` | heritability | genotypes nested in a fixed factor |
| `anova_empty_cells` | anova | empty cells: no type tests a main effect over all cells |
| `anova_covariate_zero` | anova | Type III with a covariate in an interaction (tests at covariate = 0) |
| `check_<name>` | pm_check | a check worth a look (the note of its row) |
| `variogram_spatial` | residvariogram | neighbouring plots more alike than the model allows |
| `variogram_neighbours` | residvariogram | neighbouring plots less alike than the model allows (competition) |
| `variogram_trend` | residvariogram | the variogram leaves the envelope at longer lags only |
| `variogram_incomplete` | residvariogram | more than 10% of the grid positions empty |
| `ordinal_latent_average` | predictmeans | class probabilities at latent means averaged over other variables |
| `ordinal_mean_class` | predictmeans | the mean class scores the classes as equally spaced |
| `few_permutations` | permmodels, bootmodels | fewer than 999 resamples and a p-value at the floor or near 0.05 |
| `r2_latent_scale` | R2_model | a non-gaussian model: R2 on the link scale |
| `variance_at_zero` | varcomp, se_ranef | a variance estimated as zero |
| `few_re_levels` | varcomp, se_ranef | a random effect with fewer than 5 levels |
| `variance_test_boundary` | varcomp | how to test a variance against zero (`test = TRUE`) |
| `ranef_conditional_fixed` | se_ranef | conditional SEs of a random term that contains a fixed factor |
| `smooth_linear` | semireg | a smooth estimated as a straight line (variance 0) |
| `extrapolation` | semipred, covariatemeans, pm_slopes | predictions beyond the range of the data |
| `pointwise_regions` | semipred, covariatemeans | regions of difference read from pointwise limits |
| `smooth_other_covariate` | semipred | a smooth in the prediction does not involve `covariate` (adds a constant) |
| `nonestimable_slopes` | pm_slopes | the covariate is aliased with a factor: slopes NA, their differences kept |
| `slopes_equal_by_model` | pm_slopes | no covariate x term interaction: the slopes are equal by assumption |
| `slopes_at_value` | pm_slopes | a curved term or a link: the slopes depend on `at` |

## Examples

``` r
oats <- nlme::Oats
oats$nitro <- factor(oats$nitro)
fm <- lme4::lmer(log(yield) ~ nitro * Variety + (1 | Block), data = oats)
res <- predictmeans(fm, "nitro:Variety", pairwise = TRUE)
#> The response is log(yield): the means are also given back-transformed to the original scale (trans = I for the model scale only).
advice(res)
#>                      id level
#> 1           bias_adjust  note
#> 2 many_tests_unadjusted  note
#>                                                                                                                                                                                             text
#> 1 'response' is the median of the original scale (the geometric mean for a log), not the mean. If the mean on the original scale is wanted (totals, averages to report), use bias_adjust = TRUE.
#> 2                                      66 pairs tested without adjustment: some will be significant by chance. If they are a family, use adj = "holm" (any design) or "tukey"/"mvt" (all pairs).
```
