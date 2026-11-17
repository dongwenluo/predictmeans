# Spinal Bone Mineral Density of Girls

Repeated measurements of spinal bone mineral density of 230 girls and
young women aged 9 to 26 of four ethnic groups, from a study of bone
mineral acquisition (Bachrach et al. 1999); the data of M. P. Wand's
additive mixed model example (`gamm()`).
[`HRW::femSBMD`](https://rdrr.io/pkg/HRW/man/femSBMD.html) is a larger
version (1003 measurements).

## Usage

``` r
data(SBMDfem)
```

## Format

A data frame with 547 rows and 7 columns:

- idnum:

  identification number of the subject.

- spnbmd:

  spinal bone mineral density (g/cm\\^2\\).

- age:

  age (years).

- ethnicity:

  factor: Asian, Black, Hispanic, White.

- black, hispanic, white:

  indicators (1 = yes) of the ethnic groups other than Asian.

## Source

Bachrach, L. K., Hastie, T., Wang, M.-C., Narasimhan, B. and Marcus, R.
(1999). Bone mineral acquisition in healthy Asian, Hispanic, Black, and
Caucasian youth: a longitudinal study. *Journal of Clinical
Endocrinology and Metabolism*, 84, 4702-4712.

## Examples

``` r
data(SBMDfem)
SBMDfem$idnum <- factor(SBMDfem$idnum)
# a curve over age, ethnicity as a fixed effect, a random intercept per girl
fit <- semireg(spnbmd ~ s(age, bs = "os", k = 15) + ethnicity + (1 | idnum), data = SBMDfem, prt = FALSE)
anova(fit)       # ethnicity: F test with Kenward-Roger df; the smooth: likelihood ratio test
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: spnbmd
#>                   F Df Df.res    Pr(>F)    
#> (Intercept) 823.240  1 529.55 < 2.2e-16 ***
#> ethnicity    11.297  3 226.95 6.187e-07 ***
#> age         398.080  1 531.84 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);
#> reference: exact restricted LRT for the only random effect of a gaussian fit, otherwise 0.5 chi2(0) + 0.5 chi2(1)
#>   term variance logLik logLik.without   LRT p.value                      test
#>  s_age 0.003928  700.3          532.3 335.9  <2e-16 0.5 chi2(0) + 0.5 chi2(1)
# the ethnic groups at age 15, compared in pairs
predictmeans(fit, "ethnicity", covariate = c(age = 15), pairwise = TRUE, adj = "tukey")$pairs
#>     level1   level2    estimate  std.error       df  statistic      p.value
#> 1    Asian    Black -0.10617850 0.02065830 228.7585 -5.1397500 5.896664e-07
#> 2    Asian Hispanic -0.01316308 0.02164701 228.2603 -0.6080786 5.437399e-01
#> 3    Black Hispanic  0.09301542 0.02050921 231.5557  4.5352999 9.232858e-06
#> 4    Asian    White -0.02606166 0.02163625 223.5764 -1.2045365 2.296558e-01
#> 5    Black    White  0.08011684 0.02041221 223.5667  3.9249469 1.155142e-04
#> 6 Hispanic    White -0.01289857 0.02149288 225.9151 -0.6001325 5.490194e-01
#>          p.adj p.method
#> 1 3.509783e-06    tukey
#> 2 9.294398e-01    tukey
#> 3 5.436499e-05    tukey
#> 4 6.246109e-01    tukey
#> 5 6.639323e-04    tukey
#> 6 9.319286e-01    tukey
# growth of bone density: where is the curve rising?
attr(semipred(fit, "age", "age", deriv = 1, prt = FALSE)$pred_df, "significant")
#>   level from   to direction
#> 1        8.9 16.4       > 0
```
