# Contrasts, joint tests and slopes

[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
gives the means of a model term and all their pairwise comparisons. This
vignette shows what to use when the questions are different: comparisons
with a control, trends, tests of several contrasts at once, slopes of a
covariate, and other adjustments for multiple testing. The results match
emmeans, pbkrtest and lmerTest (Kenward–Roger and Satterthwaite) and SAS
`LSMEANS`, `SLICE` and `ESTIMATE` statements.

We use the oats split plot again: three varieties on whole plots within
six blocks, four nitrogen levels on subplots.

``` r

library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$N <- factor(Oats$nitro)
fm <- lmer(yield ~ N * Variety + (1 | Block/Variety), data = Oats)
```

## Adjusting for multiple comparisons

`adj` sets how the p-values of the pairwise comparisons are adjusted.
Besides the [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
methods (`"holm"`, `"BH"`, `"bonferroni"`, …) and `"tukey"`, there are
`"sidak"`, `"scheffe"` and `"mvt"`. `"mvt"` is the single-step
multivariate t adjustment: it uses the correlations among the
comparisons, so it is exact for any family of contrasts. It is what SAS
`ADJUST=SIMULATE` estimates by simulation; here it is computed by
package mvtnorm.

``` r

p <- sapply(c("none", "tukey", "mvt", "sidak", "scheffe", "holm"), function(a)
  predictmeans(fm, "N", adj = a, plot = FALSE, prtnum = FALSE)$pairs$p.adj)
pm <- predictmeans(fm, "N", plot = FALSE, prtnum = FALSE)
data.frame(pair = paste(pm$pairs$level1, "-", pm$pairs$level2), signif(p, 3))
#>        pair     none    tukey      mvt    sidak  scheffe     holm
#> 1   0 - 0.2 6.66e-05 3.76e-04 3.66e-04 3.99e-04 1.01e-03 2.00e-04
#> 2   0 - 0.4 5.65e-10 3.36e-09 1.32e-09 3.39e-09 1.55e-08 2.82e-09
#> 3 0.2 - 0.4 1.21e-03 6.39e-03 6.38e-03 7.21e-03 1.34e-02 2.41e-03
#> 4   0 - 0.6 6.70e-13 4.26e-12 1.73e-13 4.02e-12 2.15e-11 4.02e-12
#> 5 0.2 - 0.6 1.58e-06 9.24e-06 5.43e-06 9.50e-06 3.12e-05 6.33e-06
#> 6 0.4 - 0.6 4.46e-02 1.80e-01 1.80e-01 2.39e-01 2.48e-01 4.46e-02
```

For all pairs of means Tukey and `"mvt"` agree, up to the accuracy of
the mvtnorm computation (about 1e-4, which shows only for very small
p-values); `"mvt"` also covers other families, as below. Scheffé
protects every possible contrast and is the most conservative.

## Families of contrasts

[`contrast_family()`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md)
builds standard sets of contrasts among the means of a term, with the
same coefficients and names as emmeans: `"pairwise"`, `"revpairwise"`,
`"trt.vs.ctrl"`, `"consec"`, `"poly"`, `"eff"` and `"del.eff"`. The
result goes to
[`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md).

Each nitrogen level against no nitrogen (Dunnett’s comparisons), with
the exact Dunnett adjustment:

``` r

K <- Kmatrix(fm, "N")
dunnett <- contrast_family(K, "trt.vs.ctrl")
dunnett
#>          0 0.2 0.4 0.6
#> 0.2 - 0 -1   1   0   0
#> 0.4 - 0 -1   0   1   0
#> 0.6 - 0 -1   0   0   1
```

``` r

contrastmeans(fm, "N", ctrmatrix = dunnett, adj = "mvt")[[1]]
#>         Estimate Std. Error t value df Pr(>|t|)
#> 0.2 - 0  19.5000     4.4358  4.3961 45    2e-04
#> 0.4 - 0  34.8333     4.4358  7.8529 45    0e+00
#> 0.6 - 0  44.0000     4.4358  9.9194 45    0e+00
#> attr(,"Note")
#> [1] "The p-value is adjusted by 'mvt' method, if p-value = 0 means p-value < 0.0001."
```

Linear and quadratic trends in nitrogen within each variety (`by`), with
integer orthogonal polynomial coefficients (the nitrogen levels are
equally spaced):

``` r

K2 <- Kmatrix(fm, "N:Variety")
trends <- contrast_family(K2, "poly", by = "Variety", max.degree = 2)
contrastmeans(fm, "N:Variety", ctrmatrix = trends)[[1]]
#>                        Estimate Std. Error t value df Pr(>|t|)
#> Golden Rain: linear    150.6667    24.2956  6.2014 45   0.0000
#> Golden Rain: quadratic  -8.3333    10.8653 -0.7670 45   0.4471
#> Marvellous: linear     129.1667    24.2956  5.3165 45   0.0000
#> Marvellous: quadratic  -12.1667    10.8653 -1.1198 45   0.2688
#> Victory: linear        162.1667    24.2956  6.6747 45   0.0000
#> Victory: quadratic     -10.5000    10.8653 -0.9664 45   0.3390
#> attr(,"Note")
#> [1] "The p-value is adjusted by 'none' method, if p-value = 0 means p-value < 0.0001."
```

[`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md)
(see the getting-started vignette) writes any other contrast by the
names of the means.

## Joint tests and slices

[`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md)
tests several contrasts at once with an F test. With `by`, it compares
the means of the term within each level of the `by` factor, which is SAS
`LSMEANS N*Variety / SLICE=Variety`:

``` r

pm_joint(fm, "N:Variety", by = "Variety")
#>       Variety df1 df2   F.ratio      p.value
#> 1 Golden Rain   3  45 13.022765 2.969047e-06
#> 2  Marvellous   3  45  9.952878 3.769265e-05
#> 3     Victory   3  45 15.315709 5.256928e-07
pm_joint(fm, "N:Variety", by = "N")
#>     N df1      df2   F.ratio   p.value
#> 1   0   2 30.23076 1.2245400 0.3080589
#> 2 0.2   2 30.23076 1.8814476 0.1697812
#> 3 0.4   2 30.23076 0.2156336 0.8072596
#> 4 0.6   2 30.23076 0.4010511 0.6731259
```

The first table asks whether nitrogen matters within each variety, the
second whether the varieties differ at each nitrogen level. With
`ctrmatrix`, the rows of any contrast matrix are tested together,
e.g. the departure from a straight line (quadratic and cubic terms):

``` r

curv <- contrast_family(K, "poly")[c("quadratic", "cubic"), ]
pm_joint(fm, "N", ctrmatrix = curv)
#>   df1 df2  F.ratio   p.value
#> 1   2  45 1.366873 0.2652819
```

For REML `lmer` fits the test is the Kenward–Roger F test (as
[`pbkrtest::KRmodcomp()`](https://rdrr.io/pkg/pbkrtest/man/kr__modcomp.html));
`df_method = "Satterthwaite"` gives lmerTest’s test, and `"containment"`
the SAS containment df.

## Slopes of a covariate

When nitrogen is used as a number, the question is how fast yield rises
with it in each variety.
[`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md)
gives the slopes with standard errors, df, confidence limits and tests,
and compares them (as
[`emmeans::emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.html)):

``` r

fm2 <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
pm_slopes(fm2, "nitro", "Variety", adj = "tukey")
#> Slopes of 'nitro' at 0.3 
#> 
#>      Variety slope std.error df conf.low conf.high statistic   p.value
#>  Golden Rain 75.33     11.86 51    51.53     99.14     6.353 5.726e-08
#>   Marvellous 64.58     11.86 51    40.78     88.39     5.446 1.490e-06
#>      Victory 81.08     11.86 51    57.28    104.89     6.838 9.828e-09
#> 
#> Pairwise differences of slopes (p.adj: tukey)
#> 
#>       level1     level2 estimate std.error df statistic p.value  p.adj
#>  Golden Rain Marvellous    10.75     16.77 51    0.6410  0.5244 0.7982
#>  Golden Rain    Victory    -5.75     16.77 51   -0.3429  0.7331 0.9373
#>   Marvellous    Victory   -16.50     16.77 51   -0.9839  0.3298 0.5903
```

The slopes are in yield per unit of nitrogen. With
[`poly()`](https://rdrr.io/r/stats/poly.html) or spline terms the slope
changes along the covariate; `at` sets where it is taken (the default is
the covariate’s mean).

## Any contrast from `coef()` and `vcov()`

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) give the means of a
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
result and their covariance matrix (Kenward–Roger adjusted for REML
`lmer` fits), named by the levels. Any contrast can be computed from
them:

``` r

pmv <- predictmeans(fm, "N:Variety", plot = FALSE, prtnum = FALSE)
b <- coef(pmv)
L <- (names(b) == "0.6:Victory") - (names(b) == "0.6:Golden Rain")
c(estimate = sum(L * b), SE = sqrt(drop(L %*% vcov(pmv) %*% L)))
#>  estimate        SE 
#> -6.333333  9.715020
```

They also let other packages test contrasts of the means.
marginaleffects needs the class listed first; `b1`, `b2`, … are the
means in the order of [`coef()`](https://rdrr.io/r/stats/coef.html):

``` r

options(marginaleffects_model_classes = "predictmeans")
i <- match(c("0.6:Victory", "0.6:Golden Rain"), names(b))
marginaleffects::hypotheses(pmv, hypothesis = sprintf("b%d - b%d = 0", i[1], i[2]))
#> 
#>  Hypothesis Estimate Std. Error      z Pr(>|z|)   S 2.5 % 97.5 %
#>   b12-b10=0    -6.33       9.72 -0.652    0.514 1.0 -25.4   12.7
```

marginaleffects uses normal-based tests; the standard error is the same
as above.

## More

- [`vignette("predictmeans", package = "predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/predictmeans.md):
  getting started.
- [`?contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
  [`?specify_contrast`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
  [`?contrast_family`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md),
  [`?pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md),
  [`?pm_slopes`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md).
