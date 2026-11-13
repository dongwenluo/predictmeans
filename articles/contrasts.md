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

The confidence limits that go with the adjusted p-values are in the
columns `adj.conf.low` and `adj.conf.high` (SAS `AdjLower` and
`AdjUpper`): 0 lies outside them exactly when `p.adj` is below the
level. `conf.low` and `conf.high` are the limits of each comparison on
its own. The step-down and false discovery rate methods (`"holm"`,
`"BH"`, …) have no simultaneous limits, so their pairs have no such
columns.
[`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md)
gives the same limits as `adj.lower.CL` and `adj.upper.CL` (below, for
the comparisons with a control).

``` r

pt <- predictmeans(fm, "N", adj = "tukey", plot = FALSE, prtnum = FALSE)$pairs
pt[c("level1", "level2", "estimate", "conf.low", "conf.high", "adj.conf.low",
     "adj.conf.high", "p.adj")]
#>   level1 level2   estimate  conf.low   conf.high adj.conf.low adj.conf.high
#> 1      0    0.2 -19.500000 -28.43406 -10.5659375    -31.33325     -7.666748
#> 2      0    0.4 -34.833333 -43.76740 -25.8992708    -46.66659    -23.000081
#> 3    0.2    0.4 -15.333333 -24.26740  -6.3992708    -27.16659     -3.500081
#> 4      0    0.6 -44.000000 -52.93406 -35.0659375    -55.83325    -32.166748
#> 5    0.2    0.6 -24.500000 -33.43406 -15.5659375    -36.33325    -12.666748
#> 6    0.4    0.6  -9.166667 -18.10073  -0.2326042    -20.99992      2.666585
#>          p.adj
#> 1 3.764251e-04
#> 2 3.362210e-09
#> 3 6.390152e-03
#> 4 4.263701e-12
#> 5 9.244650e-06
#> 6 1.797189e-01
```

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
#>         Estimate Std. Error t value      df Pr(>|t|) lower.CL upper.CL
#> 0.2 - 0  19.5000     4.4358  4.3961 45.0001    2e-04  10.5659  28.4341
#> 0.4 - 0  34.8333     4.4358  7.8529 45.0001    0e+00  25.8993  43.7674
#> 0.6 - 0  44.0000     4.4358  9.9194 45.0001    0e+00  35.0659  52.9341
#>         adj.lower.CL adj.upper.CL
#> 0.2 - 0       8.7182      30.2818
#> 0.4 - 0      24.0516      45.6151
#> 0.6 - 0      33.2182      54.7818
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
#>                        Estimate Std. Error t value      df Pr(>|t|) lower.CL
#> Golden Rain: linear    150.6667    24.2956  6.2014 45.0001   0.0000 101.7328
#> Golden Rain: quadratic  -8.3333    10.8653 -0.7670 45.0001   0.4471 -30.2172
#> Marvellous: linear     129.1667    24.2956  5.3165 45.0001   0.0000  80.2328
#> Marvellous: quadratic  -12.1667    10.8653 -1.1198 45.0001   0.2688 -34.0506
#> Victory: linear        162.1667    24.2956  6.6747 45.0001   0.0000 113.2328
#> Victory: quadratic     -10.5000    10.8653 -0.9664 45.0001   0.3390 -32.3839
#>                        upper.CL
#> Golden Rain: linear    199.6005
#> Golden Rain: quadratic  13.5506
#> Marvellous: linear     178.1005
#> Marvellous: quadratic    9.7172
#> Victory: linear        211.1005
#> Victory: quadratic      11.3839
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
#>       Variety df1      df2   F.ratio      p.value
#> 1 Golden Rain   3 45.00014 13.022765 2.969017e-06
#> 2  Marvellous   3 45.00014  9.952878 3.769239e-05
#> 3     Victory   3 45.00014 15.315709 5.256861e-07
pm_joint(fm, "N:Variety", by = "N")
#>     N df1     df2   F.ratio   p.value
#> 1   0   2 30.2308 1.2245400 0.3080589
#> 2 0.2   2 30.2308 1.8814476 0.1697812
#> 3 0.4   2 30.2308 0.2156336 0.8072596
#> 4 0.6   2 30.2308 0.4010511 0.6731259
```

The first table asks whether nitrogen matters within each variety, the
second whether the varieties differ at each nitrogen level. With
`ctrmatrix`, the rows of any contrast matrix are tested together,
e.g. the departure from a straight line (quadratic and cubic terms):

``` r

curv <- contrast_family(K, "poly")[c("quadratic", "cubic"), ]
pm_joint(fm, "N", ctrmatrix = curv)
#>   df1      df2  F.ratio   p.value
#> 1   2 45.00014 1.366873 0.2652818
```

For REML `lmer` fits the test is the Kenward–Roger F test (as
[`pbkrtest::KRmodcomp()`](https://rdrr.io/pkg/pbkrtest/man/kr__modcomp.html));
`df_method = "Satterthwaite"` gives lmerTest’s test, and `"containment"`
the SAS containment df.

`contrastmeans(..., joint = TRUE)` adds the same joint test to the t
tests of the individual contrasts (and, with a `permlist`, a permutation
p-value of the joint F):

``` r

contrastmeans(fm, "N", ctrmatrix = curv, joint = TRUE)[["Joint F test of the contrasts"]]
#>               df1     df2 F.ratio p.value
#> all contrasts   2 45.0001  1.3669  0.2653
```

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

## Equivalence tests

To show that two treatments do not differ by more than a margin that
matters in practice, give the margin with `equivalence` (a number `d`
for the region `-d` to `d`, or `c(lower, upper)`). The `pairs` table
then has `p.equiv`, the p-value of two one-sided t tests (TOST): below
`level`, the difference lies within the margin.

``` r

pe <- predictmeans(fm, "N", equivalence = 15, plot = FALSE, prtnum = FALSE)
pe$pairs[, c("level1", "level2", "estimate", "std.error", "p.equiv")]
#>   level1 level2   estimate std.error    p.equiv
#> 1      0    0.2 -19.500000  4.435752 0.84211224
#> 2      0    0.4 -34.833333  4.435752 0.99997388
#> 3    0.2    0.4 -15.333333  4.435752 0.52978457
#> 4      0    0.6 -44.000000  4.435752 0.99999998
#> 5    0.2    0.6 -24.500000  4.435752 0.98116610
#> 6    0.4    0.6  -9.166667  4.435752 0.09757394
```

`contrastmeans(..., equivalence = )` does the same for contrasts.

## Robust and clustered standard errors

When the error variance differs between groups and the model assumes it
is constant, `vcov = "HC3"` (or `"HC0"` to `"HC2"`) gives
heteroscedasticity-consistent standard errors for `lm`, `aov` and `glm`
models, and `vcov = ~ Block` cluster-robust ones. They work in
[`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md),
[`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
[`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md)
and
[`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md);
the df stay the model’s:

``` r

fl <- lm(yield ~ Block + N + Variety, data = Oats)
pmh <- predictmeans(fl, "N", vcov = "HC3", plot = FALSE, prtnum = FALSE)
pmh$pairs[, c("level1", "level2", "estimate", "std.error", "p.value")]
#>   level1 level2   estimate std.error      p.value
#> 1      0    0.2 -19.500000  4.864548 1.688713e-04
#> 2      0    0.4 -34.833333  5.500868 3.207893e-08
#> 3    0.2    0.4 -15.333333  5.859926 1.117925e-02
#> 4      0    0.6 -44.000000  5.212147 7.666905e-12
#> 5    0.2    0.6 -24.500000  5.589782 4.694785e-05
#> 6    0.4    0.6  -9.166667  6.151535 1.413395e-01
pmh$settings$vcov
#> [1] "HC3"
contrastmeans(fl, "N", ctrmatrix = rbind(c(-1, 1, 0, 0), c(0, -1, 1, 0)), vcov = ~ Block)[[1]]
#>      Estimate Std. Error t value df Pr(>|t|) lower.CL upper.CL
#> [1,]  19.5000     4.3162  4.5179 61        0  10.8693  28.1307
#> [2,]  15.3333     2.5716  5.9627 61        0  10.1912  20.4755
#> attr(,"Note")
#> [1] "The p-value is adjusted by 'none' method, if p-value = 0 means p-value < 0.0001. Standard errors from the covariance given by 'vcov' (CR1, clusters: Block)."
```

For any model, including mixed models, `vcov` also takes a covariance
matrix of the coefficients or a function that returns one (e.g. from the
clubSandwich package).

## More

- [`vignette("predictmeans", package = "predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/predictmeans.md):
  getting started.
- [`?contrastmeans`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md),
  [`?specify_contrast`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
  [`?contrast_family`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md),
  [`?pm_joint`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md),
  [`?pm_slopes`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md).
