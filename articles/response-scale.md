# Link and response scale

For generalized linear (mixed) models, and for models of a transformed
response such as `log(y)`, the means are estimated on the scale of the
model, the *link scale*, and then back-transformed to the scale of the
data, the *response scale*. predictmeans does the tests on the link
scale, where the model is linear, and reports both scales:

- `means`: `estimate`, `std.error`, `conf.low`, `conf.high` on the link
  scale, and `response`, `response.std.error`, `response.conf.low`,
  `response.conf.high` on the response scale;
- `pairs`: differences on the link scale and, for a log link (or
  `trans = exp`), the **ratio** of the two means with its confidence
  limits; for a logit link, the **odds ratio**.

The response-scale standard error comes from the delta method, as
emmeans gives with `type = "response"`. The confidence limits are the
back-transformed link-scale limits, so they are not symmetric around the
mean.

``` r

library(predictmeans)
```

## Binomial GLM: probabilities and odds ratios

The `Drug` data give the number of mice responding (`R`) out of `N` for
four drugs at several doses.

``` r

data(Drug)
fm_drug <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
pm_drug <- predictmeans(fm_drug, "Drug", pairwise = TRUE)
pm_drug$means[, c("Drug", "estimate", "response", "response.std.error",
                  "response.conf.low", "response.conf.high")]
#>          Drug   estimate  response response.std.error response.conf.low
#> 1    Morphine  0.1739917 0.5433885         0.03034579         0.4835759
#> 2     Amidone  0.5660805 0.6378583         0.03166332         0.5738120
#> 3 Phenadoxone  2.4185555 0.9182314         0.01536127         0.8826288
#> 4   Pethidine -1.7574624 0.1471084         0.02125403         0.1101238
#>   response.conf.high
#> 1          0.6019756
#> 2          0.6973562
#> 3          0.9437231
#> 4          0.1938091
```

`estimate` is the log odds of a response at the mean log dose;
`response` is the probability. The comparisons are tested on the
log-odds scale, and each difference is the log of an odds ratio:

``` r

pm_drug$pairs[, c("level1", "level2", "estimate", "p.value",
                  "odds.ratio", "odds.ratio.conf.low", "odds.ratio.conf.high")]
#>        level1      level2   estimate      p.value odds.ratio
#> 1    Morphine     Amidone -0.3920888 2.951890e-02  0.6756441
#> 2    Morphine Phenadoxone -2.2445638 9.302021e-24  0.1059738
#> 3     Amidone Phenadoxone -1.8524750 6.799382e-16  0.1568485
#> 4    Morphine   Pethidine  1.9314541 4.069093e-18  6.8995356
#> 5     Amidone   Pethidine  2.3235429 2.042589e-23 10.2117898
#> 6 Phenadoxone   Pethidine  4.1760179 2.573901e-37 65.1060753
#>   odds.ratio.conf.low odds.ratio.conf.high
#> 1          0.47465285            0.9617449
#> 2          0.06840188            0.1641832
#> 3          0.10004256            0.2459098
#> 4          4.46008846           10.6732394
#> 5          6.46763307           16.1234642
#> 6         34.28884136          123.6204222
```

``` r

ggplot2::autoplot(pm_drug, type = "ci")
```

![](response-scale_files/figure-html/unnamed-chunk-5-1.png)

## Poisson GLMM: counts and rate ratios

Egg counts of individual females, in containers within populations; the
innermost random effect (one level per female) models overdispersion.

``` r

data(Poisson_ind)
fm_egg <- glmer(Egg ~ Treatment + Habitat + (1 | Population/Container/individual),
                family = poisson, data = Poisson_ind)
pm_egg <- predictmeans(fm_egg, "Treatment", pairwise = TRUE)
pm_egg$means[, c("Treatment", "estimate", "response", "response.std.error",
                 "response.conf.low", "response.conf.high")]
#>   Treatment estimate response response.std.error response.conf.low
#> 1       Exp 1.790814 5.994327          0.9834990          4.341871
#> 2      Cont 1.282368 3.605166          0.5965019          2.604225
#>   response.conf.high
#> 1           8.275686
#> 2           4.990820
pm_egg$pairs[, c("level1", "level2", "p.value", "ratio", "ratio.conf.low", "ratio.conf.high")]
#>   level1 level2     p.value    ratio ratio.conf.low ratio.conf.high
#> 1    Exp   Cont 1.59407e-21 1.662705        1.50571        1.836069
```

The `ratio` is the expected count under `Exp` divided by that under
`Cont`.

For a GLMM the back-transformed mean is the response of a *typical*
unit, one whose random effects are all zero. With a log link this is
smaller than the average count over all units, which is
`exp(estimate + v/2)`, `v` being the total variance of the random
effects. Ratios are the same on both definitions, which is one reason to
report them.

`bias_adjust = TRUE` gives that population mean instead (to second
order, as emmeans’ `bias.adjust`), with the random-effect variance as
default `sigma`; `average = "response"` averages the back-transformed
cells over the other factors (here Habitat), rather than
back-transforming the average on the link scale:

``` r

pm_pop <- predictmeans(fm_egg, "Treatment", bias_adjust = TRUE, average = "response")
pm_pop$means[, c("Treatment", "response", "response.std.error")]
#>   Treatment response response.std.error
#> 1       Exp 7.428952          1.2188717
#> 2      Cont 4.467991          0.7392531
pm_pop$settings[c("bias_adjust", "average")]
#> $bias_adjust
#> [1] 0.6532427
#> 
#> $average
#> [1] "response"
```

## A transformed response: `trans = exp`

When the response is transformed before modelling, give the
back-transformation with `trans` (and `transOff` for, e.g.,
`log(y + 1)`). Differences of log means are logs of ratios of geometric
means:

``` r

data(Body_length)
fm_len <- lmer(log(BodyL) ~ Treatment + Habitat + Sex + (1 | Population/Container),
               data = Body_length)
pm_len <- predictmeans(fm_len, "Treatment", pairwise = TRUE, trans = exp)
pm_len$means[, c("Treatment", "estimate", "response", "response.std.error")]
#>   Treatment estimate response response.std.error
#> 1      Cont 2.621319 13.75386          0.3655397
#> 2       Exp 2.644563 14.07730          0.3741359
pm_len$pairs[, c("level1", "level2", "p.value", "ratio", "ratio.conf.low", "ratio.conf.high")]
#>   level1 level2      p.value    ratio ratio.conf.low ratio.conf.high
#> 1   Cont    Exp 9.949753e-06 0.977024      0.9670482       0.9871028
```

## Binomial GLMM: odds ratio and R²

The `Clinical` data give the number of favourable responses out of the
patients in each treatment group at each of eight clinics.

``` r

data(Clinical)
fm_clin <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic),
                 family = binomial, data = Clinical)
pm_clin <- predictmeans(fm_clin, "Treatment", pairwise = TRUE)
pm_clin$means[, c("Treatment", "response", "response.conf.low", "response.conf.high")]
#>   Treatment  response response.conf.low response.conf.high
#> 1   control 0.2321034        0.07559056          0.5276915
#> 2      drug 0.3873968        0.14775968          0.6975684
pm_clin$pairs[, c("level1", "level2", "p.value", "odds.ratio",
                  "odds.ratio.conf.low", "odds.ratio.conf.high")]
#>    level1 level2    p.value odds.ratio odds.ratio.conf.low odds.ratio.conf.high
#> 1 control   drug 0.04362228  0.4779715            0.234923            0.9724753
```

[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
gives the share of the variance explained, as numbers: Piepho’s R² (the
default) or Nakagawa and Schielzeth’s marginal and conditional R².

``` r

R2_model(fm_clin)
#> # Adjusted R2 for mixed models (Piepho)
#>  Total  Fixed Random Clinic 
#> 75.71%  4.77% 70.94% 70.94% 
#> 
#> Advice:
#> i For a non-gaussian model the R2 is on the link (latent) scale: not
#>   comparable with an lm R2 or across families and links; use it to compare
#>   models of the same family.
R2_model(fm_clin, method = "nakagawa")
#> # R2 for mixed models (Nakagawa and Schielzeth)
#>    Marginal Conditional 
#>       2.71%      38.70% 
#> 
#> Advice:
#> i For a non-gaussian model the R2 is on the link (latent) scale: not
#>   comparable with an lm R2 or across families and links; use it to compare
#>   models of the same family.
```

The two need not agree. For counts out of totals, as here, Piepho’s R²
takes the residual variance from the model’s working weights: with up to
37 patients per group a proportion is measured precisely, so the random
clinic effects explain most of the rest. Nakagawa’s R² uses the variance
of the latent scale of a single yes/no response, π²/3 for the logit
link, whatever the group sizes. For binary (0/1) data both use the same
latent variance.

## Comparing two lines along a covariate

`covariatemeans(contr = )` gives the difference between two levels of a
factor along a covariate. Three points decide how it is reported.

**The test is on the model scale.** On the link (or transformed) scale
the model is additive, the estimates are closest to normal and the
standard error does not depend on the mean. With a log the difference
back-transforms exactly to a ratio, with a logit to an odds ratio: the
limits leave out 1 exactly when those of the difference leave out 0.
This is what is plotted by default.

**Back-transforming is for reading the size, not for testing.** The
difference in the units of the response (`contr_scale = "response"`) is
another quantity: its band is an approximation (the delta method), it
depends on the values of the other variables, and lines that are
parallel on the model scale are not parallel on the response scale. The
values where the lines differ are still taken from the model scale, and
marked.

**“Significant at 0.05” has three meanings.**

| Question | Answer |
|----|----|
| Do the lines differ at one covariate value chosen beforehand? | the pointwise band (`band = "pointwise"`) |
| Over which range of the covariate do they differ? | the simultaneous band (`band = "simultaneous"`, Scheffé; for straight lines the Working–Hotelling band) |
| Do they differ anywhere? | the joint F test (`settings$test` in the result) |

``` r

set.seed(1)
dose <- data.frame(trt = gl(2, 40, labels = c("control", "treated")), x = runif(80, 0, 2))
dose$count <- rpois(80, exp(1 + 0.2*dose$x + 0.6*dose$x*(dose$trt == "treated")))
fm_dose <- glm(count ~ trt*x, family = poisson, data = dose)
cmp <- covariatemeans(fm_dose, "trt", covariate = "x", contr = c("treated", "control"),
                      band = "simultaneous")
cmp
#> Difference treated - control of 'trt' along x
#> 95% limits; adj.conf.*: 95% simultaneous band; asymptotic (z) df
#> Curves equal: F(2, Inf) = 28.353, p = 4.86e-13
#> They differ (simultaneous band) for x: 0.5859 to 1.984 (>)
#> response: back-transformed by the inverse link; mean
#>        x estimate std.error  df statistic p.value  p.adj conf.low conf.high
#>  0.02678  0.05425    0.2532 Inf    0.2142  0.8304 0.9773  -0.4421    0.5506
#>  0.04666  0.06624    0.2496 Inf    0.2654  0.7907 0.9654  -0.4230    0.5555
#>  0.06672  0.07834    0.2460 Inf    0.3185  0.7501 0.9505  -0.4037    0.5604
#>  0.10666  0.10242    0.2387 Inf    0.4290  0.6679 0.9121  -0.3655    0.5703
#>  0.12357  0.11262    0.2357 Inf    0.4778  0.6328 0.8921  -0.3493    0.5746
#>  0.14136  0.12335    0.2325 Inf    0.5305  0.5957 0.8687  -0.3323    0.5790
#>  0.14660  0.12651    0.2316 Inf    0.5463  0.5848 0.8614  -0.3273    0.5804
#>  0.16849  0.13971    0.2277 Inf    0.6137  0.5394 0.8284  -0.3065    0.5859
#>  0.18654  0.15059    0.2245 Inf    0.6709  0.5023 0.7985  -0.2893    0.5905
#>  0.19893  0.15807    0.2223 Inf    0.7112  0.4770 0.7765  -0.2776    0.5937
#>  adj.conf.low adj.conf.high ratio ratio.conf.low ratio.conf.high response
#>       -0.5656        0.6741 1.056         0.6427           1.734   0.1462
#>       -0.5447        0.6772 1.068         0.6551           1.743   0.1804
#>       -0.5237        0.6804 1.081         0.6678           1.751   0.2155
#>       -0.4819        0.6868 1.108         0.6939           1.769   0.2877
#>       -0.4643        0.6895 1.119         0.7052           1.776   0.3192
#>       -0.4458        0.6925 1.131         0.7172           1.784   0.3528
#>       -0.4403        0.6933 1.135         0.7208           1.787   0.3628
#>       -0.4175        0.6970 1.150         0.7360           1.797   0.4053
#>       -0.3988        0.7000 1.163         0.7488           1.805   0.4410
#>       -0.3860        0.7021 1.171         0.7576           1.811   0.4659
#>  response.std.error response.conf.low response.conf.high response.adj.conf.low
#>              0.6803           -1.1872              1.480                -1.519
#>              0.6772           -1.1470              1.508                -1.477
#>              0.6740           -1.1056              1.537                -1.434
#>              0.6676           -1.0208              1.596                -1.346
#>              0.6648           -0.9839              1.622                -1.308
#>              0.6619           -0.9444              1.650                -1.267
#>              0.6610           -0.9327              1.658                -1.255
#>              0.6573           -0.8829              1.694                -1.204
#>              0.6542           -0.8411              1.723                -1.160
#>              0.6520           -0.8120              1.744                -1.130
#>  response.adj.conf.high
#>                   1.811
#>                   1.838
#>                   1.865
#>                   1.922
#>                   1.946
#>                   1.973
#>                   1.981
#>                   2.014
#>                   2.042
#>                   2.062
#> ... 118 more rows: as.data.frame(x) gives the table.
#> autoplot(x) draws the plot (layouts: ?autoplot.pm_curves).
ggplot2::autoplot(cmp)
```

![](response-scale_files/figure-html/unnamed-chunk-11-1.png)

The table has the columns of `predictmeans()$pairs`: `estimate` (the log
ratio), `conf.low`/`conf.high` (pointwise), `p.adj` with
`adj.conf.low`/`adj.conf.high` (the simultaneous band), `ratio` with its
limits, and the response-scale difference (`response`, …);
`as.data.frame(cmp)` returns it.

The pointwise band is narrower and would call a longer range
significant; it is right only for a single value of the covariate chosen
before seeing the data.

## More

- [`vignette("R2_glmm", package = "predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/R2_glmm.md):
  adjusted R² for mixed models.
- [`vignette("contrasts", package = "predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/contrasts.md):
  contrasts, joint tests and slopes; they are on the link scale too.
- [`?predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  (arguments `trans`, `transOff`, `count`),
  [`?R2_model`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md),
  [`?covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
  (section “Comparing two lines”).
