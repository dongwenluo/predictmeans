# Adjusted R² for generalized linear mixed models

Code

``` r

library(predictmeans)
library(tidyverse)
```

## The function `R2_model()`

[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
gives an adjusted coefficient of determination (R²) for linear and
generalized linear mixed models, by the method of Piepho (2023). In
predictmeans 1.x this was
[`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md);
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
gives the same values and also works for models without random effects
(`lm`, `glm`, `gls`) and for `lme`, `glmmTMB`, `gam` and `semireg` fits.

Piepho’s R² starts from the average variance of the data: the mean of
half the expected squared difference over all pairs of observations,
which equals the expected sample variance. This is split into the
variance of the fixed effects and the variance of each random-effect
term; the shares of these in the total are the fixed, random and total
R², and the random share is split further into one share per
random-effect term. Three features matter in practice:

- **Adjusted.** The variance of the fixed effects is corrected for the
  bias from estimating them, so that for a linear model without random
  effects the result is exactly the usual adjusted R². Because of the
  correction the fixed share can be slightly negative when the fixed
  effects explain almost nothing.
- **One fit.** Only the model of interest is fitted; no null model is
  needed (hence “in one go”).
- **Covariance taken into account.** The average variance accounts for
  the covariance between observations that share a random effect, which
  Nakagawa and Schielzeth’s R² ignores.

For a generalized linear mixed model the variances are taken on the
scale of the linear predictor, where the effects are additive. The
residual variance on that scale depends on the family and link, and
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
chooses it automatically: for binary data the exact variance of the
latent distribution ($`\pi^2/3`$ for the logit link, 1 for the probit
link, $`\pi^2/6`$ for the complementary log-log link); for other data a
first-order approximation (for Poisson counts with a log link, the
inverse of the mean).

Every example below is one of Piepho (2023). The values from
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
are compared with those in the paper, which were computed with SAS PROC
GLIMMIX.

Code

``` r

args(R2_model)
```

    ## function (model, method = c("piepho", "nakagawa"), over_disp = FALSE) 
    ## NULL

- `method = "piepho"` (the default) gives the total, fixed and random
  shares and each random-effect component. `method = "nakagawa"` gives
  the marginal and conditional R² of Nakagawa and Schielzeth (2013).
- `over_disp = TRUE` treats an observation-level random effect (one
  level per observation, used to model overdispersion) as part of the
  residual variance, not as explained variance.

A caution, which the printed result repeats as advice for non-gaussian
models: an R² on the latent scale is not comparable with the R² of a
linear model, nor across families or links. Use it to compare models of
the same family and link fitted to the same data, at the same level of
aggregation (see Example 4).

## Example 1: body length (linear mixed model)

The `Body_length` data (Nakagawa and Schielzeth, 2013) are simulated
body lengths of beetles in 120 containers nested in 12 populations, with
two treatments, two habitats and both sexes. Piepho (2023, Table 2)
gives 40.09% (fixed), 33.30% (random) and 73.39% (total), and Nakagawa
and Schielzeth’s marginal and conditional R² as 39.16% and 74.09%.

Code

``` r

data(Body_length)
str(Body_length)
```

    ## 'data.frame':    960 obs. of  6 variables:
    ##  $ Population: Factor w/ 12 levels "1","10","11",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Container : Factor w/ 120 levels "1","10","100",..: 1 1 1 1 1 1 1 1 33 33 ...
    ##  $ Sex       : Factor w/ 2 levels "Female","Male": 1 1 1 1 1 1 1 1 2 2 ...
    ##  $ Habitat   : Factor w/ 2 levels "A","B": 1 1 1 1 2 2 2 2 1 1 ...
    ##  $ Treatment : Factor w/ 2 levels "Cont","Exp": 1 2 1 2 1 2 1 2 1 2 ...
    ##  $ BodyL     : num  10.8 11.2 11.6 15.6 12.4 ...

Code

``` r

mod_BodyL <- lmer(BodyL ~ Treatment + Habitat + Sex + (1 | Population/Container),
                  data = Body_length)
summary(mod_BodyL, corr = FALSE)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: BodyL ~ Treatment + Habitat + Sex + (1 | Population/Container)
    ##    Data: Body_length
    ## 
    ## REML criterion at convergence: 3059
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4339 -0.6252 -0.0314  0.6740  3.1162 
    ## 
    ## Random effects:
    ##  Groups               Name        Variance Std.Dev.
    ##  Container:Population (Intercept) 0.2345   0.4842  
    ##  Population           (Intercept) 1.3787   1.1742  
    ##  Residual                         1.1968   1.0940  
    ## Number of obs: 960, groups:  Container:Population, 120; Population, 12
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  15.21587    0.35184  43.247
    ## TreatmentExp  0.31398    0.07062   4.446
    ## HabitatB      0.08820    0.07062   1.249
    ## SexMale      -2.66858    0.11315 -23.585

Code

``` r

R2_model(mod_BodyL)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##                Total                Fixed               Random 
    ##               73.39%               40.09%               33.30% 
    ## Container:Population           Population 
    ##                5.18%               28.13%

The values agree with the paper exactly. For comparison, Nakagawa and
Schielzeth’s marginal (fixed effects) and conditional (fixed and random
effects) R², which are the counterparts of the fixed and total shares.
Their fixed share is larger because it has no correction for the
estimation of the fixed effects, so it slightly overstates the variance
explained by them:

Code

``` r

R2_model(mod_BodyL, method = "nakagawa")
```

    ## # R2 for mixed models (Nakagawa and Schielzeth)
    ##    Marginal Conditional 
    ##      39.16%      74.09%

## Example 2: egg counts (Poisson GLMM with overdispersion)

The `Poisson_ind` data are simulated egg counts of females from the same
design. Each female is one level of `individual` within her container,
so `(1 | Population/Container/individual)` includes an observation-level
random effect that models overdispersion. With `over_disp = TRUE` its
variance counts as residual, not explained, variance, as in the paper.
Piepho (2023, Table 2, Laplace fit) gives 9.13% (fixed), 41.80% (random)
and 50.93% (total).

Code

``` r

data(Poisson_ind)
str(Poisson_ind)
```

    ## 'data.frame':    480 obs. of  7 variables:
    ##  $ Population: Factor w/ 12 levels "1","10","11",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Container : Factor w/ 60 levels "1","101","103",..: 1 1 1 1 1 1 1 1 22 22 ...
    ##  $ Sex       : Factor w/ 1 level "Female": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Habitat   : Factor w/ 2 levels "B","A": 2 2 2 2 1 1 1 1 2 2 ...
    ##  $ Treatment : Factor w/ 2 levels "Exp","Cont": 2 1 2 1 2 1 2 1 2 1 ...
    ##  $ Egg       : int  6 3 10 11 8 11 8 9 10 7 ...
    ##  $ individual: Factor w/ 8 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 1 2 ...

Code

``` r

mod_poisson <- glmer(Egg ~ Treatment + Habitat + (1 | Population/Container/individual),
                     family = poisson, data = Poisson_ind)
summary(mod_poisson, corr = FALSE)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: Egg ~ Treatment + Habitat + (1 | Population/Container/individual)
    ##    Data: Poisson_ind
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##    2364.8    2389.8   -1176.4    2352.8       474 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.04030 -0.59161 -0.08938  0.44767  1.71811 
    ## 
    ## Random effects:
    ##  Groups                          Name        Variance Std.Dev.
    ##  individual:Container:Population (Intercept) 0.10000  0.3162  
    ##  Container:Population            (Intercept) 0.02262  0.1504  
    ##  Population                      (Intercept) 0.30411  0.5515  
    ## Number of obs: 480, groups:  
    ## individual:Container:Population, 480; Container:Population, 60; Population, 12
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.83998    0.16584  11.095   <2e-16 ***
    ## TreatmentCont -0.50845    0.05046 -10.077   <2e-16 ***
    ## HabitatA      -0.09833    0.04980  -1.975   0.0483 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Code

``` r

R2_model(mod_poisson, over_disp = TRUE)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##  Total  Fixed Random 
    ## 50.93%  9.13% 41.79% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

## Example 3: colour morph (binary GLMM)

The `Binary_color` data are simulated colour morphs (0 or 1) of males
from the same design. Piepho (2023, Table 2) gives 7.46% (fixed), 21.99%
(random) and 29.46% (total).

Code

``` r

load("data/Binary_color.RData")
str(Binary_color)
```

    ## 'data.frame':    480 obs. of  6 variables:
    ##  $ Population: Factor w/ 12 levels "1","10","11",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Container : Factor w/ 60 levels "10","100","102",..: 17 17 17 17 17 17 17 17 28 28 ...
    ##  $ Sex       : Factor w/ 1 level "Male": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Habitat   : Factor w/ 2 levels "B","A": 2 2 2 2 1 1 1 1 2 2 ...
    ##  $ Treatment : Factor w/ 2 levels "Exp","Cont": 2 1 2 1 2 1 2 1 2 1 ...
    ##  $ Colour    : int  0 0 0 1 0 0 0 1 0 0 ...

Code

``` r

mod_color <- glmer(Colour ~ Treatment + Habitat + (1 | Population/Container),
                   family = binomial, data = Binary_color)
summary(mod_color, corr = FALSE)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: Colour ~ Treatment + Habitat + (1 | Population/Container)
    ##    Data: Binary_color
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     573.1     594.0    -281.6     563.1       475 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3763 -0.7204 -0.3274  0.7662  3.2286 
    ## 
    ## Random effects:
    ##  Groups               Name        Variance Std.Dev.
    ##  Container:Population (Intercept) 0.005835 0.07639 
    ##  Population           (Intercept) 1.110439 1.05377 
    ## Number of obs: 480, groups:  Container:Population, 60; Population, 12
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     0.4369     0.3551   1.230  0.21862    
    ## TreatmentCont  -1.0120     0.2156  -4.694 2.68e-06 ***
    ## HabitatA       -0.6762     0.2125  -3.181  0.00147 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Code

``` r

R2_model(mod_color)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##                Total                Fixed               Random 
    ##               29.46%                7.46%               21.99% 
    ## Container:Population           Population 
    ##                0.12%               21.87% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

## Example 4: clinical trial (binomial GLMM, three links)

The `Clinical` data (Beitler and Landis, 1985) give the number of
patients with a favourable response out of the patients in each
treatment group (drug or control) at each of eight clinics. The model
has a random effect for clinics and for treatments within clinics.

### Data at clinic level

Code

``` r

data(Clinical)
Clinical$Unfavorable <- Clinical$Total - Clinical$Favorable
str(Clinical)
```

    ## 'data.frame':    16 obs. of  5 variables:
    ##  $ Clinic     : Factor w/ 8 levels "1","2","3","4",..: 1 1 2 2 3 3 4 4 5 5 ...
    ##  $ Treatment  : Factor w/ 2 levels "control","drug": 2 1 2 1 2 1 2 1 2 1 ...
    ##  $ Favorable  : int  11 10 16 22 14 7 2 1 6 0 ...
    ##  $ Total      : int  36 37 20 32 19 19 16 17 17 12 ...
    ##  $ Unfavorable: int  25 27 4 10 5 12 14 16 11 12 ...

Code

``` r

fit_clinic <- function(link) {
  glmer(cbind(Favorable, Unfavorable) ~ Treatment + (1 | Clinic/Treatment),
        family = binomial(link = link), data = Clinical)
}
mod_clinic_logit <- fit_clinic("logit")
summary(mod_clinic_logit, corr = FALSE)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: cbind(Favorable, Unfavorable) ~ Treatment + (1 | Clinic/Treatment)
    ##    Data: Clinical
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##      82.2      85.3     -37.1      74.2        12 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.4300 -0.6745 -0.2207  0.4670  1.2369 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance Std.Dev.
    ##  Treatment:Clinic (Intercept) 0.005739 0.07576 
    ##  Clinic           (Intercept) 1.931755 1.38988 
    ## Number of obs: 16, groups:  Treatment:Clinic, 16; Clinic, 8
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)    -1.1996     0.5597  -2.143   0.0321 *
    ## Treatmentdrug   0.7423     0.3215   2.308   0.0210 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Code

``` r

R2_model(mod_clinic_logit)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##            Total            Fixed           Random Treatment:Clinic 
    ##           75.73%            4.69%           71.04%            0.23% 
    ##           Clinic 
    ##           70.81% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Code

``` r

R2_model(fit_clinic("probit"))
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##            Total            Fixed           Random Treatment:Clinic 
    ##           77.90%            5.00%           72.90%            0.05% 
    ##           Clinic 
    ##           72.84% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Code

``` r

R2_model(fit_clinic("cloglog"))
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##            Total            Fixed           Random Treatment:Clinic 
    ##           72.63%            3.54%           69.09%            0.00% 
    ##           Clinic 
    ##           69.09% 
    ## 
    ## Advice:
    ## ! The fit is singular: a variance estimated as zero gives no share to its
    ##   random effect, so the random part of R2 may be understated. Consider
    ##   simplifying the random effects.
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Piepho (2023, Table 3) gives, at the clinic-by-treatment level (fixed /
random / total):

| Link    | Fixed | Random | Total |
|---------|------:|-------:|------:|
| logit   |  4.66 |  71.34 | 75.99 |
| probit  |  5.03 |  72.84 | 77.87 |
| cloglog |  3.44 |  69.21 | 72.66 |

The values from
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
are within 0.3 percentage points of these. The small differences come
from the fit: the paper used adaptive Gauss-Hermite quadrature in SAS,
while [`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) with nested
random effects can only use the Laplace approximation, which gives
slightly different variance components.

### Data at patient level

The same data with one row per patient (1 = favourable, 0 = not). The
model and its estimates are the same, but the R² values are much
smaller, and they are meant to be: they answer a different question. At
the clinic-by-treatment level the unit is a group of patients, whose
proportion of favourable responses varies much less than a single
patient’s outcome. As Piepho (2023) puts it, the patient-level analysis
is relevant when considering the value of the treatment for an
individual patient, and the clinic-by-treatment level when the
treatments are compared for entire clinics. So state the level with any
R² for binomial data, and compare R² values only at the same level.

Piepho (2023, Table 3) gives, at the patient level:

| Link    | Fixed | Random | Total |
|---------|------:|-------:|------:|
| logit   |  2.23 |  32.76 | 34.99 |
| probit  |  2.52 |  34.92 | 37.44 |
| cloglog |  1.86 |  35.87 | 37.74 |

Code

``` r

Clinical_N <- Clinical %>%
  pivot_longer(c(Favorable, Unfavorable), names_to = "outcome", values_to = "n") %>%
  uncount(n) %>%
  mutate(favour = as.integer(outcome == "Favorable")) %>%
  select(Clinic, Treatment, favour)
str(Clinical_N)
```

    ## tibble [273 × 3] (S3: tbl_df/tbl/data.frame)
    ##  $ Clinic   : Factor w/ 8 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ Treatment: Factor w/ 2 levels "control","drug": 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ favour   : int [1:273] 1 1 1 1 1 1 1 1 1 1 ...

Code

``` r

fit_patient <- function(link) {
  glmer(favour ~ Treatment + (1 | Clinic/Treatment),
        family = binomial(link = link), data = Clinical_N)
}
R2_model(fit_patient("logit"))
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##            Total            Fixed           Random Treatment:Clinic 
    ##           34.57%            2.23%           32.34%            0.11% 
    ##           Clinic 
    ##           32.24% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Code

``` r

R2_model(fit_patient("probit"))
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##            Total            Fixed           Random Treatment:Clinic 
    ##           37.48%            2.50%           34.97%            0.03% 
    ##           Clinic 
    ##           34.95% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Code

``` r

R2_model(fit_patient("cloglog"))
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##            Total            Fixed           Random Treatment:Clinic 
    ##           37.65%            1.91%           35.74%            0.00% 
    ##           Clinic 
    ##           35.74% 
    ## 
    ## Advice:
    ## ! The fit is singular: a variance estimated as zero gives no share to its
    ##   random effect, so the random part of R2 may be understated. Consider
    ##   simplifying the random effects.
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

## Example 5: foot shape of lambs (binomial GLMM with probit link)

The `Footshape_binomial` data give deformities in the feet of lambs, the
progeny of 34 sires, scored in three ordered categories (`k1`, `k2`,
`k3`; `n` lambs in total per sire), with four contrasts (`yr`, `b1`,
`b2`, `b3`) as fixed effects; 2513 lambs in all (Gilmour, Anderson and
Rae, 1987). The ordinal response is analysed here as two binary splits:
category 1 against categories 2 and 3 (`m` = `k1` lambs), and categories
1 and 2 against category 3 (`m12` = `k1` + `k2` lambs). A probit link
matches the threshold model usual for such traits, and the random share
of the R² is then the share of the latent variance due to sires.
(Piepho, 2023, also fits the ordinal response directly with a cumulative
probit model;
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
does not handle such models.)

As in Example 4 there are two levels. The sire level is relevant in
breeding programmes, where sires are selected on the proportions of
their progeny in each category; the lamb level is relevant for the
status of individual lambs. Piepho (2023, Table 5) gives (fixed / random
/ total):

| Split           | Sire level            | Lamb level           |
|-----------------|-----------------------|----------------------|
| K1 vs K2 and K3 | 50.06 / 32.72 / 82.78 | 5.54 / 4.34 / 9.88   |
| K1 and K2 vs K3 | 59.54 / 9.06 / 68.60  | 20.19 / 3.09 / 23.28 |

The values from
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
below agree with these to within 0.2 percentage points (quadrature in
SAS against the Laplace approximation in
[`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)).

### Data at sire level

Code

``` r

load("data/Footshape_binomial.RData")
str(Footshape_binomial)
```

    ## tibble [34 × 11] (S3: tbl_df/tbl/data.frame)
    ##  $ yr  : num [1:34] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ b1  : num [1:34] 1 1 1 1 1 1 1 -1 -1 -1 ...
    ##  $ b2  : num [1:34] 0 0 0 0 0 0 0 1 1 1 ...
    ##  $ b3  : num [1:34] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ k1  : num [1:34] 52 49 50 42 74 54 96 57 55 70 ...
    ##  $ k2  : num [1:34] 25 17 13 9 15 8 12 52 27 36 ...
    ##  $ k3  : num [1:34] 0 1 1 0 0 0 0 9 5 4 ...
    ##  $ sire: Factor w/ 34 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ n   : num [1:34] 77 67 64 51 89 62 108 118 87 110 ...
    ##  $ m   : num [1:34] 52 49 50 42 74 54 96 57 55 70 ...
    ##  $ m12 : num [1:34] 77 66 63 51 89 62 108 109 82 106 ...

Code

``` r

mod_sire_k1 <- glmer(cbind(m, n - m) ~ yr + b1 + b2 + b3 + (1 | sire),
                     family = binomial(link = "probit"), data = Footshape_binomial)
summary(mod_sire_k1, corr = FALSE)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( probit )
    ## Formula: cbind(m, n - m) ~ yr + b1 + b2 + b3 + (1 | sire)
    ##    Data: Footshape_binomial
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     236.0     245.1    -112.0     224.0        28 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.06404 -0.26300 -0.03548  0.43476  1.25127 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  sire   (Intercept) 0.04969  0.2229  
    ## Number of obs: 34, groups:  sire, 34
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  0.38309    0.04971   7.706 1.30e-14 ***
    ## yr           0.11199    0.04969   2.254  0.02421 *  
    ## b1           0.37284    0.07354   5.070 3.98e-07 ***
    ## b2           0.30725    0.10153   3.026  0.00248 ** 
    ## b3          -0.05590    0.06686  -0.836  0.40308    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Code

``` r

R2_model(mod_sire_k1)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##  Total  Fixed Random   sire 
    ## 82.82% 50.20% 32.62% 32.62% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Code

``` r

mod_sire_k12 <- glmer(cbind(m12, n - m12) ~ yr + b1 + b2 + b3 + (1 | sire),
                      family = binomial(link = "probit"), data = Footshape_binomial)
summary(mod_sire_k12, corr = FALSE)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( probit )
    ## Formula: cbind(m12, n - m12) ~ yr + b1 + b2 + b3 + (1 | sire)
    ##    Data: Footshape_binomial
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     147.6     156.7     -67.8     135.6        28 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6547 -0.3982  0.2360  0.5093  1.3862 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  sire   (Intercept) 0.04142  0.2035  
    ## Number of obs: 34, groups:  sire, 34
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  1.75944    0.08081  21.772  < 2e-16 ***
    ## yr           0.34671    0.07909   4.384 1.17e-05 ***
    ## b1           0.60765    0.13947   4.357 1.32e-05 ***
    ## b2           0.36198    0.11472   3.155   0.0016 ** 
    ## b3          -0.22321    0.07381  -3.024   0.0025 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Code

``` r

R2_model(mod_sire_k12)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##  Total  Fixed Random   sire 
    ## 68.47% 59.49%  8.99%  8.99% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

### Data at lamb level

With one row per lamb the R² values are much smaller than those from the
sire totals, as in Example 4.

Code

``` r

lambs <- function(data, success) {
  data %>%
    mutate(yes = .data[[success]], no = n - .data[[success]]) %>%
    select(yr, b1, b2, b3, sire, yes, no) %>%
    pivot_longer(c(yes, no), names_to = "outcome", values_to = "count") %>%
    uncount(count) %>%
    mutate(y = as.integer(outcome == "yes"))
}
Footshape_k1 <- lambs(Footshape_binomial, "m")
Footshape_k12 <- lambs(Footshape_binomial, "m12")
str(Footshape_k1)
```

    ## tibble [2,513 × 7] (S3: tbl_df/tbl/data.frame)
    ##  $ yr     : num [1:2513] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ b1     : num [1:2513] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ b2     : num [1:2513] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ b3     : num [1:2513] 0 0 0 0 0 0 0 0 0 0 ...
    ##  $ sire   : Factor w/ 34 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ outcome: chr [1:2513] "yes" "yes" "yes" "yes" ...
    ##  $ y      : int [1:2513] 1 1 1 1 1 1 1 1 1 1 ...

Code

``` r

mod_lamb_k1 <- glmer(y ~ yr + b1 + b2 + b3 + (1 | sire),
                     family = binomial(link = "probit"), data = Footshape_k1)
R2_model(mod_lamb_k1)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##  Total  Fixed Random   sire 
    ##  9.90%  5.57%  4.33%  4.33% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

Code

``` r

mod_lamb_k12 <- glmer(y ~ yr + b1 + b2 + b3 + (1 | sire),
                      family = binomial(link = "probit"), data = Footshape_k12)
R2_model(mod_lamb_k12)
```

    ## # Adjusted R2 for mixed models (Piepho)
    ##  Total  Fixed Random   sire 
    ## 23.32% 20.25%  3.07%  3.07% 
    ## 
    ## Advice:
    ## i For a non-gaussian model the R2 is on the link (latent) scale: not
    ##   comparable with an lm R2 or across families and links; use it to compare
    ##   models of the same family.

### Ordinal model at lamb level

The three categories can also be modelled together with a cumulative
probit (proportional odds) model, fitted here with
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html).
[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
does not handle such models, but the lamb-level R² is easy to compute by
hand from Piepho’s formulas: the latent residual variance is 1 (probit),
the fixed share is the bias-corrected variance of the linear predictor,
and the sire share is the sire variance times the factor for unequal
family sizes.

Code

``` r

Footshape_long <- Footshape_binomial %>%
  select(yr, b1, b2, b3, sire, k1, k2, k3) %>%
  pivot_longer(c(k1, k2, k3), names_to = "category", values_to = "count") %>%
  uncount(count) %>%
  mutate(category = factor(category, ordered = TRUE))

mod_ord <- ordinal::clmm(category ~ yr + b1 + b2 + b3 + (1 | sire), data = Footshape_long,
                         link = "probit", nAGQ = 10)

n <- nrow(Footshape_long)
X <- scale(as.matrix(Footshape_long[c("yr", "b1", "b2", "b3")]), scale = FALSE)
b <- mod_ord$beta
fixed <- (sum((X %*% b)^2) - sum(diag(crossprod(X) %*% vcov(mod_ord)[names(b), names(b)]))) / (n - 1)
sizes <- table(Footshape_long$sire)
sire <- as.numeric(ordinal::VarCorr(mod_ord)$sire) * (n - sum(sizes^2) / n) / (n - 1)
round(100 * c(Fixed = fixed, Random = sire, Total = fixed + sire) / (fixed + sire + 1), 2)
```

    ##  Fixed Random  Total 
    ##   6.75   4.18  10.93

Piepho (2023, Table 5) gives 6.75% (fixed), 4.18% (random) and 10.93%
(total) for this model.

## Summary

[`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
computes Piepho’s adjusted R² for `lmer`, `glmer` and `glmmTMB` fits
(and for models without random effects), choosing the residual variance
for the family and link automatically. It splits the explained variance
into the fixed effects and each random-effect term, and with
`over_disp = TRUE` handles an observation-level random effect for
overdispersion. On every example of Piepho (2023) it reproduces the
published values: exactly for the normal, Poisson and binary examples
fitted by the Laplace approximation, and to within 0.4 percentage points
where the paper used quadrature. R² values of non-gaussian models are on
the latent scale, and for binomial data they depend on the level of
aggregation (groups or individuals), which answers a different question
at each level: state the level, and compare R² values only between
models of the same family and link at the same level.

## References

Bates, D., Maechler, M., Bolker, B. and Walker, S. (2015). Fitting
linear mixed-effects models using lme4. *Journal of Statistical
Software*, 67(1), 1-48. <https://doi.org/10.18637/jss.v067.i01>

Beitler, P. J. and Landis, J. R. (1985). A mixed-effects model for
categorical data. *Biometrics*, 41, 991-1000.

Gilmour, A. R., Anderson, R. D. and Rae, A. L. (1987). Variance
components on an underlying scale for ordered multiple threshold
categorical data using a generalized linear mixed model. *Journal of
Animal Breeding and Genetics*, 104, 149-155.

Luo, D., Ganesh, S. and Koolaard, J. predictmeans: predicted means for
linear and semiparametric models. R package.
<https://CRAN.R-project.org/package=predictmeans>

Nakagawa, S. and Schielzeth, H. (2013). A general and simple method for
obtaining R² from generalized linear mixed-effects models. *Methods in
Ecology and Evolution*, 4, 133-142.

Piepho, H.-P. (2023). An adjusted coefficient of determination (R²) for
generalized linear mixed models in one go. *Biometrical Journal*, 65(7),
e2200290. <https://doi.org/10.1002/bimj.202200290>

R Core Team (2025). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>
