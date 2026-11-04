# Adjusted R² for generalized linear mixed models

[`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md)
computes an adjusted coefficient of determination (R²) for linear and
generalized linear mixed models with Piepho’s (2023) method. It works
for `lmer`, `glmer` and `glmmTMB` fits and picks the residual
(distribution-specific) variance for the model’s family and link
automatically.

``` r

library(predictmeans)
args(R2_glmm)
#> function (model, over_disp = FALSE) 
#> NULL
```

With `over_disp = TRUE`, an observation-level random effect (one level
per observation) is treated as overdispersion and counted as residual
variance.

## Example 1: body length (linear mixed model)

``` r

data(Body_length)
str(Body_length)
#> 'data.frame':    960 obs. of  6 variables:
#>  $ Population: Factor w/ 12 levels "1","10","11",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Container : Factor w/ 120 levels "1","10","100",..: 1 1 1 1 1 1 1 1 33 33 ...
#>  $ Sex       : Factor w/ 2 levels "Female","Male": 1 1 1 1 1 1 1 1 2 2 ...
#>  $ Habitat   : Factor w/ 2 levels "A","B": 1 1 1 1 2 2 2 2 1 1 ...
#>  $ Treatment : Factor w/ 2 levels "Cont","Exp": 1 2 1 2 1 2 1 2 1 2 ...
#>  $ BodyL     : num  10.8 11.2 11.6 15.6 12.4 ...
mod_BodyL <- lmer(BodyL ~ Treatment + Habitat + Sex + (1|Population/Container),
                  data = Body_length)
R2_glmm(mod_BodyL)
#> # Adjusted R2 for Mixed Models
#>                Total                Fixed               Random 
#>             "73.39%"             "40.09%"              "33.3%" 
#> Container:Population           Population 
#>              "5.18%"             "28.13%"
```

## Example 2: egg counts (Poisson GLMM with overdispersion)

Each female has her own level of `individual` within `Container` and
`Population`, so the innermost random effect is observation-level and
models overdispersion.

``` r

data(Poisson_ind)
mod_poisson <- glmer(Egg ~ Treatment + Habitat + (1|Population/Container/individual),
                     family = "poisson", data = Poisson_ind)
R2_glmm(mod_poisson, over_disp = TRUE)
#> # Adjusted R2 for Mixed Models
#>    Total    Fixed   Random 
#> "50.93%"  "9.13%" "41.79%"
```

## Example 3: clinical trial (binomial GLMM, three links)

The `Clinical` data (Beitler & Landis 1985) give the number of
favourable responses out of the patients in each treatment group at each
of eight clinics.

### Clinic-level data

``` r

data(Clinical)
Clinical$Unfavorable <- Clinical$Total - Clinical$Favorable
fit_grp <- function(link)
  glmer(cbind(Favorable, Unfavorable) ~ Treatment + (1|Clinic/Treatment),
        family = binomial(link = link), data = Clinical)
R2_glmm(fit_grp("logit"))
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random Treatment:Clinic 
#>         "75.73%"          "4.69%"         "71.04%"          "0.23%" 
#>           Clinic 
#>         "70.81%"
R2_glmm(fit_grp("probit"))
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random Treatment:Clinic 
#>          "77.9%"             "5%"          "72.9%"          "0.05%" 
#>           Clinic 
#>         "72.84%"
R2_glmm(fit_grp("cloglog"))
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random Treatment:Clinic 
#>         "72.63%"          "3.54%"         "69.09%"             "0%" 
#>           Clinic 
#>         "69.09%"
```

### Patient-level (binary) data

The same data with one row per patient give the same fixed effects, but
the residual variance is now that of a Bernoulli outcome.

``` r

idx <- rep(seq_len(nrow(Clinical)), Clinical$Total)
Clinical_N <- Clinical[idx, c("Clinic", "Treatment")]
Clinical_N$favour <- unlist(lapply(seq_len(nrow(Clinical)), function(i)
  rep(c(1, 0), c(Clinical$Favorable[i], Clinical$Unfavorable[i]))))
str(Clinical_N)
#> 'data.frame':    273 obs. of  3 variables:
#>  $ Clinic   : Factor w/ 8 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Treatment: Factor w/ 2 levels "control","drug": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ favour   : num  1 1 1 1 1 1 1 1 1 1 ...
fit_pat <- function(link)
  glmer(favour ~ Treatment + (1|Clinic/Treatment),
        family = binomial(link = link), data = Clinical_N)
R2_glmm(fit_pat("logit"))
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random Treatment:Clinic 
#>         "34.57%"          "2.23%"         "32.34%"          "0.11%" 
#>           Clinic 
#>         "32.24%"
R2_glmm(fit_pat("probit"))
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random Treatment:Clinic 
#>         "37.48%"           "2.5%"         "34.97%"          "0.03%" 
#>           Clinic 
#>         "34.95%"
R2_glmm(fit_pat("cloglog"))
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random Treatment:Clinic 
#>         "37.65%"          "1.91%"         "35.74%"             "0%" 
#>           Clinic 
#>         "35.74%"
```

## References

Piepho, H.-P. (2023). An adjusted coefficient of determination (R²) for
generalized linear mixed models in one go. *Biometrical Journal*, 65(7),
e2200290. <doi:10.1002/bimj.202200290>

Nakagawa, S. and Schielzeth, H. (2013). A general and simple method for
obtaining R² from generalized linear mixed-effects models. *Methods in
Ecology and Evolution*, 4, 133–142.

Beitler, P. J. and Landis, J. R. (1985). A mixed-effects model for
categorical data. *Biometrics*, 41, 991–1000.
