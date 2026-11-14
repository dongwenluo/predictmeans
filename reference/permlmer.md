# Permutation Test of random or fixed effects for `lmer`, `glmer` and `glmmTMB` models.

This function provides permutation tests for the terms in a linear mixed
model of `lmer`, and in a generalized linear mixed model of `glmer` or
`glmmTMB` (also as `permglmer()` and `permglmmTMB()`).

## Usage

``` r
permlmer(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
         method=c("auto", "working", "quantile"))
permglmer(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
          method=c("auto", "working", "quantile"))
permglmmTMB(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
            method=c("auto", "working", "quantile"))
```

## Arguments

- lmer0:

  `lmer` (or `glmer`, `glmmTMB`) model under H0, note that `lmer0` model
  must nest within `lmer1` model.

- lmer1:

  `lmer` (or `glmer`, `glmmTMB`) model under H1, note that `lmer0` model
  must nest within `lmer1` model.

- nperm:

  Number of permutation, the default value is 999.

- ncore:

  Number of cores for parallel computing, the default value is 1 (before
  predictmeans 2.0: 3).

- plot:

  Plot permutation distribution or not, the default value is FALSE.

- seed:

  Specify a random number generator seed, for reproducible results. The
  user's random number stream is left as it was.

- perms:

  Optional matrix of permutations to use instead of free ones, one row
  per observation and one column per permutation (e.g. from
  [`permindex`](https://dongwenluo.github.io/predictmeans/reference/permindex.md));
  `nperm` is then its number of columns.

- method:

  For GLMMs only: `"working"` permutes on the working (link) scale,
  `"quantile"` permutes quantile residuals of a 0/1 response and refits
  the GLMMs; `"auto"` (the default) uses `"quantile"` for 0/1 responses
  and `"working"` otherwise. See Details.

## Value

The [`anova()`](https://rdrr.io/r/stats/anova.html) table of the two
models with the permutation p-value in column `Perm-p`; for GLMMs,
attribute `"method"` says which method was used, and attribute
`"failures"` gives the number of refits that failed and were left out
(if any).

## Details

The residuals of `lmer0` are decorrelated with the symmetric square root
of their covariance under `lmer0`, permuted, correlated again and added
to the fitted values of `lmer0` (Lee and Braun, 2012); both models are
refitted to each permuted data set. The likelihood-ratio statistics are
compared after rounding to 4 decimals, so an observed statistic of 0
(the tested variance estimated at 0) gives a p-value of 1. Since
predictmeans 2.0 the square root is symmetric, so the results do not
depend on the order of the data (1.x used the Cholesky factor).

For GLMMs (`glmer` or `glmmTMB`; since predictmeans 2.0),
`method = "working"` does the same on the working (link) scale: `lmer0`
is linearised at its fit (working response \\z = \eta + (y - \mu)
d\eta/d\mu\\, working variance \\\phi V(\mu)/(a (d\mu/d\eta)^2)\\; the
covariance of z also has the random effects, as in
[`R2_glmm`](https://dongwenluo.github.io/predictmeans/reference/R2_glmm.md)),
the decorrelated working residuals are permuted, and the statistic is
the likelihood-ratio statistic of weighted linear mixed models fitted to
the working responses (REML when the two models have the same fixed
effects, ML otherwise). In simulations it kept the 5% level for counts
and for binomial data with several trials, and was much more robust than
the usual likelihood-ratio test when counts were overdispersed. For 0/1
responses the working scale fails when events are rare or clusters tiny
(36% rejections at the 5% level), so `method = "quantile"` is used:
randomised quantile residuals (Dunn and Smyth, 1996) of the response
given the conditional means of `lmer0` are permuted and mapped back to
0/1 responses, and both models are refitted (`glmer` fits with
`nAGQ = 0`, for the observed data too, so the statistic can differ from
the `Chisq` column); it kept the level (5.6-6.2%). A random effect with
one level per observation (often used for overdispersion) cannot be used
on the working scale;
[`bootmodels`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md)
can test it. `glmmTMB` models must have no zero inflation and the
default dispersion model; families poisson, binomial, nbinom1, nbinom2,
Gamma, gaussian and beta are handled.

## References

Oliver E. Lee and Thomas M. Braun (2012), *Permutation Tests for Random
Effects in Linear Mixed Models. Biometrics*, Journal 68(2).

Peter K. Dunn and Gordon K. Smyth (1996), *Randomized Quantile
Residuals*. Journal of Computational and Graphical Statistics, 5(3),
236-244.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
## Test random effects (99 permutations to keep the example quick;
## use 999 or more in practice)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
anova(fm1, fm2, fm3)
#> refitting model(s) with ML (instead of REML)
#> Data: sleepstudy
#> Models:
#> fm3: Reaction ~ Days + (1 | Subject)
#> fm2: Reaction ~ Days + ((1 | Subject) + (0 + Days | Subject))
#> fm1: Reaction ~ Days + (Days | Subject)
#>     npar    AIC    BIC  logLik -2*log(L)   Chisq Df Pr(>Chisq)    
#> fm3    4 1802.1 1814.8 -897.04    1794.1                          
#> fm2    5 1762.0 1778.0 -876.00    1752.0 42.0754  1  8.782e-11 ***
#> fm1    6 1763.9 1783.1 -875.97    1751.9  0.0639  1     0.8004    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# \donttest{
permlmer(fm3, fm2, nperm=99, ncore=1)
#> Data: sleepstudy
#> Models:
#> lmer0: Reaction ~ Days + (1 | Subject)
#> lmer1: Reaction ~ Days + ((1 | Subject) + (0 + Days | Subject))
#>       npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> lmer0    4 1794.5 1807.2 -893.23    1786.5                            
#> lmer1    5 1753.7 1769.6 -871.83    1743.7 42.796  1 6.0763e-11   0.01
permlmer(fm2, fm1, nperm=99, ncore=1)
#> Data: sleepstudy
#> Models:
#> lmer0: Reaction ~ Days + ((1 | Subject) + (0 + Days | Subject))
#> lmer1: Reaction ~ Days + (Days | Subject)
#>       npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq) Perm-p
#> lmer0    5 1753.7 1769.6 -871.83    1743.7                           
#> lmer1    6 1755.6 1774.8 -871.81    1743.6 0.041  1     0.8395   0.91
# }

## Test fixed effects
Oats$nitro <- factor(Oats$nitro)
fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# \donttest{
permlmer(fm0, fm, nperm=99, ncore=1)
#> refitting model(s) with ML (instead of REML)
#> Data: Oats
#> Models:
#> lmer0: yield ~ nitro + Variety + (1 | Block/Variety)
#> lmer1: yield ~ nitro * Variety + (1 | Block/Variety)
#>       npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> lmer0    9 616.04 636.53 -299.02    598.04                            
#> lmer1   15 625.91 660.06 -297.95    595.91 2.1375  6    0.90663   0.96

## A GLMM: do tick counts vary between locations, beyond broods?
gm1 <- glmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD), family = poisson, data = lme4::grouseticks)
gm2 <- glmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION), family = poisson,
             data = lme4::grouseticks)
permglmer(gm1, gm2, nperm = 99, ncore = 1, seed = 1)
#> Data: lme4::grouseticks
#> Models:
#> m0: TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD)
#> m1: TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION)
#>    npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq) Perm-p
#> m0    5 1988.1 2008.1 -989.04    1978.1                            
#> m1    6 1987.9 2011.9 -987.94    1975.9 2.1992  1    0.13809   0.08
# }
```
