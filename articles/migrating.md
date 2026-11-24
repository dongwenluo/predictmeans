# Migrating from predictmeans 1.x

predictmeans 2.0 removes options that did the same thing under two
names, the printing and plot settings in the calls of the analysis
functions, five functions, and the 1.x smooth interface
([`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md),
`smoothZ`). The means functions share one vocabulary: `conf_level`,
`df_method`, `at`, `weights`, `vcov`, `plot`.

**To update a 1.1.1 script, run it.** A removed argument or function
stops with its 2.0 spelling, the value you gave filled in:

``` r

library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
try(predictmeans(fm, "nitro", level = 0.1, Df = 30))
#> Error : predictmeans(): 'level' was removed in predictmeans 2.0: use conf_level = 0.9.
#>   'Df' was removed in predictmeans 2.0: use df_method = 30.
#> See vignette("migrating") for every change from 1.1.1.
```

Nothing is silently ignored. A rewritten call gives the same numbers as
1.1.1, apart from the changes listed under [Defaults that
change](#defaults-that-change) and [Results that
change](#results-that-change). The package sources check this for every
rewrite (`dev/snapshot-*.R`, `dev/check-smz-migration.R`).

## Code side by side

The models: `fm` above; `fm0`, the same model without the interaction;
`fmc`, with `nitro` as a covariate; `fs`, a smooth model; `fit` stands
for any model. Every 2.0 line below is run by `dev/check-migrating.R` in
the package sources, and every 1.1.1 line checked to stop with the 2.0
spelling.

``` r

fm0 <- lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Oats)
fmc <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = nlme::Oats)
fs <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
              data = nlme::Dialyzer)
```

### Means and comparisons

1.1.1

``` r

predictmeans(fm, "nitro:Variety",
  pairwise = TRUE, level = 0.1,
  plotord = c(2, 1), plottitle = "Yield",
  prtnum = FALSE)


predictmeans(fm, "nitro", Df = 30,
  mplot = FALSE, barplot = TRUE,
  pplot = FALSE, bkplot = FALSE)

predictmeans(fit, "A", trans = exp,
  transOff = 1)

predictmeans(fit_pois, "A", count = TRUE)

contrastmeans(fm, "nitro", ctrmatrix = K,
  Df = 30)

df_term(fm, "nitro",
  type = "Satterthwaite")
```

2.0

``` r

pm <- predictmeans(fm, "nitro:Variety",
  pairwise = TRUE, conf_level = 0.9)
autoplot(pm, x = "Variety", group = "nitro") +
  labs(title = "Yield")
# plot(pm, ...) prints it

pm <- predictmeans(fm, "nitro", df_method = 30)
autoplot(pm, type = "bar")


predictmeans(fit, "A",
  trans = function(x) exp(x) - 1)

predictmeans(fit_pois, "A")

contrastmeans(fm, "nitro", ctrmatrix = K,
  df_method = 30)

df_term(fm, "nitro",
  df_method = "Satterthwaite")
```

- Nothing is printed and no plot is drawn during the call: `print(pm)`
  (or `pm` at the prompt) shows the tables, `autoplot(pm)` (or
  `plot(pm)`) draws.
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  has no `plot` or `plot_args`: the 1.x plots are `autoplot(pm)` (means
  with one LSD bar), `type = "bar"`, `type = "ci"` (the back-transformed
  plot, with `data = TRUE` the observations) and `type = "pvalue"`; `x`,
  `group`, `panel` replace `plotord`, `lines = FALSE`
  `lineplot = FALSE`, and `+ labs()` / `+ theme_bw(14)` the titles and
  `basesz`.
- `trans` is needed only for a transformation outside the formula:
  `log(y + 1)` in the formula is detected and back-transformed, constant
  included. Poisson and negative binomial means need no `count = TRUE`.

### Covariates and curves

1.1.1

``` r

covariatemeans(fmc, "Variety",
  covariate = "nitro",
  covariateV = seq(0, 0.6, 0.1),
  level = 0.1, trellis = FALSE,
  mtitle = "Yield by nitrogen")


semipred(semireg = fs, modelterm = "QB",
  covariate = "pressure",
  covariateV = c(0.5, 1, 1.5),
  trellis = FALSE, prt = FALSE)

Kmatrix(fmc, "Variety", covariate = "nitro",
  covariateV = 0.3, prtnum = TRUE)
```

2.0

``` r

cm <- covariatemeans(fmc, "Variety",
  covariate = "nitro",
  at = seq(0, 0.6, 0.1),
  conf_level = 0.9)
autoplot(cm, trellis = FALSE,
  mtitle = "Yield by nitrogen")

sp <- semipred(fs, "QB", "pressure",
  at = c(0.5, 1, 1.5))

autoplot(sp, trellis = FALSE)

Kmatrix(fmc, "Variety", covariate = "nitro",
  at = 0.3)
```

### Permutation tests

1.1.1

``` r

permanova.lmer(fm, nperm = 999, type = 3)

permlmer(lmer0 = fm0, lmer1 = fm)
permglmer(glmer0, glmer1)
permglmmTMB(tmb0, tmb1)

permmodels(fit, fo = f, prt = FALSE)

permstats(..., level = 0.05)
```

2.0

``` r

permmodels(fm, nperm = 999, type = 3)

permlmer(fm0, fm)
permlmer(glmer0, glmer1)
permlmer(tmb0, tmb1)

res <- permmodels(fit); res

permstats(..., conf_level = 0.95)
```

[`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
uses Freedman–Lane permutations, so its p-values are not those of
[`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)’s
refits ([Results that change](#results-that-change)).

### Diagnostics, R² and letters

1.1.1

``` r

residplot(fm, level = 2, id = TRUE,
  newwd = FALSE)
CookD(fm, idn = 5, newwd = FALSE)

PMplot(pm$p_valueMatrix, level = 0.01,
  legendx = 0.7, newwd = FALSE)

varcomp(fm, level = 0.9)

R2_glmm(fm)

# LL, UL: limits of the means
ci_mcp(LL, UL)
```

2.0

``` r

residplot(fm, re_level = 2, id = 3)

CookD(fm, id = 5)

PMplot(pm$p_valueMatrix,
  conf_level = 0.99)

varcomp(fm, conf_level = 0.9)

R2_model(fm)

predictmeans(fm, "nitro", pairwise = TRUE,
  letter_ci = TRUE)$mean_table
```

- `varcomp(level = 0.9)` was a confidence level already: the value
  stays.
- [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md)
  gives the values of
  [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  as numbers, not as text with a `%` sign.
- [`ci_mcp()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  gave letters from the overlap of the intervals it was given. Overlap
  of 95% (or 83.4%) confidence intervals does not answer the pairwise
  tests; `letter_ci = TRUE` gives letters from comparison intervals,
  which overlap exactly when a pair is not significant
  ([`vignette("predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/predictmeans.md)).

## Every argument and function

| Function | 1.1.1 | 2.0 |
|----|----|----|
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md), [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md), [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md), [`permstats()`](https://dongwenluo.github.io/predictmeans/reference/permstats.md), [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md) | `level = 0.05` (a significance level) | `conf_level = 0.95` |
| [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md) | `level = 0.9` (a confidence level) | `conf_level = 0.9` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) | `Df = 10` | `df_method = 10` |
| [`df_term()`](https://dongwenluo.github.io/predictmeans/reference/df_term.md) | `type =` | `df_method =` |
| [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md), [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md), [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md) | `covariateV =` | `at =` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) | `trans = exp, transOff = 1` | `trans = function(x) exp(x) - 1` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) | `responsen =` | removed: the response is read from the model |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md) | `count = TRUE` | removed: Poisson and negative binomial means are detected |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md) | `mplot`, `barplot`, `pplot`, `bkplot`, and `plot` of earlier 2.0 versions | `autoplot(pm, type = "means")`, `"bar"`, `"pvalue"`, `"ci"`; `plot = FALSE` is still accepted |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md) | `lsd_bar`, `plotord`, `lineplot`, `plottitle`, `plotxlab`, `plotylab`, `jitterv`, `basesz`, and `plot_args` | `autoplot(pm, interval = , x = , group = , panel = , lines = , data = TRUE, jitter = ) + labs(title, x, y) + theme_bw(size)` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md) | `prtnum`, `ndecimal`, `prtplt` | removed: `print(pm, digits = )`; `autoplot(pm)` |
| [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) | `trellis`, `plotord`, `mtitle`, `ci`, `point`, `jitterv` | `autoplot(cm, panel = , group = , title = , interval = , data = , jitter = )`; `plot_args = list(...)` still works (superseded) |
| [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) | `cm$pltdf` | `cm$pred_df` (as [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)); `cm$pltdf` still gives it, with a message once per session |
| [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md) | `semipred(semireg = fit)` | `semipred(fit)` (`model =`, as in every function) |
| [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md) | `trellis`, `scales`, `plotord`, `ci`, `point`, `jitterv`, `threeD` | `autoplot(sp, panel = , scales = , group = , interval = , data = , jitter = , type = "3d")`; `plot_args = list(...)` still works (superseded); `threeD` gives one plotly widget, not a list |
| [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md), [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md), [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md), [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md), [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md) | `prt` | removed: `print(x)` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md), [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md), [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md), [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md) | `newwd` | removed: ggplot2 draws in the current device |
| [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md), [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md) | `smoothZ = list(sm = smZ(x))` | smooths in the formula: `s(x, bs = "os")` ([below](#smooth-models-smz-and-smoothz)) |
| [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md) | `ncenter`, `nscale`, `resp_scale` | removed: variables are used as they are |
| [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md) | `ngbinomial = TRUE` | `semireg_tmb(family = glmmTMB::nbinom2)` |
| [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md), [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md) | `predict_info` | removed: always kept |
| [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md) | `lmer0 =`, `lmer1 =` | `model0 =`, `model1 =` (any class: lmer, glmer, glmmTMB) |
| [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) | `fo = f` (the formula in a variable) | removed: read from the model |
| [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md) | `level = 2` (the random-effect term); `newwd`, `ask` | `re_level = 2`; removed |
| [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md) | `idn = 5` | `id = 5`, as [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md) |
| [`PMplot()`](https://dongwenluo.github.io/predictmeans/reference/PMplot.md) | `margin`, `legendx` | removed: ggplot2 sizes the margins and the legend |
| every function | dotted argument names: `conf.level`, `letterCI`, `as.is`, `test.statistic`, `max.degree`, `drop.unused.levels` | `conf_level`, `letter_ci`, `as_is`, `test_statistic`, `max_degree`, `drop_unused_levels` (one style: snake_case) |
| [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) |  | removed: smooths in the formula |
| [`permglmer()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md), [`permglmmTMB()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) |  | removed: [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md) (they were the same function) |
| [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) |  | removed: `permmodels(model, type = )` |
| [`R2_glmm()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) |  | removed: `R2_model(model)` |
| [`ci_mcp()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md) |  | removed: `predictmeans(..., letter_ci = TRUE)` |

Every other argument of 1.1.1 is kept, with its meaning.

## Smooth models: `smZ()` and `smoothZ`

Smooths are written in the formula. Their fixed terms (the covariate,
and with `by = f` the factor and its slopes) are added automatically, so
they need not be written, and cannot be forgotten. Each form below was
fitted both ways on the last build that had
[`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
(`dev/check-smz-migration.R` in the package sources): the REML
log-likelihoods are equal and the fitted values and curves agree to 1e-5
or better, unless a note says otherwise.

1.1.1

``` r

semireg(y ~ x, data = d,
  smoothZ = list(sm = smZ(x)))
```

2.0

``` r

semireg(y ~ s(x, bs = "os"), data = d)
```

1.1.1

``` r

semireg(y ~ x, data = d,
  smoothZ = list(sm = smZ(x, k = 10)))
```

2.0

``` r

semireg(y ~ s(x, bs = "os", k = 10),
  data = d)
```

1.1.1

``` r

semireg(y ~ x*f, data = d,
  smoothZ = list(sm = smZ(x, by = f,
                          group = TRUE)))
```

2.0

``` r

semireg(y ~ s(x, by = f, bs = "os"),
  data = d)
```

1.1.1

``` r

semireg(y ~ x*f, data = d,
  smoothZ = list(sm = smZ(x, by = f)))
```

2.0

``` r

# one variance for all levels
semireg(y ~ s(x, by = f, bs = "os", id = 1),
  data = d)
```

1.1.1

``` r

semireg(y ~ x + z, data = d,
  smoothZ = list(sx = smZ(x), sz = smZ(z)))
```

2.0

``` r

semireg(y ~ s(x, bs = "os") +
  s(z, bs = "os"), data = d)
```

1.1.1

``` r

semireg(y ~ x + (1 | id), data = d,
  smoothZ = list(sm = smZ(x)))
```

2.0

``` r

semireg(y ~ s(x, bs = "os") + (1 | id),
  data = d)
```

1.1.1

``` r

semireg(y ~ x + z, data = d,
  smoothZ = list(sm = smZ(cbind(x, z),
                          type = "Ztps")))
```

2.0

``` r

semireg(y ~ s(x, z, bs = "ztps", k = 6),
  data = d)
```

1.1.1

``` r

semireg(cnt ~ x*f, family = poisson,
  data = d,
  smoothZ = list(sm = smZ(x, by = f,
                          group = TRUE)))
```

2.0

``` r

semireg(cnt ~ s(x, by = f, bs = "os"),
  family = poisson, data = d)
```

1.1.1

``` r

semireg_tmb(cnt ~ x*f, family = poisson,
  data = d,
  smoothZ = list(sm = smZ(x, by = f,
                          group = TRUE)))
```

2.0

``` r

semireg_tmb(cnt ~ s(x, by = f, bs = "os"),
  family = poisson, data = d)
```

1.1.1

``` r

semireg(y ~ x, data = d,
  smoothZ = list(sm = smZ(x, intKnots = kn,
                          range.x = r)))
```

2.0

``` r

semireg(y ~ s(x, bs = "os",
  xt = list(knots = kn, range = r)),
  data = d)
```

Notes:

- `bs = "os"` is the O’Sullivan spline of
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  (`k` interior knots, 6 by default as
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md));
  `bs = "ztps"` is its thin-plate basis (`k = 6` knots per axis as
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md);
  the default of `bs = "ztps"` is 10).
- Knots of your own (the last pair above): here the fits differ, and 2.0
  is right: 1.x centred `x` before
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  saw it (`ncenter = TRUE`, the default), so knots given in the units of
  the data were in the wrong place (in the check: REML log-likelihood
  -398.53 against -389.12; 1.x with `ncenter = FALSE` gives -389.12).
- The Poisson fit through lme4 agrees to the optimizer’s tolerance
  (log-likelihood within 1e-6, fitted counts within 0.003); through
  glmmTMB to 1e-6.
- `smZ(x, drv = 1)` (the derivative of the basis) is
  `semipred(fit, ..., deriv = 1)`.
- Any other basis: `s(x, bs = "user", xt = list(basis = f))`, with
  `f(x)` returning the basis matrix; its columns must not be in the span
  of the constant and `x` (the fixed part). The other
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  types (`"bs"`, `"ns"`, the splines2 types, `"smspline"`) gave their
  B-spline or natural-spline columns iid random coefficients, which
  shrinks the curve towards zero rather than towards a straight line;
  they have no formula form that gives the same fit. `bs = "os"` is the
  penalised B-spline that does the same job properly.
- `fit$knots_lst`, `fit$range_lst`, `fit$type_lst` and
  `fit$smoothZ_call` are no longer in the result; a fit made with
  [`smZ()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  and saved before 2.0 cannot be used by
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  (it says so): refit it.

### The models of the overview article

The smooth models of the overview article
(`vignettes/articles/overview.Rmd`), in its 1.1.1 version and in 2.0:

1.1.1

``` r

mod1 <- semireg(
  apartment.price ~ construction.date,
  smoothZ = list(grp = smZ(
    construction.date, k = 25)),
  data = WarsawApts)
```

2.0

``` r

mod1 <- semireg(
  apartment.price ~
    s(construction.date, bs = "os", k = 25),
  data = WarsawApts)
```

1.1.1

``` r

mod2 <- semireg(
  apartment.price ~
    construction.date*district,
  smoothZ = list(grp = smZ(
    construction.date, k = 10,
    by = district, group = TRUE)),
  data = WarsawApts)
```

2.0

``` r

mod2 <- semireg(
  apartment.price ~
    s(construction.date, by = district,
      bs = "os", k = 10),
  data = WarsawApts)
```

1.1.1

``` r

mod <- semireg_tmb(
  mathscr ~ calwpct + log.avginc +
    compstu + (1|county) + (1|grspan),
  smoothZ = list(
    cal_grp = smZ(calwpct, k = 12),
    avg_grp = smZ(log.avginc, k = 5),
    comp_grp = smZ(compstu, k = 12)),
  data = Caschool)
```

2.0

``` r

mod <- semireg_tmb(
  mathscr ~ s(calwpct, bs = "os", k = 12) +
    s(log.avginc, bs = "os", k = 5) +
    s(compstu, bs = "os", k = 12) +
    (1 | county) + (1 | grspan),
  data = Caschool)
```

1.1.1

``` r

mod1 <- semireg_tmb(
  respirInfec ~ age + vitAdefic + sex +
    height + stunted + visit2 + visit3 +
    visit4 + visit5 + visit6 + (1|idnum),
  smoothZ = list(grp = smZ(age, k = 25)),
  family = binomial, data = indonRespir)
```

2.0

``` r

mod1 <- semireg_tmb(
  respirInfec ~ s(age, bs = "os", k = 25) +
    vitAdefic + sex + height + stunted +
    visit2 + visit3 + visit4 + visit5 +
    visit6 + (1 | idnum),
  family = binomial, data = indonRespir)
```

1.1.1

``` r

mod2 <- semireg_tmb(
  respirInfec ~ age + vitAdefic + sex +
    height + stunted + visit2 + visit3 +
    visit4 + visit5 + visit6 + (1|idnum),
  smoothZ = list(grp = smZ(age, k = 25,
    by = sex, group = TRUE)),
  family = binomial, data = indonRespir)
```

2.0

``` r

mod2 <- semireg_tmb(
  respirInfec ~
    s(age, by = sex, bs = "os", k = 25) +
    vitAdefic + height + stunted +
    visit2 + visit3 + visit4 + visit5 +
    visit6 + (1 | idnum),
  family = binomial, data = indonRespir)
```

1.1.1

``` r

mod <- semireg_tmb(
  pollenCount ~ rain + year +
    temperatureResidual + dayInSeason +
    windSpeed,
  smoothZ = list(
    s_day = smZ(dayInSeason, by = year,
                group = TRUE),
    s_temp = smZ(temperatureResidual),
    s_wind = smZ(windSpeed)),
  family = nbinom1, data = ragweed)
```

2.0

``` r

mod <- semireg_tmb(
  pollenCount ~ rain +
    s(dayInSeason, by = year, bs = "os") +
    s(temperatureResidual, bs = "os") +
    s(windSpeed, bs = "os"),
  family = nbinom1, data = ragweed)
```

1.1.1

``` r

mod <- semireg(
  log(grade) ~ (xcoord + ycoord + zcoord)^2,
  smoothZ = list(
    grp_z = smZ(zcoord),
    grp_geo = smZ(cbind(xcoord, ycoord),
                  type = "Ztps")),
  data = copper)
```

2.0

``` r

mod <- semireg(
  log(grade) ~ s(zcoord, bs = "os") +
    s(xcoord, ycoord, bs = "ztps", k = 6) +
    xcoord:ycoord + xcoord:zcoord +
    ycoord:zcoord,
  data = copper)
```

1.1.1

``` r

pred <- semipred(mod,
  covariate = c("xcoord", "ycoord"),
  threeD = TRUE, point = FALSE,
  trans = exp, prt = FALSE)
```

2.0

``` r

pred <- semipred(mod,
  covariate = c("xcoord", "ycoord"),
  trans = exp,
  plot_args = list(threeD = TRUE,
                   point = FALSE))
```

| Model | REML / ML log-likelihood, 1.1.1 code | 2.0 code |  |
|----|----|----|----|
| Warsaw, one curve | -1772.452 | -1772.452 | the same fit |
| Warsaw, a curve per district | -1751.525 | -1751.525 | the same fit |
| California schools | -1601.803 | -1601.803 | the same fit |
| Respiratory infection, one curve | -323.509 | -323.509 | the same fit |
| Respiratory infection, a curve per sex | -322.806 | -322.788 | 2.0 adds `sex:age`, the slope of each sex, which the 1.1.1 formula left out |
| Pollen, a curve per year | -1168.703 | -1164.093 | 2.0 adds `year:dayInSeason`, the slope of each year, which the 1.1.1 formula left out |
| Copper | -159.428 | -159.203 | 1.1.1 depended on the centring of the coordinates (-163.948 with `ncenter = FALSE`) |

In the two models with a curve per group, the 1.1.1 formulas had the
factor but not its interaction with the covariate. A smooth by a factor
then has no unpenalised slope per group, so the groups’ straight-line
trends are shrunk towards each other. The formula form adds the slopes
itself, which is why the 2.0 fits are better.

## Defaults that change

| Function | Setting | 1.1.1 | 2.0 |
|----|----|----|----|
| [`library(predictmeans)`](https://dongwenluo.github.io/predictmeans/) | packages attached | glmmTMB, lme4, nlme, lmerTest | lme4, nlme: add [`library(glmmTMB)`](https://github.com/glmmTMB/glmmTMB) to call `glmmTMB()` or `nbinom1`, and [`library(lmerTest)`](https://github.com/runehaubo/lmerTestR) for its [`summary()`](https://rdrr.io/r/base/summary.html) p-values and `ranova()` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md) | plots | drawn | none: `autoplot(pm)` or `plot(pm)` |
| [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md), [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md) | plots | drawn | `plot = FALSE`: `autoplot(x)`, `plot = TRUE`, or `options(predictmeans.plot = TRUE)` for a session |
| every function | printing | tables printed during the call | nothing printed: `print(x)` |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md) | standard errors of REML `lmer` fits | unadjusted, with Kenward–Roger df | Kenward–Roger adjusted (below); `vcov_adjust = FALSE` for 1.x |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`contrastmeans()`](https://dongwenluo.github.io/predictmeans/reference/contrastmeans.md), [`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md) | df | per model class | `df_method = "auto"` (table in [`vignette("predictmeans")`](https://dongwenluo.github.io/predictmeans/articles/predictmeans.md)) |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md) | `pairwise` | `FALSE`: no comparisons unless `TRUE` | `NULL`: comparisons are made when `adj`, `atvar`, `permlist`, `letter_ci` or `pairwise = TRUE` asks for them; `pairwise` only sets what is printed |
| [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md), [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md), [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md) | `ncore` | 3 | 1 |
| [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md) | `control` | `lmerControl(optimizer = "bobyqa")` | `lmerControl(optimizer = "nloptwrap")`: the same fit to the optimizer’s tolerance |
| [`semireg()`](https://dongwenluo.github.io/predictmeans/reference/semireg.md) | centring | `ncenter = TRUE` centred every covariate | variables used as they are |
| [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md) | `id` | `FALSE` | 3: the three largest residuals are labelled |

## Results that change

### Kenward–Roger adjusted standard errors

For `lmer` models fitted by REML, standard errors, confidence limits,
LSDs and p-values use the Kenward–Roger adjusted covariance, as SAS
`DDFM=KR` and emmeans do; 1.x used the unadjusted covariance with
Kenward–Roger df. For balanced designs the two agree; for unbalanced
ones the adjusted SEs are slightly larger. `vcov_adjust = FALSE` gives
the 1.x numbers:

``` r

Om <- Oats[-c(1, 2, 30), ]   # unbalanced
fm_u <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Om)
new <- predictmeans(fm_u, "Variety")
old <- predictmeans(fm_u, "Variety", vcov_adjust = FALSE)
cbind(SE_2.0 = new$means$std.error, SE_1.x = old$means$std.error)
#>        SE_2.0   SE_1.x
#> [1,] 8.086340 8.085432
#> [2,] 8.057580 8.057580
#> [3,] 8.138192 8.132226
```

### Means that cannot be estimated are NA

A mean that depends on aliased coefficients used to be a number that
depended on which coefficient R dropped. It is now `NA`, with a message;
estimable pairwise differences are still shown. Means over empty cells
keep the 1.x rule: the average of the observed cells, with a warning.

### Permutation tests

- [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  uses Freedman–Lane permutations (`method = "auto"`): each term,
  coefficient or pair is tested against the model without it, so the
  other terms keep their effects. It works for `glm`, `glmer` and
  `glmmTMB` too, without refits. `method = "raw"` is the 1.x scheme.
- Mixed-model residuals are decorrelated with the symmetric square root
  of their covariance; the Cholesky factor of 1.x gave results that
  depended on the order of the data. So `method = "raw"` p-values for
  mixed models also differ slightly from 1.x.
- [`permanova.lmer()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans-defunct.md)
  refitted each model with and without the term through
  [`permlmer()`](https://dongwenluo.github.io/predictmeans/reference/permlmer.md);
  [`permmodels()`](https://dongwenluo.github.io/predictmeans/reference/permmodels.md)
  tests the same hypotheses, with other permutation p-values.

### Other changes in numbers

- `trans = function(x) exp(x) - 1`: the pairs have no `ratio` columns.
  With an offset, 1.x gave ratios of `y + 1` labelled as ratios of `y`.
- `semireg(ngbinomial = TRUE)` took theta from a Poisson fit, which is
  not the maximum-likelihood fit (on a test data set log-likelihood
  -1133.60 against -1132.98 for
  `semireg_tmb(family = glmmTMB::nbinom2)`; theta 3.31 against 2.80,
  true value 3).
- Smooth models with knots given in the data’s units, or whose 1.x
  formula left out the slope of each group: see the [table
  above](#the-models-of-the-overview-article).

### Bugs fixed

Some 1.x results were wrong and are now right: `permlist` with `atvar`,
a user-supplied `Df` for `lmer` models, factor levels containing `":"`,
several
[`Kmatrix()`](https://dongwenluo.github.io/predictmeans/reference/Kmatrix.md)
cases with transformed covariates, and others listed in `NEWS.md`.

## The result of `predictmeans()`

The 1.x list elements are still there, so code that reads
`pm$"Predicted Means"` or `pm$"Pairwise p-value"` keeps working. The
tables `means`, `pairs` and `settings` hold the same numbers in one
layout for every function
([`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
`tidy()`, [`summary()`](https://rdrr.io/r/base/summary.html),
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)):

``` r

pm <- predictmeans(fm, "nitro", adj = "BH")
names(pm)
#>  [1] "Predicted Means"               "Standard Error of Means"      
#>  [3] "Standard Error of Differences" "LSD"                          
#>  [5] "Pairwise p-value"              "mean_table"                   
#>  [7] "p_valueMatrix"                 "means"                        
#>  [9] "pairs"                         "settings"
pm$`Predicted Means`
#> nitro
#>        0      0.2      0.4      0.6 
#>  79.3889  98.8889 114.2222 123.3889
pm$means
#>   nitro  estimate std.error       df  conf.low conf.high comparison.low
#> 1     0  79.38889  7.174754 6.791824  62.31730  96.46047       74.92186
#> 2   0.2  98.88889  7.174754 6.791824  81.81730 115.96047       94.42186
#> 3   0.4 114.22222  7.174754 6.791824  97.15064 131.29381      109.75519
#> 4   0.6 123.38889  7.174754 6.791824 106.31730 140.46047      118.92186
#>   comparison.high group  n estimable
#> 1        83.85592  A    18      TRUE
#> 2       103.35592   B   18      TRUE
#> 3       118.68925    C  18      TRUE
#> 4       127.85592     D 18      TRUE
```

## New in 2.0

- **Degrees of freedom**: `df_method = "Kenward-Roger"`,
  `"Satterthwaite"`, `"containment"`, `"between-within"`,
  `"inner-outer"`, or a number.
- **Averaging weights**: `weights = "proportional"` (SAS `OM`),
  `"outer"`, `"cells"`, `"flat"`.
- **Comparisons**: `adj = "sidak"`, `"scheffe"`, `"mvt"` (exact
  single-step);
  [`contrast_family()`](https://dongwenluo.github.io/predictmeans/reference/contrast_family.md),
  [`specify_contrast()`](https://dongwenluo.github.io/predictmeans/reference/specify_contrast.md),
  [`pm_joint()`](https://dongwenluo.github.io/predictmeans/reference/pm_joint.md),
  [`pm_slopes()`](https://dongwenluo.github.io/predictmeans/reference/pm_slopes.md);
  comparison intervals and `letter_ci`. See
  [`vignette("contrasts")`](https://dongwenluo.github.io/predictmeans/articles/contrasts.md).
- **Response scale**: back-transformation detected from the formula or
  the link, delta-method SEs, ratios and odds ratios.
- **Smooth models**: smooths in the formula (`s()`, `t2()`, `by =`),
  [`semipred()`](https://dongwenluo.github.io/predictmeans/reference/semipred.md)
  with simultaneous bands, derivatives and curve differences,
  [`anova()`](https://rdrr.io/r/stats/anova.html) of smooth terms.
- **Resampling**: Freedman–Lane permutations for GLMs and GLMMs,
  restricted permutations
  ([`permindex()`](https://dongwenluo.github.io/predictmeans/reference/permindex.md)),
  [`bootmodels()`](https://dongwenluo.github.io/predictmeans/reference/bootmodels.md).
  See
  [`vignette("resampling")`](https://dongwenluo.github.io/predictmeans/articles/resampling.md).
- **Checks**:
  [`pm_check()`](https://dongwenluo.github.io/predictmeans/reference/pm_check.md),
  [`residplot()`](https://dongwenluo.github.io/predictmeans/reference/residplot.md)
  and
  [`CookD()`](https://dongwenluo.github.io/predictmeans/reference/CookD.md)
  for every model class,
  [`varcomp()`](https://dongwenluo.github.io/predictmeans/reference/varcomp.md)
  with intervals and tests,
  [`R2_model()`](https://dongwenluo.github.io/predictmeans/reference/R2_model.md).
- **Other model classes** with a few methods
  ([`?pm_coef`](https://dongwenluo.github.io/predictmeans/reference/pm_coef.md));
  fewer dependencies; faster Kenward–Roger and Satterthwaite.
