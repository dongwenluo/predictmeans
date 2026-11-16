# Predicted Means of a Semi Paramatric Model with Covariate Variable(s)

This function produces predicted means with graph for a semi paramatric
model with new set of covariate values.

## Usage

``` r
semipred(semireg, modelterm=NULL, covariate, sm_term=NULL, contr=NULL,
           covariateV=NULL, boundary=NULL, level=0.05, trans=NULL, trellis=TRUE, 
       scales=c("fixed", "free", "free_x", "free_y"),
       plotord=NULL, ci=TRUE, point=TRUE, jitterv=0, threeD=FALSE, prt=TRUE,
       simultaneous=FALSE, contr_scale=c("model", "response"),
       band=c("pointwise", "simultaneous"),
       unconditional=getOption("predictmeans.semi_unconditional"))
```

## Arguments

- semireg:

  A list object returned by `semireg`.

- modelterm:

  Name (in "quotes") for indicating which factor term's predicted mean
  to be calculated. The `modelterm` must be given exactly as it appears
  in `semireg` model, e.g. "A" or "A:B". In case `modelterm` is the same
  as `covariate` or `NULL`, them `semipred` will product predictmeans
  with CI based on `covariate` only with out any grouping. Factors that
  are not in `modelterm`, including `by` factors of the smooth terms,
  are averaged over with equal weights, in both the parametric and the
  smooth parts of the model.

- covariate:

  Name (in "quotes") of one or two (for `Ztps` smooth) the covariate
  variables in the `semireg` shuch as "x1" or c("x1", "x2").

- sm_term:

  Names (in "quotes") of smooth terms (from `smoothZ` list in `semireg`
  model) used in the prediction such as "sm1_grp" or c("sm1_grp",
  "sm2_grp"). The default is using all smooth terms which is
  `sm_term=NULL`.

- contr:

  A numeric (or character) vector with length of two (e.g. c(4, 1) or
  c("d", "a")) which indicates to produce predicted mean with CI for
  difference between `modelterm` level 4 ("d") vs level 1 ("a") along
  `covariate`.

- covariateV:

  A numeric vector (one covariate) or a matrix with one named column per
  covariate (two covariates, e.g. for a `Ztps` smooth), then semipred
  will produce the result for `covariate` at value of `covariateV`. Each
  row of the matrix is one prediction point.

- boundary:

  A matrix or data frame of two columns, used to specify boundary of
  longitude and latitude, it is functional when the length of covariate
  is two.

- level:

  Significance level of the two-sided confidence limits: the default
  0.05 gives 95% limits (before version 2.0 the limits used
  `qnorm(1 - level)`, i.e. 90%).

- trans:

  A function object for calculating the back transformed means, e.g.
  `trans=exp`.

- trellis:

  A logical scalar. If set to TRUE (default), a trellis plots of
  predicted means with CI will be drawn.

- scales:

  Should scales be fixed ("fixed", the default), free ("free"), or free
  in one dimension ("free_x", "free_y") in a trellis graph?

- plotord:

  A numeric vector specifying the order of plotting for two way
  interaction (e.g. `plotord = c(2, 1)` will put the second variable in
  `modelterm` on the `X` axis, the first variable as the grouping
  variable, and the third one as the panel variable). The defaults are
  `c(1, 2)` for two way interactions.

- ci:

  A logical scalar to indicate whether to print confidence interval. The
  default value is TRUE.

- point:

  A logical scalar to indicate whether to print raw data points. The
  default value is TRUE.

- jitterv:

  A degree of jitter in x and y direction in the graph. The default is
  zero.

- threeD:

  A logical scalar to indicate whether to produce a 3-D plot or not. The
  default value is FALSE. 3-D plots need package plotly.

- prt:

  A logical scalar to indicate whether to produce plots on the screen or
  not. The default value is TRUE.

- band:

  `"pointwise"` (the default) or `"simultaneous"`: simultaneous
  confidence bands over the whole curve(s) instead of pointwise limits;
  the critical value is the `1 - level` quantile of the largest
  standardised deviation over the predicted points, from 10000 draws of
  the estimates' distribution (Ruppert, Wand and Carroll, 2003). The
  same argument as in
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
  where lines of a parametric model use Scheffe's critical value
  instead; for smooths that would be far too wide.

- unconditional:

  If `TRUE`, the covariance of the coefficients includes the uncertainty
  of the smoothing parameters (the variances of the smooth terms), as
  mgcv's `unconditional = TRUE`: the first-order correction \\J
  V\_\theta J'\\ is added, \\J\\ the derivative of the estimated
  coefficients in the variance parameters and \\V\_\theta\\ the
  covariance of their estimates (from the Hessian of the deviance; for
  [`semireg_tmb()`](https://dongwenluo.github.io/predictmeans/reference/semireg_tmb.md)
  fits TMB's joint precision of all parameters). The default (`NULL`) is
  `TRUE` for generalized models and `FALSE` for gaussian ones: in a
  simulation of the difference between two curves
  (`dev/sim-curve-difference.R`, 150 data sets, nominal 0.95) the bands
  without it covered 0.93 pointwise and 0.89-0.91 over the whole curve
  for Poisson and binary data, and 0.955-0.957 and 0.94-0.95 with it;
  for gaussian data 0.97 without it (0.98 with it).
  `options(predictmeans.semi_unconditional = TRUE)` (or `FALSE`) sets it
  for `semipred()` and
  [`predictmeans()`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)
  on `semireg` fits.

- simultaneous:

  `TRUE` is the same as `band = "simultaneous"` (the name before 2.0).

- contr_scale:

  With `contr` and a generalized model: `"model"` (the default) plots
  the difference on the link scale, as the ratio (log link) or the odds
  ratio (logit link) of the two curves; `"response"` the difference of
  the back-transformed curves with a delta-method band. See ‘Comparing
  two lines’ in
  [`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md):
  the test is on the model scale.

## Value

- plt:

  A ggplot object.

- pred_df:

  A data.frame with predcted data.

## Details

With `contr` the difference between two curves is computed as in
[`covariatemeans`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md),
from the contrast of the two curves and the covariance of the fixed and
smooth coefficients. For a linear model `pred_df` has the difference (in
the column of the response), `LL`, `UL`, `SE_contr` and the indicators
`includes0` (the limits include zero) and `significant`. For a
generalized model it has the columns of `covariatemeans(contr = )`:
`diff`, `ses`, `LL`, `UL` and `includes0` on the link scale; the `ratio`
or `odds.ratio` with its limits and `..includes1`; `response.diff`,
`response.ses`, `response.LL`, `response.UL` and `response.includes0` on
the response scale; and `significant` (the link scale). Before 2.0 the
band of a generalized model's difference used the standard error of the
link-scale difference multiplied by a weight evaluated at the difference
of the means, which is not a standard error of any difference. The
ranges of the covariate where the curves differ are printed and kept in
`attr(pred_df, "significant")`; with `simultaneous = TRUE` they hold for
the whole curve, and the curves differ somewhere at level `level` when
there is such a range. (A joint F test, as
[`covariatemeans()`](https://dongwenluo.github.io/predictmeans/reference/covariatemeans.md)
gives, is not given: Wald tests of penalised smooth coefficients are not
reliable.)

## References

Ruppert, D., Wand, M. P. and Carroll, R. J. (2003), *Semiparametric
Regression*, Cambridge University Press.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
data(Dialyzer, package="nlme")
str(Dialyzer)
#> Classes ‘nfnGroupedData’, ‘nfGroupedData’, ‘groupedData’ and 'data.frame':   140 obs. of  5 variables:
#>  $ Subject : Ord.factor w/ 20 levels "10"<"8"<"2"<"6"<..: 9 9 9 9 9 9 9 3 3 3 ...
#>  $ QB      : Factor w/ 2 levels "200","300": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ pressure: num  0.24 0.505 0.995 1.485 2.02 ...
#>  $ rate    : num  0.645 20.115 38.46 44.985 51.765 ...
#>  $ index   : num  1 2 3 4 5 6 7 1 2 3 ...
#>  - attr(*, "formula")=Class 'formula'  language rate ~ pressure | Subject
#>   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  - attr(*, "labels")=List of 2
#>   ..$ x: chr "Transmembrane pressure"
#>   ..$ y: chr "Ultrafiltration rate"
#>  - attr(*, "units")=List of 2
#>   ..$ x: chr "(dmHg)"
#>   ..$ y: chr "(ml/hr)"
#>  - attr(*, "outer")=Class 'formula'  language ~QB
#>   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  - attr(*, "FUN")=function (x)  
#>   ..- attr(*, "source")= chr "function (x) max(x, na.rm = TRUE)"
#>  - attr(*, "order.groups")= logi TRUE

library(ggplot2)
ggplot(Dialyzer, aes(x=pressure, y=rate, col=QB)) +
  geom_line() +
  facet_wrap(vars(Subject))


fm <- semireg(rate ~ pressure*QB+(pressure|Subject),
              smoothZ=list(
                qb_grp=smZ(pressure, k=4, by=QB, group=TRUE)
              ),
              data=Dialyzer)
#> Linear mixed model fit by REML ['semerMod']
#> Formula: rate ~ pressure * QB + (pressure | Subject) + (1 | qb_grp.200) +  
#>     (1 | qb_grp.300)
#>    Data: Dialyzer
#> REML criterion at convergence: 711.2453
#> Random effects:
#>  Groups     Name        Std.Dev. Corr 
#>  Subject    (Intercept)  2.583        
#>             pressure     1.554   0.56 
#>  qb_grp.200 (Intercept) 27.493        
#>  qb_grp.300 (Intercept) 25.483        
#>  Residual                2.197        
#> Number of obs: 140, groups:  Subject, 20; qb_grp.200, 6; qb_grp.300, 6
#> Fixed Effects:
#>    (Intercept)        pressure           QB300  pressure:QB300  
#>         30.908          12.829           7.424           7.131  
str(fm$data)
#> 'data.frame':    140 obs. of  6 variables:
#>  $ rate      : num  0.645 20.115 38.46 44.985 51.765 ...
#>  $ pressure  : num  -1.2975 -1.0325 -0.5425 -0.0525 0.4825 ...
#>  $ QB        : Factor w/ 2 levels "200","300": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Subject   : Ord.factor w/ 20 levels "10"<"8"<"2"<"6"<..: 9 9 9 9 9 9 9 3 3 3 ...
#>  $ qb_grp.200: Factor w/ 6 levels "1","2","3","4",..: 1 2 3 4 5 6 1 2 3 4 ...
#>  $ qb_grp.300: Factor w/ 6 levels "1","2","3","4",..: 1 2 3 4 5 6 1 2 3 4 ...
#>  - attr(*, "pressure_mean")= Named num 1.54
#>   ..- attr(*, "names")= chr "pressure"
#>  - attr(*, "numeric_var")= chr "pressure"
summary(fm$semer)
#> Linear mixed model fit by REML ['semerMod']
#> Formula: rate ~ pressure * QB + (pressure | Subject) + (1 | qb_grp.200) +  
#>     (1 | qb_grp.300)
#>    Data: Dialyzer
#> 
#> REML criterion at convergence: 711.2
#> 
#> Scaled residuals: 
#>      Min       1Q   Median       3Q      Max 
#> -2.51549 -0.47841 -0.00357  0.46800  2.60038 
#> 
#> Random effects:
#>  Groups     Name        Variance Std.Dev. Corr 
#>  Subject    (Intercept)   6.670   2.583        
#>             pressure      2.414   1.554   0.56 
#>  qb_grp.200 (Intercept) 755.888  27.493        
#>  qb_grp.300 (Intercept) 649.389  25.483        
#>  Residual                 4.828   2.197        
#> Number of obs: 140, groups:  Subject, 20; qb_grp.200, 6; qb_grp.300, 6
#> 
#> Fixed effects:
#>                Estimate Std. Error t value
#> (Intercept)     30.9082     0.8654  35.717
#> pressure        12.8292     0.5721  22.425
#> QB300            7.4237     1.2252   6.059
#> pressure:QB300   7.1312     0.8117   8.786
#> 
#> Correlation of Fixed Effects:
#>             (Intr) pressr QB300 
#> pressure     0.462              
#> QB300       -0.706 -0.326       
#> prssr:QB300 -0.325 -0.705  0.455
residplot(fm$semer, group="QB")

anova(fm)   # fixed effects (Kenward-Roger) and a test of each smooth term
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: rate
#>                    F Df Df.res    Pr(>F)    
#> (Intercept) 3190.734  1 18.708 < 2.2e-16 ***
#> pressure    1628.834  1 19.142 < 2.2e-16 ***
#> QB            36.679  1 18.708 8.514e-06 ***
#> pressure:QB   77.043  1 19.142 3.873e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);
#> p-values halved for the boundary (chi-square 0:1)
#>        term variance logLik logLik.without   LRT p.value
#>  qb_grp.200    755.9 -355.6         -476.4 241.5  <2e-16
#>  qb_grp.300    649.4 -355.6         -468.3 225.3  <2e-16
R2_glmm(fm$semer)
#> # Adjusted R2 for Mixed Models
#>            Total            Fixed           Random          Subject 
#>         "98.54%"         "80.31%"         "18.23%"          "1.93%" 
#> Subject.pressure       qb_grp.200       qb_grp.300 
#>          "0.65%"          "8.47%"          "7.17%" 
ap_out1 <- semipred(fm, "QB", "pressure")

str(ap_out1$pred_df)
#> 'data.frame':    230 obs. of  6 variables:
#>  $ QB      : Factor w/ 2 levels "200","300": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ pressure: num  0.235 0.24 0.245 0.25 0.255 ...
#>  $ SE      : num  1.038 1.022 1.008 0.995 0.983 ...
#>  $ rate    : num  2.72 3.03 3.34 3.65 3.96 ...
#>  $ LL      : num  0.685 1.026 1.364 1.699 2.032 ...
#>  $ UL      : num  4.76 5.03 5.32 5.6 5.89 ...
ap_out2 <- semipred(fm, "QB", "pressure", contr=c(1,2))
#> 200 differs from 300 (pointwise band, level 0.05) for pressure: 1.01 to 3.03 (<).

str(ap_out2$pred_df)
#> 'data.frame':    115 obs. of  7 variables:
#>  $ pressure   : num  0.235 0.24 0.245 0.25 0.255 ...
#>  $ rate       : num  1.045 1.03 1.015 1 0.984 ...
#>  $ LL         : num  -1.95 -1.91 -1.88 -1.85 -1.83 ...
#>  $ UL         : num  4.04 3.97 3.91 3.85 3.8 ...
#>  $ SE_contr   : num  1.53 1.5 1.48 1.46 1.44 ...
#>  $ includes0  : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ significant: logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#>  - attr(*, "significant")='data.frame':  1 obs. of  3 variables:
#>   ..$ from     : num 1.01
#>   ..$ to       : num 3.03
#>   ..$ direction: chr "<"

# \donttest{
str(sleepstudy)
#> 'data.frame':    180 obs. of  3 variables:
#>  $ Reaction: num  250 259 251 321 357 ...
#>  $ Days    : num  0 1 2 3 4 5 6 7 8 9 ...
#>  $ Subject : Factor w/ 18 levels "308","309","310",..: 1 1 1 1 1 1 1 1 1 1 ...
ggplot(sleepstudy, aes(x=Days, y=Reaction)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~x, se=FALSE) +
  facet_wrap(vars(Subject), nrow=2) +
  labs(x="Days of sleep deprivation", y="Average reaction time (ms)")


sleep.semi <- semireg(Reaction ~ Days*Subject,
                      smoothZ=list(
                        sub_grp=smZ(Days, by=Subject, group=TRUE)
                      ),
                      data=sleepstudy)
#> Linear mixed model fit by REML ['semerMod']
#> Formula: Reaction ~ Days * Subject + (1 | sub_grp.308) + (1 | sub_grp.309) +  
#>     (1 | sub_grp.310) + (1 | sub_grp.330) + (1 | sub_grp.331) +  
#>     (1 | sub_grp.332) + (1 | sub_grp.333) + (1 | sub_grp.334) +  
#>     (1 | sub_grp.335) + (1 | sub_grp.337) + (1 | sub_grp.349) +  
#>     (1 | sub_grp.350) + (1 | sub_grp.351) + (1 | sub_grp.352) +  
#>     (1 | sub_grp.369) + (1 | sub_grp.370) + (1 | sub_grp.371) +  
#>     (1 | sub_grp.372)
#>    Data: sleepstudy
#> REML criterion at convergence: 1402.094
#> Random effects:
#>  Groups      Name        Std.Dev. 
#>  sub_grp.308 (Intercept) 8.770e+01
#>  sub_grp.309 (Intercept) 1.953e+00
#>  sub_grp.310 (Intercept) 0.000e+00
#>  sub_grp.330 (Intercept) 8.867e+00
#>  sub_grp.331 (Intercept) 8.600e+00
#>  sub_grp.332 (Intercept) 5.321e+01
#>  sub_grp.333 (Intercept) 0.000e+00
#>  sub_grp.334 (Intercept) 8.331e+00
#>  sub_grp.335 (Intercept) 0.000e+00
#>  sub_grp.337 (Intercept) 0.000e+00
#>  sub_grp.349 (Intercept) 5.746e+00
#>  sub_grp.350 (Intercept) 1.773e+01
#>  sub_grp.351 (Intercept) 7.290e+00
#>  sub_grp.352 (Intercept) 1.655e+01
#>  sub_grp.369 (Intercept) 6.125e-04
#>  sub_grp.370 (Intercept) 0.000e+00
#>  sub_grp.371 (Intercept) 1.266e+01
#>  sub_grp.372 (Intercept) 0.000e+00
#>  Residual                1.708e+01
#> Number of obs: 180, groups:  
#> sub_grp.308, 8; sub_grp.309, 8; sub_grp.310, 8; sub_grp.330, 8; sub_grp.331, 8; sub_grp.332, 8; sub_grp.333, 8; sub_grp.334, 8; sub_grp.335, 8; sub_grp.337, 8; sub_grp.349, 8; sub_grp.350, 8; sub_grp.351, 8; sub_grp.352, 8; sub_grp.369, 8; sub_grp.370, 8; sub_grp.371, 8; sub_grp.372, 8
#> Fixed Effects:
#>     (Intercept)             Days       Subject309       Subject310  
#>         343.817           22.748         -128.226         -112.815  
#>      Subject330       Subject331       Subject332       Subject333  
#>         -36.741          -33.231          -47.494          -27.658  
#>      Subject334       Subject335       Subject337       Subject349  
#>         -45.052          -93.747           31.904          -65.837  
#>      Subject350       Subject351       Subject352       Subject369  
#>         -28.624          -52.067          -10.449          -37.782  
#>      Subject370       Subject371       Subject372  Days:Subject309  
#>         -52.115          -44.565          -25.931          -20.490  
#> Days:Subject310  Days:Subject330  Days:Subject331  Days:Subject332  
#>         -16.633          -19.741          -17.181          -14.881  
#> Days:Subject333  Days:Subject334  Days:Subject335  Days:Subject337  
#>         -13.606          -10.598          -25.629           -3.722  
#> Days:Subject349  Days:Subject350  Days:Subject351  Days:Subject352  
#>          -9.247           -4.093          -16.119           -8.505  
#> Days:Subject369  Days:Subject370  Days:Subject371  Days:Subject372  
#>         -11.400           -4.692          -13.129          -11.450  
residplot(sleep.semi$semer)

summary(sleep.semi$semer)
#> Linear mixed model fit by REML ['semerMod']
#> Formula: Reaction ~ Days * Subject + (1 | sub_grp.308) + (1 | sub_grp.309) +  
#>     (1 | sub_grp.310) + (1 | sub_grp.330) + (1 | sub_grp.331) +  
#>     (1 | sub_grp.332) + (1 | sub_grp.333) + (1 | sub_grp.334) +  
#>     (1 | sub_grp.335) + (1 | sub_grp.337) + (1 | sub_grp.349) +  
#>     (1 | sub_grp.350) + (1 | sub_grp.351) + (1 | sub_grp.352) +  
#>     (1 | sub_grp.369) + (1 | sub_grp.370) + (1 | sub_grp.371) +  
#>     (1 | sub_grp.372)
#>    Data: sleepstudy
#> 
#> REML criterion at convergence: 1402.1
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -2.5307 -0.3655  0.0058  0.4992  3.3329 
#> 
#> Random effects:
#>  Groups      Name        Variance  Std.Dev. 
#>  sub_grp.308 (Intercept) 7.691e+03 8.770e+01
#>  sub_grp.309 (Intercept) 3.815e+00 1.953e+00
#>  sub_grp.310 (Intercept) 0.000e+00 0.000e+00
#>  sub_grp.330 (Intercept) 7.862e+01 8.867e+00
#>  sub_grp.331 (Intercept) 7.395e+01 8.600e+00
#>  sub_grp.332 (Intercept) 2.831e+03 5.321e+01
#>  sub_grp.333 (Intercept) 0.000e+00 0.000e+00
#>  sub_grp.334 (Intercept) 6.941e+01 8.331e+00
#>  sub_grp.335 (Intercept) 0.000e+00 0.000e+00
#>  sub_grp.337 (Intercept) 0.000e+00 0.000e+00
#>  sub_grp.349 (Intercept) 3.301e+01 5.746e+00
#>  sub_grp.350 (Intercept) 3.142e+02 1.773e+01
#>  sub_grp.351 (Intercept) 5.315e+01 7.290e+00
#>  sub_grp.352 (Intercept) 2.739e+02 1.655e+01
#>  sub_grp.369 (Intercept) 3.751e-07 6.125e-04
#>  sub_grp.370 (Intercept) 0.000e+00 0.000e+00
#>  sub_grp.371 (Intercept) 1.602e+02 1.266e+01
#>  sub_grp.372 (Intercept) 0.000e+00 0.000e+00
#>  Residual                2.918e+02 1.708e+01
#> Number of obs: 180, groups:  
#> sub_grp.308, 8; sub_grp.309, 8; sub_grp.310, 8; sub_grp.330, 8; sub_grp.331, 8; sub_grp.332, 8; sub_grp.333, 8; sub_grp.334, 8; sub_grp.335, 8; sub_grp.337, 8; sub_grp.349, 8; sub_grp.350, 8; sub_grp.351, 8; sub_grp.352, 8; sub_grp.369, 8; sub_grp.370, 8; sub_grp.371, 8; sub_grp.372, 8
#> 
#> Fixed effects:
#>                 Estimate Std. Error t value
#> (Intercept)      343.817      5.689  60.440
#> Days              22.748      1.958  11.617
#> Subject309      -128.226      7.871 -16.290
#> Subject310      -112.815      7.844 -14.382
#> Subject330       -36.741      7.958  -4.617
#> Subject331       -33.231      7.957  -4.177
#> Subject332       -47.494      8.012  -5.928
#> Subject333       -27.658      7.844  -3.526
#> Subject334       -45.052      7.955  -5.663
#> Subject335       -93.747      7.844 -11.951
#> Subject337        31.904      7.844   4.067
#> Subject349       -65.837      7.937  -8.295
#> Subject350       -28.624      7.977  -3.588
#> Subject351       -52.067      7.949  -6.550
#> Subject352       -10.449      7.976  -1.310
#> Subject369       -37.782      7.844  -4.816
#> Subject370       -52.115      7.844  -6.644
#> Subject371       -44.565      7.969  -5.592
#> Subject372       -25.931      7.844  -3.306
#> Days:Subject309  -20.490      2.716  -7.545
#> Days:Subject310  -16.633      2.715  -6.127
#> Days:Subject330  -19.741      2.725  -7.245
#> Days:Subject331  -17.181      2.725  -6.306
#> Days:Subject332  -14.881      2.752  -5.407
#> Days:Subject333  -13.606      2.715  -5.012
#> Days:Subject334  -10.598      2.724  -3.890
#> Days:Subject335  -25.629      2.715  -9.440
#> Days:Subject337   -3.722      2.715  -1.371
#> Days:Subject349   -9.247      2.720  -3.399
#> Days:Subject350   -4.093      2.734  -1.497
#> Days:Subject351  -16.119      2.723  -5.920
#> Days:Subject352   -8.505      2.733  -3.112
#> Days:Subject369  -11.400      2.715  -4.199
#> Days:Subject370   -4.692      2.715  -1.728
#> Days:Subject371  -13.129      2.730  -4.810
#> Days:Subject372  -11.450      2.715  -4.217
#> 
#> Correlation matrix not shown by default, as p = 36 > 12.
#> Use print(x, correlation=TRUE)  or
#>     vcov(x)        if you need it
anova(sleep.semi)
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: Reaction
#>                      F Df Df.res    Pr(>F)    
#> (Intercept)  51325.067  1 97.685 < 2.2e-16 ***
#> Days           543.003  1 97.468 < 2.2e-16 ***
#> Subject         46.915 17 96.801 < 2.2e-16 ***
#> Days:Subject    12.033 17 97.065 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);
#> p-values halved for the boundary (chi-square 0:1)
#>         term  variance logLik logLik.without       LRT  p.value
#>  sub_grp.308 7.691e+03   -701         -710.8 1.942e+01 5.24e-06
#>  sub_grp.309 3.815e+00   -701         -701.1 3.592e-02  0.42484
#>  sub_grp.310 0.000e+00   -701         -701.0 1.906e-06  0.49945
#>  sub_grp.330 7.862e+01   -701         -703.8 5.544e+00  0.00927
#>  sub_grp.331 7.395e+01   -701         -701.2 2.815e-01  0.29787
#>  sub_grp.332 2.831e+03   -701         -723.5 4.489e+01 1.04e-11
#>  sub_grp.333 0.000e+00   -701         -701.0 1.550e-08  0.49995
#>  sub_grp.334 6.941e+01   -701         -703.3 4.463e+00  0.01732
#>  sub_grp.335 0.000e+00   -701         -701.0 1.317e-07  0.49986
#>  sub_grp.337 0.000e+00   -701         -701.0 7.133e-07  0.49966
#>  sub_grp.349 3.301e+01   -701         -701.9 1.639e+00  0.10024
#>  sub_grp.350 3.142e+02   -701         -702.9 3.632e+00  0.02834
#>  sub_grp.351 5.315e+01   -701         -701.4 7.358e-01  0.19551
#>  sub_grp.352 2.739e+02   -701         -705.1 8.012e+00  0.00232
#>  sub_grp.369 3.751e-07   -701         -701.0 0.000e+00  0.50000
#>  sub_grp.370 0.000e+00   -701         -701.0 3.327e-07  0.49977
#>  sub_grp.371 1.602e+02   -701         -704.7 7.277e+00  0.00349
#>  sub_grp.372 0.000e+00   -701         -701.0 2.723e-07  0.49979
R2_glmm(sleep.semi$semer)
#> # Adjusted R2 for Mixed Models
#>       Total       Fixed      Random sub_grp.308 sub_grp.309 sub_grp.310 
#>    "93.43%"    "58.06%"    "35.36%"    "23.48%"     "0.01%"        "0%" 
#> sub_grp.330 sub_grp.331 sub_grp.332 sub_grp.333 sub_grp.334 sub_grp.335 
#>     "0.24%"     "0.23%"     "8.64%"        "0%"     "0.21%"        "0%" 
#> sub_grp.337 sub_grp.349 sub_grp.350 sub_grp.351 sub_grp.352 sub_grp.369 
#>        "0%"      "0.1%"     "0.96%"     "0.16%"     "0.84%"        "0%" 
#> sub_grp.370 sub_grp.371 sub_grp.372 
#>        "0%"     "0.49%"        "0%" 

predout1 <- semipred(sleep.semi, "Subject", "Days")

str(predout1$pred_df)
#> 'data.frame':    1044 obs. of  6 variables:
#>  $ Subject : Factor w/ 18 levels "308","309","310",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Days    : num  0 0.184 0.367 0.551 0.735 ...
#>  $ SE      : num  16.7 15 14.8 14.7 14.4 ...
#>  $ Reaction: num  251 252 252 251 251 ...
#>  $ LL      : num  219 223 223 223 222 ...
#>  $ UL      : num  284 282 281 280 279 ...
predout2 <- semipred(sleep.semi, "Subject", "Days", contr = c(6,1))
#> 332 differs from 308 (pointwise band, level 0.05) for Days: 3.857 to 5.327 (<); 6.000 to 7.000 (>); 7.714 to 9.000 (<).

str(predout2$pred_df)
#> 'data.frame':    58 obs. of  7 variables:
#>  $ Days       : num  0 0.184 0.367 0.551 0.735 ...
#>  $ Reaction   : num  -18.48 -17.43 -15.32 -12.42 -8.98 ...
#>  $ LL         : num  -64.2 -58 -54.4 -50.8 -46.9 ...
#>  $ UL         : num  27.2 23.1 23.8 26 28.9 ...
#>  $ SE_contr   : num  23.3 20.7 20 19.6 19.3 ...
#>  $ includes0  : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ significant: logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#>  - attr(*, "significant")='data.frame':  3 obs. of  3 variables:
#>   ..$ from     : num [1:3] 3.86 6 7.71
#>   ..$ to       : num [1:3] 5.33 7 9
#>   ..$ direction: chr [1:3] "<" ">" "<"
# }
```
