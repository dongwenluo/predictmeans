# `copper` data

Copper grade at sampling points with their X, Y and Z coordinates within
a designated area. Used in the overview article to fit a semiparametric
model with a smooth in Z and a thin-plate spline in X and Y.

## Usage

``` r
data(copper)
```

## Format

A data frame with 439 rows and 4 columns:

- xcoord:

  X coordinate of the sampling point.

- ycoord:

  Y coordinate of the sampling point.

- zcoord:

  Z coordinate of the sampling point.

- grade:

  copper grade.

## Examples

``` r
data(copper)
# \donttest{
mod <- semireg(log(grade) ~ (xcoord + ycoord + zcoord)^2,
               smoothZ = list(grp_z = smZ(zcoord),
                              grp_geo = smZ(cbind(xcoord, ycoord), type = "Ztps")),
               data = copper, prt = FALSE)
#> Warning: Some predictor variables are on very different scales: consider rescaling. 
#> You may also use (g)lmerControl(autoscale = TRUE) to improve numerical stability.
#> boundary (singular) fit: see help('isSingular')
lmerTest::ranova(mod$semer)
#> Warning: Using formula(x) is deprecated when x is a character vector of length > 1.
#>   Consider formula(paste(x, collapse = " ")) instead.
#> ANOVA-like table for random-effects: Single term deletions
#> 
#> Model:
#> log ~ (xcoord + ycoord + zcoord)^2 + (1 | grp_z) + (1 | grp_geo)
#>               npar  logLik    AIC    LRT Df Pr(>Chisq)    
#> <none>          10 -159.20 338.41                         
#> (1 | grp_z)      9 -199.65 417.30 80.899  1  < 2.2e-16 ***
#> (1 | grp_geo)    9 -207.68 433.36 96.951  1  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
