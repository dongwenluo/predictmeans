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
anova(mod)   # fixed effects (Kenward-Roger) and a test of each smooth term
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
#> 
#> Response: log(grade)
#>                     F Df  Df.res    Pr(>F)    
#> (Intercept)    0.0000  1    12.6   0.99676    
#> xcoord         0.1862  1    12.6   0.67340    
#> ycoord         0.0000  1    12.6   0.99523    
#> zcoord         5.8288  1 -3789.6       NaN    
#> xcoord:ycoord  0.2840  1    15.4   0.60167    
#> xcoord:zcoord  4.0104  1   428.9   0.04585 *  
#> ycoord:zcoord 26.3669  1   424.5 4.309e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Smooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);
#> p-values halved for the boundary (chi-square 0:1)
#>     term  variance logLik logLik.without   LRT p.value
#>    grp_z 4.170e-03 -163.9         -199.7 71.41  <2e-16
#>  grp_geo 1.295e-06 -163.9         -207.7 87.46  <2e-16
# }
```
