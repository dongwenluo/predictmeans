# Multicentre Clinical Trial of a Skin Treatment

Favourable responses in a multicentre randomized clinical trial of a
topical drug against a control (Beitler & Landis 1985, Biometrics).

## Usage

``` r
data(Clinical)
```

## Format

A data frame with 16 rows and 4 columns:

- Clinic:

  factor, 8 centres.

- Treatment:

  factor, `control` or `drug`.

- Favorable:

  number of patients with a favourable response.

- Total:

  number of patients in the clinic and treatment group.

## Source

Beitler, P. J. and Landis, J. R. (1985). A mixed-effects model for
categorical data. *Biometrics*, 41, 991-1000.

## Examples

``` r
data(Clinical)
gm <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1|Clinic),
            family = binomial, data = Clinical)
predictmeans(gm, "Treatment")$means   # probability scale
#>   Treatment   estimate std.error df  conf.low conf.high comparison.low
#> 1   control -1.1964721 0.5528792  7 -2.503824 0.1108795     -1.5534398
#> 2      drug -0.4582679 0.5472385  7 -1.752281 0.8357456     -0.8115937
#>   comparison.high  response response.std.error response.conf.low
#> 1      -0.8395044 0.2321034         0.09854045        0.07559056
#> 2      -0.1049421 0.3873968         0.12987093        0.14775968
#>   response.conf.high n estimable
#> 1          0.5276915 8      TRUE
#> 2          0.6975684 8      TRUE
```
