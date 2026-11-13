# Calculate degree of freedom of a modelterm (contrast) for a lmer model

Calculate the degree of freedom of a modelterm (contrast) for a `lmer`
model using "Kenward-Roger" or "Satterthwaite" method.

## Usage

``` r
df_term(model, modelterm, covariate=NULL, ctrmatrix=NULL, ctrnames=NULL, 
          type=c("Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average"), vcov_adj=NULL)
```

## Arguments

- model:

  Model object returned by `lmer`.

- modelterm:

  Name (in "quotes") for indicating which factor term's degree of
  freedom to be calculated. The `modelterm` must be given exactly as it
  appears in the model formlar, e.g. "A" or "A:B".

- covariate:

  Name (in "quotes") of one the covariate variables in the `model`.

- ctrmatrix:

  A specified contrast matrix. If `ctrmatrix` isn't NULL, the programe
  will ignore modelterm and calculate degree of freedom for the
  `ctrmatrix`.

- ctrnames:

  Names of the specified contrasts, e.g. c("A vs D", "C vs B", ...)

- type:

  Name (in "quote") for indicating a method for claculating degree of
  freedom. The choices are "Kenward-Roger", "Satterthwaite" and
  "containment". The default method is "Kenward-Roger". Kenward-Roger
  needs a REML fit; for models fitted by ML, Satterthwaite is used with
  a warning. "containment" gives the SAS `DDFM=CONTAIN` df of
  `modelterm` (needed then) for every row, and "between-within" the SAS
  `DDFM=BETWITHIN` df (see
  [`predictmeans`](https://dongwenluo.github.io/predictmeans/reference/predictmeans.md)).
  "Kenward-Roger" uses the observed information of the variance
  parameters, as SAS `DDFM=KR`; "KR-expected" the expected information,
  as pbkrtest and emmeans; "KR-average" the average information, as
  Genstat.

- vcov_adj:

  Optional Kenward-Roger adjusted covariance matrix of the fixed
  effects, as returned by `pbkrtest::vcovAdj(model)`, to avoid computing
  it again when `df_term` is called several times for the same model.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
library(predictmeans)
# ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
df_term(fm, "nitro:Variety")
#>   0:Golden Rain    0:Marvellous       0:Victory 0.2:Golden Rain  0.2:Marvellous 
#>        16.08144        16.08144        16.08144        16.08144        16.08144 
#>     0.2:Victory 0.4:Golden Rain  0.4:Marvellous     0.4:Victory 0.6:Golden Rain 
#>        16.08144        16.08144        16.08144        16.08144        16.08144 
#>  0.6:Marvellous     0.6:Victory 
#>        16.08144        16.08144 
## The contrast has a contrast matrix as follows:
#     0:Golden Rain 0:Marvellous 0:Victory 
#[1,]            -1            0         1 
#[2,]             0            0         1 
#     0.2:Golden Rain 0.2:Marvellous 0.2:Victory 
#[1,]               0              0           0 
#[2,]               0              0           0 
#     0.4:Golden Rain  0.4:Marvellous 0.4:Victory
#[1,]               0               0           0
#[2,]               0              -1           0
#      0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#[1,]                0              0           0
#[2,]                0              0           0

# 1. Enter above contrast matrix into a pop up window, then close the window
if (interactive()) df_term(fm, "nitro:Variety")
 
# 2. Construct the contrast matrix directly
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
df_term(fm, ctrmatrix=cm, type="Satterthwaite")
#> [1] 47.32750 45.00004
```
