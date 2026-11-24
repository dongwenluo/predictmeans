# Potency of Analgesic Drugs in Mice

Data for the comparison of the effectiveness of three analgesic drugs to
a standard drug, morphine (Finney, Probit analysis, 3rd Edition 1971,
p.103). 14 groups of mice were tested for response to the drugs at a
range of doses.

## Usage

``` r
data(Drug)
```

## Format

A data frame with 14 rows and 5 columns:

- Drug:

  factor, 4 drugs (Morphine, Amidone, ...).

- Dose:

  dose (units as in Finney, 1971).

- N:

  number of mice in the group.

- R:

  number of mice responding.

- log10Dose:

  `log10(Dose)`.

## Source

Finney, D. J. (1971). *Probit Analysis*, 3rd edition. Cambridge
University Press, p. 103.

## Examples

``` r
data(Drug)
gm <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial(link = "probit"), data = Drug)
predictmeans(gm, "Drug", covariate = c(log10Dose = 0.5))$means
#>          Drug log10Dose   estimate  std.error  df    conf.low    conf.high
#> 1    Morphine       0.5 -0.1450208 0.07208002 Inf -0.28629509 -0.003746607
#> 2     Amidone       0.5  0.0928568 0.08097719 Inf -0.06585557  0.251569167
#> 3 Phenadoxone       0.5  1.2145314 0.10639704 Inf  1.00599705  1.423065795
#> 4   Pethidine       0.5 -1.3249360 0.11367167 Inf -1.54772836 -1.102143606
#>   comparison.low comparison.high   response response.std.error
#> 1     -0.2495048     -0.04053689 0.44234721         0.02845497
#> 2     -0.0245241      0.21023770 0.53699134         0.03216625
#> 3      1.0603030      1.36875980 0.88772757         0.02030152
#> 4     -1.4897093     -1.16016263 0.09259619         0.01885255
#>   response.conf.low response.conf.high n estimable
#> 1        0.38732606          0.4985053 3      TRUE
#> 2        0.47374641          0.5993130 3      TRUE
#> 3        0.84279151          0.9226415 3      TRUE
#> 4        0.06084386          0.1351996 5      TRUE
```
