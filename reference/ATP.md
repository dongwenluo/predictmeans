# ATP Content of Preserved Dog Hearts

Data from an experiment to study the effects of preserving liquids on
the enzyme content of dog hearts. There were 23 hearts and two treatment
factors, A and B, each at two levels. Measurements were made of ATP as a
percentage of total enzyme in the heart, at one and two hourly intervals
during a twelve hour period following initial preservation.

## Usage

``` r
data(ATP)
```

## Format

A data frame with 230 rows (23 hearts, 10 times) and 5 columns:

- heart:

  factor, 23 dog hearts.

- time:

  factor, time of measurement (hours after preservation; 10 levels).

- A:

  factor, treatment A (2 levels).

- B:

  factor, treatment B (2 levels).

- ATP:

  ATP as a percentage of total enzyme in the heart (%).

## Examples

``` r
data(ATP)
fm <- lmer(ATP ~ A*B*time + (1|heart), data = ATP)
predictmeans(fm, "A:B", pairwise = TRUE)$means
#>   A B estimate std.error df conf.low conf.high comparison.low comparison.high
#> 1 1 1 66.72767  2.284438 19 61.94628  71.50905       63.34497        70.11037
#> 2 1 2 66.70860  2.502476 19 61.47086  71.94634       63.00304        70.41416
#> 3 2 1 72.33483  2.284438 19 67.55345  77.11622       68.95213        75.71753
#> 4 2 2 66.81367  2.284438 19 62.03228  71.59505       63.43097        70.19637
#>   group  n estimable
#> 1     A 60      TRUE
#> 2     A 50      TRUE
#> 3     A 60      TRUE
#> 4     A 60      TRUE
```
