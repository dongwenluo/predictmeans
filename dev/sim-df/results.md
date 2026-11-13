# Results: sim-df

Percent of true null hypotheses rejected at the 5% level (two-sided t tests of single contrasts,
F tests of terms). With R data sets a true 5% lies in 5 +- 2.58 sqrt(0.05 x 0.95 / R) x 100 (99%).
Bold: outside that range.

## rm_ar1

repeated measures, 2 groups x 6 subjects x 4 times, AR(1) rho 0.6, var 1; gls corAR1

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW | z | resid |
|---|---|---|---|---|---|---|
| tG | 3.9 | 4.0 | 4.0 | 3.7 | 6.3 | 5.9 |
| tT | 5.4 | 5.4 | 5.4 | 4.9 | 6.4 | 5.4 |
| tGatT | 5.4 | 5.3 | 5.3 | 5.4 | 6.5 | 5.7 |
| tTinG | 5.0 | 5.0 | 5.0 | 4.8 | 6.0 | 5.0 |
| FT | 5.1 | 5.1 | 5.0 | 5.1 | 6.6 | 5.1 |
| FGT | 3.6 | 3.6 | 3.7 | 3.8 | 5.5 | 4.7 |

## rm_cs_lme

2 groups x 6 subjects x 4 times, subject variance 1, residual 1; lme random intercept, 10% missing

997 data sets (3 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW | z | nlme | nlmeC |
|---|---|---|---|---|---|---|---|
| tG | 5.0 | 5.0 | 5.0 | 5.0 | **7.4** | 5.0 | 5.0 |
| tT | 5.3 | 5.3 | 5.3 | 5.3 | **6.9** | 5.3 | 5.3 |
| tGatT | 4.2 | 4.2 | 4.2 | 4.5 | 5.9 | 4.5 | **3.1** |
| tTinG | 5.0 | 5.0 | 5.0 | 5.0 | 6.0 | 5.0 | 5.0 |
| FT | 5.5 | 5.4 | 5.6 | 5.6 | **8.3** | 5.6 | 5.6 |
| FGT | 5.6 | 5.6 | 5.6 | 5.6 | **7.9** | 5.6 | 5.6 |

Errors: system is computationally singular: reciprocal condition number = 1.75235e-17 1; system is computationally singular: reciprocal condition number = 1.88007e-17 1; system is computationally singular: reciprocal condition number = 8.39212e-18 1

## rm_un

as rm_ar1 with SDs 1, 1.3, 1.6, 2 and correlation 0.7^lag; gls unstructured (corSymm + varIdent)

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW | z | resid |
|---|---|---|---|---|---|---|
| tG | 4.9 | 4.9 | 4.9 | 4.9 | **7.2** | 6.2 |
| tT | 5.0 | 5.0 | 5.0 | 5.0 | **8.2** | **7.3** |
| tGatT | 5.7 | 5.7 | 5.7 | 5.7 | **8.4** | **7.1** |
| tTinG | 4.5 | 4.5 | 4.5 | 4.5 | **6.8** | 6.3 |
| FT | 4.2 | 4.2 | **9.9** | **9.9** | **18.4** | **15.9** |
| FGT | 5.2 | 5.2 | **8.9** | **8.9** | **17.2** | **15.5** |

## rm_un_12

as rm_un with 12 subjects per group

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW | z | resid |
|---|---|---|---|---|---|---|
| tG | 4.7 | 4.7 | 4.7 | 4.7 | 6.0 | 5.5 |
| tT | 5.7 | 5.7 | 5.7 | 5.7 | **6.8** | 6.4 |
| tGatT | 4.7 | 4.7 | 4.7 | 4.7 | 5.5 | 5.4 |
| tTinG | 5.0 | 5.0 | 5.0 | 5.0 | 6.2 | 6.0 |
| FT | 4.1 | 4.1 | 5.8 | 5.8 | **8.9** | **8.6** |
| FGT | 6.1 | 6.1 | **7.6** | **7.6** | **11.2** | **10.1** |

## rm_un_miss

as rm_un, 10% of the observations missing at random

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW | z | resid |
|---|---|---|---|---|---|---|
| tG | 4.6 | 4.7 | 4.7 | 4.8 | **8.6** | **7.6** |
| tT | 6.1 | 5.5 | 6.2 | 6.4 | **10.2** | **8.5** |
| tGatT | 5.6 | 4.8 | 5.2 | 5.8 | **9.6** | **8.8** |
| tTinG | **7.3** | 6.6 | **7.8** | **8.2** | **11.2** | **10.0** |
| FT | **8.7** | 5.9 | **15.7** | **18.2** | **29.7** | **26.5** |
| FGT | **8.5** | 6.0 | **15.4** | **18.1** | **28.3** | **24.1** |

## sp_4blk_miss

split plot 4 blocks x 4 A x 3 B, 15% missing at random

998 data sets (2 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | cont | BW | z | nlme | nlmeC |
|---|---|---|---|---|---|---|---|---|
| tA | 5.7 | 5.4 | 6.4 | 5.9 | 5.9 | **8.9** | 5.9 | 5.9 |
| tB | 5.6 | 5.7 | 6.1 | 6.1 | 3.6 | **7.6** | 6.1 | 6.1 |
| tAinB | 6.0 | 5.9 | 6.3 | 5.7 | 4.6 | **7.7** | 5.7 | 4.1 |
| tBinA | 4.8 | 4.9 | 4.9 | 4.9 | 4.4 | 5.6 | 4.9 | 4.9 |
| FA | 5.8 | 5.3 | 6.6 | 6.2 | 6.2 | **13.3** | 6.2 | 6.2 |
| FAB | 5.7 | 5.7 | 6.0 | 5.8 | **2.9** | **10.4** | 5.8 | 5.8 |

Errors: Singularity in backsolve at level 0, block 1 2

## sp_bal

split plot 6 blocks x 4 A (whole plots) x 3 B; var block 1, whole plot 1, residual 1; lmer and lme

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | cont | BW | z | nlme | nlmeC |
|---|---|---|---|---|---|---|---|---|
| tA | 5.9 | 5.9 | 5.9 | 5.9 | 5.9 | **7.2** | 5.9 | 5.9 |
| tB | 5.8 | 5.8 | 5.8 | 5.8 | 4.0 | 6.4 | 5.8 | 5.8 |
| tAinB | 5.3 | 5.3 | 5.3 | 5.4 | 5.3 | 5.7 | 5.4 | 4.5 |
| tBinA | 5.3 | 5.3 | 5.3 | 5.3 | 5.1 | 6.0 | 5.3 | 5.3 |
| FA | 4.0 | 4.0 | 4.0 | 4.0 | 4.0 | **8.3** | 4.0 | 4.0 |
| FAB | 5.0 | 5.0 | 5.0 | 5.0 | 4.7 | **6.9** | 5.0 | 5.0 |

## sp_miss_15

as sp_bal, 15% of the observations missing at random

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | cont | BW | z | nlme | nlmeC |
|---|---|---|---|---|---|---|---|---|
| tA | 5.6 | 5.3 | 5.5 | 5.6 | 5.6 | **7.2** | 5.6 | 5.6 |
| tB | 4.6 | 4.6 | 4.7 | 4.7 | **3.1** | 6.0 | 4.7 | 4.7 |
| tAinB | 3.8 | 3.8 | 4.0 | 3.9 | 3.4 | 4.8 | 3.9 | **2.7** |
| tBinA | 5.1 | 5.2 | 5.3 | 5.1 | 4.8 | 6.1 | 5.1 | 5.1 |
| FA | 5.6 | 5.6 | 6.0 | 5.6 | 5.6 | **10.6** | 5.6 | 5.6 |
| FAB | 3.5 | 3.5 | 3.8 | 3.5 | **2.6** | **6.9** | 3.5 | 3.5 |

## sp_miss_wp

as sp_bal, the whole plot (block 1, A1) missing

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | cont | BW | z | nlme | nlmeC |
|---|---|---|---|---|---|---|---|---|
| tA | 5.7 | 5.8 | 5.7 | 5.7 | 5.7 | **7.7** | 5.7 | 5.7 |
| tB | 5.2 | 5.2 | 5.2 | 5.2 | 3.3 | 5.6 | 5.2 | 5.2 |
| tAinB | 6.1 | 6.2 | 6.3 | 6.6 | 6.0 | **7.7** | 6.6 | 4.8 |
| tBinA | 5.4 | 5.4 | 5.4 | 5.4 | 5.3 | 5.9 | 5.4 | 5.4 |
| FA | 5.1 | 5.1 | 5.2 | 5.1 | 5.1 | **9.6** | 5.1 | 5.1 |
| FAB | 4.9 | 4.9 | 4.9 | 4.9 | 4.4 | 6.5 | 4.9 | 4.9 |

## sp_wp_small

as sp_bal, whole-plot variance 0.1

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | cont | BW | z | nlme | nlmeC |
|---|---|---|---|---|---|---|---|---|
| tA | 3.3 | **3.1** | 4.0 | 3.3 | 3.3 | 4.7 | 3.3 | 3.3 |
| tB | 5.9 | 5.9 | 5.9 | 5.9 | 3.6 | 6.7 | 5.9 | 5.9 |
| tAinB | 4.2 | 4.2 | 4.2 | 4.0 | 4.0 | 4.9 | 4.0 | 3.4 |
| tBinA | 5.8 | 5.9 | 6.1 | 5.8 | 5.6 | 6.5 | 5.8 | 5.8 |
| FA | **2.9** | **2.2** | 3.6 | **2.9** | **2.9** | 6.6 | **2.9** | **2.9** |
| FAB | 5.2 | 5.4 | 5.6 | 5.2 | 4.3 | **8.1** | 5.2 | 5.2 |

