# Results: sim-df-kr

Percent of true null hypotheses rejected at the 5% level (two-sided t tests of single contrasts,
F tests of terms). With R data sets a true 5% lies in 5 +- 2.58 sqrt(0.05 x 0.95 / R) x 100 (99%).
Bold: outside that range.

## cs_lmer_miss

2 groups x 6 subjects x 4 times; subject var 1, residual 1; 20% missing; lmer (1 | Subj)

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW |
|---|---|---|---|---|
| tG | 4.2 | 4.2 | 4.3 | 4.2 |
| tT | 4.2 | 4.2 | 4.6 | 4.4 |
| tGatT | 4.2 | 4.2 | 4.3 | 4.3 |
| tTinG | 5.4 | 5.4 | 5.5 | 5.4 |
| FT | 4.9 | 4.9 | 5.1 | 5.0 |
| FGT | 4.9 | 4.8 | 5.0 | 4.9 |

## rm_ar1_miss20

2 groups x 6 subjects x 6 times; AR(1) rho 0.6, SDs 1 to 2; 20% missing; gls corAR1 + varIdent

999 data sets (1 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW |
|---|---|---|---|---|
| tG | 4.3 | 4.3 | 4.3 | 4.1 |
| tT | 6.3 | 6.5 | 6.6 | **8.1** |
| tGatT | 5.1 | 5.0 | 5.4 | **8.8** |
| tTinG | 6.0 | 5.5 | 6.0 | **7.5** |
| FT | 5.7 | 5.3 | **8.4** | **12.8** |
| FGT | 4.5 | 4.4 | 6.6 | **11.8** |

Errors: computed "gls" fit is singular, rank 12 1

## rm_un_miss20

2 groups x 12 subjects x 4 times; SDs 1, 1.3, 1.6, 2, correlation 0.7^lag; 20% missing; gls unstructured

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW |
|---|---|---|---|---|
| tG | **7.0** | **6.9** | **7.1** | **7.1** |
| tT | 6.0 | 5.9 | 6.5 | 6.6 |
| tGatT | 6.4 | 6.1 | 6.6 | 6.7 |
| tTinG | 5.2 | 5.2 | 6.0 | 6.0 |
| FT | **6.9** | 5.7 | **11.9** | **12.4** |
| FGT | 6.2 | 5.6 | **9.7** | **11.7** |

## rs_lmer_miss

2 groups x 8 subjects x 4 times; random intercept (var 1) and slope in time (var 0.2), residual 1; 20% missing; lmer (t | Subj)

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | BW |
|---|---|---|---|---|
| tG | 6.1 | 6.2 | 6.3 | 6.2 |
| tT | 5.5 | 4.4 | 6.2 | **6.9** |
| tGatT | 5.8 | 5.9 | **6.8** | **7.4** |
| tTinG | 5.2 | 4.2 | 5.7 | 6.5 |
| FT | 5.4 | 5.0 | 6.3 | **7.3** |
| FGT | 5.6 | 4.5 | **6.8** | **7.9** |

## sp_miss_30

split plot 6 blocks x 4 A x 3 B (variances 1, 1, 1), 30% of the observations missing; lmer

996 data sets (4 failed); a true 5% lies in 3.2-6.8%.

| | KR | KRSAS | Satt | cont |
|---|---|---|---|---|
| tA | 5.3 | 5.2 | 5.8 | 5.8 |
| tB | 5.2 | 5.3 | 5.7 | 5.4 |
| tAinB | 4.6 | 4.5 | 5.1 | 4.5 |
| tBinA | 4.7 | 4.6 | 5.1 | 4.9 |
| FA | **3.2** | 3.4 | 3.8 | 3.9 |
| FAB | 4.9 | 4.9 | 5.6 | 4.9 |

Errors: NA/NaN/Inf in foreign function call (arg 1) 2; missing value where TRUE/FALSE needed 2

