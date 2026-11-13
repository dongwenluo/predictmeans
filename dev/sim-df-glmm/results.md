# Results: sim-df-glmm

Percent of true null hypotheses rejected at the 5% level (two-sided t tests of single contrasts,
F tests of terms). With R data sets a true 5% lies in 5 +- 2.58 sqrt(0.05 x 0.95 / R) x 100 (99%).
Bold: outside that range.

## gl_bern_rcbd

binary, 8 blocks x 3 treatments x 4 units; logit 0; block SD 0.5

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | resid | cont | BW | nlme | nlmeC | z |
|---|---|---|---|---|---|---|
| t12 | 3.9 | 3.9 | 3.9 | 3.9 | 3.9 | 4.0 |
| t13 | 5.4 | 5.4 | 5.4 | 5.4 | 5.4 | 5.6 |
| FTrt | 4.4 | 4.4 | 4.4 | 4.4 | 4.4 | 5.2 |

Mean df: resid_t12 92.0, resid_t13 92.0, cont_t12 86.0, cont_t13 86.0, BW_t12 86.0, BW_t13 86.0, nlme_t12 86.0, nlme_t13 86.0, nlmeC_t12 86.0, nlmeC_t13 86.0, z_t12 Inf, z_t13 Inf

## gl_binom_clus

binomial (10 trials), 2 treatments x 4 clusters x 6 units, within-cluster factor W (2 levels); logit -0.5; cluster SD 0.5

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | resid | cont | BW | nlme | nlmeC | z |
|---|---|---|---|---|---|---|
| tTrt | **12.8** | **12.7** | **7.7** | **7.7** | **7.7** | **13.7** |
| tW | 3.6 | 3.6 | **1.5** | 3.6 | 3.6 | 4.4 |
| tTrtinW | **11.0** | **11.0** | 5.0 | **11.0** | 5.0 | **12.0** |
| FTW | 4.5 | 4.5 | **1.5** | 4.5 | 4.5 | 5.1 |

Mean df: resid_tTrt 43.0, resid_tW 43.0, cont_tTrt 38.0, cont_tW 38.0, BW_tTrt 6.0, BW_tW 6.0, nlme_tTrt 6.0, nlme_tW 38.0, nlmeC_tTrt 6.0, nlmeC_tW 38.0, z_tTrt Inf, z_tW Inf

## gl_binom_clus_8

as gl_binom_clus with 8 clusters per treatment

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | resid | cont | BW | nlme | nlmeC | z |
|---|---|---|---|---|---|---|
| tTrt | **8.4** | **8.2** | 6.4 | 6.4 | 6.4 | **8.7** |
| tW | 4.5 | 4.4 | **2.9** | 4.4 | 4.4 | 4.8 |
| tTrtinW | 6.3 | 6.3 | 4.1 | 6.3 | 4.1 | 6.6 |
| FTW | 3.4 | 3.4 | **2.1** | 3.4 | 3.4 | 3.9 |

Mean df: resid_tTrt 91.0, resid_tW 91.0, cont_tTrt 78.0, cont_tW 78.0, BW_tTrt 14.0, BW_tW 14.0, nlme_tTrt 14.0, nlme_tW 78.0, nlmeC_tTrt 14.0, nlmeC_tW 78.0, z_tTrt Inf, z_tW Inf

## gl_binom_crossed

binomial (10 trials), 6 rows x 8 columns crossed, T (2 levels) on columns, W within cells; logit -0.5; SD row 0.5, column 0.5

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | io | resid | cont | z |
|---|---|---|---|---|
| tT | **7.2** | **12.0** | **12.0** | **12.3** |
| tW | 3.7 | 3.7 | 3.7 | 4.2 |
| tTinW | 5.4 | **9.1** | **9.1** | **9.5** |
| tWinT | 4.3 | 4.4 | 4.3 | 4.9 |
| FT | **7.2** | **12.0** | **12.0** | **12.3** |

Mean df: io_tT 6.0, io_tW 81.0

## gl_pois_crossed

Poisson, 4 rows x 12 columns crossed (random), treatment T (3 levels) on columns, W (2 levels) within each row x column; log mean 1.5; SD row 0.3, column 0.3; df methods of the package (inner-outer extended to crossed factors)

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | io | resid | cont | z |
|---|---|---|---|---|
| tT | **7.0** | **9.8** | **9.8** | **10.1** |
| tW | 4.8 | 4.8 | 4.8 | 5.1 |
| tTinW | 5.2 | **8.0** | **8.0** | **8.5** |
| tWinT | 3.9 | 3.9 | 3.9 | 4.6 |
| FT | **7.1** | **13.0** | **12.9** | **13.9** |

Mean df: io_tT 9.0, io_tW 78.0

## gl_pois_sp

Poisson split plot, 6 blocks x 4 A (whole plots) x 3 B; log mean 2; SD block 0.3, whole plot 0.3

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | resid | cont | BW | nlme | nlmeC | z |
|---|---|---|---|---|---|---|
| tA | **9.0** | **7.3** | **7.3** | **7.3** | **7.3** | **9.5** |
| tB | 4.6 | 4.3 | **2.6** | 4.3 | 4.3 | 5.2 |
| tAinB | 5.2 | 5.2 | 5.1 | 5.2 | 3.8 | 5.4 |
| tBinA | **2.9** | **2.5** | **2.5** | **2.5** | **2.5** | 3.5 |
| FA | **12.1** | **7.2** | **7.2** | **7.2** | **7.2** | **14.0** |
| FAB | **2.5** | **2.0** | **1.6** | **2.0** | **2.0** | 4.0 |

Mean df: resid_tA 58.0, resid_tB 58.0, cont_tA 15.0, cont_tB 40.0, BW_tA 15.0, BW_tB 10.0, nlme_tA 15.0, nlme_tB 40.0, nlmeC_tA 15.0, nlmeC_tB 40.0, z_tA Inf, z_tB Inf

## gl_pois_sp_4blk

as gl_pois_sp with 4 blocks

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | resid | cont | BW | nlme | nlmeC | z |
|---|---|---|---|---|---|---|
| tA | **10.7** | **8.1** | **8.1** | **8.1** | **8.1** | **11.8** |
| tB | 4.9 | 4.4 | **1.5** | 4.4 | 4.4 | 5.5 |
| tAinB | **6.9** | 6.6 | 5.8 | 6.6 | 4.6 | **7.1** |
| tBinA | 4.5 | 4.3 | 4.1 | 4.3 | 4.3 | 5.3 |
| FA | **15.2** | **8.0** | **8.0** | **8.0** | **8.0** | **17.9** |
| FAB | **2.4** | **1.6** | **1.2** | **1.6** | **1.6** | 4.3 |

Mean df: resid_tA 34.0, resid_tB 34.0, cont_tA 9.0, cont_tB 24.0, BW_tA 9.0, BW_tB 6.0, nlme_tA 9.0, nlme_tB 24.0, nlmeC_tA 9.0, nlmeC_tB 24.0, z_tA Inf, z_tB Inf

## gl_pois_sp_low

as gl_pois_sp with log mean 0.5 (mean count 1.6)

1000 data sets (0 failed); a true 5% lies in 3.2-6.8%.

| | resid | cont | BW | nlme | nlmeC | z |
|---|---|---|---|---|---|---|
| tA | **8.7** | 6.4 | 6.4 | 6.4 | 6.4 | **9.3** |
| tB | 4.7 | 4.5 | **2.5** | 4.5 | 4.5 | 5.1 |
| tAinB | 5.2 | 5.0 | 5.0 | 5.0 | 3.6 | 6.0 |
| tBinA | 4.1 | 4.0 | 3.9 | 4.0 | 4.0 | 4.7 |
| FA | **8.4** | 5.5 | 5.5 | 5.5 | 5.5 | **9.2** |
| FAB | **2.2** | **1.9** | **1.7** | **1.9** | **1.9** | **3.1** |

Mean df: resid_tA 58.0, resid_tB 58.0, cont_tA 15.0, cont_tB 40.0, BW_tA 15.0, BW_tB 10.0, nlme_tA 15.0, nlme_tB 40.0, nlmeC_tA 15.0, nlmeC_tB 40.0, z_tA Inf, z_tB Inf

