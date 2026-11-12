# Study 1b results

Rejection rates at the 5% level. Monte Carlo range for a true rate of 5%: 0.036-0.064 (1000 data sets),
0.031-0.069 (500).

## Pairwise comparisons

| cell | method | per pair (mean) | per pair (max) | FWER Holm | FWER max-t / Tukey | power (mean) |
|---|---|---|---|---|---|---|
| pr_lev_global | joint | 0.057 | 0.061 | 0.061 | 0.067 |  |
| pr_lev_global | flpair | 0.054 | 0.060 | 0.045 |  |  |
| pr_lev_global | flterm | 0.053 | 0.058 | 0.044 | 0.056 |  |
| pr_lev_global | param | 0.052 | 0.055 | 0.045 | 0.051 |  |
| pr_rcbd_global | joint | 0.054 | 0.068 | 0.028 | 0.058 |  |
| pr_rcbd_global | flpair | 0.058 | 0.071 | 0.033 |  |  |
| pr_rcbd_global | flterm | 0.057 | 0.072 | 0.037 | 0.065 |  |
| pr_rcbd_global | param | 0.056 | 0.069 | 0.047 | 0.059 |  |
| pr_rcbd_partial | joint | 0.044 | 0.047 | 0.040 | 0.032 | 1.000 |
| pr_rcbd_partial | flpair | 0.052 | 0.057 | 0.045 |  | 1.000 |
| pr_rcbd_partial | flterm | 0.051 | 0.052 | 0.050 | 0.042 | 1.000 |
| pr_rcbd_partial | param | 0.048 | 0.050 | 0.048 | 0.037 | 1.000 |
| pr_rcbd_skew_partial | joint | 0.046 | 0.052 | 0.027 | 0.024 | 0.998 |
| pr_rcbd_skew_partial | flpair | 0.051 | 0.057 | 0.040 |  | 0.999 |
| pr_rcbd_skew_partial | flterm | 0.052 | 0.058 | 0.036 | 0.030 | 0.999 |
| pr_rcbd_skew_partial | param | 0.048 | 0.055 | 0.034 | 0.029 | 0.999 |
| pr_sp_A_partial | joint | 0.055 | 0.059 | 0.036 | 0.034 | 0.807 |
| pr_sp_A_partial | flpair | 0.062 | 0.065 | 0.055 |  | 0.832 |
| pr_sp_A_partial | flterm | 0.062 | 0.068 | 0.054 | 0.047 | 0.831 |
| pr_sp_A_partial | param | 0.050 | 0.055 | 0.038 | 0.028 | 0.809 |
| pr_sp_B_partial | joint | 0.048 | 0.048 | 0.048 | 0.021 | 1.000 |
| pr_sp_B_partial | flpair | 0.054 | 0.054 | 0.054 |  | 1.000 |
| pr_sp_B_partial | flterm | 0.049 | 0.049 | 0.049 | 0.018 | 1.000 |
| pr_sp_B_partial | param | 0.052 | 0.052 | 0.052 | 0.017 | 1.000 |

