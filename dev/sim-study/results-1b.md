# Study 1b results

Rejection rates at the 5% level. Monte Carlo range for a true rate of 5%: 0.036-0.064 (1000 data sets),
0.031-0.069 (500).

## Pairwise comparisons

| cell | method | per pair (mean) | per pair (max) | FWER Holm | FWER max-t / Tukey | power (mean) |
|---|---|---|---|---|---|---|
| pr_rcbd_global | joint | 0.054 | 0.068 | 0.028 | 0.058 |  |
| pr_rcbd_global | flpair | 0.058 | 0.071 | 0.033 |  |  |
| pr_rcbd_global | flterm | 0.057 | 0.072 | 0.037 | 0.065 |  |
| pr_rcbd_global | param | 0.056 | 0.069 | 0.047 | 0.059 |  |
| pr_rcbd_partial | joint | 0.044 | 0.047 | 0.040 | 0.032 | 1.000 |
| pr_rcbd_partial | flpair | 0.052 | 0.057 | 0.045 |  | 1.000 |
| pr_rcbd_partial | flterm | 0.051 | 0.052 | 0.050 | 0.042 | 1.000 |
| pr_rcbd_partial | param | 0.048 | 0.050 | 0.048 | 0.037 | 1.000 |
| pr_sp_A_partial | joint | 0.055 | 0.059 | 0.036 | 0.034 | 0.807 |
| pr_sp_A_partial | flpair | 0.062 | 0.065 | 0.055 |  | 0.832 |
| pr_sp_A_partial | flterm | 0.062 | 0.068 | 0.054 | 0.047 | 0.831 |
| pr_sp_A_partial | param | 0.050 | 0.055 | 0.038 | 0.028 | 0.809 |

