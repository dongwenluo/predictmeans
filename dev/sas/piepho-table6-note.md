Subject: Possible error in Table 6 of Piepho (2023), Biometrical Journal 65(7), e2200290

Dear Professor Piepho,

While writing an R implementation of your adjusted R² for GLMMs (function `R2_model()` in the R package
predictmeans), we reproduced the examples in your 2023 paper. Tables 2, 3 and 5 agree with our results.
In Table 6, however, 12 of the 16 rows appear to have been computed with the wrong columns of the design matrix.

**What we found.** For a covariate subset with k covariates, the published fixed share is obtained when the k
estimated slopes are multiplied by the *first k columns* of the full design (YR, B1, B2, B3), not by the columns of
the covariates in the subset. The two coincide only for the subsets YR; YR, B1; YR, B1, B2; and the full model, which
is why those four rows (and Table 5) are correct.

**A likely cause.** In the program for Tables 5 and 6 (`Tables 5 & 6 foot_mv in one go.sas`), the design matrix
passed to `%omega` for the lamb-level R² comes from an auxiliary GLIMMIX call whose model statement is fixed at
`model Shape=yr b1 b2 b3;` (with `outdesign(X Z novar)=my_XZ`). The macro then keeps the first `nrows_beta`
columns of that design (`keep _X1-_X&nrows_beta;` in `Omega.sas`), where `nrows_beta` is the length of the
estimated beta after the second threshold is dropped. For a subset such as B1 alone, beta is (intercept, B1), so
the design kept is (intercept, YR): the B1 slope is multiplied by the YR column. The program we have covers only
the full model, so we cannot see the loop over the 16 subsets, but the "first k" pattern below is exactly what this
pairing gives. Using the subset's own model statement in the auxiliary GLIMMIX call (or keeping the design
columns by name) would avoid it.

**Evidence.** We refitted all 16 proportional odds models with PROC GLIMMIX (`method=quad`, `link=cumprobit`,
random sire intercept; your data and code) and computed the lamb-level R² in two ways: with the subset's own
covariates ("correct") and with the first k columns of YR, B1, B2, B3 ("first k"). The "first k" values reproduce
your Table 6 to within 0.01 in every row; the AICs agree with your table. Independent fits in R
(`ordinal::clmm`, 10 quadrature points) agree with the "correct" SAS values to within 0.01.

| Covariates | Published Ω_β / Ω_u / Ω_βu | "First k" (SAS) | Correct (SAS) |
|---|---|---|---|
| none | 0.00 / 11.58 / 11.58 | 0.00 / 11.58 / 11.58 | 0.00 / 11.58 / 11.58 |
| YR | 1.23 / 9.99 / 11.21 | 1.23 / 9.99 / 11.21 | 1.23 / 9.99 / 11.21 |
| B1 | 6.47 / 8.25 / 14.72 | 6.47 / 8.25 / 14.72 | **3.30 / 8.53 / 11.83** |
| B2 | 3.23 / 10.25 / 13.49 | 3.24 / 10.26 / 13.49 | **0.85 / 10.51 / 11.36** |
| B3 | −0.54 / 11.50 / 10.97 | −0.54 / 11.51 / 10.97 | **−0.27 / 11.48 / 11.21** |
| YR, B1 | 5.25 / 6.27 / 11.53 | 5.26 / 6.27 / 11.53 | 5.26 / 6.27 / 11.53 |
| YR, B2 | 1.19 / 9.49 / 10.68 | 1.19 / 9.50 / 10.69 | **1.97 / 9.42 / 11.39** |
| YR, B3 | 1.89 / 9.53 / 11.42 | 1.89 / 9.53 / 11.42 | **1.54 / 9.57 / 11.11** |
| B1, B2 | 13.57 / 5.35 / 18.92 | 13.58 / 5.35 / 18.92 | **4.21 / 5.93 / 10.14** |
| B1, B3 | 6.45 / 8.14 / 14.59 | 6.45 / 8.14 / 14.59 | **3.09 / 8.43 / 11.53** |
| B2, B3 | 3.45 / 10.06 / 13.51 | 3.45 / 10.06 / 13.51 | **0.68 / 10.35 / 11.02** |
| YR, B1, B2 | 6.26 / 4.58 / 10.84 | 6.26 / 4.58 / 10.84 | 6.26 / 4.58 / 10.84 |
| YR, B1, B3 | 6.10 / 5.83 / 11.93 | 6.10 / 5.83 / 11.93 | **5.66 / 5.86 / 11.52** |
| YR, B2, B3 | 1.70 / 9.06 / 10.76 | 1.70 / 9.07 / 10.77 | **2.32 / 9.01 / 11.33** |
| B1, B2, B3 | 13.60 / 5.17 / 18.76 | 13.60 / 5.17 / 18.77 | **4.16 / 5.73 / 9.89** |
| YR, B1, B2, B3 | 6.75 / 4.18 / 10.92 | 6.75 / 4.18 / 10.93 | 6.75 / 4.18 / 10.93 |

(All values in %; corrected values that differ from the published ones in bold.)

**Consequences for the text.** With the corrected values the subset with the largest fixed share is the full model
(6.75%), followed by YR, B1, B2 (6.26%), which are also the two best models by AIC and BIC; the subsets with B1 and B2
but without YR no longer stand out (4.21% and 4.16% instead of 13.57% and 13.60%). The example of a negative
coefficient of determination still holds: B3 alone gives −0.27% (published −0.54%).

We attach the SAS program and its output (`piepho-table6.sas`, `piepho-table6-results.csv`). We would be glad to
be told if we have misread the method. Thank you for the paper and for making your SAS code available; it made this
check straightforward.

Kind regards,

Dongwen Luo
Bioeconomy Science Institute, New Zealand
