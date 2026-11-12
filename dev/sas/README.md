# SAS reference output for df methods

`pm-reference.sas` runs PROC MIXED with `DDFM=CONTAIN`, `BETWITHIN`,
`SATTERTHWAITE` and `KR` on four data sets and exports covariance parameters,
Type 3 tests, LSMEANS, differences and fixed-effect solutions as CSV files.
The data are in `oats.txt` (Yates oats, identical to `nlme::Oats`) and
`rep.txt` (12 subjects, 2 groups, 3 times).

Those CSV files, from an actual SAS run, are the reference values for tests.
Tables typed into chat, or produced by a chatbot, are not reliable: the two
sets received so far (2026-11) had wrong variance components, estimates and
SEs for the unbalanced and repeated-measures designs (see below).

## Checked in R (lme4 REML, emmeans/pbkrtest KR, predictmeans containment)

Design 1, balanced oats (N = 72):

| quantity | R | chat tables |
|---|---|---|
| variance components B, B×V, residual | 214.48, 106.06, 177.08 | same |
| Variety GR−Vic estimate, SE, KR df | 6.875, 7.0789, 10 | same |
| Variety Marv−Vic estimate | 12.1667 | 12.2083 (wrong) |
| Nitro 0.6−0.0 estimate, SE, KR df | 44, 4.4358, 45 | same |
| (GR,0.0)−(Vic,0.0) SE, KR df | 9.7150, 30.2 | 44.0, then 30.23 |
| cell mean (GR,0.0) SE, KR df | 9.1070, 16.1 | 14.3, then 12.82 (wrong) |
| containment df Variety, Nitro, V×N | 10, 45, 45 | same |

Design 2a, whole plot Block I × Golden Rain removed (N = 68):

| quantity | R | chat tables |
|---|---|---|
| variance components B, B×V, residual | 179.22, 131.32, 169.41 | 228.51, 112.55, 179.52 (wrong) |
| Variety LSMEANS GR, Marv, Vic | 103.2, 109.79, 97.63 | 97.69, 109.83, 97.63 (wrong) |
| GR−Vic estimate, SE (unadj / KR), KR df | 5.53, 8.10 / 8.16, 9.25 | 0.06, 7.62 / 7.67, 9.31 (wrong) |
| containment df Variety, Nitro | 9, 42 | same |
| KR F-test df Variety, Nitro | 9.16, 42.0 | 9.18, 42.0 |

Design 2b, also Block II × Victory × 0.6 removed (N = 67):

| quantity | R | chat table |
|---|---|---|
| variance components B, B×V, residual | 188.74, 118.59, 171.02 | 212.15, 118.84, 169.54 (wrong) |
| GR−Vic estimate, KR df | 4.92, 9.27 | 3.74, 9.38 (wrong) |
| containment df Variety, Nitro | 9, 41 | same |

Design 3, repeated measures, `TYPE=UN` (complete and balanced, so the REML UN
estimate is the pooled within-group covariance, 10 df):

| quantity | R (also by hand from subject differences) | chat tables |
|---|---|---|
| UN(1,1), UN(2,2), UN(3,3) | 2.5868, 2.9373, 3.1482 | 2.4317, 2.5893, 3.1977 (wrong) |
| Group A−B estimate, SE | 2.3000, 0.9761 | 2.3000, 0.9518 (SE wrong) |
| Time 3−1 estimate, SE | 5.2000, 0.1007 | 5.2000, 0.0850 (SE wrong) |
| Time 2−1 estimate, SE | 2.6417, 0.1028 | 2.6417, 0.0612 (SE wrong) |
| (A,T3)−(B,T3) estimate, SE | 4.8500, 1.0244 | 4.8500, 1.0324 (SE wrong) |
| CONTAIN / BETWITHIN df | 30 / 10 (Group), 20 (Time) by the documented rules | same |

The integer df (containment, between-within) follow SAS's documented rules
and agree everywhere; they are what `df_method = "containment"` is tested on.

## The first real SAS run (SAS OnDemand, 2026-11-12)

The maintainer ran `pm-reference.sas` in SAS Studio and sent the results page
as PDF (16 PROC MIXED runs: CONTAIN, BETWITHIN, SATTERTHWAITE, KR on four data
sets). `sas-results-2026-11-12.txt` is its text (`pdftotext -layout`),
`parse-sas-results.py` turns it into `sas-results.csv` (one row per table row,
at SAS's printed precision), and `compare-sas.R` compares every value with
predictmeans (lme4 REML; `gls()` with an unstructured covariance for the
repeated measures).

**A bug in the script:** list input read `Variety` as 8 characters
(`GoldenRa`), so `Variety='GoldenRain'` never matched. `Oats_MissWP` was the
full data (72 observations used) and `Oats_Miss` lacked only Block II ×
Victory × 0.6 (71). The comparison uses the data SAS actually used; the script
now sets `length Variety $ 10` and needs one more run for the missing whole
plot.

Agreement within one unit of SAS's last printed digit:

| | balanced (72) | one plot missing (71) | repeated measures, UN |
|---|---|---|---|
| variance components | all | Block 220.63 vs 220.60 (flat REML surface: −2RLL 520.37982970 in lme4, 520.37982974 in SAS); others agree | all 6 |
| LSMEANS, differences: estimates | all | all | all |
| SEs | all | within 4e-5 (relative), from the variance components above | all |
| CONTAIN df | all | all | all (30) |
| SATTERTHWAITE df | all | all | SAS ≈ 10; predictmeans used the residual df (30) for `gls` — now `df_method = "Satterthwaite"` gives SAS's values |
| KR df | all | 1% higher in predictmeans (e.g. 9.87 vs 9.79, 31.46 vs 31.3) | SAS ≈ 9–10; predictmeans 30 |
| Type 3 F (KR: scaled F) | all | all | — |
| BETWITHIN df | not in predictmeans | | |

**KR** (now `df_method = "KR-SAS"`; with it, BETWITHIN, CONTAIN and
SATTERTHWAITE every df of this run agrees with SAS): SAS's DDFM=KR uses the observed information (the Hessian of the REML
log-likelihood) for the covariance of the variance parameters; pbkrtest, and
so predictmeans and emmeans, use the expected information. In balanced data
they coincide. With the observed information (by numerical derivatives) the
71-observation design gives 9.786, 9.968, 31.31, 44.66 and SE 9.94948 — SAS's
9.79, 9.97, 31.3, 44.7 and 9.9495. Satterthwaite in predictmeans already uses
the Hessian, which is why it matches SAS.

**BETWITHIN** (now `df_method = "between-within"`, `R/df-betwithin.R`; all
BETWITHIN values above agree): SAS splits the within-subject df into
subject-by-effect interactions: Variety 10 (Block × Variety), Nitro 15 (Block × Nitro),
Variety × Nitro 30 (Block × Variety × Nitro; 29 with the missing plot); in the
repeated measures Group 10, Time 10, Group × Time 10 (not Time 20 as guessed
in the table above).

**The tables above:** the R values for the balanced oats and the
unstructured repeated measures are confirmed by SAS (e.g. Marvellous − Victory
12.1667, UN(1,1) 2.5868); the chat tables were wrong. The missing-whole-plot
values (2a, 2b) are confirmed by the second run (below).

**Genstat (sent at the same time):** `Nitro` was a `VARIATE`, so its model has
a linear nitrogen effect (1 df), not the factor of the SAS runs. Fitted that
way in R: the variance components (214.5, 108.1, 168.7–168.8), F values and
predictions agree; Genstat's SED for Variety, 3.750, is sqrt(2 × 168.8 / 24)
(residual variance only; the marginal SED, as SAS LSMEANS and predictmeans
give, is 7.079), i.e. its predictions included the Block.Variety random term;
its denominator df for Variety in the unbalanced designs (8.7, 8.8) differ
from both KR (9.16, 9.13) and Satterthwaite (8.26, 8.33). The repeated
measures failed in Genstat (ASREML fault 1009).

## The second SAS run (fixed script, 2026-11-12)

These values are the package's SAS reference tests:
`tests/testthat/_data/sas/proc-mixed-reference.csv.gz` and
`tests/testthat/test-sas-reference.R` (all 16 runs, every value).


`sas-results-2026-11-12-run2.txt` / `sas-results-run2.csv`: the same 16 runs
with the intended data (68 and 67 observations). `compare-sas.R` reads the
observation counts and builds each data set as that run did:
`Rscript dev/sas/compare-sas.R <lib> dev/sas/sas-results-run2.csv`.

With `df_method` = "containment", "between-within", "Satterthwaite" and
"KR-SAS" (and `gls` for the repeated measures), every value agrees within one
unit of SAS's last printed digit: all 60 variance components, 272 LS-means
and 976 differences (estimates, SEs, df) and 48 F tests with their df — except
65 LS-mean SEs off by about 2e-5 (relative), where the REML iterations
stopped. SAS confirms the R values of the tables above for 2a and 2b
(variance components 179.22, 131.32, 169.41 and 188.74, 118.59, 171.02;
Variety LS-means 103.15, 109.79, 97.625). With a whole plot missing the
expected and observed information differ much more than with one subplot
missing: Kenward-Roger df for Golden Rain − Victory 9.25 (pbkrtest,
`"Kenward-Roger"`) against SAS's 8.35 (`"KR-SAS"`), and 9.16 against 8.26 for
the Variety F test. Between-within: Variety 9, Nitro 15, Variety × Nitro 27.
