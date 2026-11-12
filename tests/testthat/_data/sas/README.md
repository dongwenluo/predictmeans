# SAS PROC MIXED reference values

`proc-mixed-reference.csv.gz`: the tables printed by SAS PROC MIXED for
`dev/sas/pm-reference.sas` (in the package sources), run by the maintainer in
SAS OnDemand for Academics (SAS Studio), 2026-11-12 ("Results: Program 2",
the run with the fixed script). The results page was saved as PDF, its text
extracted with `pdftotext -layout` (`dev/sas/sas-results-2026-11-12-run2.txt`)
and parsed by `dev/sas/parse-sas-results.py` (`dev/sas/sas-results-run2.csv`);
this file keeps the columns the tests use. Values have the precision SAS
printed.

16 runs: four data sets (Yates' oats split plot, `nlme::Oats`; the same
without Block I x Golden Rain, 68 observations; also without Block II x
Victory x 0.6, 67; a repeated-measures trial, 12 subjects x 3 times, with an
unstructured covariance) x `DDFM=CONTAIN`, `BETWITHIN`, `SATTERTHWAITE`, `KR`.
Tables: covariance parameters (`covparms`), number of observations (`nobs`),
fixed-effect solutions (`solution`), Type 3 tests (`tests3`), LS-means
(`lsmeans`) and their differences (`diffs`).
