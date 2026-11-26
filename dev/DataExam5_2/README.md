# DataExam5.2: Finlay-Wilkinson regression (tree height, 6 sites x 37 seedlots)

Source: Williams, Harwood & Matheson (2024), *Experimental Design and Analysis for
Tree Improvement*, Table 5.2 (data), Tables 5.7 and 5.9 (results).
Model: `lm(ht ~ site + sitemean * seedlot)`; `sitemean` is one value per site,
so it is aliased with `site` (its main effect is NA).

## Bug fixed (2026-11-22)
`Kmatrix()` set `sitemean` at its overall mean in every grid row; no site has
that value, so every seedlot mean was NA. Now each row uses its own site's
value (`pm_level_covariates()`, `R/Kmatrix.R`, with a message).
CRAN 1.1.1 gave numbers, but they changed with the reference level of
`seedlot` (aliased coefficient set to 0), so they were not estimates.

## Checks done
| check | result |
|---|---|
| Book Table 5.7 (4 sites, 27 complete seedlots) | means, SE 43.82, df 52, slopes match |
| Means of seedlots seen at all 6 sites | = raw means (e.g. 13877 291.0) |
| Reference level of `seedlot` changed | means and SEs unchanged |
| `weights = "cells"` | = average of fitted values per seedlot |
| Book Table 5.9 (6 sites, raw `sitemean`) | complete seedlots match; slopes ~0.05-0.1 low |
| Table 5.9 with site means adjusted for seedlot | slopes within 0.15; means within ~5 except 13686, 14623 (2 sites; SE ~107) |

Tests: `tests/testthat/test-level-covariate.R` (the real-data part skips if
the csv is not at the author's local path).

## Open
0. Genstat 23.1 run (2026-11-22): it drops `sitemean` (aliased) and predicts
   at its overall mean, so its seedlot means are those of CRAN 1.1.1 and
   depend on the reference level (12013: Genstat 535.6, raw mean and
   predictmeans 529.0, book method 529). The book (13877 = 291, 13689 = 463)
   agrees with predictmeans, not with this Genstat setup. SAS not yet run.
1. Run `DataExam5_2.sas` and `DataExam5_2.gen`; compare with
   `predictmeans_DataExam5_2.csv` (equal and cell weights).
   Expected: SAS LSMEANS (equal) = predictmeans `weights = "equal"`;
   SAS may mark some as Non-est (it uses the overall covariate mean, like 1.1.1).
   Genstat `ADJUST=marginal`/`COMBINATIONS=present` ~ `weights = "cells"`.
2. Decide how the book got 13686 = 510 and 14623 = 592.
3. Ship DataExam5.2 in `data/` so the book test always runs.
4. Run the full test suite.
