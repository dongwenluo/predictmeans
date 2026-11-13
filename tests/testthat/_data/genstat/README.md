# Genstat REML reference values

`genstat-reml.csv.gz`: the values printed by Genstat 23.1 for
`dev/genstat/pm-reference.gen` and `dev/genstat/pm-reference-2.gen` (in the
package sources), run by the maintainer on 2026-11-14; the Output windows,
saved as RTF, as text in `dev/genstat/genstat-results-2026-11-14-run*.txt`,
parsed by `dev/genstat/parse-genstat.py`. Columns: analysis (run and response
variate), table (vc, wald_seq, wald_drop, mean, sed, sed_same, vp_pred, vp_se,
vp_sed, aov_*), term, key, value (as printed).
