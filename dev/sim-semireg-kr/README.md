# Kenward-Roger for pairs of semireg means (2026-11-26)

`sim.R`: coverage of the 95% limits of all pairwise differences (true value 0)
of semireg means, Kenward-Roger adjusted SEs and df (the default since
b05cb33) against unadjusted SEs with Satterthwaite df (`vcov_adjust = FALSE`).

| design | data sets | KR | unadjusted |
|---|---|---|---|
| A: Slate Hall layout, 25 varieties x 6 reps, `t2(row, col)` trend + row and column effects | 200 | 94.9% | 94.5% |
| B: 4 groups x 6 blocks, `s(x)` trend + block effects | 500 | 94.9% | 94.2% |

Both are close to 95%; KR is the closer, more so in the small design, as for
lmer fits. The default stays: Kenward-Roger for the pairs.
