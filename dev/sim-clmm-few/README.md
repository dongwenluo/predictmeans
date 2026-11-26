# Tests of a between-cluster treatment in a clmm with few clusters (2026-11-26, B4 G3)

`sim.R`: `clmm(y ~ trt * W + (1 | clus))`, clusters split half/half between two
treatments, 6 rows per cluster (W within, 2 x 3), 4 classes, cluster SD 1, W
effect 0.5. Rejection of the trt test at 5%, 300 data sets per row, 296-298
usable (a true 5% lies in 2.5-7.5%). The t and z columns are the Type II Wald
chi-square of trt, with the inner-outer df of trt (nclus - 2, as predictmeans
gives it) or none; not predictmeans' pair contrast, which averages over W. Run on the code of b59ded5 (the allocations of this design
are those of the final pm_perm_clmm(); the near-tie rule was relative 1e-6,
now 1e-4).

| clusters | trt effect | inner-outer t (default) | z | bootstrap (99) | cluster permutation (99) |
|---|---|---|---|---|---|
| 6 | 0 | 3.4% | 10.5% | 4.7% | 0% (20 allocations: p >= 0.1) |
| 8 | 0 | 4.4% | 10.8% | 5.1% | 3.4% (70, exact) |
| 12 | 0 | 5.4% | 8.7% | 3.7% | 2.7% (924; 99 random) |
| 8 | 1.5 | 34% | 47% | 31% | 11% |

- z tests reject 9-11% of true null hypotheses: too liberal with few clusters.
- The inner-outer df (predictmeans' default for clmm) hold 3.4-5.4%, and so does
  the parametric bootstrap (3.7-5.1%), with similar power.
- The cluster permutation is valid but conservative and weak with 12 clusters
  or fewer (its p-values are coarse: 6 clusters cannot reject at 5%); it is
  the test that needs no model for the clusters (power at 12 clusters was not
  simulated).
