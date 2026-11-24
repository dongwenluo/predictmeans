# Speed / clean-up review: work list

The rules: one item per step. Results must stay identical, at tolerance 0 or 1e-12 where only the floating-point order changes. Each change is timed, and recorded in NEWS.md and plan row O16.

| # | Where | Change | Status |
|---|---|---|---|
| 1 | `R/R2_glmm.R` (Nakagawa) | `as.matrix(Z)` and a dense `kronecker(diag(nl), vcov)` became sparse `Matrix::kronecker(Diagonal())`. | **Done**: n = 10,000, 2,000 groups with slopes, 136 s -> 0.02 s |
| 2 | `R/satt-nlme.R:213` | `Sinv` was checked: it stays sparse, with 10,000 non-zeros against n² = 4e6 for AR(1) within groups. | Checked; no change |
| 3 | `R/varcomp.R:732`, `p-adjust.R:212`, `satt-nlme.R:304` | `X %*% diag(v) %*% t(X)` becomes a column scaling. | **Done**: `tcrossprod(X * rep(sqrt(v), each = nrow(X)))` (all v >= 0 there); differences <= 3e-15; matrices are small (k×k), so this is clean-up rather than speed |
| 4 | `R/semireg-smooth.R:365` | Dense `kronecker(diag(nl), B)` becomes sparse. | **Done** (setup only, once per smooth) |
| 5 | `R/perm-engine.R:223,385`, `perm-robust.R:28` | `solve(crossprod(X))` is hoisted out of loops, or replaced by `chol2inv(chol())`. | **Done**: they were already outside the loops. One helper, `pm_xtx_inv(qrX)` (from the existing QR, so the condition number is not squared; handles pivoting), is used in the 3 places; differences ~1e-17 |
| 6 | `covariatemeans.R:276`, `semipred.R:337`, `predictmeans-stages.R:252` | `diag(K V K')` and the per-row `apply()` become `rowSums((K %*% V) * K)`. | **Done**: 3,000 means, 0.02 s -> < 0.005 s; no m×m matrix; differences ~3e-14 |
| 7 | Wald quadratic forms (`anova.R`, `inference.R`, `kr.R`, `perm-engine.R`) | Moved into one helper. | **Done**: `pm_wald_quad(L, b, V, Lb)` in `perm-engine.R`; `pm_waldF()`, `pm_wald_chisq()`, `pm_type1_F()`, `pm_inference()$ftest`, `pm_kr_ftest()` use it; no explicit inverse (the old `anova.R` form inverted); `kr.R` Theta via `crossprod`. 28 test files pass |
| 8 | `R/matrix-functions.R` / `build-variance-matrices.R` | Dense block building becomes `Matrix::bdiag()`; unused helpers are removed. | **Done**: in the non-nested case, the dense n×n `unblock()` sum is replaced by sparse `pm_place_blocks()` and split into connected blocks (`pm_sparse_components()`), so downstream whitening works per block rather than on one n×n matrix. Removed: `unblock`, `sum_blockblock`, `matrix_minus_block`, `block_minus_matrix`, `add_diag`, `add_diag_bdiag`, `sub_f` (folded into `matrix_list`). Tests: `test-matrix-functions.R` (exact vs dense build) |
| 9 | Duplicated permutation setup (`perm-fl*.R`, `permmodels.R`, `permlmer.R`) | One shared helper. | **Done**: `pm_mer_whiten()`, `pm_mer_perm_parts()` and `pm_nlme_whiten()` in `R/whiten.R`; used by `pm_fl_design()`, `pm_lmm_parts()`, `permmodels()` (raw) and `permlmer()`. Full suite: 3,044 expectations, 0 failures |

Keep the existing style: `Matrix::` sparse for anything n×n or n×q; only small p×p or q×q matrices are dense; `solve(A, B)` and `crossprod` throughout.
