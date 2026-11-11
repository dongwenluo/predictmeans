## Confirming run for the order-invariant whitening (plan step 2b), with the
## package's own permlmer() and permmodels() (install the new version first).
## Uses dev/sim-perm-lib.R. Results: dev/sim-perm-confirm/.
##
##   cf_permlmer_sp       permlmer(), split plot, test of the whole-plot variance
##                        (as lmm_sp_norm; Cholesky: 0.059, 0.056 without ties)
##   cf_permlmer_crossed  permlmer(), rows x columns (as re_crossed; Cholesky 0.064 / 0.060)
##   cf_pm_lmer_global    permmodels() on lmer, y ~ x + A + (1 | block), global null,
##                        permutation p of A (maintainer's scheme; Cholesky 0.062-0.071)
##   cf_pm_lme_global     the same with nlme::lme (the lme/gls code path)
## Reference columns: chibar (random effects) or the model's own F test p-value.
##
## Run: Rscript dev/sim-perm-confirm.R [R] [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
nperm <- if (length(args) >= 2) as.integer(args[2]) else 99
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")

pm_global <- function(d, cls) {
  if (cls == "lmer") {
    fm <- lmer(y ~ x + A + (1 | block), data = d, control = ctl)
  } else {
    fm <- nlme::lme(y ~ x + A, random = ~ 1 | block, data = d)
  }
  r <- permmodels(fm, nperm = nperm, prt = FALSE, ncore = 1)
  aT <- r$ANOVA
  pcol <- grep("^Pr|p-value", colnames(aT), value = TRUE)[1]
  c(perm = as.numeric(aT["A", "Perm_p_value"]), F_test = as.numeric(aT["A", pcol]))
}

description <- list(
  cf_permlmer_sp = "permlmer(), split plot (6 blocks x 4 whole plots x 3), normal; whole-plot variance 0",
  cf_permlmer_crossed = "permlmer(), 8 x 8 rows x columns, A (4 levels); H0 (1|row), H1 + (1|col); column variance 0",
  cf_pm_lmer_global = "permmodels() on lmer, 8 blocks x 3 A x 2, normal x, no fixed effects; p of A",
  cf_pm_lme_global = "permmodels() on lme, as cf_pm_lmer_global")

one_set <- function(cell) {
  switch(cell,
    cf_permlmer_sp = {
      d <- sp_design()
      d$y <- 0.5 * as.numeric(d$A) + 0.3 * as.numeric(d$B) + rnorm(6)[d$block] + rnorm(nrow(d))
      lmm_test(d, y ~ A * B + (1 | block), y ~ A * B + (1 | block) + (1 | wpid))
    },
    cf_permlmer_crossed = {
      d <- expand.grid(row = factor(1:8), col = factor(1:8))
      d$A <- factor(sample(rep(1:4, 16)))
      d$y <- 0.3 * as.numeric(d$A) + rnorm(8)[d$row] + rnorm(nrow(d))
      lmm_test(d, y ~ A + (1 | row), y ~ A + (1 | row) + (1 | col))
    },
    cf_pm_lmer_global = , cf_pm_lme_global = {
      d <- lev_design(8); d$x <- rnorm(nrow(d))
      d$y <- rnorm(8)[d$block] + rnorm(nrow(d))
      pm_global(d, if (cell == "cf_pm_lmer_global") "lmer" else "lme")
    },
    stop("unknown cell ", cell))
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-perm-confirm", description)
cat("SIM DONE\n")
