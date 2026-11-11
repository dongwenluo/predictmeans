## Edge cases for the permutation tests of mixed models: where are the schemes
## expected to struggle? Uses the functions in dev/sim-perm-lib.R. Type I error
## at 5% unless the cell name ends in "_power". Results: dev/sim-perm-edge/.
## Research log and summary: dev/sim-study/.
##
## Run (from the package root):
##   Rscript dev/sim-perm-edge.R [R] [nperm] [cells|all] [cores]
##   e.g. Rscript dev/sim-perm-edge.R 1000 99 all 4
## Cells already done are skipped, so an interrupted run can be restarted.
## On Windows, cores > 1 uses a socket cluster.

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
nperm <- if (length(args) >= 2) as.integer(args[2]) else 99
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")

description <- list(
  ## random effects, LMM (permlmer, Lee-Braun)
  re_few_blocks = "LMM split plot, 3 blocks (n = 36); test whole-plot variance; normal",
  re_heavy = "LMM split plot, 6 blocks; t3 errors and t3 block effects; test whole-plot variance",
  re_hetero = "LMM split plot; error SD 0.5, 1, 2 by subplot factor B (model assumes equal); test whole-plot variance",
  re_singular = "LMM split plot; block variance 0 (H0 fit often singular); test whole-plot variance",
  re_tiny = "LMM, 6 blocks x 6 plots of 1-4 units (unequal); test plot variance",
  re_crossed = "LMM, 8 x 8 row-column field trial, A (4 levels); H0 (1|row), H1 + (1|col); column variance 0",
  re_slope_corr = "LMM repeated measures, 12 subjects x 6 times; H0 (1|subj), H1 (1+time|subj) (2 parameters, 50:50 chi2(1):chi2(2))",
  re_few_blocks_power = "as re_few_blocks, whole-plot SD 0.7 (power)",
  ## fixed effects, LMM (permmodels form)
  fx_few_blocks = "LMM, 3 blocks x 3 A levels x 4 reps, covariate with extreme values related to A; test A",
  fx_hetero = "LMM, 6 blocks; A levels with 1, 2, 3 units per block and error SD 3, 1, 1 (Behrens-Fisher-like); test A",
  fx_missing = "LMM, 6 blocks, A x B (3 x 3) additive with two empty A:B combinations, large B effect; test A",
  fx_global = "LMM, 8 blocks, y ~ x + A, all fixed effects 0 (global null); test A",
  fx_assoc_skew = "LMM, 8 blocks, A tested with B associated and large; skewed errors (1000 data sets)",
  fx_lev_power = "LMM, 8 blocks, covariate with extreme values; A effect 0, 0.5, 1 (power)",
  ## random effects, GLMM (working scale)
  g_re_bin_pairs = "binary GLMM split plot, whole plots of 2 units; test whole-plot variance",
  g_re_rare = "binary GLMM split plot, about 5% events; test whole-plot variance",
  g_re_pois_low = "Poisson GLMM split plot, mean about 0.4 (mostly zeros); test whole-plot variance",
  g_re_overdisp = "counts negative binomial (size 2) fitted as Poisson GLMM; test whole-plot variance (truth 0)",
  g_re_few_blocks = "Poisson GLMM split plot, 3 blocks; test whole-plot variance",
  ## fixed effects, GLMM (working scale)
  g_fx_overdisp = "counts negative binomial (size 2) fitted as Poisson GLMM, covariate with extreme values; test A",
  g_fx_rare = "binary GLMM, about 5-10% events, A tested with B associated; test A",
  g_fx_few_blocks = "binomial (10 trials) GLMM, 3 blocks, covariate with extreme values; test A",
  g_fx_global = "Poisson GLMM, y ~ x + A, all fixed effects 0 (global null); test A")

one_set <- function(cell) {
  switch(cell,
    re_few_blocks = , re_few_blocks_power = {
      d <- sp_design(3, 3)
      sw <- if (cell == "re_few_blocks_power") 0.7 else 0
      d$y <- 0.5 * as.numeric(d$A) + rnorm(3)[d$block] + sw * rnorm(12)[d$wpid] + rnorm(nrow(d))
      lmm_test(d, y ~ A * B + (1 | block), y ~ A * B + (1 | block) + (1 | wpid))
    },
    re_heavy = , re_hetero = , re_singular = {
      d <- sp_design()
      eb <- switch(cell, re_heavy = heavy(6), re_singular = rep(0, 6), rnorm(6))
      ee <- switch(cell, re_heavy = heavy(nrow(d)), re_hetero = c(0.5, 1, 2)[d$B] * rnorm(nrow(d)), rnorm(nrow(d)))
      d$y <- 0.5 * as.numeric(d$A) + 0.3 * as.numeric(d$B) + eb[d$block] + ee
      lmm_test(d, y ~ A * B + (1 | block), y ~ A * B + (1 | block) + (1 | wpid))
    },
    re_tiny = {
      sizes <- sample(1:4, 36, replace = TRUE)
      d <- data.frame(plot = factor(rep(seq_along(sizes), sizes)))
      d$block <- factor((as.integer(d$plot) - 1) %/% 6 + 1)
      d$y <- rnorm(6)[d$block] + rnorm(nrow(d))
      lmm_test(d, y ~ 1 + (1 | block), y ~ 1 + (1 | block) + (1 | plot))
    },
    re_crossed = {
      d <- expand.grid(row = factor(1:8), col = factor(1:8))
      d$A <- factor(sample(rep(1:4, 16)))
      d$y <- 0.3 * as.numeric(d$A) + rnorm(8)[d$row] + rnorm(nrow(d))
      lmm_test(d, y ~ A + (1 | row), y ~ A + (1 | row) + (1 | col))
    },
    re_slope_corr = {
      d <- expand.grid(time = 0:5, subj = factor(1:12))
      d$y <- 1 + 0.5 * d$time + rnorm(12)[d$subj] + rnorm(nrow(d))
      lmm_test(d, y ~ time + (1 | subj), y ~ time + (1 + time | subj), k = 2)
    },
    fx_few_blocks = , fx_lev_power = {
      d <- lev_design(if (cell == "fx_few_blocks") 3 else 8)
      if (cell == "fx_few_blocks") d <- d[rep(seq_len(nrow(d)), each = 2), ]
      eff <- if (cell == "fx_lev_power") c(0, 0.5, 1)[d$A] else 0
      nb <- nlevels(d$block)
      d$y <- 3 * (d$x - mean(d$x)) + eff + rnorm(nb)[d$block] + rnorm(nrow(d))
      lmm_fixed_F(d, y ~ x + (1 | block), y ~ x + A + (1 | block), "A")
    },
    fx_hetero = {
      d <- do.call(rbind, lapply(1:6, function(b) data.frame(block = b, A = rep(1:3, 1:3))))
      d$block <- factor(d$block); d$A <- factor(d$A)
      d$y <- rnorm(6)[d$block] + c(3, 1, 1)[d$A] * rnorm(nrow(d))
      lmm_fixed_F(d, y ~ 1 + (1 | block), y ~ A + (1 | block), "A")
    },
    fx_missing = {
      d <- expand.grid(rep = 1:2, A = factor(1:3), B = factor(1:3), block = factor(1:6))
      d <- d[!((d$A == "3" & d$B == "1") | (d$A == "1" & d$B == "3")), ]
      d$y <- c(0, 2, 4)[d$B] + rnorm(6)[d$block] + rnorm(nrow(d))
      lmm_fixed_F(d, y ~ B + (1 | block), y ~ B + A + (1 | block), "A")
    },
    fx_global = {
      d <- lev_design(8); d$x <- rnorm(nrow(d))
      d$y <- rnorm(8)[d$block] + rnorm(nrow(d))
      lmm_fixed_F(d, y ~ x + (1 | block), y ~ x + A + (1 | block), "A")
    },
    fx_assoc_skew = {
      d <- assoc_design(8)
      d$y <- c(0, 3, 6)[d$B] + rnorm(8)[d$block] + skew(nrow(d))
      lmm_fixed_F(d, y ~ B + (1 | block), y ~ B + A + (1 | block), "A")
    },
    g_re_bin_pairs = , g_re_rare = {
      d <- sp_design(if (cell == "g_re_bin_pairs") 2 else 5)
      a <- if (cell == "g_re_rare") -3.2 else -0.3
      d$y <- rbinom(nrow(d), 1, plogis(a + 0.3 * (as.numeric(d$A) - 2.5) + 0.7 * rnorm(6)[d$block]))
      glmm_test(d, y ~ A + B + (1 | block), y ~ A + B + (1 | block) + (1 | wpid), binomial)
    },
    g_re_pois_low = , g_re_overdisp = , g_re_few_blocks = {
      nb <- if (cell == "g_re_few_blocks") 3 else 6
      d <- sp_design(3, nb)
      a <- if (cell == "g_re_pois_low") -1 else 1
      mu <- exp(a + 0.3 * (as.numeric(d$A) - 2.5) + 0.5 * rnorm(nb)[d$block])
      d$y <- if (cell == "g_re_overdisp") rnbinom(nrow(d), size = 2, mu = mu) else rpois(nrow(d), mu)
      glmm_test(d, y ~ A + B + (1 | block), y ~ A + B + (1 | block) + (1 | wpid), poisson)
    },
    g_fx_overdisp = , g_fx_global = {
      d <- lev_design(8)
      if (cell == "g_fx_global") d$x <- rnorm(nrow(d))
      xc <- d$x - mean(d$x)
      bl <- 0.5 * rnorm(8)[d$block]
      d$y <- if (cell == "g_fx_overdisp") rnbinom(nrow(d), size = 2, mu = exp(1 + 0.15 * xc + bl)) else
        rpois(nrow(d), exp(1 + bl))
      glmm_fixed(d, y ~ x + (1 | block), y ~ x + A + (1 | block), poisson, "A")
    },
    g_fx_rare = {
      d <- assoc_design(8, reps = 3)
      d$y <- rbinom(nrow(d), 1, plogis(c(-3.5, -2.8, -2.1)[d$B] + 0.5 * rnorm(8)[d$block]))
      glmm_fixed(d, y ~ B + (1 | block), y ~ B + A + (1 | block), binomial, "A")
    },
    g_fx_few_blocks = {
      d <- lev_design(3); d <- d[rep(seq_len(nrow(d)), each = 2), ]
      xc <- d$x - mean(d$x)
      d$n <- 10; d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.3 * xc + 0.5 * rnorm(3)[d$block]))
      glmm_fixed(d, cbind(r, n - r) ~ x + (1 | block), cbind(r, n - r) ~ x + A + (1 | block), binomial, "A")
    },
    stop("unknown cell ", cell))
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-perm-edge", description)
cat("SIM DONE\n")
