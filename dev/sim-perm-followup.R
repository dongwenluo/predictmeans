## Follow-ups to the permutation study (dev/permutation-plan.md, step 1a).
## Uses dev/sim-perm-lib.R. Results: dev/sim-perm-followup/. Log: dev/sim-study/.
##
## 1. Why is the maintainer's fixed-effect scheme liberal under the LMM global
##    null (7.1% in dev/sim-perm-edge/fx_global)? Same design, more data sets,
##    and variants of the scheme:
##      joint_full   y* = xb1[pi] + U1 e1[pi], full model (as permmodels())
##      joint_red    y* = xb0[pi] + U0 e0[pi], reduced model
##      joint_oracle as joint_full but whitened with the true V (no estimation)
##      fl           Freedman-Lane: y* = xb0 + U0 e0[pi]
##      within       y permuted within blocks (exact under the global null here)
##      satt         Satterthwaite F
##    Statistic: F of A from refits of the full model.
##    fu_lmm_global (8 blocks, n = 48) and fu_lmm_global16 (16 blocks, n = 96).
## 2. Binary GLMMs where the working scale failed (rare events; whole plots of
##    2): test of the whole-plot variance by
##      routeB  quantile back-transform: randomised PIT residuals q of the H0
##              fit (conditional on its BLUPs) permuted, mapped back to 0/1 with
##              the H0 conditional means, both GLMMs refitted, LRT
##      boot    parametric bootstrap from the H0 fit (simulate(), refit), LRT
##      (routeB and boot: all fits with nAGQ = 0, observed and resampled alike)
##      glmm_chibar  Laplace LRT with the 50:50 chi-bar reference
##
## Run: Rscript dev/sim-perm-followup.R [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
nperm <- if (length(args) >= 1) as.integer(args[1]) else 99
cores <- if (length(args) >= 3) as.integer(args[3]) else 1
source("dev/sim-perm-lib.R")

## ---- 1. LMM global null --------------------------------------------------------
lmm_global_diag <- function(d, f0, f1, term, sb = 1, se = 1) {
  m1 <- lme4::lmer(f1, data = d, REML = TRUE, control = ctl)
  m0 <- lme4::lmer(f0, data = d, REML = TRUE, control = ctl)
  Fs <- function(y) stats::anova(refit(m1, y))[term, "F value"]
  obs <- stats::anova(m1)[term, "F value"]
  n <- nrow(d)
  one <- function(m) {
    U <- lower_chol(getME(m, "Z") %*% getME(m, "Lambda"), rep(1, n))
    xb <- as.vector(getME(m, "X") %*% fixef(m))
    list(U = U, xb = xb, e = as.vector(Matrix::solve(U, d$y - xb)))
  }
  F1 <- one(m1); F0 <- one(m0)
  Zb <- getME(m1, "Z")                                    # block indicator matrix
  Uo <- lower_chol(Zb * (sb / se), rep(1, n))             # true V / se^2 = (sb/se)^2 ZZ' + I
  eo <- as.vector(Matrix::solve(Uo, d$y - F1$xb))
  blocks <- split(seq_len(n), d$block)
  s <- matrix(NA_real_, nperm, 5, dimnames = list(NULL, c("joint_full", "joint_red", "joint_oracle", "fl", "within")))
  for (k in seq_len(nperm)) {
    i <- sample.int(n)
    w <- seq_len(n); for (b in blocks) w[b] <- b[sample.int(length(b))]
    s[k, ] <- c(
      tryCatch(Fs(F1$xb[i] + as.vector(F1$U %*% F1$e[i])), error = function(e) NA),
      tryCatch(Fs(F0$xb[i] + as.vector(F0$U %*% F0$e[i])), error = function(e) NA),
      tryCatch(Fs(F1$xb[i] + as.vector(Uo %*% eo[i])), error = function(e) NA),
      tryCatch(Fs(F0$xb + as.vector(F0$U %*% F0$e[i])), error = function(e) NA),
      tryCatch(Fs(d$y[w]), error = function(e) NA))
  }
  satt <- stats::anova(lmerTest::as_lmerModLmerTest(m1))[term, "Pr(>F)"]
  c(apply(s, 2, pval_perm, obs = obs), satt = satt)
}

## ---- 2. binary GLMM: route B and parametric bootstrap ------------------------------
pit_q <- function(y, mu) {
  lo <- stats::pbinom(y - 1, 1, mu); hi <- stats::pbinom(y, 1, mu)
  u <- lo + stats::runif(length(y)) * (hi - lo)
  stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}
## Resampled and observed statistics use nAGQ = 0 fits (about 150 times
## faster than Laplace refits on rare-event data, which also often fail);
## the Laplace LRT is kept only for the chi-bar reference.
fctl <- glmerControl(optimizer = "bobyqa", calc.derivs = FALSE, check.conv.singular = "ignore")
## LRTs are rounded to 4 decimals (values below 1e-4 are optimizer noise at
## the boundary, i.e. 0), so ties at 0 are not counted as "more extreme".
glmm_re_routeB <- function(d, f0, f1) {
  lr <- function(m1, m0) round(max(2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0))), 0), 4)
  lap <- tryCatch(lr(glmer(f1, data = d, family = binomial, control = fctl),
                     glmer(f0, data = d, family = binomial, control = fctl)), error = function(e) NA)
  g0 <- glmer(f0, data = d, family = binomial, nAGQ = 0, control = fctl)
  g1 <- glmer(f1, data = d, family = binomial, nAGQ = 0, control = fctl)
  obs <- lr(g1, g0)
  lr_y <- function(y) lr(refit(g1, y), refit(g0, y))
  mu0 <- stats::fitted(g0)                                 # conditional on the H0 BLUPs
  q <- pit_q(d$y, mu0)
  sB <- vapply(seq_len(nperm), function(k) {
    ys <- stats::qbinom(stats::pnorm(q[sample.int(nrow(d))]), 1, mu0)
    tryCatch(lr_y(ys), error = function(e) NA)
  }, numeric(1))
  ysim <- stats::simulate(g0, nsim = nperm)
  sP <- vapply(ysim, function(ys) tryCatch(lr_y(ys), error = function(e) NA), numeric(1))
  c(routeB = pval_perm(sB, obs), boot = pval_perm(sP, obs), glmm_chibar = chibar(lap))
}

description <- list(
  fu_lmm_global = "LMM global null as fx_global (8 blocks x 3 A x 2, normal x), variants of the joint scheme; 2000 data sets",
  fu_lmm_global16 = "as fu_lmm_global with 16 blocks (n = 96); 1000 data sets",
  fu_bin_rare = "binary GLMM split plot, about 5% events (as g_re_rare): route B, parametric bootstrap; 500 data sets",
  fu_bin_pairs = "binary GLMM split plot, whole plots of 2 (as g_re_bin_pairs): route B, parametric bootstrap; 500 data sets")
R_of <- c(fu_lmm_global = 2000, fu_lmm_global16 = 1000, fu_bin_rare = 500, fu_bin_pairs = 500)

one_set <- function(cell) {
  switch(cell,
    fu_lmm_global = , fu_lmm_global16 = {
      d <- lev_design(if (cell == "fu_lmm_global") 8 else 16); d$x <- rnorm(nrow(d))
      d$y <- rnorm(nlevels(d$block))[d$block] + rnorm(nrow(d))
      lmm_global_diag(d, y ~ x + (1 | block), y ~ x + A + (1 | block), "A")
    },
    fu_bin_rare = , fu_bin_pairs = {
      d <- sp_design(if (cell == "fu_bin_pairs") 2 else 5)
      a <- if (cell == "fu_bin_rare") -3.2 else -0.3
      d$y <- rbinom(nrow(d), 1, plogis(a + 0.3 * (as.numeric(d$A) - 2.5) + 0.7 * rnorm(6)[d$block]))
      glmm_re_routeB(d, y ~ A + B + (1 | block), y ~ A + B + (1 | block) + (1 | wpid))
    },
    stop("unknown cell ", cell))
}

cells <- if (length(args) >= 2 && args[2] != "all") strsplit(args[2], ",")[[1]] else names(description)
for (cell in cells) run_cells(cell, one_set, R_of[[cell]], nperm, cores, "dev/sim-perm-followup", description)
cat("SIM DONE\n")
