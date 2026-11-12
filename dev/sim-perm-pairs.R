## Study 1b, part 1: pairwise comparisons by permutation (plan §10, step 1b).
## Uses dev/sim-perm-lib.R. Results: dev/sim-perm-pairs/.
##
## For each data set, the p-value of every pairwise comparison of the means
## of one factor, by four methods:
##   joint   the maintainer's scheme as in permmodels(): full-model fitted
##           values and decorrelated residuals permuted together, model
##           refitted, |t| of every pair from the refit (one permlist for all
##           pairs, as predictmeans(permlist = ) uses it);
##   flpair  Freedman-Lane for each pair: the model with that pair's two means
##           equal (linear constraint), residuals decorrelated with the full
##           model's V (fixed), permuted, added back; |t| of the pair by GLS
##           with V fixed (batched, no refits);
##   flterm  the same with all means of the factor equal (one reduced model,
##           one set of permutations for all pairs, like a permlist);
##   param   predictmeans() t tests (Kenward-Roger for lmer; nlme's df for
##           lme/gls) and Tukey adjustment.
## Familywise error over the true null pairs: Holm for every method, max-|t|
## (single step) for joint and flterm, Tukey for param. Each is saved as the
## smallest adjusted p-value among the true null pairs (a family error when
## it is <= alpha).
##
## Run: Rscript dev/sim-perm-pairs.R [R] [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
nperm <- if (length(args) >= 2) as.integer(args[2]) else 199
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")
suppressPackageStartupMessages(library(nlme))
pm_whiten <- predictmeans:::pm_whiten
pm_whiten_blocks <- predictmeans:::pm_whiten_blocks
pm_perm_p <- predictmeans:::pm_perm_p
pm_perm_tstats <- predictmeans:::pm_perm_tstats

## rows: differences of the means of 'term' (averaged with equal weights over
## the other factors, covariates at their means), pairs (i, j) with i < j
pair_K <- function(m, d, term) {
  tt <- stats::delete.response(stats::terms(m))
  vars <- all.vars(tt)
  grid <- lapply(vars, function(v) if (is.factor(d[[v]])) factor(levels(d[[v]]), levels = levels(d[[v]])) else mean(d[[v]]))
  grid <- expand.grid(stats::setNames(grid, vars))
  X <- stats::model.matrix(tt, grid)
  lev <- levels(d[[term]])
  M <- rowsum(X, grid[[term]])[lev, , drop = FALSE] / as.vector(table(grid[[term]])[lev])
  pr <- utils::combn(length(lev), 2)
  K <- M[pr[2, ], , drop = FALSE] - M[pr[1, ], , drop = FALSE]
  rownames(K) <- paste0(lev[pr[1, ]], "v", lev[pr[2, ]])
  K
}

## basis of the null space of the rows of K (columns of N: K N = 0)
null_basis <- function(K) {
  q <- qr(t(K))
  Q <- qr.Q(q, complete = TRUE)
  Q[, -seq_len(q$rank), drop = FALSE]
}

## |t| of the rows of K by OLS on the whitened data, for each column of Y
ols_absT <- function(qrX, XtXi, K, Y) {
  n <- nrow(Y); p <- qrX$rank
  B <- qr.coef(qrX, Y)
  s2 <- colSums(qr.resid(qrX, Y)^2) / (n - p)
  se1 <- sqrt(rowSums((K %*% XtXi) * K))
  abs(K %*% B) / outer(se1, sqrt(s2))
}

## Freedman-Lane with V fixed: per pair, and for all pairs at once
fl_pairs <- function(Xw, yw, K, idx) {
  n <- nrow(Xw)
  qrX <- qr(Xw); XtXi <- chol2inv(qr.R(qrX))
  obs <- as.vector(ols_absT(qrX, XtXi, K, matrix(yw)))
  perm_from <- function(N) {
    X0 <- Xw %*% N
    f0 <- qr.fitted(qr(X0), yw); e0 <- yw - f0
    f0 + matrix(e0[idx], n)
  }
  p_pair <- vapply(seq_len(nrow(K)), function(j) {
    Tj <- ols_absT(qrX, XtXi, K[j, , drop = FALSE], perm_from(null_basis(K[j, , drop = FALSE])))
    pm_perm_p(obs[j], Tj)
  }, numeric(1))
  Tall <- ols_absT(qrX, XtXi, K, perm_from(null_basis(K)))
  list(pair = p_pair, term = pm_perm_p(obs, Tall), term_max = pm_perm_p(obs, matrix(rep(apply(Tall, 2, max), each = length(obs)), length(obs))))
}

## whitened X and y with V fixed at the fitted model's
whitened <- function(m, d) {
  X <- stats::model.matrix(stats::delete.response(stats::terms(m)), d)
  if (inherits(m, "merMod")) {
    W <- pm_whiten(getME(m, "Z") %*% getME(m, "Lambda"))
    y <- getME(m, "y")
  } else {
    # the methods are internal (not registered), so call them directly
    bsm <- if (inherits(m, "lme")) predictmeans:::build_Sigma_mats.lme else predictmeans:::build_Sigma_mats.gls
    V_list <- bsm(m, sigma_scale = FALSE)
    ord <- predictmeans:::get_order(A = V_list, B = cbind(X, y = 0))[[2]]
    if (!identical(as.integer(ord), seq_len(nrow(d)))) stop("data not in group order")
    W <- pm_whiten_blocks(V_list)
    y <- nlme::getResponse(m)
  }
  list(Xw = W$solve(X), yw = W$solve(as.numeric(y)))
}

## all p-values for one fitted model; 'null' = names of the true null pairs
pair_methods <- function(m, d, term, null) {
  K <- pair_K(m, d, term)
  k <- nrow(K)
  # joint: permmodels() refits, |t| from each refit's coefficients and covariance
  pmod <- permmodels(m, nperm = nperm, prt = FALSE, ncore = 1)
  Tj <- pm_perm_tstats(pmod$permlist, K)
  b <- if (inherits(m, "merMod")) fixef(m) else if (inherits(m, "lme")) fixef(m) else coef(m)
  V <- as.matrix(vcov(m))
  obs <- abs(as.vector(K %*% b)) / sqrt(diag(K %*% V %*% t(K)))
  p_joint <- pm_perm_p(obs, Tj)
  p_joint_max <- pm_perm_p(obs, matrix(rep(apply(Tj, 2, max), each = k), k))
  # Freedman-Lane with V fixed
  w <- whitened(m, d)
  idx <- vapply(seq_len(nperm), function(i) sample.int(nrow(d)), integer(nrow(d)))
  fl <- fl_pairs(w$Xw, w$yw, K, idx)
  # parametric t tests and Tukey from predictmeans()
  pr <- as.data.frame(predictmeans(m, term, adj = "tukey", plot = FALSE, prtnum = FALSE)$pairs)
  key <- paste0(pmin(as.character(pr$level1), as.character(pr$level2)), "v",
                pmax(as.character(pr$level1), as.character(pr$level2)))
  p_par <- pr$p.value[match(rownames(K), key)]
  p_tuk <- pr$p.adj[match(rownames(K), key)]
  nul <- rownames(K) %in% null
  fam <- function(padj) if (any(nul)) min(padj[nul]) else NA_real_
  out <- c(stats::setNames(p_joint, paste0("joint.", rownames(K))),
           stats::setNames(fl$pair, paste0("flpair.", rownames(K))),
           stats::setNames(fl$term, paste0("flterm.", rownames(K))),
           stats::setNames(p_par, paste0("param.", rownames(K))),
           joint.holm = fam(p.adjust(p_joint, "holm")), joint.maxT = fam(p_joint_max),
           flpair.holm = fam(p.adjust(fl$pair, "holm")),
           flterm.holm = fam(p.adjust(fl$term, "holm")), flterm.maxT = fam(fl$term_max),
           param.holm = fam(p.adjust(p_par, "holm")), param.tukey = fam(p_tuk))
  out
}

rcbd <- function(nblock = 8, reps = 1) {
  d <- expand.grid(rep = seq_len(reps), A = factor(1:4), block = factor(seq_len(nblock)))
  d[order(d$block, d$A), ]
}

description <- list(
  pr_rcbd_global = "lmer y ~ A + (1|block), 8 blocks x 4 levels, normal, all A means equal; pairs of A",
  pr_rcbd_partial = "as pr_rcbd_global, A means (0,0,0,3): null pairs 1v2 1v3 2v3",
  pr_rcbd_skew_partial = "as pr_rcbd_partial with skewed errors and block effects (exponential)",
  pr_sp_A_partial = "split plot, 6 blocks, A (whole plot) means (0,0,0,2), B (0,1,2); pairs of A (null 1v2 1v3 2v3)",
  pr_sp_B_partial = "split plot, 6 blocks, B means (0,0,2), A (0,1,2,3); pairs of B (null 1v2)",
  pr_lev_global = "lmer y ~ x + A + (1|block), 8 blocks x 3 A x 2, x lognormal larger for A3 with slope 1, A null; pairs of A",
  pr_lme_varident = "lme y ~ A, random 1|block, varIdent by A (sd 1,1,2,3), 8 blocks x 4 x 2, A null; pairs of A",
  pr_gls_ar1 = "gls y ~ G + time, corAR1 within subject (rho 0.6), 3 groups x 6 subjects x 6 times, G null; pairs of G",
  pr_rcbd_power = "as pr_rcbd_global with A means (0,0.4,0.8,1.2); power of each pair")

one_set <- function(cell) {
  switch(cell,
    pr_rcbd_global = , pr_rcbd_partial = , pr_rcbd_skew_partial = , pr_rcbd_power = {
      d <- rcbd()
      mu <- switch(cell, pr_rcbd_global = c(0, 0, 0, 0), pr_rcbd_power = c(0, 0.4, 0.8, 1.2), c(0, 0, 0, 3))
      if (cell == "pr_rcbd_skew_partial") {
        d$y <- mu[d$A] + skew(8)[d$block] + skew(nrow(d))
      } else d$y <- mu[d$A] + rnorm(8)[d$block] + rnorm(nrow(d))
      m <- lmer(y ~ A + (1 | block), data = d, control = ctl)
      null <- switch(cell, pr_rcbd_global = c("1v2", "1v3", "2v3", "1v4", "2v4", "3v4"), pr_rcbd_power = character(),
                     c("1v2", "1v3", "2v3"))
      pair_methods(m, d, "A", null)
    },
    pr_sp_A_partial = , pr_sp_B_partial = {
      d <- sp_design()
      d <- d[order(d$block, d$wp, d$sub), ]
      if (cell == "pr_sp_A_partial") {
        muA <- c(0, 0, 0, 2); muB <- c(0, 1, 2)
      } else {
        muA <- c(0, 1, 2, 3); muB <- c(0, 0, 2)
      }
      d$y <- muA[d$A] + muB[d$B] + rnorm(6)[d$block] + rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
      m <- lmer(y ~ A * B + (1 | block) + (1 | block:A), data = d, control = ctl)
      if (cell == "pr_sp_A_partial") pair_methods(m, d, "A", c("1v2", "1v3", "2v3")) else pair_methods(m, d, "B", "1v2")
    },
    pr_lev_global = {
      d <- lev_design(8)
      d <- d[order(d$block, d$A, d$rep), ]
      d$y <- d$x + rnorm(8)[d$block] + rnorm(nrow(d))
      m <- lmer(y ~ x + A + (1 | block), data = d, control = ctl)
      pair_methods(m, d, "A", c("1v2", "1v3", "2v3"))
    },
    pr_lme_varident = {
      d <- rcbd(8, 2)
      d$y <- rnorm(8)[d$block] + c(1, 1, 2, 3)[d$A] * rnorm(nrow(d))
      m <- lme(y ~ A, random = ~ 1 | block, weights = varIdent(form = ~ 1 | A), data = d)
      pair_methods(m, d, "A", c("1v2", "1v3", "2v3", "1v4", "2v4", "3v4"))
    },
    pr_gls_ar1 = {
      d <- expand.grid(t = 1:6, subject = factor(1:18))
      d$G <- factor(rep(1:3, each = 6))[d$subject]
      d$time <- factor(d$t)
      e <- unlist(lapply(1:18, function(s) as.numeric(stats::arima.sim(list(ar = 0.6), n = 6, sd = sqrt(1 - 0.36)))))
      d$y <- c(0, 0.5, 1, 1.5, 2, 2.5)[d$t] + e
      m <- gls(y ~ G + time, correlation = corAR1(form = ~ t | subject), data = d)
      pair_methods(m, d, "G", c("1v2", "1v3", "2v3"))
    },
    stop("unknown cell ", cell))
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-perm-pairs", description)
cat("SIM DONE\n")
