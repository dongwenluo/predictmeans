## Generalized semireg() fits: normal-based tests and bands (before) vs t
## with the inner-outer df of the design's strata (pm_sr_ddf(), the GLMM
## default of predictmeans()). Poisson counts in a split-plot design with
## curves over time: 6 blocks, 3 treatments on whole plots, 11 times; no
## treatment effect; log mean 1 + sin(t/2)/2 + block + whole plot (SD 0.3
## and SD 0.2 or 0.4). The covariance includes the smoothing-parameter
## correction (unconditional, the default for generalized fits).
## Per data set: the comparison a1 vs a2 at t = 5 (rejection at 0.05), and
## over 21 points the difference a2 - a1 (true value 0) and the curve of a1
## on the link scale (true value 1 + sin(t/2)/2): pointwise coverage (mean
## over the points) and simultaneous coverage.
## Result (2026-11-17, 1000 data sets per setting, MC SE 0.007; t on the 10
## whole-plot df):
##   plot SD  a1 vs a2 at t = 5    difference: pointwise / simultaneous    curve: pointwise / simultaneous
##   0.2      z 5.8%, t 3.1%        z 0.946 / 0.957, t 0.973 / 0.996        z 0.922 / 0.909, t 0.955 / 0.959
##   0.4      z 9.2%, t 5.9%        z 0.923 / 0.933, t 0.954 / 0.981        z 0.927 / 0.924, t 0.957 / 0.970
## Decision: the inner-outer df (the GLMM default): z is liberal; t is
## conservative for differences when the whole-plot variance is small, as the
## same rule is for glmer fits.
##   Rscript dev/sim-semireg-glmm-df.R [nsim] [out.rds]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 500
q <- function(e) suppressWarnings(suppressMessages(e))
tt <- seq(0, 10, by = 0.5)
one <- function(seed, sd_plot) {
  set.seed(seed)
  d <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  eta <- 1 + sin(d$t/2)/2 + rnorm(6, 0, 0.3)[d$Block] + rnorm(18, 0, sd_plot)[interaction(d$Block, d$A)]
  d$y <- rpois(nrow(d), exp(eta))
  f <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), family = poisson, data = d,
                 prt = FALSE))
  pm <- q(predictmeans(f, "A", covariate = c(t = 5), pairwise = TRUE))
  pr <- pm$pairs[1, ]
  out <- c(rej_z = 2*pnorm(-abs(pr$estimate/pr$std.error)) < 0.05, rej_t = pr$p.value < 0.05, df_pair = pr$df)
  des <- q(pm_semi_design(f, "A", "t", tt, NULL, f$data))
  D <- as.matrix(des$D_mx)
  g <- des$pred_df$A
  Da1 <- D[g == "a1", , drop = FALSE]
  for (what in c("diff", "curve")) {
    L <- if (what == "diff") D[g == "a2", , drop = FALSE] - Da1 else Da1
    truth <- if (what == "diff") 0 else 1 + sin(tt/2)/2
    est <- as.numeric(L %*% des$beta_u)
    se <- des$sig.epsHat*sqrt(rowSums((L %*% as.matrix(des$sm_Cov)) * L))
    df <- as.numeric(pm_sr_ddf(f, L, des$cov_idx, if (what == "curve") "A"))
    crit <- list(z_pw = qnorm(0.975), t_pw = qt(0.975, df),
                 z_sim = pm_sim_crit(L, des$sm_Cov, 0.05), t_sim = pm_sim_crit(L, des$sm_Cov, 0.05, df = df))
    inside <- lapply(crit, function(cr) abs(est - truth) <= cr*se)
    cov <- vapply(names(inside), function(nm) if (grepl("_pw", nm)) mean(inside[[nm]]) else all(inside[[nm]]),
                  numeric(1))
    out <- c(out, stats::setNames(c(cov, median(df)), paste(what, c(names(crit), "df"), sep = ".")))
  }
  out
}
res <- list()
for (sd_plot in c(0.2, 0.4)) {
  r <- parallel::mclapply(seq_len(nsim), function(i) tryCatch(one(3000 + i, sd_plot), error = function(e) NULL),
                          mc.cores = parallel::detectCores())
  r <- do.call(rbind, r[!vapply(r, is.null, logical(1))])
  res[[paste0("sd_plot_", sd_plot)]] <- r
  m <- colMeans(r)
  cat(sprintf("sd_plot %.1f (%d data sets, MC SE %.3f):\n", sd_plot, nrow(r), sqrt(0.95*0.05/nrow(r))))
  cat(sprintf("  a1 vs a2 at t = 5: reject z %.3f, t %.3f (df %.0f)\n", m[["rej_z"]], m[["rej_t"]],
              median(r[, "df_pair"])))
  for (what in c("diff", "curve")) {
    v <- function(k) m[[paste0(what, ".", k)]]
    cat(sprintf("  %-5s pointwise z %.3f t %.3f | simultaneous z %.3f t %.3f | df %.0f\n", what,
                v("z_pw"), v("t_pw"), v("z_sim"), v("t_sim"), median(r[, paste0(what, ".df")])))
  }
}
if (length(args) > 1) saveRDS(res, args[2])
