## Coverage of semipred() bands of a gaussian semireg() fit in a split-plot
## design with curves over time: normal limits (predictmeans 1.x, and gam-
## style) vs t limits with Satterthwaite df at each point (pm_sr_satt_ddf();
## simultaneous: pm_sim_crit(df = )). 6 blocks, 3 treatments on whole plots,
## 11 times; no treatment effect, true curve 2 + sin(t/2).
## Per data set, over 21 points: the difference a2 - a1 (true value 0) and the
## curve of a1; pointwise coverage (mean over the points) and simultaneous
## coverage (the band covers at every point).
##   Rscript dev/sim-semipred-df.R [nsim] [out.rds]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 500
q <- function(e) suppressWarnings(suppressMessages(e))
tt <- seq(0, 10, by = 0.5)
one <- function(seed, sd_plot) {
  set.seed(seed)
  d <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  d$y <- 2 + sin(d$t/2) + rnorm(6, 0, 0.5)[d$Block] + rnorm(18, 0, sd_plot)[interaction(d$Block, d$A)] +
    rnorm(nrow(d), 0, 0.3)
  f <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), data = d, prt = FALSE))
  des <- q(pm_semi_design(f, "A", "t", tt, NULL, f$data))
  D <- as.matrix(des$D_mx)
  g <- des$pred_df$A
  Da1 <- D[g == "a1", , drop = FALSE]
  out <- NULL
  for (what in c("diff", "curve")) {
    L <- if (what == "diff") D[g == "a2", , drop = FALSE] - Da1 else Da1
    truth <- if (what == "diff") 0 else 2 + sin(tt/2)
    est <- as.numeric(L %*% des$beta_u)
    se <- des$sig.epsHat*sqrt(rowSums((L %*% as.matrix(des$sm_Cov)) * L))
    df <- pm_sr_satt_ddf(f, L, des$cov_idx)
    df[is.na(df)] <- Inf
    crit <- list(z_pw = qnorm(0.975), t_pw = qt(0.975, df),
                 z_sim = pm_sim_crit(L, des$sm_Cov, 0.05), t_sim = pm_sim_crit(L, des$sm_Cov, 0.05, df = df))
    inside <- lapply(crit, function(cr) abs(est - truth) <= cr*se)
    cov <- c(z_pw = mean(inside$z_pw), t_pw = mean(inside$t_pw), z_sim = all(inside$z_sim), t_sim = all(inside$t_sim))
    out <- c(out, stats::setNames(c(cov, median(df)), paste(what, c(names(crit), "df"), sep = ".")))
  }
  out
}
res <- list()
for (sd_plot in c(0.2, 0.4)) {
  r <- parallel::mclapply(seq_len(nsim), function(i) tryCatch(one(2000 + i, sd_plot), error = function(e) NULL),
                          mc.cores = parallel::detectCores())
  r <- do.call(rbind, r[!vapply(r, is.null, logical(1))])
  res[[paste0("sd_plot_", sd_plot)]] <- r
  m <- colMeans(r)
  cat(sprintf("sd_plot %.1f (%d data sets, MC SE ~%.3f for simultaneous):\n", sd_plot, nrow(r), sqrt(0.95*0.05/nrow(r))))
  for (what in c("diff", "curve")) {
    cat(sprintf("  %-5s pointwise z %.3f t %.3f | simultaneous z %.3f t %.3f | median df %.1f\n", what,
                m[paste0(what, ".z_pw")], m[paste0(what, ".t_pw")], m[paste0(what, ".z_sim")], m[paste0(what, ".t_sim")],
                median(r[, paste0(what, ".df")])))
  }
}
if (length(args) > 1) saveRDS(res, args[2])
