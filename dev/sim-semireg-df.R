## Type I error of comparisons of treatment means at a covariate value in a
## split-plot design with smooth curves over time (semireg()): normal-based
## tests (predictmeans 1.x and gam-style z) vs Satterthwaite df for the
## combination of fixed and smooth coefficients (pm_sr_satt_ddf()).
## No treatment effect; 6 blocks, 3 treatments on whole plots, 11 times.
##   Rscript dev/sim-semireg-df.R [nsim] [out.rds]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 500
q <- function(e) suppressWarnings(suppressMessages(e))
one <- function(seed, sd_plot) {
  set.seed(seed)
  d <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  d$y <- 2 + sin(d$t/2) + rnorm(6, 0, 0.5)[d$Block] + rnorm(18, 0, sd_plot)[interaction(d$Block, d$A)] +
    rnorm(nrow(d), 0, 0.3)
  f <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), data = d, prt = FALSE))
  pm <- q(predictmeans(f, "A", covariate = c(t = 5), pairwise = TRUE))
  pr <- pm$pairs
  z <- pr$estimate/pr$std.error
  c(p_z = 2*pnorm(-abs(z[1])), p_t = pr$p.value[1], df = pr$df[1])
}
res <- list()
for (sd_plot in c(0.2, 0.4)) {
  r <- t(vapply(seq_len(nsim), function(i) one(1000 + i, sd_plot), numeric(3)))
  res[[paste0("sd_plot_", sd_plot)]] <- r
  cat(sprintf("sd_plot %.1f: reject at 0.05 z %.3f, Satterthwaite t %.3f (median df %.1f), MC SE %.3f\n", sd_plot,
              mean(r[, "p_z"] < 0.05), mean(r[, "p_t"] < 0.05), median(r[, "df"]), sqrt(0.05*0.95/nsim)))
}
if (length(args) > 1) saveRDS(res, args[2])
