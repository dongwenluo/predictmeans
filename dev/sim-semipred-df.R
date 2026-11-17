## Coverage of semipred() bands of a gaussian semireg() fit: normal limits
## (predictmeans 1.x, gam-style) vs t limits with Satterthwaite df at each
## point (pm_sr_satt_ddf()), the df either from all variance parameters
## ("t_all") or with the smoothing variances held at their estimates
## ("t_design": only the design's variance components and the residual
## count); simultaneous bands: pm_sim_crit(df = ).
## Designs: a split-plot with curves over time (6 blocks, 3 treatments on
## whole plots, 11 times; whole-plot SD 0.2 or 0.4), and curves by group
## without design strata (3 groups of 40, x uniform); no treatment effect,
## true curve 2 + sin(t/2).
## Per data set, over 21 points: the difference a2 - a1 (true value 0) and the
## curve of a1; pointwise coverage (mean over the points) and simultaneous
## coverage (the band covers at every point).
##   Rscript dev/sim-semipred-df.R [nsim] [out.rds]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 500
q <- function(e) suppressWarnings(suppressMessages(e))
tt <- seq(0, 10, by = 0.5)
fit_of <- function(seed, design) {
  set.seed(seed)
  if (design == "smooth") {
    d <- data.frame(t = runif(120, 0, 10), A = gl(3, 40, labels = c("a1", "a2", "a3")))
    d$y <- 2 + sin(d$t/2) + rnorm(120, 0, 0.3)
    return(q(semireg(y ~ s(t, by = A, bs = "os", k = 8), data = d, prt = FALSE)))
  }
  sd_plot <- if (design == "split0.2") 0.2 else 0.4
  d <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  d$y <- 2 + sin(d$t/2) + rnorm(6, 0, 0.5)[d$Block] + rnorm(18, 0, sd_plot)[interaction(d$Block, d$A)] +
    rnorm(nrow(d), 0, 0.3)
  q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), data = d, prt = FALSE))
}
one <- function(seed, design) {
  f <- fit_of(seed, design)
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
    dfs <- list(all = pm_sr_satt_ddf(f, L, des$cov_idx),
                design = pm_sr_satt_ddf(f, L, des$cov_idx, hold = names(f$refit$smoothZ)))
    dfs <- lapply(dfs, function(x) { x[is.na(x)] <- Inf; x })
    crit <- list(z_pw = qnorm(0.975), tall_pw = qt(0.975, dfs$all), tdes_pw = qt(0.975, dfs$design),
                 z_sim = pm_sim_crit(L, des$sm_Cov, 0.05),
                 tall_sim = pm_sim_crit(L, des$sm_Cov, 0.05, df = dfs$all),
                 tdes_sim = pm_sim_crit(L, des$sm_Cov, 0.05, df = dfs$design))
    inside <- lapply(crit, function(cr) abs(est - truth) <= cr*se)
    cov <- vapply(names(inside), function(nm) if (grepl("_pw", nm)) mean(inside[[nm]]) else all(inside[[nm]]),
                  numeric(1))
    out <- c(out, stats::setNames(c(cov, min(dfs$all), median(dfs$all), min(dfs$design), median(dfs$design)),
                                  paste(what, c(names(crit), "dfall_min", "dfall_med", "dfdes_min", "dfdes_med"),
                                        sep = ".")))
  }
  out
}
res <- list()
for (design in c("split0.2", "split0.4", "smooth")) {
  r <- parallel::mclapply(seq_len(nsim), function(i) tryCatch(one(2000 + i, design), error = function(e) NULL),
                          mc.cores = parallel::detectCores())
  r <- do.call(rbind, r[!vapply(r, is.null, logical(1))])
  res[[design]] <- r
  m <- colMeans(r)
  cat(sprintf("%s (%d data sets, MC SE ~%.3f for simultaneous):\n", design, nrow(r), sqrt(0.95*0.05/nrow(r))))
  for (what in c("diff", "curve")) {
    v <- function(k) m[[paste0(what, ".", k)]]
    cat(sprintf("  %-5s pointwise z %.3f t_all %.3f t_design %.3f | simultaneous z %.3f t_all %.3f t_design %.3f\n",
                what, v("z_pw"), v("tall_pw"), v("tdes_pw"), v("z_sim"), v("tall_sim"), v("tdes_sim")))
    cat(sprintf("        median df: all %.1f (smallest point %.1f), design %.1f (smallest point %.1f)\n",
                median(r[, paste0(what, ".dfall_med")]), median(r[, paste0(what, ".dfall_min")]),
                median(r[, paste0(what, ".dfdes_med")]), median(r[, paste0(what, ".dfdes_min")])))
  }
}
if (length(args) > 1) saveRDS(res, args[2])
