## Type I error of comparisons of treatment means at a time point in a
## repeated-measures design with curves over time and AR(1) residuals,
## fitted with semireg_tmb(): normal-based tests (z) vs Satterthwaite df for
## the combination of fixed and smooth coefficients (pm_sr_satt_ddf_tmb()).
## 3 treatments x 6 subjects, 10 times, random subject intercept, AR(1)
## residuals (rho = 0.6), no treatment effect; a1 vs a2 at t = 5.
##   Rscript dev/sim-semireg-tmb-df.R [nsim] [out.rds] [log]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 300
out <- if (length(args) > 1) args[2] else file.path(tempdir(), "sim-semireg-tmb-df.rds")
logf <- if (length(args) > 2) args[3] else file.path(tempdir(), "sim-semireg-tmb-df.log")
q <- function(e) suppressWarnings(suppressMessages(e))
one <- function(seed, ns = 6, nt = 10, rho = 0.6) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(nt), id = factor(seq_len(3*ns)))
  d$A <- factor(c("a1", "a2", "a3")[(as.integer(d$id) - 1) %/% ns + 1])
  e <- unlist(lapply(seq_len(3*ns), function(i) as.numeric(stats::arima.sim(list(ar = rho), nt, sd = 0.4))))
  d$y <- 2 + sin(d$t/3) + rnorm(3*ns, 0, 0.4)[d$id] + e
  d$tf <- factor(d$t)
  r <- tryCatch({
    f <- q(semireg_tmb(y ~ s(t, by = A, bs = "os", k = 5) + (1 | id) + ar1(tf + 0 | id), dispformula = ~0,
                       data = d, prt = FALSE))
    pr <- q(predictmeans(f, "A", covariate = c(t = 5), pairwise = TRUE))$pairs
    z <- pr$estimate[1]/pr$std.error[1]
    c(p_z = 2*pnorm(-abs(z)), p_t = pr$p.value[1], df = pr$df[1])
  }, error = function(e) c(p_z = NA, p_t = NA, df = NA))
  cat(seed, r, "\n", file = logf, append = TRUE)
  r
}
res <- do.call(rbind, parallel::mclapply(seq_len(nsim), function(i) one(3000 + i), mc.cores = parallel::detectCores()))
saveRDS(res, out)
ok <- stats::complete.cases(res)
cat(sprintf("AR(1) rho 0.6, %d of %d fits: reject at 0.05 z %.3f, Satterthwaite t %.3f (median df %.1f), MC SE %.3f\n",
            sum(ok), nsim, mean(res[ok, "p_z"] < 0.05), mean(res[ok, "p_t"] < 0.05), median(res[ok, "df"]),
            sqrt(0.05*0.95/sum(ok))))
