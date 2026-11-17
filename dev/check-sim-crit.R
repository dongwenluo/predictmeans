## The critical value of simultaneous bands (pm_sim_crit(), R/p-adjust.R)
## against mvtnorm: qmvnorm() / qmvt() of the band's correlation, for 37
## points along a fitted curve and 5 points, df Inf, 30 and 8; and, for df
## varying along the band (no mvtnorm counterpart), against the exact
## comonotone construction with qchisq().
## Result (2026-11-17, after replacing the Wilson-Hilferty approximation of
## the chi-square, which made the value 0.77% too large at df 8, by the exact
## quantile transform):
##   df Inf: +0.21% (37 points), -0.22% (5 points); df 30: +0.20%, -0.23%;
##   df 8: +0.55% (5 points; mvtnorm fails for 37 nearly collinear points).
##   With 2e5 draws: 0.15-0.26% from mvtnorm, about the Monte Carlo error of a
##   95% quantile; the default 1e4 draws: about 0.3%.
##   df 5 to 60 along the band: 0.12% from the exact construction.
##   Rscript dev/check-sim-crit.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(11)
d <- data.frame(x = c(0, 10, runif(118, 0, 10))); d$y <- sin(d$x / 1.5) + rnorm(120, 0, 0.4)
f <- q(semireg(y ~ s(x, bs = "os", k = 10), data = d, prt = FALSE))
xs <- seq(0.5, 9.5, by = 0.25)
des <- q(pm_semi_design(f, "x", "x", xs, NULL, f$data))
D <- as.matrix(des$D_mx); C <- as.matrix(des$sm_Cov)
cases <- list("37 points of a curve" = seq_len(nrow(D)), "5 points" = round(seq(1, nrow(D), length = 5)))
for (nm in names(cases)) {
  i <- cases[[nm]]
  Di <- D[i, , drop = FALSE]
  R <- cov2cor(Di %*% C %*% t(Di))
  for (df in c(Inf, 30, 8)) {
    ours <- pm_sim_crit(Di, C, 0.05, df = df)[1]
    alg <- mvtnorm::GenzBretz(maxpts = 2e6, abseps = 1e-5)
    ref <- tryCatch(if (is.finite(df)) mvtnorm::qmvt(0.95, tail = "both.tails", corr = R, df = df, algorithm = alg)$quantile else
      mvtnorm::qmvnorm(0.95, tail = "both.tails", corr = R, algorithm = alg)$quantile, error = function(e) NA_real_)
    big <- pm_sim_crit(Di, C, 0.05, df = df, nsim = 2e5)[1]
    cat(sprintf("%-22s df %-4s: pm_sim_crit %.4f (nsim 2e5: %.4f) | mvtnorm %.4f | difference %.4f (%.2f%%)\n",
                nm, format(df), ours, big, ref, ours - ref, 100 * (ours / ref - 1)))
  }
}
## df varying along the band: Wilson-Hilferty vs the exact quantile transform of one chi-square
i <- cases[[1]]; Di <- D[i, , drop = FALSE]
dfv <- seq(5, 60, length = length(i))
wh <- pm_sim_crit(Di, C, 0.05, df = dfv)
se <- sqrt(rowSums((Di %*% C) * Di)); e <- eigen(C, symmetric = TRUE); A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
set.seed(20240403); nsim <- 1e5
Z <- matrix(rnorm(ncol(A) * nsim), ncol(A)); u <- pnorm(rnorm(nsim))
T0 <- abs(Di %*% A %*% Z) / se
M <- sapply(seq_len(nsim), function(j) max(pt(T0[, j] / sqrt(qchisq(u[j], dfv) / dfv), dfv)))
exact <- qt(quantile(M, 0.95), dfv)
cat(sprintf("df 5 to 60 along the curve: pm_sim_crit vs exact chi-square quantiles: max difference %.4f (%.2f%%)\n",
            max(abs(wh - exact)), 100 * max(abs(wh / exact - 1))))
