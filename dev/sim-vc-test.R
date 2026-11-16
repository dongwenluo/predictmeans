## Type I error of the tests of single variance components (pm_vc_test()):
## the exact restricted LRT (default for gaussian lme4 REML fits) vs the
## 0.5 chi2(0) + 0.5 chi2(1) mixture, data with the tested variance 0.
## (a) semireg s(x, f1*f2): the f2 and f1:f2 deviation curves (the shapes
##     differ by f1 only); the mixture rejected 1-3% (dev/sim-semireg-lrt.R)
## (b) lmer, crossed random effects (1 | a) + (1 | b), b's variance 0
##   Rscript dev/sim-vc-test.R [nsim] [out.rds] [log]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 200
out <- if (length(args) > 1) args[2] else file.path(tempdir(), "sim-vc-test.rds")
logf <- if (length(args) > 2) args[3] else file.path(tempdir(), "sim-vc-test.log")
q <- function(e) suppressWarnings(suppressMessages(e))
pv <- function(fit, rows, method) {
  t <- q(pm_vc_test(fit, rows = rows, method = method, nsim = 4000))
  stats::setNames(t$p.value, t$row)[rows]
}
case_a <- function(seed) {
  set.seed(seed)
  d <- expand.grid(x = seq(0, 1, length = 25), f1 = factor(c("p", "q")), f2 = factor(c("u", "v", "w")), rep = 1:4)
  d$id <- interaction(d$f1, d$f2, d$rep)
  d$y <- sin(2*pi*d$x) + 1.5*(d$f1 == "q")*sin(pi*d$x)^2 + rnorm(nlevels(d$id), 0, 0.2)[d$id] +
    rnorm(nrow(d), 0, 0.25)
  fit <- q(semireg(y ~ s(x, f1*f2, bs = "os", k = 8) + (1 | id), data = d, prt = FALSE))
  rows <- c("s_x_f2.(Intercept)", "s_x_f1_f2.(Intercept)")
  c(exact = pv(fit, rows, "exact"), mixture = pv(fit, rows, "mixture"))
}
case_b <- function(seed) {
  set.seed(seed)
  d <- expand.grid(a = factor(1:8), b = factor(1:6), r = 1:2)
  d$y <- rnorm(8, 0, 0.7)[d$a] + rnorm(nrow(d))
  fit <- q(lme4::lmer(y ~ 1 + (1 | a) + (1 | b), d))
  c(exact = pv(fit, "b.(Intercept)", "exact"), mixture = pv(fit, "b.(Intercept)", "mixture"))
}
run <- function(fun, tag) {
  r <- parallel::mclapply(seq_len(nsim), function(i) {
    v <- tryCatch(fun(5000 + i), error = function(e) NULL)
    cat(tag, i, if (is.null(v)) "error" else round(v, 4), "\n", file = logf, append = TRUE)
    v
  }, mc.cores = parallel::detectCores())
  do.call(rbind, r[!vapply(r, is.null, logical(1))])
}
res <- list(a = run(case_a, "a"), b = run(case_b, "b"))
saveRDS(res, out)
for (nm in names(res)) {
  r <- res[[nm]]
  cat(sprintf("case %s (%d fits): rejection at 0.05: %s (MC SE %.3f)\n", nm, nrow(r),
              paste(sprintf("%s %.3f", colnames(r), colMeans(r < 0.05)), collapse = ", "), sqrt(0.05*0.95/nrow(r))))
}
