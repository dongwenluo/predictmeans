## Type I error of anova()'s likelihood ratio tests of the smooth variance
## components (halved p-values, chi-square 0:1) in s(x, f1*f2): data whose
## curves have the same shape for every level of f2 and f1:f2.
## Scenario "f1": the shape differs by f1; scenario "none": one shape.
##   Rscript dev/sim-semireg-lrt.R [nsim] [out.rds]
suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 200
q <- function(e) suppressWarnings(suppressMessages(e))
one <- function(seed, f1_shape) {
  set.seed(seed)
  d <- expand.grid(x = seq(0, 1, length = 25), f1 = factor(c("p", "q")), f2 = factor(c("u", "v", "w")), rep = 1:4)
  d$id <- interaction(d$f1, d$f2, d$rep)
  d$y <- sin(2*pi*d$x) + f1_shape*(d$f1 == "q")*sin(pi*d$x)^2 + rnorm(nlevels(d$id), 0, 0.2)[d$id] +
    rnorm(nrow(d), 0, 0.25)
  fit <- q(semireg(y ~ s(x, f1*f2, bs = "os", k = 8) + (1 | id), data = d, prt = FALSE))
  sm <- anova(fit)$smooth
  stats::setNames(sm$p.value, sm$term)[c("s_x_f1", "s_x_f2", "s_x_f1_f2")]
}
res <- list()
for (sc in c(none = 0, f1 = 1.5)) {
  r <- t(vapply(seq_len(nsim), function(i) one(2000 + i, sc), numeric(3)))
  res[[as.character(sc)]] <- r
  cat(sprintf("f1 shape %.1f: rejection at 0.05: f1 %.3f, f2 %.3f, f1:f2 %.3f (MC SE %.3f)\n", sc,
              mean(r[, 1] < 0.05), mean(r[, 2] < 0.05), mean(r[, 3] < 0.05), sqrt(0.05*0.95/nsim)))
}
if (length(args) > 1) saveRDS(res, args[2])
