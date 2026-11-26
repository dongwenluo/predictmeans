## Tests of a between-cluster treatment in a cumulative logit mixed model with
## few clusters (2026-11-26, B4 G3, dev/review-b4.md). clmm(y ~ trt * W +
## (1 | clus)), clusters split half/half between two treatments, 6 rows per
## cluster (W within clusters, 2 levels x 3), 4 classes, cluster SD 1.
## Rejection at 5% of the trt test: inner-outer t (predictmeans' df, nclus - 2),
## z (Wald chi-square), parametric bootstrap (bootmodels(), 99 samples) and
## whole-cluster permutation (permmodels(perms = permindex(group =)), 99).
## Run from the package root: Rscript dev/sim-clmm-few/sim.R  (about 40 min, 4 cores)
suppressMessages({library(ordinal); pkgload::load_all(".", quiet = TRUE)})
one_set <- function(nclus, effect, seed) {
  set.seed(seed)
  d <- expand.grid(rep = 1:3, W = factor(1:2), clus = factor(seq_len(nclus)))
  d$trt <- factor(ifelse(as.integer(d$clus) <= nclus / 2, "a", "b"))
  eta <- rnorm(nclus)[d$clus] + effect * (d$trt == "b") + 0.5 * (d$W == "2")
  d$y <- cut(eta + rlogis(nrow(d)), c(-Inf, -1, 0.5, 2, Inf), labels = 1:4, ordered_result = TRUE)
  m <- tryCatch(suppressWarnings(clmm(y ~ trt * W + (1 | clus), data = d)), error = function(e) NULL)
  if (is.null(m)) return(rep(NA, 4))
  a <- tryCatch(suppressMessages(pm_anova(m, type = "II", test_statistic = "Chisq")), error = function(e) NULL)
  if (is.null(a)) return(rep(NA, 4))
  chi <- a["trt", "Chisq"]
  p_z <- stats::pchisq(chi, 1, lower.tail = FALSE)
  p_t <- 2 * stats::pt(-sqrt(chi), nclus - 2)                     # inner-outer df of trt
  b <- tryCatch(suppressMessages(bootmodels(m, nsim = 99, seed = seed)), error = function(e) NULL)
  p_b <- if (is.null(b)) NA else as.numeric(b$ANOVA["trt", "Boot_p_value"])
  pr <- tryCatch(suppressMessages(permmodels(m, perms = permindex(d, group = "clus", nsim = 99, seed = seed))),
                 error = function(e) NULL)
  p_p <- if (is.null(pr)) NA else as.numeric(pr$ANOVA["trt", "Perm_p_value"])
  c(t = p_t, z = p_z, boot = p_b, perm = p_p)
}
cases <- expand.grid(nclus = c(6, 8, 12), effect = 0)
cases <- rbind(cases, data.frame(nclus = 8, effect = 1.5))
reps <- 300
out <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i) {
  P <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) one_set(cases$nclus[i], cases$effect[i],
                                                                            1e5 * i + r), mc.cores = 4))
  ok <- stats::complete.cases(P)
  data.frame(nclus = cases$nclus[i], effect = cases$effect[i], sets = sum(ok),
             t(colMeans(P[ok, , drop = FALSE] < 0.05)))
}))
print(out, digits = 3)
saveRDS(out, "dev/sim-clmm-few/results.rds")
