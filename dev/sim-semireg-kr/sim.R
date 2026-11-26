## Coverage of 95% limits of pairwise differences of semireg means
## (2026-11-26): Kenward-Roger adjusted (the default since b05cb33) against
## unadjusted SEs with Satterthwaite df (vcov_adjust = FALSE). All true
## differences are 0, so a pair's limits cover when they contain 0.
## A: the Slate Hall layout (25 varieties x 6 reps, 10 x 15 field),
##    yield ~ variety + t2(row, col) + (1 | rowf) + (1 | colf), a smooth field
##    trend, row and column effects and plot error.
## B: 4 groups x 6 blocks, y ~ g + s(x) + (1 | blk), a sine trend in x.
## Run from the package root: Rscript dev/sim-semireg-kr/sim.R  (about 20 min, 4 cores)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
sh <- read.csv(system.file("extdata", "slatehall.csv", package = "predictmeans"))
one_A <- function(seed) {
  set.seed(seed)
  d <- sh
  d$variety <- factor(d$variety); d$rowf <- factor(d$fieldrow); d$colf <- factor(d$fieldcolumn)
  r <- d$fieldrow / 10; cl <- d$fieldcolumn / 15
  d$yield <- 2 * sin(2 * r) * cos(3 * cl) + rnorm(10, 0, 0.4)[d$fieldrow] + rnorm(15, 0, 0.4)[d$fieldcolumn] +
    rnorm(nrow(d), 0, 0.6)
  f <- tryCatch(suppressMessages(semireg(yield ~ variety + t2(fieldrow, fieldcolumn) + (1 | rowf) + (1 | colf),
                                         data = d)), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA))
  cov <- function(p) mean(p$pairs$conf.low <= 0 & p$pairs$conf.high >= 0)
  c(KR = cov(suppressMessages(predictmeans(f, "variety"))),
    unadj = cov(suppressMessages(predictmeans(f, "variety", vcov_adjust = FALSE))))
}
one_B <- function(seed) {
  set.seed(seed)
  d <- expand.grid(g = factor(1:4), blk = factor(1:6))
  d$x <- runif(nrow(d))
  d$y <- sin(3 * d$x) + rnorm(6, 0, 0.5)[d$blk] + rnorm(nrow(d), 0, 0.4)
  f <- tryCatch(suppressMessages(semireg(y ~ g + s(x) + (1 | blk), data = d)), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA))
  cov <- function(p) mean(p$pairs$conf.low <= 0 & p$pairs$conf.high >= 0)
  c(KR = cov(suppressMessages(predictmeans(f, "g"))),
    unadj = cov(suppressMessages(predictmeans(f, "g", vcov_adjust = FALSE))))
}
run <- function(fun, reps, base) {
  P <- do.call(rbind, parallel::mclapply(seq_len(reps), function(r) fun(base + r), mc.cores = 4))
  P <- P[stats::complete.cases(P), , drop = FALSE]
  c(sets = nrow(P), colMeans(P))
}
out <- rbind(A_slatehall_t2 = run(one_A, 200, 1000), B_small_s = run(one_B, 500, 5000))
print(round(out, 4))
saveRDS(out, "dev/sim-semireg-kr/results.rds")
