# Timings for Kenward-Roger and Satterthwaite df and predictmeans().
# Usage: Rscript --no-init-file dev/benchmarks.R [library]
# Compares with pbkrtest and lmerTest when they are installed.

lib <- commandArgs(TRUE)[1]
suppressPackageStartupMessages({
  if (is.na(lib)) library(predictmeans) else library(predictmeans, lib.loc = lib)
})
tm <- function(expr, reps = 3) {
  e <- substitute(expr); pf <- parent.frame()
  min(vapply(seq_len(reps), function(i) system.time(eval(e, pf))[["elapsed"]], numeric(1)))
}
sim <- function(n_block, n_rep = 2) {
  set.seed(1)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), Block = factor(seq_len(n_block)), rep = seq_len(n_rep))
  d$WP <- interaction(d$Block, d$A)
  d$y <- as.numeric(d$A) + as.numeric(d$B) + rnorm(n_block)[d$Block] + rnorm(nlevels(d$WP))[d$WP] + rnorm(nrow(d))
  d
}
has <- function(p) requireNamespace(p, quietly = TRUE)
rows <- list()
for (nb in c(10, 50, 250)) {
  d <- sim(nb)
  m <- lme4::lmer(y ~ A*B + (1|Block) + (1|WP), data = d)
  L <- Kmatrix(m, "A:B")$K
  r <- data.frame(n = nrow(d),
    kr_vcov = tm(predictmeans:::pm_vcov_kr(m)),
    kr_df = tm(df_term(m, ctrmatrix = L, type = "Kenward-Roger")),
    satt_df = tm(df_term(m, ctrmatrix = L, type = "Satterthwaite")),
    predictmeans = tm(predictmeans(m, "A:B", adj = "BH", plot = FALSE, prtnum = FALSE)),
    containment = tm(predictmeans(m, "A:B", df_method = "containment", plot = FALSE, prtnum = FALSE)),
    pbkrtest_vcov = if (has("pbkrtest") && nrow(d) <= 3000) tm(pbkrtest::vcovAdj(m), 1) else NA,
    lmerTest_satt = if (has("lmerTest")) tm({mt <- lmerTest::as_lmerModLmerTest(m); for (i in seq_len(nrow(L))) lmerTest::contest(mt, L[i, ], ddf = "Satterthwaite")}, 1) else NA)
  rows[[length(rows) + 1]] <- r
  print(r, row.names = FALSE)
}
res <- do.call(rbind, rows)
cat("\nSeconds (minimum of 3 runs; pbkrtest/lmerTest 1 run):\n"); print(res, row.names = FALSE)
