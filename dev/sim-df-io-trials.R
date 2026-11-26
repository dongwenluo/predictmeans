## Inner-outer df of aggregated binomial GLMMs: count rows or trials? (2026-11-26)
## Size of the t test of a within-block treatment contrast at the 5% level,
## true null, glmer(cbind(y, m - y) ~ trt + (1 | blk), nAGQ = 9), block SD 1,
## p = 0.5. The residual stratum has rows - rank df (SAS GLIMMIX with
## events/trials) or trials - rank df (the same data as one row per trial).
## Run: Rscript dev/sim-df-io-trials.R  (about 5 minutes)
suppressMessages(library(lme4))
one_case <- function(ntrt, nblk, m, reps = 1000, seed = 20261126) {
  set.seed(seed)
  d <- expand.grid(trt = factor(seq_len(ntrt)), blk = factor(seq_len(nblk))); d$m <- m
  rk <- ntrt + nblk - 1
  df_rows <- nrow(d) - rk; df_trials <- m * nrow(d) - rk
  p <- t(replicate(reps, {
    b <- rnorm(nblk)
    d$y <- rbinom(nrow(d), m, plogis(b[d$blk]))
    f <- tryCatch(suppressWarnings(suppressMessages(
      glmer(cbind(y, m - y) ~ trt + (1 | blk), data = d, family = binomial, nAGQ = 9))), error = function(e) NULL)
    if (is.null(f)) c(NA, NA, NA) else {
      z <- coef(summary(f))["trt2", "z value"]
      c(2 * pt(-abs(z), df_rows), 2 * pt(-abs(z), df_trials), 2 * pnorm(-abs(z)))
    }
  }))
  p <- p[stats::complete.cases(p), , drop = FALSE]
  data.frame(design = sprintf("%d trt x %d blocks, %d trials/row", ntrt, nblk, m), df_rows, df_trials,
             fits = nrow(p), rows = mean(p[, 1] < 0.05), trials = mean(p[, 2] < 0.05), z = mean(p[, 3] < 0.05))
}
res <- rbind(one_case(4, 3, 25), one_case(4, 3, 5), one_case(4, 6, 2), one_case(2, 3, 10))
print(res, digits = 3)
