## Inner-outer df of aggregated binomial GLMMs: count rows or trials? (2026-11-26)
## Size of the t test of a within-block treatment contrast at the 5% level,
## true null, glmer(cbind(y, m - y) ~ trt + (1 | blk), nAGQ = 9), block SD 1,
## p = 0.5. The residual stratum has rows - rank df (predictmeans before
## 2026-11-26) or trials - rank df (the same data as one row per trial; the
## df pm_inner_outer_df() gives, checked below). 'od': SD of an extra
## row-level logit effect (overdispersion); 'olre': the fit has (1 | obs), so
## the treatment is in the rows' stratum and its df count rows.
## Run: Rscript dev/sim-df-io-trials.R  (about 10 minutes)
suppressMessages(library(lme4))
pkgload::load_all(".", quiet = TRUE)
one_case <- function(ntrt, nblk, m, od = 0, olre = FALSE, reps = 1000, seed = 20261126) {
  set.seed(seed)
  d <- expand.grid(trt = factor(seq_len(ntrt)), blk = factor(seq_len(nblk))); d$m <- m
  d$obs <- factor(seq_len(nrow(d)))
  f <- if (olre) cbind(y, m - y) ~ trt + (1 | blk) + (1 | obs) else cbind(y, m - y) ~ trt + (1 | blk)
  rk <- ntrt + nblk - 1
  df_rows <- nrow(d) - rk
  df_trials <- if (olre) df_rows else m * nrow(d) - rk
  p <- t(replicate(reps, {
    b <- rnorm(nblk)
    d$y <- rbinom(nrow(d), m, plogis(b[d$blk] + rnorm(nrow(d), 0, od)))
    g <- tryCatch(suppressWarnings(suppressMessages(
      glmer(f, data = d, family = binomial, nAGQ = if (olre) 1 else 9))), error = function(e) NULL)
    if (is.null(g)) c(NA, NA, NA) else {
      if (!exists("checked", inherits = FALSE)) {
        stopifnot(pm_inner_outer_df(g)[["trt2"]] == df_trials)   # the package's df
        checked <<- TRUE
      }
      z <- coef(summary(g))["trt2", "z value"]
      c(2 * pt(-abs(z), df_rows), 2 * pt(-abs(z), df_trials), 2 * pnorm(-abs(z)))
    }
  }))
  rm(list = intersect("checked", ls(globalenv())), envir = globalenv())
  p <- p[stats::complete.cases(p), , drop = FALSE]
  data.frame(design = sprintf("%d trt x %d blocks, %d trials/row%s%s", ntrt, nblk, m,
                              if (od > 0) sprintf(", row SD %.1f", od) else "", if (olre) ", fit (1 | obs)" else ""),
             df_rows, df_package = df_trials, fits = nrow(p), rows = mean(p[, 1] < 0.05),
             package = mean(p[, 2] < 0.05), z = mean(p[, 3] < 0.05))
}
res <- rbind(one_case(4, 3, 25), one_case(4, 3, 5), one_case(4, 6, 2), one_case(2, 3, 10),
             one_case(4, 3, 25, od = 0.3, reps = 400), one_case(4, 3, 25, od = 0.3, olre = TRUE, reps = 400))
print(res, digits = 3)
