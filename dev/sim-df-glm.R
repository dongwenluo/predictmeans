## glm with a fixed dispersion (Poisson, binomial): which reference
## distribution for Wald tests? (dev/review-defaults.md, step 5)
## predictmeans 1.x and 2.0 use t with the residual df (as for lm);
## emmeans, SAS GLIMMIX and Stata use z (chi-square for several df).
## Type I error at 5% of the Wald test of treatment 1 vs 2 (t / z) and of
## all treatments together (F / chi-square), no treatment effects.
## 4000 data sets per design; a true 5% lies in 4.1-5.9% (99%).
## Rscript dev/sim-df-glm.R   (writes dev/sim-df-glm/results.md)
set.seed(20261114)
nsim <- 4000
designs <- list(
  list(name = "Poisson, 3 trt x 4, mean 3", fam = "poisson", k = 3, n = 4, mu = 3, block = FALSE),
  list(name = "Poisson, 3 trt x 4, mean 10", fam = "poisson", k = 3, n = 4, mu = 10, block = FALSE),
  list(name = "Poisson, 3 trt x 10, mean 3", fam = "poisson", k = 3, n = 10, mu = 3, block = FALSE),
  list(name = "Poisson RCBD, 3 trt x 4 blocks, mean 5", fam = "poisson", k = 3, n = 4, mu = 5, block = TRUE),
  list(name = "binomial (10 trials), 3 trt x 4, p 0.3", fam = "binomial", k = 3, n = 4, mu = 0.3, m = 10, block = FALSE),
  list(name = "binomial (10 trials), 3 trt x 10, p 0.3", fam = "binomial", k = 3, n = 10, mu = 0.3, m = 10, block = FALSE),
  list(name = "binary, 3 trt x 20, p 0.3", fam = "binomial", k = 3, n = 20, mu = 0.3, m = 1, block = FALSE),
  list(name = "binary RCBD, 3 trt x 8 blocks (5 units), p 0.3", fam = "binomial", k = 3, n = 8, mu = 0.3, m = 5,
       block = TRUE))
res <- lapply(designs, function(ds) {
  d <- data.frame(trt = factor(rep(seq_len(ds$k), each = ds$n)), blk = factor(rep(seq_len(ds$n), ds$k)))
  form <- if (ds$block) ~ trt + blk else ~ trt
  X <- model.matrix(form, d)
  fam <- if (ds$fam == "poisson") poisson() else binomial()
  bmu <- if (ds$block) rnorm(ds$n, 0, 0.3)[d$blk] else 0      # fixed block effects for the design
  eta <- fam$linkfun(ds$mu) + bmu
  p <- t(replicate(nsim, {
    y <- if (ds$fam == "poisson") rpois(nrow(d), fam$linkinv(eta)) else rbinom(nrow(d), ds$m, fam$linkinv(eta)) / ds$m
    w <- if (ds$fam == "poisson") rep(1, nrow(d)) else rep(ds$m, nrow(d))
    f <- suppressWarnings(glm.fit(X, y, weights = w, family = fam))
    V <- tryCatch(chol2inv(qr.R(f$qr))[order(f$qr$pivot), order(f$qr$pivot)], error = function(e) NULL)
    if (is.null(V) || f$rank < ncol(X)) return(rep(NA, 4))
    b <- f$coefficients
    l <- replace(numeric(ncol(X)), 2, 1)                      # trt2 - trt1
    tv <- sum(l * b) / sqrt(drop(l %*% V %*% l))
    L <- diag(ncol(X))[2:ds$k, , drop = FALSE]
    W <- drop(t(L %*% b) %*% solve(L %*% V %*% t(L), L %*% b))
    dfr <- f$df.residual
    c(2 * pt(-abs(tv), dfr), 2 * pnorm(-abs(tv)), pf(W / (ds$k - 1), ds$k - 1, dfr, lower.tail = FALSE),
      pchisq(W, ds$k - 1, lower.tail = FALSE))
  }))
  r <- round(100 * colMeans(p < 0.05, na.rm = TRUE), 1)
  data.frame(design = ds$name, resid.df = nrow(d) - ncol(X), t.resid = r[1], z = r[2], F.resid = r[3], chisq = r[4],
             fits = sum(!is.na(p[, 1])))
})
out <- do.call(rbind, res)
dir.create("dev/sim-df-glm", showWarnings = FALSE)
writeLines(c("# glm, fixed dispersion: % of true nulls rejected at 5% (4000 data sets per design)", "",
             "| design | residual df | t, residual df (now) | z | F, residual df (now) | chi-square | fits |",
             "|---|---|---|---|---|---|---|",
             apply(out, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |"))),
           "dev/sim-df-glm/results.md")
print(out, row.names = FALSE)
