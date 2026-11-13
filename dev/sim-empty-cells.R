## Empty cells: how should margins be compared when some cells of an
## interaction have no data? (dev/review-defaults.md)
##
## A 3 x 4 factorial A x B, 3 replicates per cell, y ~ A * B (lm), sigma = 1.
## Empty cells: (A1, B4), or (A1, B4) and (A2, B3). True: no A effect. B4
## differs from B1-B3 by delta (a site or year effect). The A1 vs A2
## comparison by
##   present  - predictmeans' default (1.x): each margin averages its own
##              observed cells (A1 over B1-B3, A2 over B1-B4);
##   common   - both margins over the B levels observed for both (Type IV-like);
##   additive - the means of y ~ A + B (right only without interaction).
## Scenario "inter": A1 and A2 differ only at B4, where A1 has no data, so the
## common-cells null still holds but the additive fit is biased.
## 2000 data sets per cell; a true 5% lies in 3.9-6.1% (99%).
## Rscript dev/sim-empty-cells.R   (writes dev/sim-empty-cells/results.md)

suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(20261113)
nsim <- 2000
d0 <- expand.grid(rep = 1:3, B = factor(paste0("B", 1:4)), A = factor(paste0("A", 1:3)))
patterns <- list(one = c("A1:B4"), two = c("A1:B4", "A2:B3"))
res <- list()
for (pat in names(patterns)) {
  d <- d0[!paste(d0$A, d0$B, sep = ":") %in% patterns[[pat]], ]
  d$y <- rnorm(nrow(d))
  fit <- lm(y ~ A * B, data = d)
  fa <- lm(y ~ A + B, data = d)
  # the present rule as predictmeans applies it, from its K matrix
  Kp <- suppressWarnings(Kmatrix(fit, "A")$K)
  Lp <- Kp["A1", ] - Kp["A2", ]
  # common cells: A1 - A2 averaged over the B levels observed for both
  Kc <- Kmatrix(fit, "A:B")$K
  cells <- unique(paste(d$A, d$B, sep = ":"))
  com <- levels(d$B)[paste0("A1:", levels(d$B)) %in% cells & paste0("A2:", levels(d$B)) %in% cells]
  Lc <- colMeans(Kc[paste0("A1:", com), , drop = FALSE] - Kc[paste0("A2:", com), , drop = FALSE])
  Ka <- Kmatrix(fa, "A")$K
  La <- Ka["A1", ] - Ka["A2", ]
  X <- model.matrix(fit); keep <- !is.na(coef(fit)); X <- X[, keep]
  Xa <- model.matrix(fa)
  Lp <- Lp[colnames(X)]; Lc <- Lc[colnames(X)]
  # check against predictmeans() once
  pmp <- suppressWarnings(predictmeans(fit, "A", pairwise = TRUE, plot = FALSE, prtnum = FALSE))$pairs
  r1 <- pmp[pmp$level1 == "A1" & pmp$level2 == "A2", ]
  b <- qr.coef(qr(X), d$y); V <- summary(fit)$sigma^2 * solve(crossprod(X))
  stopifnot(abs(sum(Lp * b) - r1$estimate) < 1e-8, abs(sqrt(drop(Lp %*% V %*% Lp)) - r1$std.error) < 1e-8)
  XtXi <- solve(crossprod(X)); XatXai <- solve(crossprod(Xa))
  dfr <- nrow(X) - ncol(X); dfa <- nrow(Xa) - ncol(Xa)
  B4 <- d$B == "B4"; inter <- d$A == "A2" & d$B == "B4"
  for (scen in c("B4 effect", "inter")) for (delta in c(0, 1, 2, 4)) {
    mu <- if (scen == "B4 effect") delta * B4 else delta * inter
    p <- matrix(NA, nsim, 3, dimnames = list(NULL, c("present", "common", "additive")))
    for (s in seq_len(nsim)) {
      y <- mu + rnorm(nrow(d))
      b <- XtXi %*% crossprod(X, y); s2 <- sum((y - X %*% b)^2) / dfr
      ba <- XatXai %*% crossprod(Xa, y); s2a <- sum((y - Xa %*% ba)^2) / dfa
      tt <- c(sum(Lp * b) / sqrt(s2 * drop(Lp %*% XtXi %*% Lp)),
              sum(Lc * b) / sqrt(s2 * drop(Lc %*% XtXi %*% Lc)),
              sum(La * ba) / sqrt(s2a * drop(La %*% XatXai %*% La)))
      p[s, ] <- 2 * pt(-abs(tt), c(dfr, dfr, dfa))
    }
    res[[length(res) + 1]] <- data.frame(empty = pat, scenario = scen, delta = delta,
                                         t(round(100 * colMeans(p < 0.05), 1)))
  }
}
out <- do.call(rbind, res)
dir.create("dev/sim-empty-cells", showWarnings = FALSE)
writeLines(c("# Empty cells: % of true nulls rejected at 5% (A1 vs A2; 2000 data sets each)", "",
             "| empty cells | scenario | delta | present (1.x default) | common cells | additive model |",
             "|---|---|---|---|---|---|",
             apply(out, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |"))),
           "dev/sim-empty-cells/results.md")
print(out, row.names = FALSE)
