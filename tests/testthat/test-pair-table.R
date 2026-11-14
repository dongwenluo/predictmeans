## One pair table (pm_pair_table(), R/predictmeans-stages.R) feeds the pairs
## table, the printed p-value matrices and the letters, so they agree.

conf_fit <- function() {
  # drop "fine" only and always with v3: the v3 means cannot be estimated
  set.seed(9)
  dc <- expand.grid(C = factor(c("c1", "c2")), vol = factor(c("v1", "v2", "v3", "v4")), rep = 1:3)
  dc$drop <- factor(ifelse(dc$vol == "v3", "fine", "coarse"))
  dc$y <- as.numeric(dc$vol) + as.numeric(dc$C) + rnorm(nrow(dc))
  lm(y ~ C*vol + drop, data = dc)
}

test_that("the printed p-value matrix and the pairs table hold the same adjusted p-values", {
  fm <- lm(yield ~ Block + nitro*Variety, data = transform(as.data.frame(nlme::Oats), nitro = factor(nitro)))
  for (a in c("tukey", "holm", "scheffe")) {
    p <- predictmeans(fm, "nitro:Variety", adj = a, plot = FALSE, prtnum = FALSE)
    P <- p$`Pairwise p-value`
    pr <- p$pairs
    expect_equal(P[cbind(pr$level2, pr$level1)], round(pr$p.adj, 4), label = a)
  }
})

test_that("with means that cannot be estimated, the matrices, pairs and limits agree", {
  fm <- conf_fit()
  pa <- suppressWarnings(predictmeans(fm, "C:vol", atvar = "C", adj = "tukey", plot = FALSE, prtnum = FALSE))
  pr <- pa$pairs
  # within each level of C three means take part in the estimable pairs (v3 cannot be estimated): k = 3
  w <- which(pr$by == "c1" & !is.na(pr$estimate))
  expect_equal(pr$p.adj[w], ptukey(sqrt(2)*abs(pr$statistic[w]), 3, pr$df[w], lower.tail = FALSE))
  M <- pa$p_valueMatrix[[1]]
  l1 <- sub("^c1:", "", pr$level1[w]); l2 <- sub("^c1:", "", pr$level2[w])
  expect_equal(M[cbind(l2, l1)], pr$p.adj[w])
  for (a in c("tukey", "scheffe")) {
    q <- suppressWarnings(predictmeans(fm, "C:vol", atvar = "C", adj = a, plot = FALSE, prtnum = FALSE))$pairs
    ok <- !is.na(q$p.adj)
    expect_identical(q$adj.conf.low[ok] > 0 | q$adj.conf.high[ok] < 0, q$p.adj[ok] < 0.05, label = a)
  }
})

test_that("family sizes: k means in the estimable pairs, Scheffe rank of their contrasts", {
  # 4 means, the pairs with mean 3 not estimable: k = 3, rank 2
  ij <- cbind(c(1, 1, 2, 1, 2, 3), c(2, 3, 3, 4, 4, 4))
  L <- t(apply(ij, 1, function(x) { v <- numeric(4); v[x[1]] <- 1; v[x[2]] <- -1; v }))
  est <- !(ij[, 1] == 3 | ij[, 2] == 3)
  f <- pm_families(rep("all", 6), est, ij = ij, L = L)[[1]]
  expect_equal(c(f$m, f$k, f$r), c(3, 3, 2))
  # pairs of 12 means whose contrasts have rank 5 (as emmeans' "scheffe method with rank 5")
  X <- cbind(diag(4)[rep(1:4, 3), ], diag(3)[rep(1:3, each = 4), ])
  pp <- pm_pair_index(12)
  f2 <- pm_families(rep("all", 66), rep(TRUE, 66), ij = cbind(pp$j, pp$i), L = pm_pair_matrix(12) %*% X)[[1]]
  expect_equal(c(f2$k, f2$r), c(12, 5))
})

test_that("with a permlist the letters use 'level', as without", {
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fm <- lm(yield ~ Block + nitro*Variety, data = O)
  pl <- suppressWarnings(permmodels(fm, nperm = 99, seed = 1, prt = FALSE))
  p <- predictmeans(fm, "nitro", permlist = pl, adj = "holm", level = 0.1, plot = FALSE, prtnum = FALSE)
  expect_equal(p$settings$level, 0.1)
  expect_match(p$mean_table$Note, "level '0.1'")
  expect_true(all(p$pairs$p.adj < 0.1) && any(p$pairs$p.adj >= 0.05))   # the level matters here
  expect_identical(anyDuplicated(p$means$group), 0L)                     # every mean its own letter at 0.1
})
