# Kmatrix() on covariates that enter the model in complex ways. emmeans is an
# independent implementation of the same reference-grid idea (equal weights
# over factor levels, covariates at their means), so its linear functions
# must agree with K row by row.

sim_data <- function() {
  set.seed(1)
  d <- expand.grid(A = factor(c("a1","a2","a3")), B = factor(c("b1","b2")),
                   blk = factor(1:4), rep = 1:2)
  n <- nrow(d)
  d$x <- runif(n, 1, 10)
  d$z <- rnorm(n)
  d$y <- 10 + as.numeric(d$A) + 0.5*as.numeric(d$B) + sin(d$x) + 0.3*d$z +
    rnorm(4)[d$blk] + rnorm(n)
  d$ch <- as.character(d$B)
  d$cnt <- rpois(n, 3)
  d$nf <- rep(c(0, 10, 20), length = n)
  d
}

expect_K_matches_emmeans <- function(model, modelterm, specs, ...) {
  K <- suppressWarnings(Kmatrix(model, modelterm)$K)
  em <- emmeans::emmeans(model, specs, ...)
  L <- em@linfct
  rownames(L) <- do.call(paste, c(lapply(em@grid[all.vars(specs)], as.character), sep = ":"))
  expect_equal(unclass(K)[, colnames(K)], L[rownames(K), colnames(K)],
               ignore_attr = TRUE, tolerance = 1e-12)
}

test_that("transformed covariates match emmeans", {
  skip_if_not_installed("emmeans")
  d <- sim_data()
  expect_K_matches_emmeans(lm(y ~ A + poly(x, 2), d), "A", ~A)
  expect_K_matches_emmeans(lm(y ~ A + splines::ns(x, 3), d), "A", ~A)
  expect_K_matches_emmeans(lm(y ~ A + log(x), d), "A", ~A)
  expect_K_matches_emmeans(lm(y ~ A + scale(x), d), "A", ~A)
  expect_K_matches_emmeans(lm(y ~ A + I(x^2) + x, d), "A", ~A)
})

test_that("covariates that only appear in interactions match emmeans", {
  skip_if_not_installed("emmeans")
  d <- sim_data()
  expect_K_matches_emmeans(lm(y ~ A + A:x, d), "A", ~A)
  expect_K_matches_emmeans(lm(y ~ A*B + x:z, d), "A:B", ~A:B)
  expect_K_matches_emmeans(lm(y ~ A*B + x:z, d), "B", ~B)
})

test_that("character predictors, factor() terms and offsets match emmeans", {
  skip_if_not_installed("emmeans")
  d <- sim_data()
  expect_K_matches_emmeans(lm(y ~ A + ch + x, d), "A", ~A)
  expect_K_matches_emmeans(lm(y ~ factor(nf) + x, d), "nf", ~nf)
  expect_K_matches_emmeans(glm(cnt ~ A + x + offset(log(x)), family = poisson, d), "A", ~A)
  m <- lm(y ~ factor(nf) + x, d)
  expect_equal(Kmatrix(m, "factor(nf)")$K, Kmatrix(m, "nf")$K, ignore_attr = TRUE)
})

test_that("mixed models with transformed covariates match emmeans", {
  skip_if_not_installed("emmeans")
  d <- sim_data()
  expect_K_matches_emmeans(suppressMessages(lmer(y ~ A + poly(x, 2) + (1|blk), d)), "A", ~A)
  expect_K_matches_emmeans(suppressMessages(lmer(y ~ A + A:x + (1|blk), d)), "A", ~A)
  expect_K_matches_emmeans(lme(y ~ A + poly(x, 2), random = ~1|blk, data = d), "A", ~A)
  expect_K_matches_emmeans(glmmTMB(y ~ A + poly(x, 2) + (1|blk), data = d), "A", ~A)
})

test_that("varied covariates and cov_reduce match emmeans", {
  skip_if_not_installed("emmeans")
  d <- sim_data()
  m <- lm(y ~ A*poly(x, 2), d)
  K <- Kmatrix(m, "A", "x", covariateV = c(2, 5, 8))$K
  em <- emmeans::emmeans(m, ~A|x, at = list(x = c(2, 5, 8)))
  L <- em@linfct
  rownames(L) <- paste(em@grid$A, em@grid$x, sep = ":")
  expect_equal(unclass(K), L[rownames(K), colnames(K)], ignore_attr = TRUE, tolerance = 1e-12)

  m2 <- lm(y ~ A + poly(x, 2), d)
  K2 <- Kmatrix(m2, "A", cov_reduce = median)$K
  L2 <- emmeans::emmeans(m2, ~A, cov.reduce = median)@linfct
  expect_equal(unclass(K2), L2, ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("models fitted without a data argument or with a transformed response work", {
  d <- sim_data()
  m1 <- with(d, lm(y ~ A + poly(x, 2)))
  m2 <- lm(y ~ A + poly(x, 2), d)
  expect_equal(Kmatrix(m1, "A")$K, Kmatrix(m2, "A")$K)
  m3 <- lm(log(y) ~ A + x, d)
  expect_equal(Kmatrix(m3, "A", data = d)$K, Kmatrix(m3, "A")$K)
  expect_identical(Kmatrix(m3, "A")$response, "log(y)")
})

test_that("cov_reduce = 'average' averages over the observed covariate values", {
  d <- sim_data()
  m <- lm(y ~ A + poly(x, 2), d)
  K <- Kmatrix(m, "A", cov_reduce = "average")$K
  TT <- delete.response(terms(m))
  manual <- t(sapply(levels(d$A), function(a) {
    nd <- transform(d, A = factor(a, levels = levels(d$A)))
    colMeans(model.matrix(TT, model.frame(TT, nd)))
  }))
  expect_equal(unclass(K), manual, ignore_attr = TRUE, tolerance = 1e-12)
  expect_error(Kmatrix(m, "A", cov_reduce = "median"), "cov_reduce")
})

test_that("non-estimable rows are flagged", {
  d <- sim_data()
  dm <- d[!(d$A == "a3" & d$B == "b2"), ]
  m <- lm(y ~ A*B, dm)
  K <- suppressWarnings(Kmatrix(m, "A:B")$K)
  expect_identical(unname(attr(K, "estimable")), c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_true(all(attr(Kmatrix(lm(y ~ A*B, d), "A:B")$K, "estimable")))
})

test_that("missing cells are averaged over the observed cells only", {
  set.seed(4)
  d <- expand.grid(nf = c(0, 10, 20), A = factor(c("a1","a2")), rep = 1:3)
  d$y <- rnorm(nrow(d)) + d$nf/10
  d <- d[!(d$nf == 20 & d$A == "a1"), ]
  m <- lm(y ~ factor(nf)*A, d)
  K <- suppressWarnings(Kmatrix(m, "A")$K)
  cells <- unique(d[c("nf", "A")])
  TT <- delete.response(terms(m))
  X <- model.matrix(TT, model.frame(TT, cells, xlev = list(A = levels(d$A))))
  ref <- t(sapply(levels(d$A), function(a) colMeans(X[cells$A == a, , drop = FALSE])))
  expect_equal(unclass(K)[, colnames(ref)], ref, ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("a matrix covariateV gives joint prediction points", {
  d <- sim_data()
  m <- lm(y ~ A + x*z, d)
  pts <- cbind(x = c(2, 5, 8), z = c(-1, 0, 1))
  kk <- Kmatrix(m, "A", c("x", "z"), covariateV = pts)
  expect_equal(nrow(kk$K), 3 * 3)
  expect_equal(unique(kk$fctnames[c("x", "z")]), as.data.frame(pts), ignore_attr = TRUE)
  expect_equal(unname(kk$K[, "x:z"]), rep(pts[, "x"] * pts[, "z"], times = 3))
})

test_that("mgcv smooths use the lpmatrix and exclude random-effect smooths", {
  skip_if_not_installed("mgcv")
  d <- sim_data()
  g1 <- mgcv::gam(y ~ A + s(x), data = d)
  lp <- predict(g1, newdata = data.frame(A = levels(d$A), x = mean(d$x)), type = "lpmatrix")
  expect_equal(unclass(Kmatrix(g1, "A")$K), lp, ignore_attr = TRUE, tolerance = 1e-12)

  lpa <- t(sapply(levels(d$A), function(a) {
    colMeans(predict(g1, newdata = transform(d, A = factor(a, levels = levels(d$A))), type = "lpmatrix"))
  }))
  expect_equal(unclass(Kmatrix(g1, "A", cov_reduce = "average")$K), lpa,
               ignore_attr = TRUE, tolerance = 1e-12)

  g2 <- mgcv::gam(y ~ A*B + s(x, by = A) + s(blk, bs = "re"), data = d)
  K2 <- Kmatrix(g2, "A:B")$K
  expect_true(all(K2[, grep("^s\\(blk\\)", colnames(K2))] == 0))
  K3 <- Kmatrix(g2, "A:B", exclude = character(0))$K
  expect_true(any(K3[, grep("^s\\(blk\\)", colnames(K3))] != 0))
  expect_error(Kmatrix(g2, "A", exclude = "s(nothing)"), "Unknown smooth")
  expect_error(Kmatrix(g2, "blk"), "excluded")
})

test_that("predictmeans() works for gam models and matches predict.gam", {
  skip_if_not_installed("mgcv")
  d <- sim_data()
  g <- mgcv::gam(y ~ A*B + s(x, by = A) + s(blk, bs = "re"), data = d)
  nd <- expand.grid(A = levels(d$A), B = levels(d$B))
  nd$x <- mean(d$x)
  nd$blk <- d$blk[1]
  pg <- predict(g, newdata = nd, se.fit = TRUE, exclude = "s(blk)")
  pm <- predictmeans(g, "A:B", plot = FALSE, prtnum = FALSE)
  mt <- pm$mean_table
  mt <- mt[order(mt$B, mt$A), ]
  expect_equal(mt$Mean, round(as.vector(pg$fit), 4))
  expect_equal(mt$SE, round(as.vector(pg$se.fit), 5))
  expect_equal(mt$Df[1], round(g$df.residual, 2))
})

test_that("predictmeans() works with transformed covariates", {
  skip_if_not_installed("emmeans")
  d <- sim_data()
  m <- suppressMessages(lmer(y ~ A + poly(x, 2) + (1|blk), d))
  pm <- predictmeans(m, "A", plot = FALSE, prtnum = FALSE)
  em <- summary(emmeans::emmeans(m, ~A, lmer.df = "kenward-roger"))
  expect_equal(pm$mean_table$Mean, round(em$emmean, 4))
})

test_that("vectorised Kenward-Roger df equal pbkrtest::Lb_ddf() row by row", {
  data("ATP", package = "predictmeans", envir = environment())
  fm <- lmer(ATP ~ A*B*time + (1|heart), data = ATP)
  Va <- pbkrtest::vcovAdj(fm)
  K <- Kmatrix(fm, "A:time")$K
  pairs <- t(combn(nrow(K), 2, function(ij) K[ij[1], ] - K[ij[2], ]))
  L <- rbind(K, pairs)
  ref <- apply(L, 1, function(l) pbkrtest::Lb_ddf(l, V0 = vcov(fm), Vadj = Va))
  expect_equal(unname(predictmeans:::pm_kr_ddf(L, vcov(fm), Va)), unname(ref), tolerance = 1e-10)
  fm2 <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = transform(as.data.frame(nlme::Oats), nitro = factor(nitro)))
  Va2 <- pbkrtest::vcovAdj(fm2)
  K2 <- Kmatrix(fm2, "nitro:Variety")$K
  ref2 <- apply(K2, 1, function(l) pbkrtest::Lb_ddf(l, V0 = vcov(fm2), Vadj = Va2))
  expect_equal(unname(predictmeans:::pm_kr_ddf(K2, vcov(fm2), Va2)), unname(ref2), tolerance = 1e-10)
})
