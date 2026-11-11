## permmodels(method = "freedman-lane") for lm and aov models

fl_data <- function() {
  set.seed(1)
  d <- data.frame(A = factor(rep(1:3, each = 8)), B = factor(rep(1:2, 12)))
  d$x <- exp(rnorm(24) + 0.8 * (d$A == "3"))
  d$w <- runif(24, 0.5, 2)
  d$y <- 3 * d$x + rnorm(24)
  d
}

## the permutations permmodels() uses after set.seed(seed)
fl_idx <- function(seed, nperm, n) {
  set.seed(seed)
  vapply(seq_len(nperm), function(i) sample.int(n), integer(n))
}

test_that("Freedman-Lane p-values equal a direct computation", {
  d <- fl_data()
  fm <- lm(y ~ x + A * B, data = d)
  r <- permmodels(fm, nperm = 99, type = "II", seed = 3, prt = FALSE, method = "freedman-lane")
  idx <- fl_idx(3, 99, 24)

  # type II test of A: reduced model without A and A:B
  red <- lm(y ~ x + B, data = d)
  f0 <- fitted(red); e0 <- resid(red)
  Fobs <- pm_anova(fm, "II", "F")["A", "F value"]
  Fs <- apply(idx, 2, function(i) {
    dd <- d; dd$y <- f0 + e0[i]
    pm_anova(lm(y ~ x + A * B, data = dd), "II", "F")["A", "F value"]
  })
  expect_equal(as.numeric(r$ANOVA["A", "Perm_p_value"]), round((sum(round(Fs, 6) >= round(Fobs, 6)) + 1) / 100, 4))

  # t test of the slope: reduced model without that column
  X <- model.matrix(fm)
  f0 <- qr.fitted(qr(X[, colnames(X) != "x"]), d$y); e0 <- d$y - f0
  ts <- apply(idx, 2, function(i) {
    dd <- d; dd$y <- f0 + e0[i]
    coef(summary(lm(y ~ x + A * B, data = dd)))["x", "t value"]
  })
  tobs <- coef(summary(fm))["x", "t value"]
  expect_equal(r$COEFFICENTS["x", "Perm_p_value"],
               round((sum(round(abs(ts), 6) >= round(abs(tobs), 6)) + 1) / 100, 4))
  expect_true(is.na(r$COEFFICENTS["(Intercept)", "Perm_p_value"]))
  expect_null(r$permlist)
})

test_that("Freedman-Lane works with weights and type I tests", {
  d <- fl_data()
  fw <- lm(y ~ x + A, data = d, weights = w)
  r <- permmodels(fw, nperm = 99, type = "I", seed = 3, prt = FALSE, method = "freedman-lane")
  idx <- fl_idx(3, 99, 24)
  red <- lm(y ~ x, data = d, weights = w)
  f0 <- fitted(red); e0 <- resid(red) * sqrt(d$w)
  Fobs <- anova(fw)["A", "F value"]
  Fs <- apply(idx, 2, function(i) {
    dd <- d; dd$y <- f0 + e0[i] / sqrt(d$w)
    anova(lm(y ~ x + A, data = dd, weights = w))["A", "F value"]
  })
  expect_equal(as.numeric(r$ANOVA["A", "Perm_p_value"]), round((sum(round(Fs, 6) >= round(Fobs, 6)) + 1) / 100, 4))
})

test_that("Freedman-Lane: large effects, aov, seeds and unsupported models", {
  d <- fl_data()
  fm <- lm(y ~ x + A, data = d)
  r1 <- permmodels(fm, nperm = 49, type = "III", seed = 5, prt = FALSE, method = "freedman-lane")
  r2 <- permmodels(fm, nperm = 49, type = "III", seed = 5, prt = FALSE, method = "freedman-lane")
  expect_identical(r1, r2)
  expect_equal(as.numeric(r1$ANOVA["x", "Perm_p_value"]), 1 / 50)
  expect_equal(r1$ANOVA["(Intercept)", "Perm_p_value"], "")

  fa <- aov(y ~ x + A, data = d)
  ra <- permmodels(fa, nperm = 49, seed = 5, prt = FALSE, method = "freedman-lane")
  expect_null(ra$COEFFICENTS)
  expect_equal(as.numeric(ra$ANOVA["x", "Perm_p_value"]), 1 / 50)

  expect_error(permmodels(glm(y ~ A, data = d), nperm = 9, method = "freedman-lane"),
               "lm and aov")
  fl2 <- lm(cbind(y, x) ~ A, data = d)
  expect_error(permmodels(fl2, nperm = 9, method = "freedman-lane"), "single response")
})
