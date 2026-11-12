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

  # glm models are permuted on the working scale (test-perm-fl-classes.R)
  expect_s3_class(permmodels(glm(y ~ A, data = d), nperm = 9, prt = FALSE, method = "freedman-lane")$perm, "pm_perm")
  fl2 <- lm(cbind(y, x) ~ A, data = d)
  expect_error(permmodels(fl2, nperm = 9, method = "freedman-lane"), "single response")
})

test_that("Freedman-Lane p-values of pairs and contrasts equal direct refits", {
  d <- fl_data()
  fm <- lm(y ~ x + A, data = d)
  r <- permmodels(fm, nperm = 99, seed = 5, prt = FALSE, method = "freedman-lane")
  expect_s3_class(r$perm, "pm_perm")
  idx <- fl_idx(5, 99, 24)
  # pair 2 v 3: the reduced model merges levels 2 and 3 of A
  dd <- d; dd$A23 <- factor(ifelse(d$A == "1", "1", "23"))
  red <- lm(y ~ x + A23, data = dd)
  f0 <- fitted(red); e0 <- resid(red)
  tobs <- abs(coef(summary(lm(y ~ x + relevel(A, "2"), data = d)))["relevel(A, \"2\")3", "t value"])
  ts <- apply(idx, 2, function(i) {
    dd$y <- f0 + e0[i]
    abs(coef(summary(lm(y ~ x + relevel(A, "2"), data = dd)))["relevel(A, \"2\")3", "t value"])
  })
  p23 <- (sum(round(ts, 6) >= round(tobs, 6)) + 1) / 100
  L <- rbind(c(0, 0, 1, 0), c(0, 0, -1, 1))       # 1 v 2 and 2 v 3 on (Intercept, x, A2, A3)
  expect_equal(pm_perm_fl_t(r$perm, L, fm)[2], p23)
  # 1 v 2 is the coefficient A2: the same test as its t test in the table
  expect_equal(pm_perm_fl_t(r$perm, L, fm)[1], unname(r$COEFFICENTS["A2", "Perm_p_value"]), tolerance = 1e-4)

  # predictmeans() and contrastmeans() use it
  pm <- predictmeans(fm, "A", permlist = r, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_identical(pm$pairs$p.method, rep("permutation", 3))
  expect_identical(pm$settings$nperm, 99L)
  key <- paste(pm$pairs$level1, pm$pairs$level2)
  expect_equal(pm$pairs$p.value[key == "2 3"], p23)
  cm <- contrastmeans(fm, "A", ctrmatrix = rbind(c(0, -1, 1)), permlist = r)
  expect_equal(unname(cm[[1]][1, "Permuted Pr(>|t|)"]), round(p23, 4))
  # a different model
  expect_error(predictmeans(lm(y ~ A, data = d), "A", permlist = r, plot = FALSE, prtnum = FALSE),
               "different model")
})

test_that("with a permlist, observed t statistics use the unadjusted covariance", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  Om <- Oats[-c(1, 5, 9, 20, 33), ]
  m <- lme4::lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Om)
  pl <- permmodels(m, nperm = 9, prt = FALSE, ncore = 1, seed = 1)
  a <- predictmeans(m, "nitro", permlist = pl, plot = FALSE, prtnum = FALSE, ncore = 1)
  b <- predictmeans(m, "nitro", vcov_adjust = FALSE, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_equal(a$pairs$std.error, b$pairs$std.error)
  expect_false(a$settings$vcov_adjust)
  expect_warning(predictmeans(m, "nitro", permlist = pl, vcov_adjust = TRUE, plot = FALSE, prtnum = FALSE, ncore = 1),
                 "unadjusted covariance")
})
