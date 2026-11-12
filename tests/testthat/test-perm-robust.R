## Heteroscedasticity-robust (HC3) Freedman-Lane permutations: permmodels(robust = TRUE)

robust_data <- function() {
  set.seed(2)
  d <- data.frame(A = factor(rep(1:3, c(6, 6, 18))))
  d$x <- rnorm(30)
  d$y <- d$x + c(3, 1, 1)[d$A] * rnorm(30)
  d
}

## HC3 covariance by hand
hc3 <- function(fm) {
  X <- model.matrix(fm)
  XtXi <- solve(crossprod(X))
  h <- hatvalues(fm)
  XtXi %*% crossprod(X * (residuals(fm) / (1 - h)), X * (residuals(fm) / (1 - h))) %*% XtXi
}

test_that("robust tables: HC3 standard errors and Wald F tests", {
  d <- robust_data()
  fm <- lm(y ~ x + A, data = d)
  V <- hc3(fm)
  r3 <- permmodels(fm, nperm = 99, type = "III", seed = 1, robust = TRUE, prt = FALSE)
  expect_equal(unname(r3$COEFFICENTS[, "Robust SE"]), unname(round(sqrt(diag(V)), 4)))
  b <- coef(fm)[c("A2", "A3")]
  FA <- drop(b %*% solve(V[c("A2", "A3"), c("A2", "A3")], b)) / 2
  expect_equal(r3$ANOVA["A", "F value"], round(FA, 4))
  # type II, as car::Anova(white.adjust = "hc3")
  r2 <- permmodels(fm, nperm = 99, type = "II", seed = 1, robust = TRUE, prt = FALSE)
  expect_equal(r2$ANOVA[, "F value"], c(8.8503, 3.8582))
  expect_match(attr(r2$ANOVA, "heading"), "robust \\(HC3\\)")
  expect_true(isTRUE(r2$perm$robust))
  # type I runs (sequential hypotheses)
  r1 <- permmodels(fm, nperm = 99, seed = 1, robust = TRUE, prt = FALSE)
  expect_equal(rownames(r1$ANOVA), c("x", "A"))
  expect_output(print(r1$perm), "robust HC3")
})

test_that("robust p-value of a coefficient by hand (given permutations)", {
  d <- robust_data()
  fm <- lm(y ~ x + A, data = d)
  set.seed(5)
  P <- replicate(19, sample(30))
  r <- permmodels(fm, perms = P, type = "III", robust = TRUE, prt = FALSE)
  X <- model.matrix(fm)
  red <- lm(d$y ~ X[, -2] - 1)                     # reduced model without x
  tobs <- abs(coef(fm)["x"]) / sqrt(hc3(fm)["x", "x"])
  tperm <- apply(P, 2, function(i) {
    ys <- fitted(red) + residuals(red)[i]
    f <- lm(ys ~ X - 1)
    abs(coef(f)[2]) / sqrt(hc3(f)[2, 2])
  })
  expect_equal(r$COEFFICENTS["x", "Perm_p_value"], round((sum(round(tperm, 6) >= round(tobs, 6)) + 1) / 20, 4))
})

test_that("pairs and contrasts with a robust permlist use robust t and robust limits", {
  d <- robust_data()
  fm <- lm(y ~ x + A, data = d)
  r <- permmodels(fm, nperm = 199, seed = 1, robust = TRUE, prt = FALSE)
  pr <- predictmeans(fm, "A", permlist = r, plot = FALSE, prtnum = FALSE, level = 0.1)$pairs
  # limits: estimate +/- critical value x robust SE of the pair
  V <- hc3(fm)
  L <- rbind(c(0, 0, -1, 0), c(0, 0, 0, -1), c(0, 0, 1, -1))
  se <- sqrt(diag(L %*% V %*% t(L)))
  half <- (pr$conf.high - pr$conf.low) / 2
  ratio <- half / se
  expect_true(all(ratio > 1 & ratio < 4))
  expect_equal(pr$conf.high - pr$estimate, pr$lsd)
  # contrastmeans() gives the same limits for the same pairs
  cn <- rbind(c(-1, 1, 0), c(-1, 0, 1), c(0, -1, 1))
  ct <- contrastmeans(fm, "A", ctrmatrix = cn, permlist = r, level = 0.1)[[1]]
  expect_equal(unname(abs(ct[, "upper.CL"] - ct[, "lower.CL"])), round(pr$conf.high - pr$conf.low, 4), tolerance = 1e-3)
  # joint F of the contrasts: the robust Wald F of A
  cj <- contrastmeans(fm, "A", ctrmatrix = cn, permlist = r, joint = TRUE)
  expect_true(cj[[3]]$Perm_p_value > 0 && cj[[3]]$Perm_p_value <= 1)
})

test_that("robust = TRUE: glm on the working scale; errors for mixed models and raw", {
  d <- robust_data()
  set.seed(3)
  d$count <- rpois(30, exp(1 + 0.3 * d$x))
  g <- glm(count ~ x + A, family = poisson, data = d)
  rg <- permmodels(g, nperm = 49, type = "III", seed = 1, robust = TRUE, prt = FALSE)
  expect_match(attr(rg$ANOVA, "heading"), "working")
  expect_true(all(rg$ANOVA$Perm_p_value > 0 & rg$ANOVA$Perm_p_value <= 1))
  fm <- lm(y ~ x + A, data = d)
  expect_error(permmodels(fm, nperm = 9, robust = TRUE, method = "raw", prt = FALSE), "freedman-lane")
  m <- lme4::lmer(y ~ x + (1 | A), data = d)
  expect_error(permmodels(m, nperm = 9, robust = TRUE, prt = FALSE), "lm, aov and glm")
})
