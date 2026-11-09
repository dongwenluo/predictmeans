test_that("R2_model(): Piepho values as numbers equal R2_glmm()", {
  fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  g <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  for (m in list(fm, g)) {
    txt <- suppressMessages(capture.output(r <- R2_glmm(m)))
    v <- R2_model(m)
    expect_s3_class(v, "pm_r2")
    expect_equal(names(v), names(r))
    expect_equal(round(100*unclass(v), 2), as.numeric(sub("%", "", r)), ignore_attr = TRUE)
  }
})

test_that("R2_model(method = 'nakagawa'): published sleepstudy values; random slopes as Johnson (2014)", {
  fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  v <- R2_model(fm, "nakagawa")
  # performance::r2_nakagawa(): marginal 0.279, conditional 0.799
  expect_equal(round(unclass(v), 3), c(Marginal = 0.279, Conditional = 0.799), ignore_attr = TRUE)
  X <- getME(fm, "X")
  S <- VarCorr(fm)$Subject
  vr <- mean(rowSums((X %*% S) * X))
  expect_equal(attr(v, "variances")[["random"]], vr, tolerance = 1e-10)
  expect_equal(attr(v, "variances")[["observation"]], sigma(fm)^2)
  # same model fitted by glmmTMB
  ft <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  expect_equal(unclass(R2_model(ft, "nakagawa")), unclass(v), tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("R2_model(method = 'nakagawa'): latent variance for binomial, error for other families", {
  g <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  v <- R2_model(g, "nakagawa")
  expect_equal(attr(v, "variances")[["observation"]], pi^2/3)
  vf <- var(drop(getME(g, "X") %*% fixef(g)))
  vr <- VarCorr(g)$herd[1, 1]
  expect_equal(unclass(v)[["Marginal"]], vf/(vf + vr + pi^2/3), tolerance = 1e-10)
  gp <- glmer(incidence ~ period + (1 | herd), data = cbpp, family = poisson)
  expect_error(R2_model(gp, "nakagawa"), "gaussian and binomial")
})
