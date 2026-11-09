# glmmTMB >= 1.1.12 stopped storing Lambdat, Lind and lower in its
# random-effect terms and renamed the dispersion parameter (betad ->
# betadisp). These tests cover the code that handles both versions.

test_that("R2_glmm() on a gaussian glmmTMB fit equals the lmer result", {
  skip_if_not_installed("glmmTMB")
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  f <- yield ~ nitro + Variety + (1 | Block/Variety)
  r_lmer <- R2_glmm(lme4::lmer(f, data = Oats))
  # gaussian glmmTMB has a dispersion parameter (betad / betadisp)
  r_tmb <- R2_glmm(glmmTMB::glmmTMB(f, data = Oats, REML = TRUE))
  pct <- function(x) as.numeric(sub("%", "", x))
  expect_equal(pct(r_tmb[1:3]), pct(r_lmer[1:3]), tolerance = 0.01)
})

test_that("glmmTMB random-effect covariance from VarCorr() equals lmer's, correlated terms too", {
  skip_if_not_installed("glmmTMB")
  f <- Reaction ~ Days + (Days | Subject) + (1 | Days)
  m1 <- suppressMessages(lme4::lmer(f, data = lme4::sleepstudy, REML = FALSE))
  m2 <- glmmTMB::glmmTMB(f, data = lme4::sleepstudy, REML = FALSE)
  G1 <- sigma(m1)^2 * as.matrix(Matrix::tcrossprod(lme4::getME(m1, "Lambda")))
  G2 <- sigma(m2)^2 * as.matrix(Matrix::tcrossprod(predictmeans:::pm_tmb_lambda(m2)))
  expect_equal(G2, G1, tolerance = 2e-3, ignore_attr = TRUE)
  # a singular covariance (zero variance) gives a factor instead of an error
  expect_equal(predictmeans:::pm_chol_lower(matrix(c(4, 2, 2, 1), 2)), matrix(c(2, 1, 0, 0), 2))
})

test_that("R2_glmm(): glmmTMB with correlated terms equals lmer; components split by column", {
  skip_if_not_installed("glmmTMB")
  pct <- function(x) as.numeric(sub("%", "", x))
  f <- Reaction ~ Days + (Days | Subject)
  m1 <- lme4::lmer(f, data = lme4::sleepstudy, REML = FALSE)
  m2 <- glmmTMB::glmmTMB(f, data = lme4::sleepstudy, REML = FALSE)
  invisible(capture.output(r1 <- R2_glmm(m1), r2 <- R2_glmm(m2)))
  expect_identical(names(r1), c("Total", "Fixed", "Random", "Subject", "Subject.Days"))
  expect_identical(names(r2), names(r1))
  expect_equal(pct(r2), pct(r1), tolerance = 1e-3)
  # components by hand: sigma_k^2 tr(Z_k' P Z_k), plus the covariance part,
  # make up the random-effect total
  d <- lme4::sleepstudy
  n <- nrow(d)
  P <- diag(n) - 1/n
  Zi <- model.matrix(~ 0 + Subject, d)
  Zs <- Zi * d$Days
  v <- as.matrix(lme4::VarCorr(m1)$Subject)
  parts <- c(v[1, 1] * sum(diag(t(Zi) %*% P %*% Zi)), v[2, 2] * sum(diag(t(Zs) %*% P %*% Zs)),
             2 * v[1, 2] * sum(diag(t(Zi) %*% P %*% Zs)))
  expect_equal(pct(r1[4:5]) / pct(r1[3]), parts[1:2] / sum(parts), tolerance = 1e-3)
})

test_that("R2_glmm(): glmmTMB with zero-inflation random effects or aliased columns", {
  skip_if_not_installed("glmmTMB")
  set.seed(1)
  d <- data.frame(g = factor(rep(1:10, each = 20)), x = rnorm(200))
  d$y <- rpois(200, exp(0.5 + 0.3 * d$x + rnorm(10)[d$g])) * rbinom(200, 1, 0.8)
  m <- glmmTMB::glmmTMB(y ~ x + (1 | g), ziformula = ~ (1 | g), family = poisson, data = d)
  invisible(capture.output(r <- R2_glmm(m)))
  expect_true(all(is.finite(as.numeric(sub("%", "", r)))))
  # B is a copy of A: glmmTMB drops its columns; same R2 as without B
  set.seed(1)
  da <- data.frame(g = factor(rep(1:6, each = 8)), A = factor(rep(1:4, 12)))
  da$y <- rnorm(48) + rnorm(6)[da$g]
  da$B <- da$A
  ma <- suppressMessages(glmmTMB::glmmTMB(y ~ A + B + (1 | g), data = da))
  m0 <- glmmTMB::glmmTMB(y ~ A + (1 | g), data = da)
  invisible(capture.output(ra <- R2_glmm(ma), r0 <- R2_glmm(m0)))
  expect_identical(ra, r0)
})

test_that("semireg_tmb() with predict_info and semipred() run, with R2_glmm()", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("HRW")
  data(WarsawApts, package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  fit <- suppressWarnings(semireg_tmb(areaPerMzloty ~ construction.date,
           smoothZ = list(grp = smZ(construction.date, k = 25)),
           data = W, prt = FALSE))
  sp <- semipred(fit, covariate = "construction.date", prt = FALSE)$pred_df
  expect_true(all(is.finite(sp$SE)) && all(sp$SE > 0))
  r2 <- R2_glmm(fit$semer)
  expect_true(all(is.finite(as.numeric(sub("%", "", r2)))))
  # with a correlated random term as well
  W$surf <- as.numeric(scale(W$surface))
  fit2 <- suppressWarnings(semireg_tmb(areaPerMzloty ~ construction.date + surf + (surf | district),
           smoothZ = list(grp = smZ(construction.date, k = 10)), data = W, prt = FALSE))
  sp2 <- semipred(fit2, covariate = "construction.date", prt = FALSE)$pred_df
  expect_true(all(is.finite(sp2$SE)) && all(sp2$SE > 0))
})
