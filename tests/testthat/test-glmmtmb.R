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

test_that("glmmTMB random-effect covariance from VarCorr(): equals lmer's; blocks in the order of Z", {
  skip_if_not_installed("glmmTMB")
  # a correlated term, not singular, so both fits reach the same optimum
  f <- Reaction ~ Days + (Days | Subject)
  m1 <- lme4::lmer(f, data = lme4::sleepstudy, REML = FALSE)
  m2 <- glmmTMB::glmmTMB(f, data = lme4::sleepstudy, REML = FALSE)
  G1 <- sigma(m1)^2 * as.matrix(Matrix::tcrossprod(lme4::getME(m1, "Lambda")))
  G2 <- sigma(m2)^2 * as.matrix(Matrix::tcrossprod(predictmeans:::pm_tmb_lambda(m2)))
  expect_equal(G2, G1, tolerance = 2e-3, ignore_attr = TRUE)
  # two terms: each term's block of G is its VarCorr() matrix, at its columns
  # of Z, and the terms do not covary (checked against glmmTMB itself, as
  # lmer and glmmTMB versions can disagree on this nearly singular fit)
  m3 <- suppressWarnings(glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject) + (1 | Days), data = lme4::sleepstudy))
  re <- predictmeans:::pm_re_terms(m3)
  G3 <- sigma(m3)^2 * as.matrix(Matrix::tcrossprod(predictmeans:::pm_tmb_lambda(m3, re)))
  expect_identical(ncol(G3), ncol(lme4::getME(m3, "Z")))
  vc <- glmmTMB::VarCorr(m3)$cond
  for (j in seq_along(re)) {
    nc <- length(re[[j]]$comps)
    lev1 <- re[[j]]$cols[seq_len(nc)]
    expect_equal(G3[lev1, lev1], as.matrix(vc[[j]]), tolerance = 1e-8, ignore_attr = TRUE)
  }
  expect_true(all(G3[re[[1]]$cols, re[[2]]$cols] == 0))
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

test_that("varcomp() for glmmTMB: lmer's variances and SEs, Wald CIs, fit left as it was", {
  skip_if_not_installed("glmmTMB")
  f <- Reaction ~ Days + (Days | Subject)
  m1 <- lme4::lmer(f, data = lme4::sleepstudy)
  m2 <- glmmTMB::glmmTMB(f, data = lme4::sleepstudy, REML = TRUE)
  v1 <- suppressMessages(varcomp(m1))
  v2 <- varcomp(m2)
  expect_identical(rownames(v2), rownames(v1))
  expect_equal(v2[, "vcov"], v1[, "vcov"], tolerance = 1e-4)
  expect_equal(v2[, "SE"], v1[, "SE"], tolerance = 1e-3)
  # variances, with ci_method = "Wald": glmmTMB's Wald CIs for the SDs, squared
  ci_sd <- confint(m2, parm = "theta_")[1:2, 1:2]
  v3 <- varcomp(m2, ci_method = "Wald")
  expect_equal(unname(as.matrix(v3[c(1, 3), 3:4])), unname(ci_sd^2), tolerance = 1e-4)
  expect_true(all(v2[, 3] < v2[, "vcov"] & v2[, "vcov"] < v2[, 4]))
  # the fitted model is unchanged
  expect_equal(as.numeric(glmmTMB::VarCorr(m2)$cond$Subject[1, 1]), v2[1, "vcov"], tolerance = 1e-6)
  expect_identical(colnames(varcomp(m2, ci = FALSE)), c("vcov", "SE"))
})
