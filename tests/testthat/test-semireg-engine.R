## The fitting engine of semireg() and semireg_tmb() (R/semireg-engine.R):
## one fit, the fit's own deviance function, covariance without inverting
## Lambda, the data kept out of the caller's frame.

test_that("semireg(): Satterthwaite df and variance components equal the same model fitted by lme", {
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  fit <- semireg(areaPerMzloty ~ construction.date, smoothZ = list(grp = smZ(construction.date, k = 25)),
                 data = WarsawApts, prt = FALSE)
  # the same mixed model in nlme, as Harezlak, Ruppert and Wand write it
  d <- fit$data
  d$Z <- as.matrix(smZ(d$construction.date, k = 25))
  d$one <- factor(1)
  fl <- nlme::lme(areaPerMzloty ~ construction.date, random = list(one = nlme::pdIdent(~ -1 + Z)), data = d,
                  method = "REML")
  expect_equal(as.numeric(logLik(fit$semer)), as.numeric(logLik(fl)), tolerance = 1e-6)
  # Satterthwaite and Kenward-Roger df with the deviance of the fitted model
  # (before 2.0 lmerTest's used that of the placeholder model: 1871 df instead of 150)
  L1 <- matrix(c(0, 1), 1)
  ddf_lme <- pm_inference(fl, df_method = "Satterthwaite", vcov_adjust = FALSE)$df(L1, "construction.date")
  expect_equal(anova(fit$semer, ddf = "Satterthwaite")["construction.date", "Df.res"], ddf_lme, tolerance = 1e-3)
  kr_lme <- pm_inference(fl, df_method = "Kenward-Roger")$df(L1, "construction.date")
  expect_equal(anova(fit$semer)["construction.date", "Df.res"], kr_lme, tolerance = 1e-3)
  # varcomp(): the SEs of the variance components (before 2.0: SE 7.7 for 1.45)
  v1 <- varcomp(fit)
  v2 <- varcomp(fl)
  expect_equal(v1["residual", "vcov"], v2["residual", "vcov"], tolerance = 1e-4)
  expect_equal(v1["residual", "SE"], v2["residual", "SE"], tolerance = 1e-3)
  expect_equal(v1[1, "SE"], v2[1, "SE"], tolerance = 1e-3)
})

test_that("semireg(): a smooth with its variance at zero (a straight line) and a negative binomial fit", {
  set.seed(3)
  n <- 300
  d <- data.frame(x = runif(n, 0, 3), g = gl(2, n/2, labels = c("a", "b")), id = gl(20, 1, n))
  eta <- 0.3*sin(2*d$x) + 0.5*(d$g == "b")*d$x/3 + rnorm(20, 0, 0.3)[d$id]
  d$yp <- rpois(n, exp(0.5 + eta))
  d$yb <- rbinom(n, 1, plogis(-1 + 2*eta))
  d$ynb <- MASS::rnegbin(n, exp(1 + eta), 2)
  # before 2.0: "Matrix seems negative semi-definite"
  fit <- suppressWarnings(semireg(ynb ~ x*g, smoothZ = list(grp = smZ(x, k = 8, by = g)), data = d, ngbinomial = TRUE,
                                  prt = FALSE))
  expect_s3_class(fit, "semireg")
  expect_true(grepl("Negative Binomial", stats::family(fit$semer)$family))
  expect_true(all(is.finite(as.matrix(fit$CovMat))))
  sp <- semipred(fit, "g", "x", prt = FALSE)$pred_df
  expect_true(all(is.finite(sp$SE_mu)))
  # a variance estimated as zero: its coefficients have zero covariance
  th <- lme4::getME(fit$semer, "theta")
  if (any(th == 0)) {
    cs <- as.matrix(fit$cov_lst$grp)
    expect_true(all(abs(cs[-(1:4), -(1:4)]) < 1e-10))
  }
})

test_that("semireg(): the covariance of the coefficients is (C'WC + D)^-1 when Lambda is invertible", {
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  fit <- semireg(areaPerMzloty ~ construction.date, smoothZ = list(grp = smZ(construction.date, k = 10)),
                 data = WarsawApts, prt = FALSE)
  s <- fit$semer
  C <- cbind(lme4::getME(s, "X"), lme4::getME(s, "Z"))
  L <- lme4::getME(s, "Lambda")
  p <- ncol(lme4::getME(s, "X"))
  D <- Matrix::.bdiag(list(Matrix::Matrix(0, p, p), solve(Matrix::tcrossprod(L))))
  expect_equal(as.matrix(fit$CovMat), as.matrix(solve(crossprod(C) + D)), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("semireg(): cbind() responses, data as an expression, nothing written into the caller's frame", {
  data(cbpp, package = "lme4", envir = environment())
  set.seed(1)
  cbpp$x <- runif(nrow(cbpp))
  # before 2.0 the second column of cbind() was centred: "Response is constant"
  b <- semireg(cbind(incidence, size - incidence) ~ period + x + (1 | herd), smoothZ = list(sx = smZ(x, k = 5)),
               family = binomial, data = cbpp, prt = FALSE)
  expect_s3_class(b, "semireg")
  expect_equal(unname(lme4::fixef(b$semer)[2:4]),
               unname(lme4::fixef(lme4::glmer(cbind(incidence, size - incidence) ~ period + x + (1 | herd),
                                              family = binomial, data = cbpp))[2:4]), tolerance = 0.05)
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  f <- function() {
    before <- ls()
    m <- semireg(areaPerMzloty ~ construction.date, smoothZ = list(grp = smZ(construction.date, k = 10)),
                 data = subset(WarsawApts, surface > 30), prt = FALSE)
    list(m = m, new = setdiff(ls(), c(before, "before", "m")))
  }
  r <- f()
  expect_s3_class(r$m, "semireg")
  expect_length(r$new, 0)
  expect_equal(nrow(r$m$data), sum(WarsawApts$surface > 30))
  t1 <- semireg_tmb(areaPerMzloty ~ construction.date, smoothZ = list(grp = smZ(construction.date, k = 10)),
                    data = subset(WarsawApts, surface > 30), prt = FALSE)
  expect_equal(as.numeric(logLik(t1$semer)), as.numeric(logLik(r$m$semer)), tolerance = 1e-3)
})

test_that("unconditional = TRUE: the smoothing-parameter correction, close to mgcv's; lme4 and glmmTMB agree", {
  skip_on_cran()
  q <- function(e) suppressWarnings(suppressMessages(e))
  set.seed(1)
  n <- 200
  d <- data.frame(x = runif(n, 0, 3), f = gl(2, n/2, labels = c("a", "b")))
  d$y <- sin(2*d$x) + 0.6*(d$f == "b")*d$x + rnorm(n, 0, 0.4)
  fit <- q(semireg(y ~ s(x, by = f), data = d, prt = FALSE))
  g <- mgcv::gam(y ~ f + s(x, by = f), data = d, method = "REML")
  xs <- c(0.2, 1.5, 2.8)
  nd <- expand.grid(x = xs, f = levels(d$f))
  a <- semipred(fit, "f", "x", covariateV = xs, prt = FALSE)$pred_df
  b <- semipred(fit, "f", "x", covariateV = xs, prt = FALSE, unconditional = TRUE)$pred_df
  a <- a[order(a$f, a$x), ]
  b <- b[order(b$f, b$x), ]
  pu <- as.numeric(predict(g, nd, se.fit = TRUE, unconditional = TRUE)$se.fit)
  expect_true(all(b$SE >= a$SE))
  # first-order correction (mgcv adds a second-order term too): within 1%
  expect_equal(b$SE, unname(pu), tolerance = 0.01)
  # the option turns it on for predictmeans() too
  op <- options(predictmeans.semi_unconditional = TRUE)
  on.exit(options(op))
  pm <- predictmeans(fit, "f", covariate = c(x = 1.5))
  expect_equal(pm$means$std.error, b$SE[c(2, 5)], tolerance = 1e-8)
  options(op)
  # Poisson: the lme4 correction (stage-1 deviance) and glmmTMB's joint precision agree
  d$cnt <- rpois(n, exp(0.5 + sin(2*d$x)))
  p <- q(semireg(cnt ~ s(x), family = poisson, data = d, prt = FALSE))
  t1 <- q(semireg_tmb(cnt ~ s(x), family = poisson, data = d, prt = FALSE))
  pb <- q(semipred(p, NULL, "x", covariateV = xs, prt = FALSE, unconditional = TRUE))$pred_df
  tb <- q(semipred(t1, NULL, "x", covariateV = xs, prt = FALSE, unconditional = TRUE))$pred_df
  pa <- semipred(p, NULL, "x", covariateV = xs, prt = FALSE)$pred_df
  expect_true(all(pb$SE_mu >= pa$SE_mu))
  expect_equal(pb$SE_mu, tb$SE_mu, tolerance = 0.01)
})

test_that("anova() of a semireg fit: Kenward-Roger with a smooth variance at zero, smooth LRTs, no lmerTest", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  data(Dialyzer, package = "nlme", envir = environment())
  fit <- q(semireg(rate ~ s(pressure, QB, k = 6) + (1 | Subject), data = Dialyzer, prt = FALSE))
  expect_false(methods::is(fit$semer, "lmerModLmerTest"))
  a <- q(anova(fit))
  expect_equal(a$smooth$variance[a$smooth$term == "s_pressure_QB"], 0, tolerance = 1e-8)
  # the deviation block is estimated at zero: KR leaves it out, so the df are those of the
  # model without it (before: 0.54 and 0.04 df for QB and pressure:QB)
  b <- q(semireg(rate ~ s(pressure, k = 6) + QB + pressure:QB + (1 | Subject), data = Dialyzer, prt = FALSE))
  ab <- anova(b$semer)
  expect_equal(a$fixed["QB", "Df.res"], ab["QB", "Df.res"], tolerance = 1e-4)
  expect_equal(a$fixed["pressure:QB", "Df.res"], ab["QB:pressure", "Df.res"], tolerance = 1e-4)
  expect_gt(a$fixed["QB", "Df.res"], 10)
  # no random term left: the REML log-likelihood of the linear model
  set.seed(7)
  d <- data.frame(x = runif(200, 0, 3), f = gl(2, 100))
  d$y <- sin(2*d$x) + rnorm(200, 0, 0.3)
  f1 <- q(semireg(y ~ f + s(x, k = 8), data = d, prt = FALSE))
  s1 <- anova(f1)$smooth
  expect_equal(s1$LRT, 2*(as.numeric(logLik(f1$semer)) - as.numeric(logLik(lm(y ~ f + x, d), REML = TRUE))),
               tolerance = 1e-6)
  expect_output(print(anova(f1)), "Smooth terms")
})

test_that("semireg_tmb(): a glmmTMB family named without glmmTMB attached (family = nbinom1), as in 1.x", {
  set.seed(4)
  d <- data.frame(x = runif(120, 0, 3), g = gl(6, 20))
  d$y <- rnbinom(120, mu = exp(1 + sin(d$x)), size = 3)
  # stats and the packages after it on the search path, not glmmTMB (load_all() puts imports earlier on it)
  e <- new.env(parent = as.environment("package:stats"))
  e$d <- d
  fit <- suppressWarnings(suppressMessages(evalq(
    predictmeans::semireg_tmb(y ~ x + (1 | g), smoothZ = list(sm = predictmeans::smZ(x, k = 5)), family = nbinom1,
                              data = d, prt = FALSE), e)))
  expect_equal(fit$semer$modelInfo$family$family, "nbinom1")
})
