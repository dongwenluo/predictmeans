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
  # anova(): lmerTest's Satterthwaite df with the deviance of the fitted model
  # (before 2.0 that of the placeholder model: 1871 df instead of 150)
  ddf_lme <- pm_inference(fl, df_method = "Satterthwaite", vcov_adjust = FALSE)$df(matrix(c(0, 1), 1), "construction.date")
  expect_equal(anova(fit$semer)["construction.date", "DenDF"], ddf_lme, tolerance = 1e-3)
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
