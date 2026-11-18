## Tensor products in semireg(): t2() (te() is fitted as t2(), ti() refused)

q <- function(e) suppressWarnings(suppressMessages(e))
tensor_data <- function() {
  set.seed(1)
  d <- data.frame(x = runif(300), z = runif(300), g = gl(10, 30))
  d$y <- sin(3*d$x)*cos(2*d$z) + rnorm(10, 0, 0.3)[d$g] + rnorm(300, 0, 0.2)
  d
}

test_that("t2(): the same model as gam() with the unconstrained t2 basis and its penalties (paraPen)", {
  d <- tensor_data()
  fs <- q(semireg(y ~ t2(x, z) + (1 | g), data = d))
  expect_equal(sort(names(fs$u_lst)), c("t2_x_z_1", "t2_x_z_2", "t2_x_z_3"))
  expect_true(all(c("x", "z", "x:z") %in% attr(terms(fs$semer, fixed.only = TRUE), "term.labels")))
  sm <- mgcv::smoothCon(mgcv::t2(x, z), data = d, absorb.cons = FALSE)[[1]]
  X <- mgcv::PredictMat(sm, d)
  gm <- mgcv::gam(y ~ X - 1 + s(g, bs = "re"), data = d, paraPen = list(X = sm$S), method = "REML")
  expect_lt(max(abs(fitted(fs$semer) - fitted(gm)))/sd(d$y), 1e-5)
  nd <- data.frame(x = c(0.2, 0.5, 0.8), z = c(0.3, 0.5, 0.7))
  pg <- as.numeric(mgcv::PredictMat(sm, nd) %*% coef(gm)[seq_len(ncol(X))])
  ps <- semipred(fs, covariate = c("x", "z"), at = cbind(x = nd$x, z = nd$z))$pred_df$y
  expect_equal(ps, pg, tolerance = 1e-4)
  # without other random effects it is also gam()'s own t2() (mgcv's constraint on the whole
  # product then makes no difference to the fit)
  f0 <- q(semireg(y ~ t2(x, z), data = d))
  g0 <- mgcv::gam(y ~ t2(x, z), data = d, method = "REML")
  expect_lt(max(abs(fitted(f0$semer) - fitted(g0)))/sd(d$y), 1e-5)
  expect_equal(semipred(f0, covariate = c("x", "z"), at = cbind(x = nd$x, z = nd$z))$pred_df$y,
               as.numeric(predict(g0, newdata = nd)), tolerance = 1e-5)
  # the glmmTMB engine
  ft <- q(semireg_tmb(y ~ t2(x, z) + (1 | g), data = d))
  expect_equal(-ft$semer$fit$objective, as.numeric(logLik(fs$semer)), tolerance = 1e-6)
})

test_that("t2(..., by = f): a surface per level, fixed terms f + f:x + f:z + f:x:z; te() and ti()", {
  set.seed(1)
  d <- data.frame(x = runif(400), z = runif(400), A = gl(2, 200, labels = c("a1", "a2")), g = gl(20, 20))
  d$y <- sin(3*d$x)*cos(2*d$z) + 0.5*(d$A == "a2")*d$x*d$z + rnorm(20, 0, 0.3)[d$g] + rnorm(400, 0, 0.2)
  fb <- q(semireg(y ~ t2(x, z, by = A) + (1 | g), data = d))
  expect_true(all(c("A", "A:x", "A:z", "A:x:z") %in% attr(terms(fb$semer, fixed.only = TRUE), "term.labels")))
  expect_length(fb$u_lst, 6)
  pm <- q(predictmeans(fb, "A", covariate = c(x = 0.5, z = 0.5), pairwise = TRUE))
  expect_true(is.finite(pm$pairs$df) && pm$pairs$df > 10)
  expect_equal(nrow(anova(fb)$smooth), 6)
  # te(): fitted as t2(), with a message
  expect_message(semireg(y ~ te(x, z) + (1 | g), data = d), "fitted as t2\\(\\)")
  ft <- q(semireg(y ~ te(x, z) + (1 | g), data = d))
  f2 <- q(semireg(y ~ t2(x, z) + (1 | g), data = d))
  expect_equal(as.numeric(logLik(ft$semer)), as.numeric(logLik(f2$semer)))
  expect_error(semireg(y ~ ti(x, z), data = d), "ti\\(\\) is not available")
  expect_error(semireg(y ~ t2(x, A), data = d), "give it as 'by'")
})
