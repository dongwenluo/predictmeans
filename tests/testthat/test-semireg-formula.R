## Smooths written in the formula of semireg() / semireg_tmb()
## (R/semireg-smooth.R): the same fits and predictions as mgcv::gam(REML),
## the fixed parts added automatically, the own bases equal to smZ().

sf_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = runif(n, 0, 3), z = runif(n), f = gl(2, n/2, labels = c("a", "b")), w = runif(n),
                  id = gl(20, 1, n))
  d$y <- sin(2*d$x) + 0.6*(d$f == "b")*d$x + cos(3*d$z) + rnorm(20, 0, 0.3)[d$id] + rnorm(n, 0, 0.3)
  d
}

test_that("formula smooths: the same fit as gam(method = 'REML'), fixed parts added automatically", {
  d <- sf_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  cases <- list(
    list(y ~ s(x), y ~ s(x), "x"),
    list(y ~ s(x, bs = "cr", k = 12), y ~ s(x, bs = "cr", k = 12), "x"),
    list(y ~ s(x, by = f), y ~ f + s(x, by = f), c("f", "f:x")),          # f added: gam needs it written
    list(y ~ s(x) + s(z), y ~ s(x) + s(z), c("x", "z")),
    list(y ~ s(x, z), y ~ s(x, z), c("x", "z")),
    list(y ~ s(x, by = w), y ~ s(x, by = w), c("w", "w:x")),             # numeric by, not centred
    list(y ~ s(x, by = f, id = 1), y ~ f + s(x, by = f, id = 1), c("f", "f:x")),
    list(y ~ s(x, bs = "ps", k = 15), y ~ s(x, bs = "ps", k = 15), "x"),
    list(y ~ f + s(x, f, bs = "fs", k = 6), y ~ f + s(x, f, bs = "fs", k = 6), "f"),
    list(y ~ s(x, bs = "cs"), y ~ s(x, bs = "cs"), character(0)),
    list(y ~ s(x, by = f) + (1 | id), y ~ f + s(x, by = f) + s(id, bs = "re"), c("f", "f:x")))
  for (cs in cases) {
    s <- q(semireg(cs[[1]], data = d, prt = FALSE))
    g <- mgcv::gam(cs[[2]], data = d, method = "REML")
    lab <- deparse(cs[[1]])
    expect_equal(unname(fitted(s$semer)), unname(fitted(g)), tolerance = 1e-4, label = lab)
    fx <- attr(stats::terms(lme4::nobars(stats::formula(s$semer))), "term.labels")
    expect_true(all(cs[[3]] %in% fx), label = lab)
  }
  # writing the fixed part or not: the same model
  a <- q(semireg(y ~ s(x, by = f), data = d, prt = FALSE))
  b <- q(semireg(y ~ f*x + s(x, by = f), data = d, prt = FALSE))
  expect_equal(as.numeric(logLik(a$semer)), as.numeric(logLik(b$semer)), tolerance = 1e-8)
})

test_that("bs = 'os' and 'ztps' are the bases of smZ(): the same fits as the smoothZ form", {
  d <- sf_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  pairs <- list(
    list(q(semireg(y ~ x, smoothZ = list(sm = smZ(x, k = 12)), data = d, prt = FALSE)),
         q(semireg(y ~ s(x, bs = "os", k = 12), data = d, prt = FALSE))),
    list(q(semireg(y ~ x*f, smoothZ = list(sm = smZ(x, k = 8, by = f, group = TRUE)), data = d, prt = FALSE)),
         q(semireg(y ~ s(x, by = f, bs = "os", k = 8), data = d, prt = FALSE))),
    list(q(semireg(y ~ x + z, smoothZ = list(sm = smZ(cbind(x, z), k = 6, type = "Ztps")), data = d, prt = FALSE)),
         q(semireg(y ~ s(x, z, bs = "ztps", k = 6), data = d, prt = FALSE))))
  for (p in pairs) {
    expect_equal(as.numeric(logLik(p[[1]]$semer)), as.numeric(logLik(p[[2]]$semer)), tolerance = 1e-8)
    expect_equal(unname(fitted(p[[1]]$semer)), unname(fitted(p[[2]]$semer)), tolerance = 1e-6)
  }
})

test_that("semipred() and predictmeans() on formula smooths equal predict.gam(): curves, SEs, differences", {
  d <- sf_data()
  fit <- suppressMessages(semireg(y ~ s(x, by = f) + s(z), data = d, prt = FALSE))
  g <- mgcv::gam(y ~ f + s(x, by = f) + s(z), data = d, method = "REML")
  xs <- c(0.5, 1.5, 2.5)
  sp <- semipred(fit, "f", "x", covariateV = xs, prt = FALSE)$pred_df
  sp <- sp[order(sp$f, sp$x), ]
  nd <- expand.grid(x = xs, f = levels(d$f))
  nd$z <- mean(d$z)
  pg <- predict(g, nd, se.fit = TRUE)
  expect_equal(sp$y, as.numeric(pg$fit), tolerance = 1e-5)
  expect_equal(sp$SE, as.numeric(pg$se.fit), tolerance = 1e-4)
  pm <- predictmeans(fit, "f", covariate = c(x = 1.5))
  expect_equal(pm$means$estimate, as.numeric(pg$fit[c(2, 5)]), tolerance = 1e-5)
  dc <- semipred(fit, "f", "x", covariateV = xs, contr = c("b", "a"), prt = FALSE)$pred_df
  Xp <- predict(g, nd, type = "lpmatrix")
  L <- Xp[4:6, ] - Xp[1:3, ]
  expect_equal(dc$y, as.numeric(L %*% coef(g)), tolerance = 1e-5)
  expect_equal(dc$SE_contr, unname(sqrt(diag(L %*% g$Vp %*% t(L)))), tolerance = 1e-4)
  # variance components and R2 are named by the smooth blocks
  expect_true(all(c("s_x_fa.(Intercept)", "s_x_fb.(Intercept)", "s_z.(Intercept)") %in% rownames(varcomp(fit))))
  expect_true(all(c("s_x_fa", "s_x_fb", "s_z") %in% names(R2_model(fit))))
})

test_that("formula smooths in GLMMs and glmmTMB; negative binomial SEs of semireg_tmb() no longer scaled by theta", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  set.seed(4)
  n <- 600
  d <- data.frame(x = runif(n, 0, 3), f = gl(2, n/2, labels = c("a", "b")), id = gl(30, 1, n))
  u <- rnorm(30, 0, 0.4)
  d$cnt <- MASS::rnegbin(n, exp(0.3 + sin(2*d$x) + 0.4*(d$f == "b") + u[d$id]), 3)
  p <- q(semireg(cnt ~ s(x) + f + (1 | id), family = poisson, data = d, prt = FALSE))
  gp <- mgcv::gam(cnt ~ s(x) + f + s(id, bs = "re"), family = poisson, data = d, method = "REML")
  expect_lt(max(abs(predict(p$semer) - predict(gp))), 0.05)      # Laplace vs LAML: close
  # glmmTMB's sigma() is theta for nbinom2: the covariance used to be scaled by theta^2
  t1 <- q(semireg_tmb(cnt ~ x + f + (1 | id), smoothZ = list(sx = smZ(x, k = 8)), family = glmmTMB::nbinom2,
                      data = d, prt = FALSE))
  l1 <- q(semireg(cnt ~ x + f + (1 | id), smoothZ = list(sx = smZ(x, k = 8)), ngbinomial = TRUE, data = d,
                  prt = FALSE))
  a <- semipred(t1, "f", "x", covariateV = c(-1, 0, 1), prt = FALSE)$pred_df
  b <- semipred(l1, "f", "x", covariateV = c(-1, 0, 1), prt = FALSE)$pred_df
  expect_equal(a$SE_mu, b$SE_mu, tolerance = 0.1)
  t2 <- q(semireg_tmb(cnt ~ s(x) + f + (1 | id), family = glmmTMB::nbinom2, data = d, prt = FALSE))
  expect_equal(names(t2$u_lst), "s_x")
  expect_true(all(is.finite(semipred(t2, "f", "x", prt = FALSE)$pred_df$SE_mu)))
})

test_that("formula smooths: refusals with a message", {
  d <- sf_data()
  expect_error(semireg(y ~ t2(x, z), data = d, prt = FALSE), "Tensor-product")
  expect_error(semireg(y ~ te(x, z), data = d, prt = FALSE), "Tensor-product")
  expect_error(semireg(y ~ s(x), smoothZ = list(sm = smZ(z)), data = d, prt = FALSE), "not both")
  expect_error(semireg(y ~ s(x):f, data = d, prt = FALSE), "term of its own")
})
