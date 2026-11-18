## Smooths written in the formula of semireg() / semireg_tmb()
## (R/semireg-smooth.R): the same fits and predictions as mgcv::gam(REML),
## the fixed parts added automatically, the own bases equal to the mixed
## model with the O'Sullivan or thin-plate basis (as 1.x's smZ()).

## The same model fitted by nlme: y ~ fixed, coefficients of the basis Z iid.
lme_basis <- function(fixed, Z, d) {
  d$Z <- as.matrix(Z)
  d$one <- factor(1)
  nlme::lme(fixed, random = list(one = nlme::pdIdent(~ -1 + Z)), data = d, method = "REML")
}
os_knots <- function(x, k) stats::quantile(unique(x), seq(0, 1, length = k + 2)[-c(1, k + 2)])
os_range <- function(x) c(1.01*min(x) - 0.01*max(x), 1.01*max(x) - 0.01*min(x))

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
    s <- q(semireg(cs[[1]], data = d))
    g <- mgcv::gam(cs[[2]], data = d, method = "REML")
    lab <- deparse(cs[[1]])
    expect_equal(unname(fitted(s$semer)), unname(fitted(g)), tolerance = 1e-4, label = lab)
    fx <- attr(stats::terms(lme4::nobars(stats::formula(s$semer))), "term.labels")
    expect_true(all(cs[[3]] %in% fx), label = lab)
  }
  # writing the fixed part or not: the same model
  a <- q(semireg(y ~ s(x, by = f), data = d))
  b <- q(semireg(y ~ f*x + s(x, by = f), data = d))
  expect_equal(as.numeric(logLik(a$semer)), as.numeric(logLik(b$semer)), tolerance = 1e-8)
})

test_that("bs = 'os' and 'ztps': the mixed model with the O'Sullivan or thin-plate basis (as 1.x's smZ())", {
  d <- sf_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  pairs <- list(
    list(q(semireg(y ~ s(x, bs = "os", k = 12), data = d)),
         lme_basis(y ~ x, ZOSull(d$x, intKnots = os_knots(d$x, 12), range.x = os_range(d$x)), d)),
    list(q(semireg(y ~ s(x, z, bs = "ztps", k = 6), data = d)),
         lme_basis(y ~ x + z, Ztps(cbind(d$x, d$z), k = 6), d)))
  for (p in pairs) {
    expect_equal(as.numeric(logLik(p[[1]]$semer)), as.numeric(logLik(p[[2]])), tolerance = 1e-6)
    expect_equal(unname(fitted(p[[1]]$semer)), unname(fitted(p[[2]])), tolerance = 1e-5, ignore_attr = TRUE)
  }
})

test_that("semipred() and predictmeans() on formula smooths equal predict.gam(): curves, SEs, differences", {
  d <- sf_data()
  fit <- suppressMessages(semireg(y ~ s(x, by = f) + s(z), data = d))
  g <- mgcv::gam(y ~ f + s(x, by = f) + s(z), data = d, method = "REML")
  xs <- c(0.5, 1.5, 2.5)
  sp <- semipred(fit, "f", "x", at = xs)$pred_df
  sp <- sp[order(sp$f, sp$x), ]
  nd <- expand.grid(x = xs, f = levels(d$f))
  nd$z <- mean(d$z)
  pg <- predict(g, nd, se.fit = TRUE)
  expect_equal(sp$y, as.numeric(pg$fit), tolerance = 1e-5)
  expect_equal(sp$SE, as.numeric(pg$se.fit), tolerance = 1e-4)
  pm <- predictmeans(fit, "f", covariate = c(x = 1.5))
  expect_equal(pm$means$estimate, as.numeric(pg$fit[c(2, 5)]), tolerance = 1e-5)
  dc <- semipred(fit, "f", "x", at = xs, contr = c("b", "a"))$pred_df
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
  p <- q(semireg(cnt ~ s(x) + f + (1 | id), family = poisson, data = d))
  gp <- mgcv::gam(cnt ~ s(x) + f + s(id, bs = "re"), family = poisson, data = d, method = "REML")
  expect_lt(max(abs(predict(p$semer) - predict(gp))), 0.05)      # Laplace vs LAML: close
  # glmmTMB's sigma() is theta for nbinom2: the covariance used to be scaled by theta^2
  t1 <- q(semireg_tmb(cnt ~ s(x, bs = "os", k = 8) + f + (1 | id), family = glmmTMB::nbinom2,
                      data = d))
  # lme4 with the negative binomial of glmmTMB's theta: the same model
  th <- glmmTMB::sigma(t1$semer)
  l1 <- q(semireg(cnt ~ s(x, bs = "os", k = 8) + f + (1 | id),
                  family = MASS::negative.binomial(theta = th), data = d))
  a <- semipred(t1, "f", "x", at = c(0.5, 1.5, 2.5))$pred_df
  b <- semipred(l1, "f", "x", at = c(0.5, 1.5, 2.5))$pred_df
  expect_equal(a$SE_mu, b$SE_mu, tolerance = 0.05)
  t2 <- q(semireg_tmb(cnt ~ s(x) + f + (1 | id), family = glmmTMB::nbinom2, data = d))
  expect_equal(names(t2$u_lst), "s_x")
  expect_true(all(is.finite(semipred(t2, "f", "x")$pred_df$SE_mu)))
})

test_that("formula smooths: refusals with a message", {
  d <- sf_data()
  # tensor products: t2() is fitted, te() as t2() with a message, ti() refused (test-semireg-tensor.R)
  expect_s3_class(suppressWarnings(semireg(y ~ t2(x, z), data = d)), "semireg")
  expect_message(suppressWarnings(semireg(y ~ te(x, z), data = d)), "fitted as t2")
  expect_error(semireg(y ~ ti(x, z), data = d), "ti\\(\\) is not available")
  expect_error(semireg(y ~ s(x), smoothZ = list(), data = d), "'smoothZ' was removed in predictmeans 2.0")
  expect_error(semireg(y ~ s(x):f, data = d), "term of its own")
})

test_that("s(x, f): common curve + deviation curves with one variance, as the same model in lme", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  set.seed(7)
  n <- 480
  d <- data.frame(x = runif(n, 0, 3), f1 = gl(2, n/2, labels = c("a", "b")), f2 = gl(3, 1, n, labels = c("p", "q", "r")))
  d$y <- sin(2*d$x) + 0.5*(d$f1 == "b")*sin(3*d$x) + 0.3*(d$f2 == "r") + rnorm(n, 0, 0.3)
  fit <- q(semireg(y ~ s(x, f1, k = 8), data = d))
  expect_setequal(names(fit$u_lst), c("s_x", "s_x_f1"))
  fx <- attr(stats::terms(lme4::nobars(stats::formula(fit$semer))), "term.labels")
  expect_setequal(fx, c("x", "f1", "x:f1"))
  # Harezlak, Ruppert and Wand's form: one pdIdent block for the common curve, one for the deviations
  dd <- fit$data
  bl <- lme4::getME(fit$semer, "Ztlist")
  dd$Zc <- as.matrix(t(bl[["s_x.(Intercept)"]]))
  dd$Zd <- as.matrix(t(bl[["s_x_f1.(Intercept)"]]))
  dd$one <- factor(1)
  fl <- nlme::lme(y ~ x*f1, random = list(one = nlme::pdIdent(~ -1 + Zc), one = nlme::pdIdent(~ -1 + Zd)), data = dd,
                  method = "REML")
  expect_equal(as.numeric(logLik(fit$semer)), as.numeric(logLik(fl)), tolerance = 1e-6)
  # each level's curve is fixed part + common curve + its deviation
  sp <- semipred(fit, "f1", "x", at = 1.5)$pred_df
  b <- lme4::fixef(fit$semer)
  Zn <- predictmeans:::pm_sm_predict(fit$smooths, data.frame(x = 1.5, f1 = factor(c("a", "b"), levels = c("a", "b"))))
  by_hand <- b[1] + b["x"]*1.5 + c(0, b["f1b"] + b["x:f1b"]*1.5) +
    as.numeric(Zn$s_x %*% fit$u_lst$s_x) + as.numeric(Zn$s_x_f1 %*% fit$u_lst$s_x_f1)
  expect_equal(sp$y[order(sp$f1)], unname(by_hand), tolerance = 1e-8)
  # s(x, f1*f2): a deviation term per factor term, each with its variance; no f1:f2 deviation in these data
  f2 <- q(semireg(y ~ s(x, f1*f2, k = 8), data = d))
  expect_setequal(names(f2$u_lst), c("s_x", "s_x_f1", "s_x_f2", "s_x_f1_f2"))
  fx2 <- attr(stats::terms(lme4::nobars(stats::formula(f2$semer))), "term.labels")
  expect_setequal(fx2, c("x", "f1", "f2", "x:f1", "x:f2", "f1:f2", "x:f1:f2"))
  v <- varcomp(f2)
  expect_lt(v["s_x_f1_f2.(Intercept)", "vcov"], 0.1*v["s_x_f1.(Intercept)", "vcov"])
  expect_equal(nrow(predictmeans(f2, "f1:f2", covariate = c(x = 1.5))$means), 6)
  # s(x, f1:f2) is s(x, cell) with the cells as one factor
  f3 <- q(semireg(y ~ s(x, f1:f2, k = 8), data = d))
  d$cell <- interaction(d$f1, d$f2, sep = ":")
  f3b <- q(semireg(y ~ s(x, cell, k = 8), data = d))
  expect_equal(as.numeric(logLik(f3$semer)), as.numeric(logLik(f3b$semer)), tolerance = 1e-8)
  expect_equal(nrow(semipred(f3, "f1:f2", "x", at = c(1, 2))$pred_df), 12)
  # surfaces: s(x, z, f)
  d$z <- runif(n)
  f4 <- q(semireg(y ~ s(x, z, f1, k = 12), data = d))
  expect_setequal(names(f4$u_lst), c("s_x_z", "s_x_z_f1"))
  expect_setequal(attr(stats::terms(lme4::nobars(stats::formula(f4$semer))), "term.labels"),
                  c("x", "z", "f1", "x:f1", "z:f1"))
  # refusals
  expect_error(semireg(y ~ s(x, f1, f2), data = d), "one term")
  expect_error(semireg(y ~ s(x, f1, by = f2), data = d), "cannot also have 'by'")
  expect_error(semireg(y ~ s(f1:f2), data = d), "needs a covariate")
})

test_that("random curves (s(x) | g): (x | g) plus a spline deviation per level with one variance", {
  d <- sf_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  ctl <- lme4::lmerControl(optimizer = "bobyqa")
  a <- q(semireg(y ~ s(x, bs = "os", k = 8) + (s(x, k = 5) | id), data = d, control = ctl))
  expect_true("s_x_id" %in% names(a$u_lst))
  # population curve: the random curves are left out (the default is the population smooth)
  pa <- q(semipred(a, "x", "x"))$pred_df
  pb <- q(semipred(a, "x", "x", sm_term = "s_x"))$pred_df
  expect_equal(pa[[3]], pb[[3]], tolerance = 1e-10)
  expect_equal(pa$SE, pb$SE, tolerance = 1e-10)
  expect_error(semipred(a, "x", "x", sm_term = "s_x_id"), "smooth blocks")
  expect_error(semireg(y ~ s(x) + (s(x, by = f) | id), data = d), "cannot have 'by'")
})

test_that("xt = list(knots, range) sets the knots and range of the own bases", {
  d <- sf_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  kn <- c(0.5, 1, 1.5, 2, 2.5)
  a <- q(semireg(y ~ s(x, bs = "os", xt = list(knots = kn, range = c(0, 3))), data = d))
  b <- lme_basis(y ~ x, ZOSull(d$x, intKnots = kn, range.x = c(0, 3)), d)
  expect_equal(as.numeric(logLik(a$semer)), as.numeric(logLik(b)), tolerance = 1e-6)
  km <- as.matrix(expand.grid(x = c(0.5, 1.5, 2.5), z = c(0.25, 0.75)))
  a <- q(semireg(y ~ s(x, z, bs = "ztps", xt = list(knots = km)), data = d))
  b <- lme_basis(y ~ x + z, Ztps(cbind(d$x, d$z), knots = km), d)
  expect_equal(as.numeric(logLik(a$semer)), as.numeric(logLik(b)), tolerance = 1e-6)
})

test_that("smooths with fewer distinct covariate values than basis columns: exact map, as gam()", {
  d <- sf_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  d$xr <- round(d$x*2)/2                       # 7 distinct values, os basis of 14 columns
  s1 <- q(semireg(y ~ s(xr, bs = "os", k = 10), data = d))
  g1 <- mgcv::gam(y ~ s(xr, bs = "os", k = 10), data = d, method = "REML")
  expect_equal(unname(fitted(s1$semer)), unname(fitted(g1)), tolerance = 1e-4)
  nd <- data.frame(xr = seq(0, 3, by = 0.25))
  p1 <- q(semipred(s1, "xr", "xr"))$pred_df
  i <- match(nd$xr, round(p1$xr, 8))
  expect_gt(sum(!is.na(i)), 5)
  expect_equal(p1$y[i][!is.na(i)], as.numeric(predict(g1, nd))[!is.na(i)], tolerance = 1e-3)
  # a factor smooth with levels of three observations (fewer than the basis columns)
  d2 <- d[!(d$id %in% c("3", "4")) | ave(seq_len(nrow(d)), d$id, FUN = seq_along) <= 3, ]
  s2 <- q(semireg(y ~ s(x, id, bs = "fs", k = 5), data = droplevels(d2)))
  g2 <- mgcv::gam(y ~ s(x, id, bs = "fs", k = 5), data = droplevels(d2), method = "REML")
  expect_equal(unname(fitted(s2$semer)), unname(fitted(g2)), tolerance = 1e-4)
})

test_that("objects of the session in xt or k of random curves and t2() are not taken for data", {
  set.seed(1)
  d <- data.frame(x = rep(seq(0, 1, length = 10), 12), z = runif(120), id = factor(rep(1:12, each = 10)),
                  g = factor(rep(1:2, each = 60)))
  d$y <- sin(3 * d$x) + d$z + rnorm(12, 0, 0.3)[d$id] + rnorm(120, 0, 0.2)
  rg <- c(-0.05, 1.05)
  kk <- 4
  fm <- suppressMessages(semireg(y ~ s(x, by = g, bs = "os", k = 6, xt = list(range = rg)) +
                                   (s(x, bs = "os", k = kk, xt = list(range = rg)) | id), data = d))
  expect_false(any(c("rg", "kk") %in% fm$sm_vars))
  expect_true(all(is.finite(semipred(fm, "g", "x")$pred_df$SE)))
  expect_true(all(is.finite(semipred(fm, "g", "x", contr = c("2", "1"))$pred_df$SE_contr)))
  expect_true(all(is.finite(predictmeans(fm, "g", covariate = c(x = 0.5))$means$std.error)))
  ft <- suppressMessages(semireg(y ~ t2(x, z, k = kk) + (1 | id), data = d))
  expect_false("kk" %in% ft$sm_vars)
  expect_true(all(is.finite(semipred(ft, covariate = c("x", "z"))$pred_df$SE)))
})

test_that("print() and summary() of a fit show the model, not the list", {
  fit <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (1 | Subject), data = nlme::Dialyzer)
  out <- capture.output(print(fit))
  expect_true(any(grepl("REML criterion", out)))
  expect_true(any(grepl("s(pressure):QB200", out, fixed = TRUE)))
  expect_false(any(grepl("^\\$data", out)))
  expect_lt(length(out), 60)
  expect_s3_class(summary(fit), "summary.merMod")
  ft <- semireg_tmb(rate ~ s(pressure, bs = "os") + (1 | Subject), data = nlme::Dialyzer)
  out <- capture.output(print(ft))
  expect_true(any(grepl("s_pressure", out)))
  expect_lt(length(out), 60)
  expect_s3_class(summary(ft), "summary.glmmTMB")
})
