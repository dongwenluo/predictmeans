## bs = "user": a basis function of the user's own in the formula of semireg()

q <- function(e) suppressWarnings(suppressMessages(e))
user_data <- function() {
  set.seed(11)
  sp <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  sp$y <- 2 + sin(sp$t/2) + 0.08*(sp$A == "a2")*sp$t + rnorm(6, 0, 0.5)[sp$Block] + rnorm(nrow(sp), 0, 0.3)
  sp
}
kn <- c(2, 4, 6, 8)
rg <- c(-0.1, 10.1)
ob <- function(x) ZOSull(x, intKnots = kn, range.x = rg)

test_that("bs = \"user\" with the O'Sullivan basis is the bs = \"os\" fit (curves, slopes, means, fixed terms)", {
  sp <- user_data()
  fo <- q(semireg(y ~ s(t, by = A, bs = "os", xt = list(knots = kn, range = rg)) + (1 | Block), data = sp))
  fu <- q(semireg(y ~ s(t, by = A, bs = "user", xt = list(basis = ob)) + (1 | Block), data = sp))
  expect_equal(as.numeric(logLik(fu$semer)), as.numeric(logLik(fo$semer)), tolerance = 1e-8)
  expect_true(all(c("A", "A:t") %in% attr(terms(fu$semer, fixed.only = TRUE), "term.labels")))
  for (dv in 0:1) {
    po <- semipred(fo, "A", "t", deriv = dv)$pred_df
    pu <- semipred(fu, "A", "t", deriv = dv)$pred_df
    expect_equal(pu$y, po$y, tolerance = 1e-8)
    expect_equal(pu$SE, po$SE, tolerance = 1e-8)
  }
  expect_equal(predictmeans(fu, "A", covariate = c(t = 5))$means, predictmeans(fo, "A", covariate = c(t = 5))$means,
               tolerance = 1e-8)
  # xt as the function itself; the crossed-factor syntax
  fx <- q(semireg(y ~ s(t, A, bs = "user", xt = ob) + (1 | Block), data = sp))
  expect_setequal(names(fx$u_lst), c("s_t", "s_t_A"))
  # the glmmTMB engine
  a <- q(semireg_tmb(y ~ s(t, bs = "user", xt = ob) + (1 | Block), data = sp))
  b <- q(semireg_tmb(y ~ s(t, bs = "os", xt = list(knots = kn, range = rg)) + (1 | Block), data = sp))
  expect_equal(a$semer$fit$objective, b$semer$fit$objective, tolerance = 1e-8)
})

test_that("bs = \"user\" with two covariates: the Ztps basis gives the bs = \"ztps\" fit", {
  set.seed(3)
  s2 <- data.frame(u = runif(200), v = runif(200), g = gl(2, 100, labels = c("a", "b")))
  s2$w <- sin(3*s2$u) + cos(3*s2$v) + 0.4*(s2$g == "b") + rnorm(200, sd = 0.2)
  kk <- as.matrix(expand.grid(u = seq(0.1, 0.9, length = 4), v = seq(0.1, 0.9, length = 4)))
  zb <- function(u, v) unclass(as.matrix(Ztps(cbind(u, v), knots = kk)))
  fz <- q(semireg(w ~ g + s(u, v, bs = "ztps", xt = list(knots = kk)), data = s2))
  fu <- q(semireg(w ~ g + s(u, v, bs = "user", xt = list(basis = zb)), data = s2))
  expect_equal(as.numeric(logLik(fu$semer)), as.numeric(logLik(fz$semer)), tolerance = 1e-8)
  pts <- cbind(u = c(0.2, 0.5, 0.8), v = c(0.3, 0.5, 0.7))
  expect_equal(semipred(fu, "g", c("u", "v"), at = pts)$pred_df$w,
               semipred(fz, "g", c("u", "v"), at = pts)$pred_df$w, tolerance = 1e-8)
})

test_that("bs = \"user\": the basis function is checked before the fit", {
  sp <- user_data()
  fit <- function(xt) semireg(y ~ s(t, bs = "user", xt = xt) + (1 | Block), data = sp)
  one_row <- function(f) function(x) f(matrix(ob(x), length(x)), x)
  expect_error(fit(list(basis = 3)), "needs the basis function")
  expect_error(fit(list(basis = ob, null = "quadratic")), "must be \"linear\", \"constant\" or \"none\"")
  expect_error(fit(function(x) stop("no basis here")), "failed on the data: no basis here")
  expect_error(fit(function(x) ob(x)[-1, ]), "a row per value \\(198 rows\\), not 197")
  expect_error(fit(function(x) { B <- ob(x); B[x == 3, 1] <- NA; B }), "NA or infinite values \\(first at t = 3\\)")
  expect_error(fit(function(x) ob(x) + stats::rnorm(length(x), 0, 1e-3)), "different values on two calls")
  # knots from the data: prediction would use other ones
  expect_error(fit(function(x) ZOSull(x, intKnots = stats::quantile(x, 1:4/5), range.x = range(x) + c(-0.1, 0.1))),
               "depends on the data it is given")
  # a single value must give one row (prediction at one covariate value)
  # (a column instead of a row for one value; before 2.0.0 ZOSull() itself gave a wrong shape here)
  expect_error(fit(function(x) { B <- ob(x); if (length(x) == 1) t(B) else B }), "for 1 value it returned")
  expect_error(fit(function(x) { if (length(x) < 20) stop("needs 20 values"); ob(x) }),
               "failed on a single value \\(as at prediction\\): needs 20 values")
  expect_error(fit(one_row(function(B, x) cbind(B, B[, 1] + B[, 2]))), "linearly dependent columns \\(rank 6 of 7\\)")
  expect_error(fit(one_row(function(B, x) cbind(1, B))), "span of the fixed part \\(the constant\\).*null = \"none\"")
  expect_error(fit(one_row(function(B, x) cbind(x, B))), "the constant and the straight line.*null = \"constant\"")
  expect_s3_class(q(fit(list(basis = one_row(function(B, x) cbind(x, B)), null = "constant"))), "semireg")
  # no fixed part: a periodic basis
  fr <- function(x) cbind(sin(2*pi*x/10), cos(2*pi*x/10), sin(4*pi*x/10), cos(4*pi*x/10))
  f0 <- q(fit(list(basis = fr, null = "none")))
  expect_false("t" %in% attr(terms(f0$semer, fixed.only = TRUE), "term.labels"))
  # failing just beyond the ends of the data: a warning (slopes there need it)
  inside <- function(x) { if (any(x < 0 | x > 10)) stop("outside [0, 10]"); ob(x) }
  expect_warning(fit(inside), "fails just beyond the ends of the data")
})
