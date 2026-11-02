# Supporting a new model class through the adapter generics (?pm_coef).
# The made-up classes below wrap lm()/glm() fits and provide only a few
# methods; predictmeans() must then give the same results as for the fits.

oats_fac <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}
reg <- function(generic, cls, fn, ns = "stats") registerS3method(generic, cls, fn, envir = asNamespace(ns))

test_that("a class with only coef, vcov, terms and df.residual works", {
  fit_toylm <- function(formula, data) {
    structure(list(fit = lm(formula, data), call = match.call()), class = "toylm")
  }
  reg("coef", "toylm", function(object, ...) coef(object$fit))
  reg("vcov", "toylm", function(object, ...) vcov(object$fit))
  reg("terms", "toylm", function(x, ...) terms(x$fit))
  reg("df.residual", "toylm", function(object, ...) object$fit$df.residual)

  Oats <- oats_fac()
  toy <- fit_toylm(yield ~ Block + nitro*Variety, data = Oats)
  ref <- lm(yield ~ Block + nitro*Variety, data = Oats)
  expect_error(model.frame(toy))
  expect_equal(Kmatrix(toy, "nitro:Variety")$K, Kmatrix(ref, "nitro:Variety")$K)
  a <- predictmeans(toy, "nitro:Variety", adj = "tukey", plot = FALSE, prtnum = FALSE)
  b <- predictmeans(ref, "nitro:Variety", adj = "tukey", plot = FALSE, prtnum = FALSE)
  expect_equal(a$means[c("estimate", "std.error", "df", "conf.low", "conf.high", "group", "n")],
               b$means[c("estimate", "std.error", "df", "conf.low", "conf.high", "group", "n")])
  expect_equal(a$pairs, b$pairs)
  expect_equal(a$mean_table, b$mean_table)
  cm <- specify_contrast(Kmatrix(toy, "nitro"), lin = `0.6` - `0`)
  expect_equal(contrastmeans(toy, "nitro", ctrmatrix = cm), contrastmeans(ref, "nitro", ctrmatrix = cm))
})

test_that("pm_data() and pm_df() methods are used", {
  dat <- oats_fac()
  fit_toyz <- function(formula, data) {
    structure(list(fit = lm(formula, data), call = quote(toyz()), stored = data), class = "toyz")
  }
  reg("coef", "toyz", function(object, ...) coef(object$fit))
  reg("vcov", "toyz", function(object, ...) vcov(object$fit))
  reg("terms", "toyz", function(x, ...) terms(x$fit))
  reg("pm_data", "toyz", function(model, ...) model$stored, ns = "predictmeans")
  reg("pm_df", "toyz", function(model, L, ...) rep(Inf, nrow(L)), ns = "predictmeans")
  toy <- fit_toyz(yield ~ nitro + Variety, dat)
  pm <- predictmeans(toy, "nitro", plot = FALSE, prtnum = FALSE)
  expect_true(all(is.infinite(pm$means$df)))
  expect_equal(pm$means$conf.high - pm$means$estimate, qnorm(0.975) * pm$means$std.error)
  expect_true(all(is.infinite(pm$pairs$df)))
  expect_identical(pm$settings$df_method, "model (pm_df)")
})

test_that("a pm_linkinv() method back-transforms the means", {
  set.seed(1)
  d <- data.frame(A = factor(rep(c("a1", "a2", "a3"), each = 20)), x = runif(60))
  d$y <- rpois(60, exp(1 + as.numeric(d$A)/3 + d$x))
  fit_toypois <- function(formula, data) {
    structure(list(fit = glm(formula, poisson, data), call = match.call()), class = "toypois")
  }
  reg("coef", "toypois", function(object, ...) coef(object$fit))
  reg("vcov", "toypois", function(object, ...) vcov(object$fit))
  reg("terms", "toypois", function(x, ...) terms(x$fit))
  reg("df.residual", "toypois", function(object, ...) object$fit$df.residual)
  reg("pm_linkinv", "toypois", function(model, ...) exp, ns = "predictmeans")
  toy <- fit_toypois(y ~ A + x, d)
  pm <- predictmeans(toy, "A", plot = FALSE, prtnum = FALSE)
  ref <- predictmeans(glm(y ~ A + x, poisson, d), "A", plot = FALSE, prtnum = FALSE)
  expect_equal(pm$means$response, exp(pm$means$estimate))
  expect_equal(pm$means$response, ref$means$response)
})

test_that("an object without estimates gives a clear message", {
  bad <- structure(list(call = quote(f(y ~ x)), terms = terms(y ~ x)), class = "nocoef")
  reg("terms", "nocoef", function(x, ...) x$terms)
  expect_error(Kmatrix(bad, "x"), "pm_coef")
})
