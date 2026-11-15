## Back-transformation (R/back-transform.R): the response transformation in
## the formula, Box-Cox and Yeo-Johnson, the mean of the back-transformed
## response, and the Box-Cox check of pm_check().

bt_oats <- function() transform(as.data.frame(nlme::Oats), nitro = factor(nitro))

test_that("the transformation of the response in the formula is found, with its inverse", {
  y <- c(0.5, 2, 7, 30)
  chk <- function(expr, fwd) {
    inv <- pm_parse_trans(expr)
    expect_false(is.null(inv), label = deparse(expr))
    expect_equal(inv(fwd(y)), y, tolerance = 1e-10, label = deparse(expr))
  }
  chk(quote(log(y)), log)
  expect_identical(pm_parse_trans(quote(log(y))), exp)               # exp itself: ratios in the pairs
  chk(quote(log(y + 1)), function(y) log(y + 1))
  chk(quote(log(1 + y)), log1p)
  chk(quote(log(y - 0.1)), function(y) log(y - 0.1))
  chk(quote(log(y, 2)), log2)
  chk(quote(log10(y)), log10)
  chk(quote(log2(y)), log2)
  chk(quote(log1p(y)), log1p)
  chk(quote(sqrt(y)), sqrt)
  chk(quote(sqrt(y + 0.5)), function(y) sqrt(y + 0.5))
  chk(quote(y^0.25), function(y) y^0.25)
  chk(quote(I(y^2)), function(y) y^2)
  chk(quote(1/y), function(y) 1/y)
  chk(quote(100/(y + 1)), function(y) 100/(y + 1))
  chk(quote(pm_boxcox(y, 0.3)), function(y) pm_boxcox(y, 0.3))
  chk(quote(pm_boxcox(y, 0)), log)
  chk(quote(predictmeans::pm_yeojohnson(y, -0.5)), function(y) pm_yeojohnson(y, -0.5))
  p <- c(0.1, 0.4, 0.9)
  expect_equal(pm_parse_trans(quote(qlogis(p)))(qlogis(p)), p)
  expect_equal(pm_parse_trans(quote(asin(sqrt(p))))(asin(sqrt(p))), p)
  # not transformations of one response
  for (e in list(quote(y), quote(cbind(a, b)), quote(log(y/x)), quote(y^x), quote(scale(y)), quote(log(y) + 1)))
    expect_null(pm_parse_trans(e), label = deparse(e))
  # Yeo-Johnson takes zero and negative values
  v <- c(-3, -0.5, 0, 0.5, 4)
  for (l in c(-1, 0, 0.5, 1, 2, 2.5)) expect_equal(pm_yeojohnson_inv(l)(pm_yeojohnson(v, l)), v, tolerance = 1e-10)
  expect_equal(pm_yeojohnson(v, 1), v)
  expect_error(pm_boxcox(c(1, 0), 0.5), "positive")
})

test_that("predictmeans() back-transforms a transformed response by default; trans = I keeps the model scale", {
  O <- bt_oats()
  m <- lm(log(yield) ~ Block + nitro, data = O)
  expect_message(r <- predictmeans(m, "nitro", pairwise = TRUE, prtnum = FALSE), "The response is log\\(yield\\)")
  expect_equal(r$means$response, exp(r$means$estimate))
  expect_true("Bk_Mean" %in% names(r$mean_table$Table))
  expect_true("ratio" %in% names(r$pairs))                              # differences of logs as ratios
  expect_identical(r$settings$back_transform$estimate, "median of the original scale (bias_adjust = TRUE gives the mean)")
  expect_match(r$settings$back_transform$by, "inverse of log\\(yield\\)")
  # the same as giving trans
  r2 <- predictmeans(m, "nitro", pairwise = TRUE, trans = exp, prtnum = FALSE)
  expect_equal(r$means, r2$means)
  expect_silent(r0 <- predictmeans(m, "nitro", trans = I, prtnum = FALSE))
  expect_null(r0$means$response)
  expect_null(r0$settings$back_transform)
  # an offset in the formula; other transformations
  r1 <- suppressMessages(predictmeans(lm(log(yield + 10) ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_equal(r1$means$response, exp(r1$means$estimate) - 10)
  rs <- suppressMessages(predictmeans(lm(sqrt(yield) ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_equal(rs$means$response, rs$means$estimate^2)
  rb <- suppressMessages(predictmeans(lm(pm_boxcox(yield, 0.3) ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_equal(rb$means$response, (0.3*rb$means$estimate + 1)^(1/0.3))
  # an untransformed response, and a link: as before
  expect_silent(ru <- predictmeans(lm(yield ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_null(ru$means$response)
  # mixed model, and covariatemeans()
  rl <- suppressMessages(predictmeans(lme4::lmer(log(yield) ~ nitro + (1 | Block), data = O), "nitro", prtnum = FALSE))
  expect_equal(rl$means$response, exp(rl$means$estimate))
  O$nit <- as.numeric(as.character(O$nitro))
  expect_message(cm <- covariatemeans(lm(log(yield) ~ nit * Variety, data = O), "Variety", covariate = "nit", plot = FALSE),
                 "back-transformed")
  # the back-transformed fitted values (1.x: exp(fitted + SE/2)): between the limits, at their geometric mean
  mfit <- lm(log(yield) ~ nit * Variety, data = O)
  nd <- data.frame(nit = cm$pltdf$xvar, Variety = cm$pltdf$Variety)
  expect_equal(cm$pltdf$Mean, unname(exp(predict(mfit, nd))), tolerance = 1e-8)
  expect_equal(cm$pltdf$Mean, sqrt(cm$pltdf$LL*cm$pltdf$UL), tolerance = 1e-8)
})

test_that("bias_adjust gives the mean of the back-transformed response: exact for log and sqrt", {
  O <- bt_oats()
  ml <- lm(log(yield) ~ Block + nitro, data = O)
  s2 <- sigma(ml)^2
  r <- suppressMessages(predictmeans(ml, "nitro", bias_adjust = TRUE, prtnum = FALSE))
  expect_equal(r$means$response, exp(r$means$estimate + s2/2), tolerance = 1e-9)
  expect_equal(r$means$response.std.error, exp(r$means$estimate + s2/2)*r$means$std.error, tolerance = 1e-9)
  expect_identical(r$settings$back_transform$estimate, "mean (bias adjusted)")
  ms <- lm(sqrt(yield) ~ Block + nitro, data = O)
  rs <- suppressMessages(predictmeans(ms, "nitro", bias_adjust = TRUE, prtnum = FALSE))
  expect_equal(rs$means$response, rs$means$estimate^2 + sigma(ms)^2, tolerance = 1e-9)
  # Gauss-Hermite: E f(Z) for a standard normal Z
  gh <- pm_gauss_hermite()
  expect_equal(sum(gh$w), 1)
  expect_equal(sum(gh$w*gh$z^4), 3, tolerance = 1e-10)
  expect_equal(sum(gh$w*exp(1.5*gh$z)), exp(1.5^2/2), tolerance = 1e-10)
  # a logit with a large random effect: the integral, where the second-order formula is off
  eta <- c(-1, 0.3); s <- 1.2
  bf <- pm_back_funs(stats::plogis, stats::plogis, NULL, FALSE)
  exact <- pm_bias_funs(bf, s^2)$h(eta)
  truth <- sapply(eta, function(e) integrate(function(z) plogis(e + s*z)*dnorm(z), -Inf, Inf)$value)
  expect_equal(exact, truth, tolerance = 1e-6)
  expect_gt(max(abs(pm_bias_funs(bf, s^2, "taylor")$h(eta) - truth)), 0.005)
  # a reciprocal near zero has no finite mean: the second-order formula, with a warning
  bfr <- pm_back_funs(function(x) 1/x, function(x) 1/x, NULL, FALSE)
  expect_warning(pm_bias_funs(bfr, 4)$h(0.5), "second-order formula")
})

test_that("summary() says what the back-transformed means are", {
  O <- bt_oats()
  out <- function(m, ...) capture.output(print(summary(suppressMessages(predictmeans(m, "nitro", prtnum = FALSE, ...)))))
  expect_true(any(grepl("'response' is the median of the original scale", out(lm(log(yield) ~ Block + nitro, data = O)))))
  O$cnt <- round(O$yield/10)
  expect_true(any(grepl("the inverse link: 'response' is the mean\\.", out(glm(cnt ~ Block + nitro, poisson, data = O)))))
  gm <- suppressMessages(lme4::glmer(cnt ~ nitro + (1 | Block), poisson, data = O))
  expect_true(any(grepl("'response' is the typical cluster", out(gm))))
  expect_true(any(grepl("'response' is the mean \\(bias adjusted\\)", out(gm, bias_adjust = TRUE))))
})

test_that("pm_check() gives the Box-Cox power the residuals suggest", {
  set.seed(2)
  d <- data.frame(A = factor(rep(1:4, each = 15)))
  d$y <- exp(1 + 0.3*as.numeric(d$A) + rnorm(60, sd = 0.5))
  m <- lm(y ~ A, d)
  bc <- pm_boxcox_profile(m)
  ref <- MASS::boxcox(y ~ A, data = d, lambda = seq(-3, 3, 0.01), plotit = FALSE)
  expect_equal(bc$lambda, ref$x[which.max(ref$y)])
  expect_equal(c(bc$lower, bc$upper), range(ref$x[ref$y > max(ref$y) - qchisq(0.95, 1)/2]), tolerance = 0.011)
  expect_identical(bc$name, "log(y)")
  pc <- pm_check(m)
  row <- pc[pc$check == "Transformation of the response (Box-Cox)", ]
  expect_true(row$look)
  expect_match(row$note, "log\\(y\\) is within the interval")
  # no transformation needed, a transformed response, a glm: no flag or no row
  O <- bt_oats()
  p0 <- pm_check(lm(yield ~ Block + nitro, data = O))
  expect_false(p0$look[p0$check == "Transformation of the response (Box-Cox)"])
  expect_null(pm_boxcox_profile(lm(log(yield) ~ Block + nitro, data = O)))
  expect_null(pm_boxcox_profile(glm(round(yield) ~ nitro, poisson, data = O)))
})
