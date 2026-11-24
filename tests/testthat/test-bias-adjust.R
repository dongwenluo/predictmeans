## Bias-adjusted back-transformed means and averaging on the response scale
## (predictmeans(bias_adjust =, sigma =, average =)), against emmeans

glmm_data <- function() {
  set.seed(3)
  d <- expand.grid(A = factor(1:3), B = factor(1:2), rep = 1:4, block = factor(1:6))
  d$y <- rpois(nrow(d), exp(0.5 + 0.4 * as.integer(d$A) + 0.8 * (d$B == "2") + rnorm(6, 0, 0.6)[d$block]))
  d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.3 * as.integer(d$A) + rnorm(6, 0, 0.8)[d$block]))
  d$w <- exp(1 + 0.2 * as.integer(d$A) + rnorm(6, 0, 0.3)[d$block] + rnorm(nrow(d), 0, 0.2))
  d
}
resp <- function(pm) pm$means[, c("response", "response.std.error")]

test_that("bias adjustment by hand: the mean of h(eta + u); 'taylor' is h(eta) + sigma^2/2 h''(eta)", {
  d <- glmm_data()
  m <- lme4::glmer(y ~ A * B + (1 | block), family = poisson, data = d)
  s2 <- as.data.frame(lme4::VarCorr(m))$vcov[1]
  pm <- predictmeans(m, "A", bias_adjust = TRUE, plot = FALSE)
  eta <- pm$means$estimate
  # log link: E exp(eta + u) = exp(eta + sigma^2/2), exactly
  expect_equal(pm$means$response, exp(eta + s2 / 2), tolerance = 1e-8)
  expect_equal(pm$means$response.std.error, exp(eta + s2 / 2) * pm$means$std.error, tolerance = 1e-8)
  # the second-order formula (emmeans' bias.adjust), lower
  pt <- predictmeans(m, "A", bias_adjust = "taylor", plot = FALSE)
  expect_equal(pt$means$response, exp(eta) * (1 + s2 / 2), tolerance = 1e-6)          # h'' = exp
  expect_equal(pt$means$response.std.error, exp(eta) * (1 + s2 / 2) * pt$means$std.error, tolerance = 1e-6)
  expect_true(all(pt$means$response < pm$means$response))
  expect_error(predictmeans(m, "A", bias_adjust = "yes", plot = FALSE), "TRUE, FALSE or \"taylor\"")
  expect_equal(pm$settings$bias_adjust, sqrt(s2))
  expect_identical(pm$settings$average, "link")
  p0 <- predictmeans(m, "A", plot = FALSE)
  expect_true(is.na(p0$settings$bias_adjust))
  expect_equal(p0$means$response, exp(eta))
  # a given sigma
  expect_equal(predictmeans(m, "A", bias_adjust = TRUE, sigma = 1, plot = FALSE)$means$response,
               exp(eta + 0.5), tolerance = 1e-8)
  expect_equal(predictmeans(m, "A", bias_adjust = "taylor", sigma = 1, plot = FALSE)$means$response,
               exp(eta) * 1.5, tolerance = 1e-6)
})

test_that("bias_adjust = 'taylor' equals emmeans (bias.adjust, regrid) for Poisson, binomial and a log response", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  d <- glmm_data()
  m <- lme4::glmer(y ~ A * B + (1 | block), family = poisson, data = d)
  s <- sqrt(as.data.frame(lme4::VarCorr(m))$vcov[1])
  e1 <- summary(emmeans::emmeans(m, ~ A, type = "response", bias.adjust = TRUE, sigma = s))
  expect_equal(unname(as.matrix(resp(predictmeans(m, "A", bias_adjust = "taylor", plot = FALSE)))),
               unname(as.matrix(e1[, c("rate", "SE")])), tolerance = 1e-6)
  e2 <- summary(emmeans::emmeans(emmeans::regrid(emmeans::ref_grid(m)), ~ A))
  expect_equal(unname(as.matrix(resp(predictmeans(m, "A", average = "response", plot = FALSE)))),
               unname(as.matrix(e2[, c("rate", "SE")])), tolerance = 1e-6)
  mb <- lme4::glmer(cbind(r, 10 - r) ~ A + B + (1 | block), family = binomial, data = d)
  sb <- sqrt(as.data.frame(lme4::VarCorr(mb))$vcov[1])
  eb <- summary(emmeans::emmeans(emmeans::regrid(emmeans::ref_grid(mb), bias.adjust = TRUE, sigma = sb), ~ A))
  pb <- predictmeans(mb, "A", bias_adjust = "taylor", average = "response", plot = FALSE)
  expect_equal(unname(as.matrix(resp(pb))), unname(as.matrix(eb[, c("prob", "SE")])), tolerance = 1e-6)
  ml <- lme4::lmer(log(w) ~ A + B + (1 | block), data = d)
  sl <- sqrt(sum(as.data.frame(lme4::VarCorr(ml))$vcov))          # default: random + residual variance
  el <- summary(emmeans::emmeans(ml, ~ A, type = "response", bias.adjust = TRUE, sigma = sl))
  pl <- predictmeans(ml, "A", trans = exp, bias_adjust = "taylor", plot = FALSE)
  expect_equal(unname(as.matrix(resp(pl))), unname(as.matrix(el[, c("response", "SE")])), tolerance = 1e-6)
})

test_that("identity scale and random slopes", {
  d <- glmm_data()
  fl <- lm(w ~ A + B, data = d)
  expect_warning(predictmeans(fl, "A", bias_adjust = TRUE, plot = FALSE), "back-transformed")
  d$x <- as.integer(d$B) - 1.5
  ms <- suppressWarnings(lme4::glmer(y ~ A + x + (x | block), family = poisson, data = d))
  expect_error(predictmeans(ms, "A", bias_adjust = TRUE, plot = FALSE), "sigma")
  pm <- predictmeans(ms, "A", bias_adjust = TRUE, sigma = 0.5, plot = FALSE)
  expect_equal(pm$means$response, exp(pm$means$estimate + 0.25 / 2), tolerance = 1e-8)
})
