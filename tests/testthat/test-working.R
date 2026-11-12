## The working (linearised) model of GLMs and GLMMs, and permlmer() for GLMMs

glmm_data <- function() {
  set.seed(3)
  d <- expand.grid(sub = 1:3, wp = 1:4, block = 1:6)
  d$A <- factor(d$wp); d$B <- factor(d$sub); d$block <- factor(d$block)
  d$wpid <- interaction(d$block, d$wp)
  d$t <- runif(nrow(d), 1, 3)
  d$y <- rpois(nrow(d), d$t * exp(0.5 + 0.3 * as.numeric(d$A) + rnorm(6, 0, .5)[d$block]))
  d$yb <- rbinom(nrow(d), 1, plogis(-1 + 0.3 * as.numeric(d$A) + rnorm(6, 0, .5)[d$block]))
  d$k <- rbinom(nrow(d), 10, plogis(-1 + 0.3 * as.numeric(d$A)))
  d
}

test_that("working response and variances of glm, glmer and glmmTMB fits agree", {
  d <- glmm_data()
  # glm: z = eta (without offset) + working residuals; resvar = phi / working weights
  m <- glm(y ~ A + offset(log(t)), family = poisson, data = d)
  W <- pm_working(m)
  expect_equal(W$z, unname(m$linear.predictors - log(d$t) + residuals(m, "working")))
  # (weights(m, "working") come from glm's last iteration, not quite at the fit)
  expect_equal(W$resvar, unname(1 / fitted(m)))
  expect_equal(W$resvar, unname(1 / weights(m, "working")), tolerance = 1e-3)
  expect_equal(W$xb, unname(as.vector(model.matrix(m) %*% coef(m))))
  # glmer with an offset, and glmmTMB with the same model
  g <- lme4::glmer(y ~ A + offset(log(t)) + (1 | block), family = poisson, data = d)
  Wg <- pm_working(g)
  eta <- g@resp$eta
  expect_equal(Wg$z, eta - log(d$t) + (d$y - g@resp$mu) / exp(eta))
  expect_equal(Wg$resvar, 1 / g@resp$sqrtWrkWt()^2)
  tm <- glmmTMB::glmmTMB(y ~ A + offset(log(t)) + (1 | block), family = poisson, data = d)
  Wt <- pm_working(tm)
  expect_equal(Wt$z, Wg$z, tolerance = 1e-4)
  expect_equal(Wt$resvar, Wg$resvar, tolerance = 1e-4)
  expect_equal(as.matrix(tcrossprod(Wt$A)), as.matrix(tcrossprod(Wg$A)), tolerance = 1e-3)
  # binomial with trials: proportions, trials as weights
  gb <- lme4::glmer(cbind(k, 10 - k) ~ A + (1 | block), family = binomial, data = d)
  tb <- glmmTMB::glmmTMB(cbind(k, 10 - k) ~ A + (1 | block), family = binomial, data = d)
  expect_equal(pm_working(gb)$y, d$k / 10)
  expect_equal(pm_working(tb)$resvar, pm_working(gb)$resvar, tolerance = 1e-3)
  expect_false(pm_working(gb)$binary)
  expect_true(pm_working(lme4::glmer(yb ~ A + (1 | block), family = binomial, data = d))$binary)
  # the working formula: response replaced, offsets dropped
  expect_identical(deparse(pm_working_formula(y ~ A + offset(log(t)) + (1 | block))), ".pm_z ~ A + (1 | block)")
  expect_identical(deparse(pm_working_formula(y ~ offset(log(t)))), ".pm_z ~ 1")
})

test_that("glmmTMB family variances match simulated responses", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  dd <- data.frame(x = runif(n))
  mu <- exp(1 + dd$x)
  dd$nb2 <- rnbinom(n, mu = mu, size = 2)
  dd$nb1 <- rnbinom(n, mu = mu, size = mu / 1.5)
  dd$gam <- rgamma(n, shape = 4, rate = 4 / mu)
  dd$bet <- rbeta(n, plogis(dd$x) * 10, (1 - plogis(dd$x)) * 10)
  fits <- list(glmmTMB::glmmTMB(nb2 ~ x, family = glmmTMB::nbinom2, data = dd),
               glmmTMB::glmmTMB(nb1 ~ x, family = glmmTMB::nbinom1, data = dd),
               glmmTMB::glmmTMB(gam ~ x, family = Gamma(link = "log"), data = dd),
               glmmTMB::glmmTMB(bet ~ x, family = glmmTMB::beta_family(), data = dd))
  for (m in fits) {
    W <- pm_working(m)
    # Var(y) = resvar x (dmu/deta)^2 on the response scale
    v_model <- W$resvar * (m$modelInfo$family$mu.eta(W$eta))^2
    v_sim <- apply(as.matrix(simulate(m, nsim = 3000, seed = 2)), 1, var)
    expect_equal(mean(v_sim / v_model), 1, tolerance = 0.03)
  }
})

test_that("permlmer() tests GLMM random and fixed effects on the working scale", {
  skip_on_cran()
  d <- glmm_data()
  g0 <- lme4::glmer(y ~ A + B + (1 | block), family = poisson, data = d)
  g1 <- lme4::glmer(y ~ A + B + (1 | block) + (1 | wpid), family = poisson, data = d)
  set.seed(5); a <- runif(1)
  set.seed(5); r <- permlmer(g0, g1, nperm = 19, ncore = 1, seed = 1); b <- runif(1)
  expect_identical(a, b)                         # the seed leaves the user's stream alone
  expect_identical(attr(r, "method"), "working")
  expect_identical(r$`Perm-p`[2], permglmer(g0, g1, nperm = 19, ncore = 1, seed = 1)$`Perm-p`[2])
  # the observed statistic: REML LRT of the weighted working LMMs
  W <- pm_working(g0)
  dz <- d; dz$z <- W$z; dz$w <- 1 / W$resvar
  l0 <- lme4::lmer(z ~ A + B + (1 | block), data = dz, weights = w)
  l1 <- lme4::lmer(z ~ A + B + (1 | block) + (1 | wpid), data = dz, weights = w)
  obs <- round(max(as.numeric(2 * (logLik(l1) - logLik(l0))), 0), 4)
  expect_true(r$`Perm-p`[2] >= 1 / 20 && r$`Perm-p`[2] <= 1)
  expect_true(obs > 0)
  # fixed effects
  g0f <- lme4::glmer(y ~ B + (1 | block), family = poisson, data = d)
  rf <- permlmer(g0f, g0, nperm = 19, ncore = 1, seed = 1)
  expect_true(rf$`Perm-p`[2] > 0)
  # 0/1 responses: the quantile route
  b0 <- lme4::glmer(yb ~ A + (1 | block), family = binomial, data = d)
  b1 <- lme4::glmer(yb ~ A + (1 | block) + (1 | wpid), family = binomial, data = d)
  rb <- permlmer(b0, b1, nperm = 9, ncore = 1, seed = 1)
  expect_identical(attr(rb, "method"), "quantile")
  expect_error(permlmer(g0, g1, nperm = 9, ncore = 1, method = "quantile"), "0/1 responses")
  expect_error(permlmer(g1, g0f, nperm = 9, ncore = 1), "nested")
  # glmmTMB
  t0 <- glmmTMB::glmmTMB(y ~ A + B + (1 | block), family = poisson, data = d)
  t1 <- glmmTMB::glmmTMB(y ~ A + B + (1 | block) + (1 | wpid), family = poisson, data = d)
  rt <- permglmmTMB(t0, t1, nperm = 19, ncore = 1, seed = 1)
  expect_equal(rt$`Perm-p`[2], r$`Perm-p`[2])
  expect_error(permlmer(g0, t1, nperm = 9, ncore = 1), "Both models")
})

test_that("R2_glmm() of glmmTMB models uses the right working variance", {
  skip_on_cran()
  set.seed(4)
  d <- expand.grid(rep = 1:5, A = factor(1:4), g = factor(1:12))
  d$t <- runif(nrow(d), 1, 4)
  mu <- exp(1 + 0.3 * as.numeric(d$A) + rnorm(12, 0, .5)[d$g])
  d$po <- rpois(nrow(d), d$t * mu)
  d$nb <- rnbinom(nrow(d), mu = mu, size = 2)
  # an offset: glmmTMB and glmer give the same R2 (glmmTMB ignored the offset)
  expect_equal(pm_r2_piepho(glmmTMB::glmmTMB(po ~ A + offset(log(t)) + (1 | g), family = poisson, data = d)),
               pm_r2_piepho(lme4::glmer(po ~ A + offset(log(t)) + (1 | g), family = poisson, data = d)),
               tolerance = 1e-3)
  # nbinom2: residual variance mu (1 + mu/theta) / mu^2 on the log scale
  m <- glmmTMB::glmmTMB(nb ~ A + (1 | g), family = glmmTMB::nbinom2, data = d)
  W <- pm_working(m)
  expect_equal(W$resvar, (1 + W$mu / sigma(m)) / W$mu)
  expect_gt(pm_r2_piepho(m)[["Total"]], 20)
})
