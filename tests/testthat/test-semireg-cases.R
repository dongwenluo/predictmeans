## What predictmeans adds to a semi-parametric fit (dev/review-semireg.md,
## section 5): means and comparisons at a covariate value with
## Satterthwaite df in a split-plot design, tests of the smooth variance
## components of s(x, f1*f2), rate ratios between count curves, R2 shares.

q <- function(e) suppressWarnings(suppressMessages(e))
split_plot <- function(seed = 11, eff = c(a1 = 0, a2 = 0.5, a3 = 1)) {
  set.seed(seed)
  d <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  d$y <- 2 + sin(d$t/2) + eff[d$A]*d$t/10 + rnorm(6, 0, 0.5)[d$Block] +
    rnorm(18, 0, 0.4)[interaction(d$Block, d$A)] + rnorm(nrow(d), 0, 0.3)
  d
}

test_that("split plot with curves over time: means at t = 5 with Satterthwaite df for fixed + smooth coefficients", {
  d <- split_plot()
  fit <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), data = d, prt = FALSE))
  pm <- q(predictmeans(fit, "A", covariate = c(t = 5), pairwise = TRUE))
  expect_equal(pm$settings$df_method, "Satterthwaite (fixed + smooth coefficients)")
  # whole-plot comparisons: far fewer df than the 180 residual df, more than the 10 whole-plot error df
  expect_true(all(pm$means$df > 5 & pm$means$df < 15))
  expect_true(all(pm$pairs$df > 10 & pm$pairs$df < 30))
  expect_equal(pm$pairs$p.value, 2*pt(-abs(pm$pairs$estimate/pm$pairs$std.error), pm$pairs$df))
  tq <- qt(0.975, pm$means$df)
  expect_equal(pm$means$conf.high - pm$means$estimate, tq*pm$means$std.error)
  # for fixed effects alone the df are those of lmer's Satterthwaite (R/satterthwaite.R)
  L <- diag(length(lme4::fixef(fit$semer)))[2:3, , drop = FALSE]
  expect_equal(pm_sr_satt_ddf(fit, L, seq_len(ncol(L))), unname(pm_satt_ddf(L, pm_satt_parts(fit$semer))),
               tolerance = 1e-2)
  # generalized fits: the inner-outer df of the design's strata, as glmer gets
  # for the same design (the smooth blocks are not strata); z without strata
  d$cnt <- rpois(nrow(d), exp(d$y/2))
  gf <- q(semireg(cnt ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), family = poisson, data = d,
                  prt = FALSE))
  gm <- q(predictmeans(gf, "A", covariate = c(t = 5), pairwise = TRUE))
  io <- pm_inner_outer_df(q(lme4::glmer(cnt ~ A + A:t + (1 | Block) + (1 | Block:A), family = poisson, data = d)))
  expect_equal(gm$means$df, rep(attr(io, "terms")[["A"]], 3))
  expect_equal(gm$means$df, rep(10, 3))                   # 18 whole plots - 6 blocks - 3 treatments + 1
  expect_equal(gm$pairs$df, rep(10, 3))
  expect_equal(gm$settings$df_method, "inner-outer (design strata; smooths within)")
  expect_equal(gm$pairs$p.value, 2*pt(-abs(gm$pairs$estimate/gm$pairs$std.error), 10))
  sp <- q(semipred(gf, "A", "t"))$pred_df
  expect_equal(sp$df, rep(10, nrow(sp)))
  g0 <- q(semireg(cnt ~ s(t, by = A, bs = "os", k = 6), family = poisson, data = d, prt = FALSE))
  g0m <- q(predictmeans(g0, "A", covariate = c(t = 5), pairwise = TRUE))
  expect_true(all(is.infinite(g0m$means$df)))
  expect_equal(g0m$settings$df_method, "normal (semireg)")
  expect_null(q(semipred(g0, "A", "t"))$pred_df$df)
})

test_that("s(x, f1*f2): a variance per deviation term; anova() tests which deviation curves are needed", {
  set.seed(1)
  d <- expand.grid(x = seq(0, 1, length = 25), f1 = factor(c("p", "q")), f2 = factor(c("u", "v", "w")), rep = 1:4)
  d$id <- interaction(d$f1, d$f2, d$rep)
  # the shape differs by f1 only
  d$y <- sin(2*pi*d$x) + (d$f1 == "q")*1.5*sin(pi*d$x)^2 + rnorm(nlevels(d$id), 0, 0.2)[d$id] +
    rnorm(nrow(d), 0, 0.25)
  fit <- q(semireg(y ~ s(x, f1*f2, bs = "os", k = 8) + (1 | id), data = d, prt = FALSE))
  sm <- anova(fit)$smooth
  expect_setequal(sm$term, c("s_x", "s_x_f1", "s_x_f2", "s_x_f1_f2"))
  expect_lt(sm$p.value[sm$term == "s_x_f1"], 0.01)
  expect_gt(sm$p.value[sm$term == "s_x_f2"], 0.05)
  expect_gt(sm$p.value[sm$term == "s_x_f1_f2"], 0.05)
  vc <- varcomp(fit)
  expect_true(all(c("s_x_f1.(Intercept)", "s_x_f2.(Intercept)", "s_x_f1_f2.(Intercept)") %in% rownames(vc)))
  r2 <- R2_model(fit)
  expect_true(all(c("s_x", "s_x_f1", "id") %in% names(r2)))
})

test_that("counts with a negative binomial and zero inflation (glmmTMB engine): rate ratio between curves", {
  set.seed(5)
  d <- expand.grid(x = seq(0, 3, length = 30), trt = factor(c("ctl", "trt")), site = factor(1:8))
  mu <- exp(1 + sin(d$x) + (d$trt == "trt")*0.5*d$x + rnorm(8, 0, 0.3)[d$site])
  d$y <- rnbinom(nrow(d), mu = mu, size = 3)
  d$y[runif(nrow(d)) < 0.1] <- 0
  fit <- q(semireg_tmb(y ~ s(x, by = trt, bs = "os", k = 6) + (1 | site), ziformula = ~ 1,
                       family = glmmTMB::nbinom2, data = d, prt = FALSE))
  out <- q(semipred(fit, "trt", "x", contr = c("trt", "ctl"), band = "simultaneous"))$pred_df
  expect_true(all(c("ratio", "ratio.LL", "ratio.UL", "ratio.includes1", "significant") %in% names(out)))
  expect_equal(out$ratio, exp(out$diff), tolerance = 1e-10)
  expect_equal(out$ratio.includes1, out$includes0)
  # the ratio grows with x (true log ratio 0.5 x) and differs from 1 at large x
  expect_gt(cor(out$x, log(out$ratio)), 0.9)
  expect_true(any(out$significant[out$x > 2]))
  expect_false(is.null(attr(out, "significant")))
})

test_that("semireg_tmb(): REML by default for gaussian fits (as documented), Satterthwaite df as the lme4 engine", {
  d <- split_plot()
  f <- y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A)
  a <- q(semireg(f, data = d, prt = FALSE))
  b <- q(semireg_tmb(f, data = d, prt = FALSE))
  expect_true(b$semer$modelInfo$REML)
  expect_equal(-b$semer$fit$objective, as.numeric(logLik(a$semer)), tolerance = 1e-6)
  pa <- q(predictmeans(a, "A", covariate = c(t = 5), pairwise = TRUE))
  pb <- q(predictmeans(b, "A", covariate = c(t = 5), pairwise = TRUE))
  expect_equal(pb$settings$df_method, "Satterthwaite (fixed + smooth coefficients)")
  expect_equal(pb$means$df, pa$means$df, tolerance = 1e-3)
  expect_equal(pb$pairs$df, pa$pairs$df, tolerance = 1e-3)
  expect_equal(pb$pairs$p.value, pa$pairs$p.value, tolerance = 1e-4)
  # other families: maximum likelihood (Laplace), and the inner-outer df of
  # the design's strata as with the lme4 engine (A varies within blocks: the
  # residual stratum)
  d$cnt <- rpois(nrow(d), exp(d$y/2))
  g <- q(semireg_tmb(cnt ~ s(t, by = A, bs = "os", k = 6) + (1 | Block), family = poisson, data = d, prt = FALSE))
  expect_false(g$semer$modelInfo$REML)
  g4 <- q(semireg(cnt ~ s(t, by = A, bs = "os", k = 6) + (1 | Block), family = poisson, data = d, prt = FALSE))
  dg <- q(predictmeans(g, "A", covariate = c(t = 5)))$means$df
  expect_equal(dg, q(predictmeans(g4, "A", covariate = c(t = 5)))$means$df)
  expect_equal(dg, rep(nrow(d) - 6 - 5, 3))                # N - rank([X, blocks])
})

test_that("semireg_tmb() with AR(1) residuals: means at a time with Satterthwaite df", {
  set.seed(9)
  d <- expand.grid(t = 1:8, id = factor(1:12))
  d$A <- factor(rep(c("a1", "a2", "a3"), each = 32))
  d$y <- sin(d$t/3) + rnorm(12, 0, 0.4)[d$id] +
    unlist(lapply(1:12, function(i) as.numeric(stats::arima.sim(list(ar = 0.6), 8, sd = 0.4))))
  d$tf <- factor(d$t)
  fit <- q(semireg_tmb(y ~ s(t, by = A, bs = "os", k = 4) + (1 | id) + ar1(tf + 0 | id), dispformula = ~0, data = d,
                       prt = FALSE))
  pm <- q(predictmeans(fit, "A", covariate = c(t = 4), pairwise = TRUE))
  expect_true(all(is.finite(pm$means$df) & pm$means$df > 2 & pm$means$df < 60))
  expect_true(all(is.finite(pm$pairs$df) & pm$pairs$df > 2 & pm$pairs$df < 60))
  expect_equal(pm$pairs$p.value, 2*pt(-abs(pm$pairs$estimate/pm$pairs$std.error), pm$pairs$df))
})
