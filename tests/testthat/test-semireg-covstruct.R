## glmmTMB covariance structures (ar1(), ou(), cs(), diag(), ...) next to
## smooths in semireg_tmb(): kept as written (findbars() returns only the
## bar inside them, and a formula rebuilt from it had turned ar1() into an
## unstructured term without a message).

q <- function(e) suppressWarnings(suppressMessages(e))
ar1_data <- function(seed = 7, ns = 20, nt = 8) {
  set.seed(seed)
  d <- expand.grid(time = 1:nt, id = factor(1:ns))
  d$trt <- factor(rep(c("a", "b"), each = nt*ns/2))
  e <- unlist(lapply(seq_len(ns), function(i) as.numeric(stats::arima.sim(list(ar = 0.7), nt, sd = 0.3))))
  d$y <- sin(d$time/3) + (d$trt == "b")*0.05*d$time + e + rnorm(nrow(d), 0, 0.2)
  d$tf <- factor(d$time)
  d$times <- glmmTMB::numFactor(d$time)
  d
}
struct_codes <- function(m) vapply(glmmTMB::VarCorr(m)$cond, function(v) names(attr(v, "blockCode")), "")

test_that("semireg_tmb() keeps ar1(), ou(), cs(), diag() terms: fit, refit for anova(), the same model as smoothZ", {
  d <- ar1_data()
  for (st in c("ar1", "cs", "diag", "ou")) {
    term <- if (st == "ou") "ou(times + 0 | id)" else paste0(st, "(tf + 0 | id)")
    f1 <- stats::as.formula(paste("y ~ s(time, by = trt, bs = 'os', k = 4) +", term))
    f2 <- stats::as.formula(paste("y ~ time*trt +", term))
    # dispformula = ~0: the structure is the residual covariance (with a
    # nugget, cs() is not identified and glmmTMB flags the Hessian)
    a <- q(semireg_tmb(f1, dispformula = ~0, data = d))
    b <- q(semireg_tmb(f2, smoothZ = list(sm = smZ(time, k = 4, by = trt, group = TRUE)), dispformula = ~0, data = d))
    expect_true(st %in% struct_codes(a$semer), label = paste(st, "formula form"))
    expect_true(st %in% struct_codes(b$semer), label = paste(st, "smoothZ form"))
    expect_equal(a$semer$fit$objective, b$semer$fit$objective, tolerance = 1e-6, label = st)
    # anova()'s refit without a smooth block keeps the structure
    red <- q(pm_sr_refit(a, drop = names(a$refit$smoothZ)[1]))
    expect_true(st %in% struct_codes(red), label = paste(st, "refit"))
  }
  # ou() with equally spaced times is ar1()
  a <- q(semireg_tmb(y ~ s(time, by = trt, bs = "os", k = 4) + ar1(tf + 0 | id), dispformula = ~0, data = d))
  o <- q(semireg_tmb(y ~ s(time, by = trt, bs = "os", k = 4) + ou(times + 0 | id), dispformula = ~0, data = d))
  expect_equal(a$semer$fit$objective, o$semer$fit$objective, tolerance = 1e-6)
})

test_that("pure AR1 residuals (dispformula = ~0): curves and means with finite SEs", {
  d <- ar1_data()
  fit <- q(semireg_tmb(y ~ s(time, by = trt, bs = "os", k = 4) + ar1(tf + 0 | id), dispformula = ~0, data = d))
  expect_true("ar1" %in% struct_codes(fit$semer))
  cur <- q(semipred(fit, "trt", "time"))$pred_df
  expect_true(all(is.finite(cur$SE) & cur$SE > 0))
  pm <- q(predictmeans(fit, "trt", covariate = c(time = 6), pairwise = TRUE))
  expect_true(all(is.finite(pm$means$std.error) & pm$means$std.error > 0))
})

test_that("semireg() (lme4) refuses covariance structures; bootmodels() keeps them in its null model", {
  d <- ar1_data()
  expect_error(semireg(y ~ s(time) + ar1(tf + 0 | id), data = d), "semireg_tmb")
  m <- q(glmmTMB::glmmTMB(y ~ trt + ar1(tf + 0 | id), data = d))
  expect_true("ar1" %in% struct_codes(pm_boot_parts(m)$null_fit))
  expect_equal(pm_random_terms(y ~ x + (1 | g) + ar1(tf + 0 | id)), c("(1 | g)", "ar1(tf + 0 | id)"))
  expect_equal(deparse(pm_nobars(y ~ x + (1 | g) + ar1(tf + 0 | id))), "y ~ x")
})
