## permlmer(): an LRT of 0 up to optimizer noise must not look extreme
test_that("permlmer() gives p = 1 when the tested variance is estimated at 0", {
  set.seed(22)
  d <- expand.grid(sub = 1:3, wp = 1:4, block = 1:6)
  d$A <- factor(d$wp); d$block <- factor(d$block); d$wpid <- interaction(d$block, d$wp)
  d$y <- rnorm(6)[d$block] + rnorm(nrow(d))
  ctl <- lmerControl(check.conv.singular = "ignore")
  m0 <- lmer(y ~ A + (1 | block), data = d, control = ctl)
  m1 <- lmer(y ~ A + (1 | block) + (1 | wpid), data = d, control = ctl)
  lrt <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  expect_lt(abs(lrt), 1e-4)
  res <- suppressMessages(permlmer(m0, m1, nperm = 19, ncore = 1, seed = 1))
  expect_equal(res$`Perm-p`[2], 1)
})

test_that("permlmer() and permmodels(method = \"raw\") use the prior weights and offset of lmer fits", {
  skip_on_cran()
  set.seed(6121)
  d <- data.frame(g = gl(12, 5), x = rnorm(60), w = runif(60, 0.3, 3), o = rnorm(60))
  d$y <- 1 + 0.4 * d$x + rnorm(12)[d$g] + rnorm(60)/sqrt(d$w) + d$o
  d$y[3] <- NA
  fw <- lmer(y ~ x + (1 | g), d, weights = w, offset = o)
  # before 2.0.0 the refit from the formula's variables stopped ("object 'w' not found")
  for (ty in c("I", "II", "III")) {
    r <- permmodels(fw, nperm = 19, seed = 1, method = "raw", type = ty)
    expect_true(all(r$ANOVA["x", "Perm_p_value"] >= 1/20), label = ty)
  }
  # the identity permutation gives the data back: the refit is the fit itself
  # (weights, offset and the NA row kept)
  r <- permmodels(fw, perms = matrix(seq_len(nobs(fw)), ncol = 1), method = "raw")
  expect_equal(r$permlist[[1]]$coef, lme4::fixef(fw), tolerance = 1e-5)
  # the decorrelation: V = Z Lambda Lambda' Z' + diag(1/w), and W W' = V
  A <- as.matrix(lme4::getME(fw, "Z") %*% lme4::getME(fw, "Lambda"))
  W <- predictmeans:::pm_whiten(A, predictmeans:::pm_mer_resvar(fw))
  V <- tcrossprod(A) + diag(1/weights(fw))
  Wm <- W$mult(diag(nrow(V)))
  expect_equal(tcrossprod(Wm), V, ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(W$solve(Wm), diag(nrow(V)), ignore_attr = TRUE, tolerance = 1e-8)
  expect_null(predictmeans:::pm_mer_resvar(lmer(y ~ x + (1 | g), d)))
  # weights and offset unchanged under permutation: the permuted data of the
  # null model, refitted, give the null model back when the permutation is the identity
  f0 <- lmer(y ~ 1 + (1 | g), d, weights = w, offset = o, REML = FALSE)
  f1 <- lmer(y ~ x + (1 | g), d, weights = w, offset = o, REML = FALSE)
  id <- matrix(seq_len(nobs(f1)), ncol = 1)
  a <- permlmer(f0, f1, perms = id[, rep(1, 3)])
  expect_equal(a$`Perm-p`[2], 1)
  a <- permlmer(f0, f1, nperm = 19, seed = 1)
  expect_gte(a$`Perm-p`[2], 1/20)
})

test_that("permlmer(): different fixed effects are tested by the ML likelihood ratio that anova() prints", {
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  m0 <- lmer(yield ~ Variety + (1 | Block/Variety), Oats)
  m1 <- lmer(yield ~ nitro + Variety + (1 | Block/Variety), Oats)
  ml <- function(m) as.numeric(logLik(refitML(m)))
  chisq_ml <- 2 * (ml(m1) - ml(m0))
  # before 2.0.0 the permuted statistic was logLik(REML = FALSE) at the REML estimates: 66.43 here
  expect_gt(abs(2 * (as.numeric(logLik(m1, REML = FALSE)) - as.numeric(logLik(m0, REML = FALSE))) - chisq_ml), 0.05)
  id <- matrix(seq_len(nobs(m1)), ncol = 1)
  a <- permlmer(m0, m1, perms = id[, rep(1, 3)])
  expect_equal(a$Chisq[2], chisq_ml, tolerance = 1e-6)
  # the identity permutation reproduces the observed statistic exactly: the same ML statistic is permuted
  expect_equal(a$`Perm-p`[2], 1)
  # ML fits given: the same result
  expect_equal(permlmer(refitML(m0), refitML(m1), perms = id[, rep(1, 3)])$Chisq[2], chisq_ml, tolerance = 1e-6)
  a <- permlmer(m0, m1, nperm = 19, seed = 1)
  expect_equal(a$`Perm-p`[2], 1/20)
})

test_that("refits of lmer fits with a missing response work with every lme4 version", {
  # lme4 1.1-35 refit() dropped the NA rows twice ("replacement has 58 rows, data has 59")
  set.seed(3)
  d <- data.frame(g = gl(10, 6), x = rnorm(60))
  d$y <- d$x + rnorm(10)[d$g] + rnorm(60)
  d$y[c(3, 17)] <- NA
  m <- lmer(y ~ x + (1 | g), d)
  r <- predictmeans:::pm_mer_refit(m, lme4::getME(m, "y"))
  expect_equal(lme4::fixef(r), lme4::fixef(m), tolerance = 1e-6)
  expect_identical(nobs(r), nobs(m))
  m0 <- lmer(y ~ 1 + (1 | g), d, REML = FALSE)
  res <- suppressMessages(permlmer(m0, update(m, REML = FALSE), nperm = 9, seed = 1))
  expect_true(is.finite(res$`Perm-p`[2]))
})
