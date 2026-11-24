# Kenward-Roger adjusted covariance and F test, ported from pbkrtest.

test_that("pm_vcov_kr(information = \"expected\") and pm_kr_ftest() equal pbkrtest", {
  skip_if_not_installed("pbkrtest")
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  set.seed(7)
  Ou <- Oats[-sample(nrow(Oats), 9), ]
  data("sleepstudy", package = "lme4", envir = environment())
  models <- list(
    nested = lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou),
    slopes = lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy))
  for (nm in names(models)) {
    m <- models[[nm]]
    ref <- pbkrtest::vcovAdj(m)
    new <- predictmeans:::pm_vcov_kr(m, "expected")
    expect_equal(unclass(new)[, ], as.matrix(ref), ignore_attr = TRUE, tolerance = 1e-10, label = nm)
    expect_equal(attr(new, "W"), as.matrix(attr(ref, "W")), ignore_attr = TRUE, tolerance = 1e-10)
    for (k in seq_along(attr(ref, "P"))) {
      expect_equal(attr(new, "P")[[k]], as.matrix(attr(ref, "P")[[k]]), ignore_attr = TRUE, tolerance = 1e-10)
    }
    L <- diag(length(lme4::fixef(m)))[-1, , drop = FALSE]
    kr_ref <- unlist(pbkrtest::KRmodcomp(m, L)$test["Ftest", c("stat", "ddf", "F.scaling", "p.value")])
    kr_new <- predictmeans:::pm_kr_ftest(m, L, PhiA = new)[c("stat", "ddf", "F.scaling", "p.value")]
    expect_equal(kr_new, kr_ref, ignore_attr = TRUE, tolerance = 1e-10)
  }
})

test_that("Kenward-Roger df and SEs in predictmeans() don't need pbkrtest", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  pm <- predictmeans(fm, "nitro", vcov_adjust = TRUE, plot = FALSE)
  expect_equal(pm$settings$df_method, "Kenward-Roger")
  expect_true(all(pm$means$df > 0))
  # balanced: the observed and the expected information give the same df
  pe <- predictmeans(fm, "nitro", df_method = "KR-expected", plot = FALSE)
  expect_equal(pe$pairs$df, pm$pairs$df, tolerance = 1e-4)
  expect_equal(pe$settings$df_method, "KR-expected")
})

test_that("Kenward-Roger uses the prior weights of lmer fits", {
  set.seed(4817)
  d <- data.frame(g = gl(10, 6), x = rnorm(60), w = runif(60, 0.3, 3))
  d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(60)/sqrt(d$w)
  d$iw <- 1/d$w
  L <- matrix(c(0, 1), 1)
  kr <- function(m, info = "expected") predictmeans:::pm_kr_ftest(m, L, predictmeans:::pm_vcov_kr(m, info))
  # constant weights only rescale sigma: same adjusted F test as no weights
  # (pbkrtest and lmerTest ignore the weights here and give other df)
  f1 <- lmer(y ~ x + (1 | g), d)
  f2 <- lmer(y ~ x + (1 | g), d, weights = rep(2, 60))
  for (info in c("expected", "observed")) {
    expect_equal(kr(f2, info), kr(f1, info), tolerance = 1e-6, label = info)
  }
  # unequal weights: as lme() with varFixed(~ 1/w), the same model
  fw <- lmer(y ~ x + (1 | g), d, weights = w)
  lw <- nlme::lme(y ~ x, random = ~ 1 | g, data = d, weights = nlme::varFixed(~ iw))
  expect_equal(unclass(predictmeans:::pm_vcov_kr(fw, "expected"))[, ],
               unclass(predictmeans:::pm_vcov_kr_nlme(lw, "expected"))[, ], ignore_attr = TRUE, tolerance = 1e-4)
  expect_equal(kr(fw)[c("stat", "ddf")],
               predictmeans:::pm_kr_ftest(lw, L, predictmeans:::pm_vcov_kr_nlme(lw, "expected"))[c("stat", "ddf")],
               tolerance = 1e-4)
  expect_gt(kr(fw, "observed")[["ddf"]], 0)
})

test_that("Kenward-Roger with an offset equals the fit to y - offset (and pbkrtest)", {
  set.seed(9133)
  d <- data.frame(g = gl(10, 6), x = rnorm(60), o = rnorm(60))
  d$y <- 1 + 0.5 * d$x + rnorm(10)[d$g] + rnorm(60) + d$o
  L <- matrix(c(0, 1), 1)
  fo <- lmer(y ~ x + (1 | g), d, offset = o)
  fy <- lmer(I(y - o) ~ x + (1 | g), d)
  for (info in c("expected", "observed", "average")) {
    expect_equal(predictmeans:::pm_kr_ftest(fo, L, predictmeans:::pm_vcov_kr(fo, info)),
                 predictmeans:::pm_kr_ftest(fy, L, predictmeans:::pm_vcov_kr(fy, info)),
                 tolerance = 1e-6, label = info)
  }
  skip_if_not_installed("pbkrtest")
  expect_equal(predictmeans:::pm_kr_ftest(fo, L, predictmeans:::pm_vcov_kr(fo, "expected"))[["ddf"]],
               pbkrtest::KRmodcomp(fo, L)$test["Ftest", "ddf"], tolerance = 1e-8)
})

test_that("an invalid Kenward-Roger df stops, and a failed adjustment falls back to Satterthwaite with a message", {
  expect_snapshot(predictmeans:::pm_kr_ddf_formula(A1 = 0.85, A2 = 0.1, q = 1), error = TRUE)
  expect_equal(predictmeans:::pm_kr_ddf_formula(c(0.1, NA), c(0.1, NA), 1)$ddf[2], NA_real_)
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Oats)
  local_mocked_bindings(pm_vcov_kr = function(...) stop("simulated failure"), .package = "predictmeans")
  expect_snapshot(dfm <- predictmeans:::pm_df_method(fm, "auto"))
  expect_identical(dfm$df_type, "Satterthwaite")
  expect_s3_class(dfm$vcov_adj, "pm_satt")
})

test_that("KR: observed information that is not positive definite falls back to the expected one, with a message", {
  expect_true(predictmeans:::pm_kr_info_check(diag(c(2, 1)), "observed"))
  expect_snapshot(ok <- predictmeans:::pm_kr_info_check(matrix(c(1, 2, 2, 1), 2), "observed"))
  expect_false(ok)
  expect_false(suppressMessages(predictmeans:::pm_kr_info_check(diag(c(1, 0)), "average")))   # singular
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  set.seed(7)
  Ou <- Oats[-sample(nrow(Oats), 9), ]        # unbalanced: observed and expected information differ
  fits <- list(
    lmer = list(m = lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Ou),
                kr = function(m, i) predictmeans:::pm_vcov_kr(m, i)),
    lme = list(m = nlme::lme(yield ~ nitro + Variety, random = ~ 1 | Block/Variety, data = Ou),
               kr = function(m, i) predictmeans:::pm_vcov_kr_nlme(m, i)))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expected <- f$kr(f$m, "expected")
    observed <- f$kr(f$m, "observed")
    expect_identical(attr(observed, "information"), "observed")
    expect_false(isTRUE(all.equal(unclass(observed)[, ], unclass(expected)[, ])))
    for (info in c("observed", "average")) {
      fb <- with_mocked_bindings(f$kr(f$m, info), pm_kr_info_check = function(I, information) FALSE,
                                 .package = "predictmeans")
      expect_identical(attr(fb, "information"), "expected", label = paste(nm, info))
      expect_equal(unclass(fb)[, ], unclass(expected)[, ], ignore_attr = TRUE, tolerance = 1e-10)
      expect_equal(attr(fb, "W"), attr(expected, "W"), ignore_attr = TRUE, tolerance = 1e-10)
    }
  }
})
