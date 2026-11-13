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
  pm <- predictmeans(fm, "nitro", vcov_adjust = TRUE, plot = FALSE, prtnum = FALSE)
  expect_equal(pm$settings$df_method, "Kenward-Roger")
  expect_true(all(pm$means$df > 0))
  # balanced: the observed and the expected information give the same df
  pe <- predictmeans(fm, "nitro", df_method = "KR-expected", plot = FALSE, prtnum = FALSE)
  expect_equal(pe$pairs$df, pm$pairs$df, tolerance = 1e-4)
  expect_equal(pe$settings$df_method, "KR-expected")
})
