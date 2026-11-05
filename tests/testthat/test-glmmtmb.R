# glmmTMB >= 1.1.12 stopped storing Lambdat, Lind and lower in its
# random-effect terms and renamed the dispersion parameter (betad ->
# betadisp). These tests cover the code that handles both versions.

test_that("R2_glmm() on a gaussian glmmTMB fit equals the lmer result", {
  skip_if_not_installed("glmmTMB")
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  f <- yield ~ nitro + Variety + (1 | Block/Variety)
  r_lmer <- R2_glmm(lme4::lmer(f, data = Oats))
  # gaussian glmmTMB has a dispersion parameter (betad / betadisp)
  r_tmb <- R2_glmm(glmmTMB::glmmTMB(f, data = Oats, REML = TRUE))
  pct <- function(x) as.numeric(sub("%", "", x))
  expect_equal(pct(r_tmb[1:3]), pct(r_lmer[1:3]), tolerance = 0.01)
})

test_that("random-effect terms rebuilt with mkReTrms match what glmmTMB stored", {
  skip_if_not_installed("glmmTMB")
  data(Clinical)
  fit <- glmmTMB::glmmTMB(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic/Treatment),
                          family = binomial, data = Clinical, doFit = FALSE)
  stored <- fit$condList$reTrms
  skip_if(is.null(stored$Lambdat), "this glmmTMB does not store Lambdat (nothing to compare)")
  mk <- if (requireNamespace("reformulas", quietly = TRUE)) reformulas::mkReTrms else lme4::mkReTrms
  rebuilt <- mk(predictmeans:::pm_findbars(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic/Treatment)),
                fit$fr, reorder.terms = FALSE)
  expect_equal(as.matrix(rebuilt$Lambdat), as.matrix(stored$Lambdat))
  expect_identical(rebuilt$Lind, stored$Lind)
  expect_equal(rebuilt$lower, stored$lower)
})

test_that("semireg_tmb() with predict_info and semipred() run, with R2_glmm()", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("HRW")
  data(WarsawApts, package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  fit <- suppressWarnings(semireg_tmb(areaPerMzloty ~ construction.date * district,
           smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
           data = W, prt = FALSE))
  sp <- semipred(fit, "district", "construction.date", prt = FALSE)$pred_df
  expect_true(all(is.finite(sp$SE)) && all(sp$SE > 0))
  expect_equal(nlevels(factor(sp$district)), nlevels(W$district))
  r2 <- R2_glmm(fit$semer)
  expect_true(all(is.finite(as.numeric(sub("%", "", r2)))))
})
