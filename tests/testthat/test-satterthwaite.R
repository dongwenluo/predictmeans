test_that("vectorised Satterthwaite df equal lmerTest", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  set.seed(7)
  Ou <- Oats[-sample(nrow(Oats), 9), ]
  data("sleepstudy", package = "lme4", envir = environment())
  models <- list(lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou),
                 lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou, REML = FALSE),
                 lme4::lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy))
  for (m in models) {
    L <- rbind(diag(length(lme4::fixef(m))), c(1, -1, rep(0, length(lme4::fixef(m)) - 2)))
    mt <- lmerTest::as_lmerModLmerTest(m)
    ref <- apply(L, 1, function(l) lmerTest::contest1D(mt, l)$df)
    new <- predictmeans:::pm_satt_ddf(L, predictmeans:::pm_satt_parts(m))
    expect_equal(unname(new), ref, tolerance = 1e-8)
    expect_equal(unname(df_term(m, ctrmatrix = L, type = "Satterthwaite")), ref, tolerance = 1e-8)
  }
})
