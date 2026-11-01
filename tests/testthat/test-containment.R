# Containment df (SAS DDFM=CONTAIN). The integer df follow the SAS rule
# (rank contributions to [X Z]); the SEs are the unadjusted REML ones.

oats_sas <- function() {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  names(Oats)[names(Oats) == "nitro"] <- "Nitro"
  Oats
}

test_that("balanced split plot: Variety 10, Nitro and Variety:Nitro 45", {
  Oats <- oats_sas()
  m <- suppressMessages(lme4::lmer(yield ~ Variety*Nitro + (1|Block) + (1|Block:Variety), data = Oats))
  m2 <- suppressMessages(lme4::lmer(yield ~ Variety*Nitro + (1|Block/Variety), data = Oats))
  cdf <- predictmeans:::pm_containment_df(m)
  expect_equal(as.vector(cdf), c(5, 10, 45, 45))
  expect_equal(cdf, predictmeans:::pm_containment_df(m2))

  pv <- predictmeans(m, "Variety", df_method = "containment", pairwise = TRUE, plot = FALSE, prtnum = FALSE)
  expect_true(all(pv$means$df == 10) && all(pv$pairs$df == 10))
  expect_equal(pv$pairs$std.error[1], 7.078902, tolerance = 1e-6)   # unadjusted
  expect_false(pv$settings$vcov_adjust)
  pc <- predictmeans(m, "Variety:Nitro", df_method = "containment", plot = FALSE, prtnum = FALSE)
  expect_true(all(pc$means$df == 45))
  expect_equal(pc$means$std.error[1], 9.10701, tolerance = 1e-6)
  expect_equal(pc$means$conf.low[1], 80 - qt(0.975, 45) * 9.10701, tolerance = 1e-6)
})

test_that("missing whole plot and subplot: Variety 9, Nitro 41", {
  Oats <- oats_sas()
  Om <- Oats[!(Oats$Block == "I" & Oats$Variety == "Golden Rain") &
             !(Oats$Block == "II" & Oats$Variety == "Victory" & Oats$Nitro == "0.6"), ]
  m <- suppressMessages(lme4::lmer(yield ~ Variety*Nitro + (1|Block/Variety), data = Om))
  expect_equal(as.vector(predictmeans:::pm_containment_df(m)), c(5, 9, 41, 41))
  cm <- contrastmeans(m, "Variety", ctrmatrix = rbind(c(1, 0, -1)), df_method = "containment")
  expect_equal(unname(cm[[1]][, "df"]), 9)
  expect_equal(unname(df_term(m, "Nitro", type = "containment")), rep(41, 4))
})

test_that("random slopes: (Days | Subject) gives Subject and Days:Subject", {
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  cdf <- predictmeans:::pm_containment_df(m)
  expect_equal(as.vector(cdf), c(17, 17))
  expect_equal(attr(cdf, "residual"), 180 - 36)
  # a fixed term no random effect contains gets the residual df
  m2 <- lme4::lmer(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  expect_equal(unname(predictmeans:::pm_containment_df(m2)["Days"]), 180 - 19)
})

test_that("df_method is ignored with a warning for other model classes", {
  fm <- lm(yield ~ Variety*Nitro, data = oats_sas())
  expect_warning(predictmeans(fm, "Variety", df_method = "containment", plot = FALSE, prtnum = FALSE),
                 "lmer models only")
})
