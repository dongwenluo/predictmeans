test_that("bootmodels(): bootstrap p-values, repeatable, usable as permlist", {
  skip_on_cran()
  data(Oats, package = "nlme", envir = environment())
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  set.seed(4); u1 <- runif(1); set.seed(4)
  b1 <- bootmodels(fm, nsim = 19, seed = 1, prt = FALSE)
  expect_identical(runif(1), u1)
  b2 <- bootmodels(fm, nsim = 19, seed = 1, prt = FALSE)
  expect_identical(b1$ANOVA, b2$ANOVA)
  expect_identical(attr(b1, "resampling"), "parametric bootstrap")
  expect_length(b1$permlist, 19)
  # nitro has a very large effect: no simulated F reaches the observed one
  expect_equal(b1$ANOVA["nitro", "Boot_p_value"], round(1/20, 4))
  expect_true(all(b1$ANOVA$Boot_p_value >= 1/20 - 1e-8 & b1$ANOVA$Boot_p_value <= 1))
  r <- predictmeans(fm, "nitro", permlist = b1, plot = FALSE, ncore = 1)
  expect_true(all(r$pairs$p.method == "bootstrap"))
  b3 <- bootmodels(fm, nsim = 9, type = "III", seed = 1, prt = FALSE)
  expect_true(is.na(b3$ANOVA["(Intercept)", "Boot_p_value"]))
})

test_that("bootmodels() works for glmer and refuses other models", {
  skip_on_cran()
  g <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  b <- bootmodels(g, nsim = 9, seed = 2, prt = FALSE)
  expect_equal(nrow(b$COEFFICENTS), 4)
})

test_that("bootmodels() for lm, aov, glm, gls, lme and glmmTMB; lm p-values approach the exact F tests", {
  skip_on_cran()
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  # lm: F is pivotal, so the bootstrap p-values approach those of the F tests
  fl <- lm(yield ~ Block + nitro*Variety, O)
  b <- bootmodels(fl, nsim = 999, seed = 1, prt = FALSE)
  pF <- anova(fl)[["Pr(>F)"]]
  ok <- !is.na(pF) & pF > 0.01
  expect_equal(b$ANOVA$Boot_p_value[ok], pF[ok], tolerance = 0.05)
  expect_true(is.na(b$ANOVA["Residuals", "Boot_p_value"]))
  expect_equal(bootmodels(aov(yield ~ Block + nitro*Variety, O), nsim = 19, seed = 1, prt = FALSE)$ANOVA,
               bootmodels(fl, nsim = 19, seed = 1, prt = FALSE)$ANOVA)
  # glm: Poisson (deviance), binomial with counts, negative binomial; usable as permlist
  gp <- glm(breaks ~ wool*tension, poisson, warpbreaks)
  bg <- bootmodels(gp, nsim = 19, seed = 1, prt = FALSE)
  expect_equal(rownames(bg$ANOVA), c("wool", "tension", "wool:tension"))
  r <- predictmeans(gp, "tension", permlist = bg, plot = FALSE)
  expect_true(all(r$pairs$p.method == "bootstrap"))
  expect_length(bootmodels(glm(cbind(incidence, size - incidence) ~ period, binomial, lme4::cbpp), nsim = 9, seed = 1,
                           prt = FALSE)$permlist, 9)
  expect_length(bootmodels(MASS::glm.nb(breaks ~ wool*tension, warpbreaks), nsim = 9, seed = 1, prt = FALSE)$permlist, 9)
  # lme and gls: refits to responses simulated from the null model's covariance
  bl <- bootmodels(nlme::lme(yield ~ nitro*Variety, random = ~1 | Block/Variety, O), nsim = 19, seed = 1, prt = FALSE)
  expect_equal(bl$ANOVA["nitro", "Boot_p_value"], round(1/20, 4))
  bs <- bootmodels(nlme::gls(distance ~ age*Sex, nlme::Orthodont, correlation = nlme::corCompSymm(form = ~1 | Subject)),
                   nsim = 9, seed = 1, prt = FALSE)
  expect_true(is.na(bs$ANOVA["(Intercept)", "Boot_p_value"]))
  # glmmTMB: no Type I tests, Type II used
  skip_if_not_installed("glmmTMB")
  tm <- glmmTMB::glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = lme4::cbpp)
  expect_message(bt <- bootmodels(tm, nsim = 9, seed = 1, prt = FALSE), "Type II")
  expect_equal(rownames(bt$ANOVA), "period")
  # refused: quasi families and aov with Error()
  expect_error(bootmodels(glm(breaks ~ wool, quasipoisson, warpbreaks), nsim = 5), "quasi family")
  expect_error(bootmodels(aov(yield ~ nitro + Error(Block), O), nsim = 5), "Error\\(\\)")
})

test_that("bootmodels(): lme and gls responses are simulated with the null model's covariance of the data", {
  g <- nlme::gls(distance ~ age*Sex, nlme::Orthodont, correlation = nlme::corCompSymm(form = ~1 | Subject))
  bp <- predictmeans:::pm_boot_parts(g)
  Y <- do.call(cbind, predictmeans:::pm_with_seed(1, bp$simulate(4000)))
  rho <- coef(bp$null_fit$modelStruct$corStruct, unconstrained = FALSE)[[1]]
  s2 <- bp$null_fit$sigma^2
  expect_equal(mean(apply(Y, 1, var)), s2, tolerance = 0.03)
  # two measurements of one subject, and of two subjects
  expect_equal(cov(Y[1, ], Y[2, ]), rho*s2, tolerance = 0.06)
  expect_equal(cov(Y[1, ], Y[5, ]), 0, tolerance = 0.06*s2)
  expect_equal(rowMeans(Y), as.numeric(fitted(bp$null_fit)), tolerance = 0.02)
})
