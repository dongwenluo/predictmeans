## permmodels(method = "freedman-lane") for mixed models, GLMs and GLMMs

test_that("Freedman-Lane for lmer/lme/gls uses the model's Wald F and t", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  Om <- Oats[-c(1, 5, 9, 20, 33), ]
  m1 <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Om)
  r1 <- permmodels(m1, nperm = 49, type = "II", seed = 1, prt = FALSE, method = "freedman-lane")
  # car's type II hypotheses (as pm_anova()), unadjusted covariance: Wald chi-square / df
  a1 <- pm_anova(m1, type = "II", test.statistic = "Chisq")
  expect_equal(r1$ANOVA[, "F value"], round(unname(a1[, "Chisq"] / a1[, "Df"]), 4))
  ct <- coef(summary(lmerTest::as_lmerModLmerTest(m1)))
  expect_equal(unname(r1$COEFFICENTS[, "t value"]), round(unname(ct[, "t value"]), 4))
  # the same model by lme: the same F (both REML)
  m2 <- nlme::lme(yield ~ nitro * Variety, random = ~ 1 | Block/Variety, data = Om)
  r2 <- permmodels(m2, nperm = 49, type = "II", seed = 1, prt = FALSE, method = "freedman-lane")
  expect_equal(r2$ANOVA[, "F value"], r1$ANOVA[, "F value"], tolerance = 1e-3)
  # a pair that is one coefficient: the same test as that coefficient's
  pm <- predictmeans(m1, "nitro", permlist = r1, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_equal(pm$pairs$p.value[1], unname(r1$COEFFICENTS["nitro0.2", "Perm_p_value"]), tolerance = 1e-4)
  # lme with shuffled rows and the matching permutations: the same p-values
  o <- sample(nrow(Om))
  m3 <- nlme::lme(yield ~ nitro * Variety, random = ~ 1 | Block/Variety, data = Om[o, ])
  i1 <- permindex(Om, block = "Block", nsim = 29, seed = 2)
  i3 <- i1[o, ]; i3[] <- match(i1[o, ], o)
  expect_identical(permmodels(m2, perms = i1, prt = FALSE, method = "freedman-lane")$ANOVA$Perm_p_value,
                   permmodels(m3, perms = i3, prt = FALSE, method = "freedman-lane")$ANOVA$Perm_p_value)
  # gls
  g <- nlme::gls(yield ~ nitro + Variety, correlation = nlme::corCompSymm(form = ~ 1 | Block), data = Om)
  rg <- permmodels(g, nperm = 29, type = "III", seed = 1, prt = FALSE, method = "freedman-lane")
  expect_equal(rg$ANOVA[, "F value"], round(unname(pm_anova(g, "III")[-1, "Chisq"] / pm_anova(g, "III")[-1, "Df"]), 4),
               tolerance = 1e-3)
})

test_that("Freedman-Lane on the working scale for glm, glmer and glmmTMB", {
  # a gaussian glm is the linear model: the same p-values as lm
  d <- data.frame(A = factor(rep(1:3, each = 8)), x = exp(rnorm(24)))
  d$y <- d$x + rnorm(24)
  a <- permmodels(lm(y ~ x + A, data = d), nperm = 49, type = "II", seed = 3, prt = FALSE, method = "freedman-lane")
  b <- permmodels(glm(y ~ x + A, data = d), nperm = 49, type = "II", seed = 3, prt = FALSE, method = "freedman-lane")
  expect_equal(as.numeric(a$ANOVA[c("x", "A"), "Perm_p_value"]), b$ANOVA[c("x", "A"), "Perm_p_value"])
  expect_equal(unname(a$COEFFICENTS[, "Perm_p_value"]), unname(b$COEFFICENTS[, "Perm_p_value"]))
  skip_on_cran()
  data(Clinical, package = "predictmeans", envir = environment())
  gm <- lme4::glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic), family = binomial, data = Clinical)
  tm <- glmmTMB::glmmTMB(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic), family = binomial,
                         data = Clinical)
  rg <- permmodels(gm, nperm = 49, seed = 1, prt = FALSE, method = "freedman-lane")
  rt <- permmodels(tm, nperm = 49, seed = 1, prt = FALSE, method = "freedman-lane")
  expect_equal(rg$ANOVA$Perm_p_value, rt$ANOVA$Perm_p_value)
  expect_match(attr(rg$ANOVA, "heading"), "working")
  expect_error(permmodels(gm, nperm = 9, prt = FALSE), "freedman-lane")
  pm <- predictmeans(gm, "Treatment", permlist = rg, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_equal(pm$pairs$p.value, unname(rg$COEFFICENTS["Treatmentdrug", "Perm_p_value"]), tolerance = 1e-4)
})

test_that("raw permutations of lmer: type I compares lmerTest's type I F", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  Om <- Oats[-c(1, 5, 9, 20, 33, 41, 50), ]
  m <- lme4::lmer(yield ~ Variety + nitro + Variety:nitro + (1 | Block), data = Om)
  r <- permmodels(m, nperm = 9, prt = FALSE, ncore = 1, seed = 1)
  a <- anova(lmerTest::as_lmerModLmerTest(m), type = "I")
  expect_equal(r$ANOVA[, "F value"], round(unname(a[, "F value"]), 4))
  X <- lme4::getME(m, "X")
  Fs <- pm_type1_F(pm_type1_L(X, attr(terms(m), "term.labels")), lme4::fixef(m), vcov(m))
  expect_equal(unname(Fs[, 1]), unname(a[, "F value"]), tolerance = 1e-8)
})
