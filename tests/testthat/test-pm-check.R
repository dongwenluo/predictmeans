## pm_check(): a table of model checks with plain notes.

test_that("a sound split plot: nothing to look at; interaction and straight-line notes", {
  o <- as.data.frame(nlme::Oats)
  o$Variety <- factor(gsub(" ", "", o$Variety))
  o$Nitro <- factor(sprintf("%.1f", o$nitro))
  d <- droplevels(o[!(o$Block == "I" & o$Variety == "GoldenRain"), ])
  fm <- lme4::lmer(yield ~ Variety * Nitro + (1 | Block/Variety), data = d)
  ck <- pm_check(fm)
  expect_s3_class(ck, "pm_check")
  expect_false(any(ck$look))
  expect_setequal(ck$check, c("Fit", "Outliers (|residual| > 3)", "Normality of residuals",
                              "Normality of random effects: Variety:Block", "Normality of random effects: Block",
                              "Equal variances: Variety", "Equal variances: Nitro", "Interaction: Variety:Nitro",
                              "Straight line in Nitro"))
  # the interaction: SAS Type 3 F 0.41 on 6 and 42 df, p = 0.8656
  expect_equal(ck$p.value[ck$check == "Interaction: Variety:Nitro"], 0.8656, tolerance = 1e-3)
  # curvature: quadratic + cubic contrasts with the actual spacing of the rates
  P <- stats::poly(c(0, 0.2, 0.4, 0.6), 3)
  C <- t(P[, 2:3]); colnames(C) <- c("0.0", "0.2", "0.4", "0.6")
  expect_equal(ck$p.value[ck$check == "Straight line in Nitro"], pm_joint(fm, "Nitro", ctrmatrix = C)$p.value)
  expect_output(print(ck), "Checks inform, they do not decide")
})

test_that("problems are flagged: overdispersion, outliers, unequal variances, aliasing", {
  ck <- pm_check(glm(breaks ~ wool * tension, family = poisson, data = warpbreaks))
  expect_true(ck$look[ck$check == "Overdispersion (Pearson chi-square / df)"])
  expect_true(ck$look[ck$check == "Outliers (|residual| > 3)"])
  set.seed(2)
  h <- data.frame(trt = factor(rep(c("a", "b", "c"), each = 20)))
  h$y <- rnorm(60, c(10, 12, 14)[h$trt], c(1, 1, 4)[h$trt])
  h$y[5] <- 25
  ck <- pm_check(lm(y ~ trt, h))
  expect_true(ck$look[ck$check == "Equal variances: trt"])
  expect_true(ck$look[ck$check == "Outliers (|residual| > 3)"])
  expect_match(ck$result[ck$check == "Outliers (|residual| > 3)"], "^1 of 60.*: 5 ")
  h$trt2 <- h$trt
  ck <- suppressWarnings(pm_check(lm(y ~ trt + trt2, h)))
  expect_identical(ck$result[ck$check == "Estimability"], "2 coefficient(s) aliased")
  # a binomial model with trials: overdispersion checked; a 0/1 response: not
  data(Drug, envir = environment())
  expect_true("Overdispersion (Pearson chi-square / df)" %in%
                pm_check(glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug))$check)
  set.seed(1)
  b <- data.frame(x = rnorm(100)); b$y <- rbinom(100, 1, plogis(b$x))
  expect_false("Overdispersion (Pearson chi-square / df)" %in% pm_check(glm(y ~ x, family = binomial, data = b))$check)
})

test_that("other classes run: lme, gls, glmer, glmmTMB, aov with Error()", {
  o <- as.data.frame(nlme::Oats); o$nitro <- factor(o$nitro)
  for (m in list(nlme::lme(yield ~ nitro * Variety, random = ~ 1 | Block/Variety, data = o),
                 nlme::gls(yield ~ nitro * Variety, correlation = nlme::corCompSymm(form = ~ 1 | Block), data = o),
                 lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial,
                             data = lme4::cbpp),
                 glmmTMB::glmmTMB(breaks ~ wool + tension + (1 | tension:wool), family = poisson, data = warpbreaks),
                 aov(yield ~ nitro * Variety + Error(Block/Variety), o))) {
    ck <- suppressWarnings(pm_check(m))
    expect_s3_class(ck, "pm_check")
    expect_true(all(c("Fit", "Outliers (|residual| > 3)") %in% ck$check))
  }
})
