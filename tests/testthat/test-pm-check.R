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
                              "Equal variances: Variety", "Equal variances: Nitro", "Spread vs fitted values",
                              "Balance of cells", "Interaction: Variety:Nitro",
                              "Straight line in Nitro"))
  # the interaction: SAS Type 3 F 0.41 on 6 and 42 df, p = 0.8656
  expect_equal(ck$p.value[ck$check == "Interaction: Variety:Nitro"], 0.8656, tolerance = 1e-3)
  # curvature: quadratic + cubic contrasts with the actual spacing of the rates
  P <- stats::poly(c(0, 0.2, 0.4, 0.6), 3)
  C <- t(P[, 2:3]); colnames(C) <- c("0.0", "0.2", "0.4", "0.6")
  expect_equal(ck$p.value[ck$check == "Straight line in Nitro"], pm_joint(fm, "Nitro", ctrmatrix = C)$p.value)
  expect_output(print(ck), "Checks inform, they do not decide")
  expect_output(print(ck), "Nothing flagged")
  expect_true(all(ck$section %in% pm_check_sections))
  expect_identical(ck$result[ck$check == "Balance of cells"], "unbalanced (5 to 6 per cell)")
  expect_identical(nrow(summary(ck)), 0L)
})

test_that("new checks: zeros, underdispersion, 0/1 dispersion, separation, GVIF, few RE levels", {
  set.seed(5307)
  z <- data.frame(x = rnorm(200))
  z$y <- ifelse(runif(200) < .3, 0, rpois(200, exp(1 + .3 * z$x)))
  ck <- pm_check(glm(y ~ x, family = poisson, data = z))
  expect_true(ck$look[ck$check == "Zero counts"])
  expect_true(all(summary(ck)$look))
  # underdispersed counts: binomial(10, 0.5) fitted as Poisson
  u <- data.frame(x = rnorm(200)); u$y <- rbinom(200, 10, 0.5)
  od <- pm_check_dispersion(glm(y ~ x, family = poisson, data = u))
  expect_lt(od$ratio, 0.7)
  ck <- pm_check(glm(y ~ x, family = poisson, data = u))
  expect_true(ck$look[ck$check == "Dispersion (Pearson chi-square / df)"])
  s <- data.frame(g = factor(rep(c("a", "b"), each = 20)), x = rnorm(40)); s$y <- as.numeric(s$g == "b")
  ck <- pm_check(suppressWarnings(glm(y ~ g + x, family = binomial, data = s)))
  expect_true(ck$look[ck$check == "Separation (binomial)"])
  expect_identical(ck$status[ck$check == "Dispersion"], "i")
  # GVIF: equals 1/(1 - R^2) for two covariates
  cv <- data.frame(a = rnorm(50)); cv$b <- cv$a + rnorm(50, sd = 0.1); cv$y <- rnorm(50)
  m <- lm(y ~ a + b, cv)
  expect_equal(unname(pm_check_gvif(m)[1]^2), 1 / (1 - cor(cv$a, cv$b)^2))
  expect_true(pm_check(m)$look[pm_check(m)$check == "Collinearity (GVIF^(1/(2 df)))"])
  o <- as.data.frame(nlme::Oats); o$nitro <- factor(o$nitro)
  fm <- lme4::lmer(yield ~ nitro + (1 | Block), data = o)
  expect_false("Random-effect levels: Block" %in% pm_check(fm)$check)   # 6 levels
  fm3 <- lme4::lmer(yield ~ nitro + (1 | Variety), data = o)
  expect_identical(pm_check(fm3)$status[pm_check(fm3)$check == "Random-effect levels: Variety"], "i")
})

test_that("problems are flagged: overdispersion, outliers, unequal variances, aliasing", {
  ck <- pm_check(glm(breaks ~ wool * tension, family = poisson, data = warpbreaks))
  expect_true(ck$look[ck$check == "Dispersion (Pearson chi-square / df)"])
  expect_false("Overdispersion (Pearson chi-square / df)" %in% ck$check)
  set.seed(1)
  b <- data.frame(x = rnorm(100)); b$y <- rbinom(100, 1, plogis(b$x))
  expect_identical(pm_check(glm(y ~ x, family = binomial, data = b))$result[
    pm_check(glm(y ~ x, family = binomial, data = b))$check == "Dispersion"], "0/1 response: not identifiable")
  expect_true(ck$look[ck$check == "Outliers (|residual| > 3)"])
  set.seed(2)
  h <- data.frame(trt = factor(rep(c("a", "b", "c"), each = 20)))
  h$y <- rnorm(60, c(10, 12, 14)[h$trt], c(1, 1, 4)[h$trt])
  h$y[5] <- 25
  ck <- pm_check(lm(y ~ trt, h))
  expect_true(ck$look[ck$check == "Equal variances: trt"])
  expect_true(ck$look[ck$check == "Outliers (|residual| > 3)"])
  expect_match(ck$result[ck$check == "Outliers (|residual| > 3)"], "^1 of 60.*: 5 ")
  # rows worth a look with a note are stored as advice; summary() keeps it
  a <- pm_advice_get(ck)
  expect_true(all(c("check_equal_variances_trt", "check_outliers_residual_3") %in% a$id))
  expect_true(all(a$level == "caution"))
  expect_identical(pm_advice_get(summary(ck)), a)
  expect_equal(nrow(a), sum(ck$look & nzchar(ck$note)))
  h$trt2 <- h$trt
  ck <- suppressWarnings(pm_check(lm(y ~ trt + trt2, h)))
  expect_identical(ck$result[ck$check == "Estimability"], "2 coefficient(s) aliased")
  # a binomial model with trials: overdispersion checked; a 0/1 response: not
  data(Drug, envir = environment())
  expect_true("Dispersion (Pearson chi-square / df)" %in%
                pm_check(glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug))$check)
  set.seed(1)
  b <- data.frame(x = rnorm(100)); b$y <- rbinom(100, 1, plogis(b$x))
  expect_false("Dispersion (Pearson chi-square / df)" %in% pm_check(glm(y ~ x, family = binomial, data = b))$check)
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
