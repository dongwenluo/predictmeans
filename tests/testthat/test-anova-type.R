test_that("pm_anova_type() takes I/II/III and 1/2/3 in any form", {
  for (t in list("I", "i", 1, "1")) expect_identical(pm_anova_type(t), "I")
  for (t in list("II", 2L, "2")) expect_identical(pm_anova_type(t), "II")
  for (t in list("III", "iii", 3, " 3 ")) expect_identical(pm_anova_type(t), "III")
  expect_error(pm_anova_type(4), "must be one of")
  expect_error(pm_anova_type(NA), "must be one of")
  expect_error(pm_anova_type("I", c("II", "III")), "not available")
})

test_that("anova.lmerMod() gives the same table for every spelling of the type", {
  dat <- Oats
  dat$nitro <- factor(dat$nitro)
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = dat)
  expect_equal(unclass(anova(fm, type = 3)), unclass(anova(fm, type = "III")))
  expect_equal(unclass(anova(fm, type = 2)), unclass(anova(fm, type = "II")))
  expect_equal(unclass(anova(fm, type = 1)), unclass(anova(fm, type = "I")))
  expect_equal(unclass(pm_anova(fm, 3, "F")), unclass(pm_anova(fm, "III", "F")))
})

test_that("type advice appears only for unbalanced or empty-cell data", {
  dat <- Oats
  dat$nitro <- factor(dat$nitro)
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block), data = dat)
  expect_identical(as.vector(pm_design_balance(fm)), "balanced")
  expect_length(pm_anova_advice(fm, "III"), 0)
  ub <- lme4::lmer(yield ~ nitro * Variety + (1 | Block), data = dat[-c(1, 5), ])
  expect_identical(as.vector(pm_design_balance(ub)), "unbalanced")
  expect_match(pm_anova_advice(ub, "I"), "order")
  expect_match(paste(attr(anova(ub, type = 2), "heading"), collapse = ""), "Type II tests")
  expect_false(any(grepl("Advice:", attr(anova(ub, advice = FALSE), "heading"))))
  em <- lme4::lmer(yield ~ nitro * Variety + (1 | Block),
                   data = dat[!(dat$nitro == "0" & dat$Variety == "Victory"), ])
  expect_identical(as.vector(pm_design_balance(em)), "empty")
  expect_match(paste(pm_anova_advice(em, "III"), collapse = " "), "Empty cells")
})

test_that("anova advice is stored for code and AI, and the option hides only the printed note", {
  dat <- Oats
  dat$nitro <- factor(dat$nitro)
  ub <- lme4::lmer(yield ~ nitro * Variety + (1 | Block), data = dat[-c(1, 5), ])
  a <- pm_advice_get(anova(ub, type = 1))
  expect_identical(a$id, "anova_unbalanced_typeI")
  expect_identical(a$level, "caution")
  expect_identical(a$text, pm_anova_advice(ub, "I"))
  op <- options(predictmeans.advice = FALSE)
  tab <- tryCatch(anova(ub, type = 1), finally = options(op))
  expect_false(any(grepl("Advice:", attr(tab, "heading"))))
  expect_identical(pm_advice_get(tab)$id, "anova_unbalanced_typeI")
  expect_equal(nrow(pm_advice_get(anova(ub, advice = FALSE))), 0)
  em <- lme4::lmer(yield ~ nitro * Variety + (1 | Block),
                   data = dat[!(dat$nitro == "0" & dat$Variety == "Victory"), ])
  expect_identical(pm_advice_get(anova(em))$id, c("anova_empty_cells", "anova_unbalanced_typeIII"))
})

test_that("Type I Wald tests of a glmer: the first term equals its Type II test when alone", {
  gm <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = lme4::cbpp)
  expect_equal(pm_anova_type1(gm)[1, "Chisq"], pm_anova(gm, "II")[1, "Chisq"])
})
