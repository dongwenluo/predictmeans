# Type II/III tests used by permmodels(), compared with car::Anova().

expect_anova_equal <- function(new, ref, label) {
  expect_identical(names(new), names(ref), label = label)
  expect_identical(rownames(new), rownames(ref), label = label)
  expect_equal(as.matrix(new), as.matrix(ref), tolerance = 1e-10, ignore_attr = TRUE, label = label)
}

test_that("pm_anova() equals car::Anova() for lm, glm, lme, gls and lmer", {
  skip_if_not_installed("car")
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  set.seed(7)
  Ou <- Oats[-sample(nrow(Oats), 9), ]
  data("Drug", package = "predictmeans", envir = environment())
  # car's model.matrix.gls() only finds data in the global environment
  assign("Ou_anova_test", Ou, envir = globalenv())
  on.exit(rm("Ou_anova_test", envir = globalenv()), add = TRUE)
  models <- list(
    lm = lm(yield ~ Block + nitro*Variety, Ou),
    glm = glm(cbind(R, N - R) ~ Drug + log10Dose, binomial, Drug),
    lme = lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = Ou),
    gls = gls(yield ~ nitro*Variety, correlation = corCompSymm(form = ~1|Block), data = Ou_anova_test),
    lmer = lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou))
  for (nm in names(models)) {
    m <- models[[nm]]
    tests <- if (nm == "glm") c("LR", "Wald", "F") else if (nm == "lmer") c("Chisq", "F") else "Chisq"
    for (ty in c("II", "III")) for (ts in tests) {
      ref <- suppressMessages(car::Anova(m, type = ty, test.statistic = ts))
      expect_anova_equal(predictmeans:::pm_anova(m, type = ty, test.statistic = ts, kr_information = "expected"), ref,
                         paste(nm, ty, ts))
    }
  }
})

test_that("permmodels() type II works with aliased coefficients", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  Oats$n2 <- Oats$nitro
  m <- lm(yield ~ nitro + n2 + Variety, Oats)
  p <- suppressMessages(permmodels(m, method = "raw", nperm = 9, type = "II", prt = FALSE, ncore = 1, seed = 1))
  expect_true("Variety" %in% rownames(p$ANOVA))
  skip_if_not_installed("car")
  ref <- suppressMessages(car::Anova(m, type = "II"))
  expect_anova_equal(predictmeans:::pm_anova(m, "II"), ref, "aliased lm")
})

test_that("Kenward-Roger F tables use the observed information, as pm_joint()", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  Oats$Variety <- factor(Oats$Variety, ordered = FALSE)
  set.seed(7)
  Ou <- Oats[-sample(nrow(Oats), 9), ]
  m <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Ou)
  a <- predictmeans:::pm_anova(m, type = "III", test.statistic = "F")
  jt <- pm_joint(m, "nitro:Variety", ctrmatrix = t(kronecker(contr.sum(4), contr.sum(3))))
  expect_equal(a["nitro:Variety", "Df.res"], jt$df2, tolerance = 1e-8)
  expect_equal(a["nitro:Variety", "F"], jt$F.ratio, tolerance = 1e-8)
  ae <- predictmeans:::pm_anova(m, type = "III", test.statistic = "F", kr_information = "expected")
  expect_false(isTRUE(all.equal(a[, "Df.res"], ae[, "Df.res"])))     # unbalanced: they differ
})
