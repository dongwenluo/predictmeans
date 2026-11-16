## anova() of lmer fits (anova.lmerMod, R/anova.R): the package's own F tests
## replace the tables of lmerTest, which predictmeans attached before 2.0.

test_that("anova() of an lmer fit: Type III/I with Satterthwaite df equal lmerTest; KR default; comparisons by lme4", {
  skip_if_not_installed("lmerTest")
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fm <- lme4::lmer(yield ~ nitro*Variety + (1 | Block/Variety), data = O)
  lt <- lmerTest::as_lmerModLmerTest(fm)
  a3 <- anova(fm, ddf = "Satterthwaite")
  r3 <- anova(lt, type = "III")
  expect_equal(a3[rownames(r3), "F"], r3[["F value"]], tolerance = 1e-6)
  expect_equal(a3[rownames(r3), "Df.res"], r3[["DenDF"]], tolerance = 1e-4)
  a1 <- anova(fm, type = "I")
  r1 <- anova(lt, type = "I")
  expect_equal(a1[["F value"]], r1[["F value"]], tolerance = 1e-6)
  expect_equal(a1[["DenDF"]], r1[["DenDF"]], tolerance = 1e-4)
  # the default: Type III with Kenward-Roger df (balanced: the same df here)
  ak <- anova(fm)
  expect_match(attr(ak, "heading")[1], "Kenward-Roger")
  expect_equal(ak[rownames(r3), "Df.res"], r3[["DenDF"]], tolerance = 1e-4)
  # several models: lme4's likelihood ratio comparison, with the models' names
  f0 <- lme4::lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = O, REML = FALSE)
  f1 <- lme4::lmer(yield ~ nitro*Variety + (1 | Block/Variety), data = O, REML = FALSE)
  cmp <- anova(f0, f1)
  expect_equal(rownames(cmp), c("f0", "f1"))
  expect_true("Chisq" %in% names(cmp))
  # lme4's own table on request
  expect_true("Sum Sq" %in% names(anova(fm, ddf = "lme4")))
})

test_that("with lmerTest loaded too: its fits keep its anova(), comparisons keep names, ranova() refits semireg smooths", {
  skip_if_not_installed("lmerTest")
  skip_on_cran()
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  lt0 <- lmerTest::lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = O, REML = FALSE)
  lt1 <- lmerTest::lmer(yield ~ nitro*Variety + (1 | Block/Variety), data = O, REML = FALSE)
  expect_match(attr(anova(lt1), "heading")[1], "Satterthwaite")
  expect_equal(rownames(anova(lt0, lt1)), c("lt0", "lt1"))   # lmerTest -> anova.lmerMod -> lme4
  data(Dialyzer, package = "nlme", envir = environment())
  fit <- suppressWarnings(semireg(rate ~ s(pressure, k = 6) + QB + (1 | Subject), data = Dialyzer, prt = FALSE))
  rt <- suppressWarnings(suppressMessages(lmerTest::ranova(fit$semer)))
  a <- anova(fit)$smooth
  expect_equal(rt[["LRT"]][rownames(rt) == "(1 | s_pressure)"], a$LRT, tolerance = 1e-6)
  expect_error(update(fit$semer), "new random part")
})
