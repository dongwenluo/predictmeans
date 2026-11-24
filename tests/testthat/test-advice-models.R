## Advice rules of permmodels()/bootmodels(), R2_model() and varcomp() (R/advice.R)

ids_of <- function(x) pm_advice_get(x)$id

oats <- function() {
  d <- Oats
  d$nitro <- factor(d$nitro)
  d
}

test_that("permmodels(): few permutations are advised, enough are not", {
  skip_on_cran()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = oats())
  p49 <- permmodels(fm, nperm = 49, seed = 1)
  expect_identical(ids_of(p49), "few_permutations")
  expect_match(pm_advice_get(p49)$text, "at the floor")
  expect_true(any(grepl("^Advice:", capture.output(print(p49)))))
  expect_length(ids_of(permmodels(fm, nperm = 999, seed = 1)), 0)
})

test_that("pm_perm_advice(): p-values near 0.05 with few permutations", {
  out <- list(ANOVA = data.frame(Df = 1, Perm_p_value = 0.06))
  expect_identical(pm_perm_advice(out, 99)$id, "few_permutations")
  expect_null(pm_perm_advice(list(ANOVA = data.frame(Df = 1, Perm_p_value = 0.5)), 99))
  expect_null(pm_perm_advice(out, 4999))
})

test_that("R2_model(): latent scale for non-gaussian fits only", {
  gm <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = lme4::cbpp)
  expect_identical(ids_of(R2_model(gm)), "r2_latent_scale")
  expect_identical(ids_of(R2_model(gm, method = "nakagawa")), "r2_latent_scale")
  fm <- lme4::lmer(yield ~ nitro + (1 | Block), data = oats())
  expect_length(ids_of(R2_model(fm)), 0)
  expect_true(any(grepl("^Advice:", capture.output(print(R2_model(gm))))))
})

test_that("varcomp(): variance at zero and few levels", {
  set.seed(5071)
  d <- data.frame(g = factor(rep(1:4, each = 10)), x = rnorm(40))
  d$y <- d$x + rnorm(40)
  fm <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d))
  skip_if_not(lme4::isSingular(fm))
  v <- varcomp(fm)
  expect_setequal(ids_of(v), c("variance_at_zero", "few_re_levels", "variance_test_boundary"))
  expect_true(any(grepl("^Advice:", capture.output(print(v)))))
  # the printed table has no advice attribute line
  expect_false(any(grepl("attr", capture.output(print(v)))))
  ok <- lme4::lmer(yield ~ nitro + (1 | Block), data = oats())
  expect_length(ids_of(varcomp(ok)), 0)
})
