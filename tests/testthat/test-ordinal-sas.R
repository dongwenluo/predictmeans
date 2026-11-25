# Ordinal models against SAS (dev/sas/ordinal.sas, run by the maintainer 2026-11-26;
# results: dev/sas/ordinal-results-2026-11-26.pdf). SAS prints estimates and SEs to 4
# decimals: tolerance 1e-4; the class probabilities to 14-16 digits: 1e-6 (convergence).
# SAS writes logit P(Y <= j) = alpha_j + x'b, so its b is minus clm's.

test_that("cumulative logit: clm = SAS Newton (observed information), vglm = SAS Fisher (expected)", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  skip_if_not_installed("VGAM")
  data(wine, package = "ordinal", envir = environment())
  f <- rating ~ temp * contact
  sas <- c(-1.4113, 1.1436, 3.3771, 4.9420, -2.3212, -1.3475, -0.3595)
  se_newton <- c(0.5454, 0.5097, 0.6382, 0.7509, 0.7009, 0.6604, 0.9238)
  se_fisher <- c(0.5386, 0.5108, 0.6406, 0.7525, 0.6992, 0.6629, 0.9208)
  cm <- ordinal::clm(f, data = wine)
  expect_lt(max(abs(c(cm$alpha, -cm$beta) - sas)), 1e-4)
  expect_lt(max(abs(sqrt(diag(vcov(cm))) - se_newton)), 1e-4)
  vm <- VGAM::vglm(f, VGAM::cumulative(parallel = TRUE), data = wine, control = VGAM::vglm.control(epsilon = 1e-12))
  expect_lt(max(abs(VGAM::coef(vm) - sas)), 1e-4)
  expect_lt(max(abs(sqrt(diag(VGAM::vcov(vm))) - se_fisher)), 1e-4)       # why vglm and clm SEs differ
  # class probabilities of the four cells (SAS OUTPUT PREDPROBS)
  P <- rbind(c(0.19603508959566, 0.5622964047613, 0.2086490205024, 0.02592966527815, 0.00708981986249),
             c(0.05959593693189, 0.38960329628917, 0.43467286832494, 0.08938866095626, 0.02673923749775),
             c(0.0233747827134, 0.21210388562606, 0.50642705738585, 0.1902823812734, 0.0678118930013),
             c(0.00432308507495, 0.04859512494076, 0.28982099044683, 0.37103561926005, 0.28622518027741))
  pp <- predictmeans(cm, "temp:contact", mean_type = "prob")$probs
  expect_lt(max(abs(pp$estimate - as.vector(t(P)))), 1e-6)
  pv <- predictmeans(vm, "temp:contact", mean_type = "prob")$probs
  expect_lt(max(abs(pv$estimate - as.vector(t(P)))), 1e-6)
})

test_that("probit, cloglog and adjacent-category links = SAS", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  skip_if_not_installed("VGAM")
  data(wine, package = "ordinal", envir = environment())
  f <- rating ~ temp * contact
  cp <- ordinal::clm(f, data = wine, link = "probit")
  expect_lt(max(abs(c(cp$alpha, -cp$beta) - c(-0.8030, 0.6951, 2.0069, 2.9100, -1.4277, -0.7954, -0.1430))), 1e-4)
  cl <- ordinal::clm(f, data = wine, link = "cloglog")
  expect_lt(max(abs(c(cl$alpha, -cl$beta) - c(-1.7156, 0.3260, 1.7705, 2.6343, -1.6697, -0.9228, 0.1202))), 1e-4)
  # SAS LINK=ALOGIT: log P(Y = j) / P(Y = j + 1) = alpha_j + x'b, minus vglm's acat
  va <- VGAM::vglm(f, VGAM::acat(parallel = TRUE), data = wine)
  expect_lt(max(abs(-VGAM::coef(va) - c(-0.9855, 0.9069, 2.5324, 2.7924, -1.7106, -0.9885, 0.0401))), 1e-4)
})

test_that("clmm = SAS GLIMMIX (Laplace); df within judges 57 as SAS", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  mm <- ordinal::clmm(rating ~ temp * contact + (1 | judge), data = wine)
  expect_lt(max(abs(c(mm$alpha, -mm$beta) - c(-1.6628, 1.4439, 4.1687, 6.0440, -2.9467, -1.7174, -0.2315))), 1e-4)
  expect_lt(max(abs(sqrt(diag(vcov(mm)))[1:7] - c(0.7003, 0.6641, 0.8438, 0.9898, 0.7611, 0.7020, 0.9540))), 1e-4)
  expect_lt(abs(mm$ST[[1]][1, 1]^2 - 1.27603083846416), 1e-4)
  r <- predictmeans(mm, "temp:contact", pairwise = TRUE)
  expect_equal(unique(r$pairs$df), 57)                       # SAS Den DF of temp, contact, temp*contact
  # a treatment between 6 clusters: SAS's containment default gives 53 (near z); predictmeans keeps the
  # inner-outer 4, which held the level in the simulation of dev/review-classes.md (z: 11% at 5%)
})
