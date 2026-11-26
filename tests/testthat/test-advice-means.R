## Advice rules of predictmeans() and contrastmeans() (R/advice.R)

ids_of <- function(x) pm_advice_get(x)$id

test_that("a covariate aliased with a factor, and means with large SEs, get advice", {
  data("DataExam5.2", package = "predictmeans", envir = environment())
  m <- lm(ht ~ site + sitemean * seedlot, DataExam5.2)
  r <- suppressMessages(predictmeans(m, "seedlot"))
  expect_setequal(ids_of(r), c("aliased_covariate", "large_se"))
  a <- pm_advice_get(r)
  # the seedlots measured at 2 sites only (SE about 107 vs 30)
  expect_match(a$text[a$id == "large_se"], "13686")
  expect_match(a$text[a$id == "large_se"], "14623")
  expect_identical(a$level[1], "caution")
})

test_that("unadjusted pairs, back-transformed GLMM means and bias_adjust", {
  gm <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = lme4::cbpp)
  r <- suppressMessages(predictmeans(gm, "period", pairwise = TRUE))
  # unequal LSDs, but no pair judged differently; cbpp is overdispersed (fitted with (1 | obs) in lme4's examples):
  # the df count trials, so the advice says when to count rows
  expect_setequal(ids_of(r), c("bias_adjust", "many_tests_unadjusted", "overdispersion_trials"))
  # bias_adjust done: no bias advice; an adjustment: no many-tests advice
  r2 <- suppressMessages(predictmeans(gm, "period", bias_adjust = TRUE, adj = "holm"))
  expect_false(any(c("bias_adjust", "many_tests_unadjusted") %in% ids_of(r2)))
  # pairs not asked for: no many-tests advice
  expect_false("many_tests_unadjusted" %in% ids_of(suppressMessages(predictmeans(gm, "period"))))
  # printed by print() and summary(), hidden by the option
  expect_true(any(grepl("^Advice:", capture.output(print(summary(r))))))
  op <- options(predictmeans.advice = FALSE)
  out <- tryCatch(capture.output(print(summary(r))), finally = options(op))
  expect_false(any(grepl("^Advice:", out)))
})

test_that("a log response: median, bias_adjust advised; a GLM mean: no bias advice", {
  dat <- Oats
  dat$nitro <- factor(dat$nitro)
  fm <- lm(log(yield) ~ nitro + Variety, data = dat)
  expect_true("bias_adjust" %in% ids_of(suppressMessages(predictmeans(fm, "nitro"))))
  g <- glm(breaks ~ wool + tension, family = poisson, data = warpbreaks)
  expect_false("bias_adjust" %in% ids_of(suppressMessages(predictmeans(g, "tension"))))
  # balanced lm, few means: nothing to say
  expect_length(ids_of(suppressMessages(predictmeans(lm(yield ~ nitro + Variety, data = dat), "nitro"))), 0)
})

test_that("empty cells and non-estimable means are advised", {
  dat <- Oats
  dat$nitro <- factor(dat$nitro)
  em <- lm(yield ~ nitro * Variety, data = dat[!(dat$nitro == "0" & dat$Variety == "Victory"), ])
  expect_true("empty_cells" %in% ids_of(suppressWarnings(suppressMessages(predictmeans(em, "nitro")))))
})

test_that("contrastmeans(): many unadjusted contrasts are advised and printed", {
  dat <- Oats
  dat$nitro <- factor(dat$nitro)
  fm <- lm(yield ~ nitro + Variety, data = dat)
  C <- rbind(diag(4)[-1, ] - c(1, 0, 0, 0), c(1, 1, -1, -1) / 2, c(1, -1, 1, -1) / 2, c(1, -1, -1, 1) / 2)
  cm <- contrastmeans(fm, "nitro", ctrmatrix = C)
  expect_identical(ids_of(cm), "many_tests_unadjusted")
  expect_true(any(grepl("^Advice:", capture.output(print(cm)))))
  expect_length(ids_of(contrastmeans(fm, "nitro", ctrmatrix = C, adj = "holm")), 0)
})
