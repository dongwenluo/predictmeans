## Singular fits are noted and recorded; adjusted p-values use each pair's df (SAS ADJDFE=ROW)

test_that("a singular lmer fit gives a message and settings$singular", {
  set.seed(4)
  d <- expand.grid(A = factor(1:3), rep = 1:4, block = factor(1:5))
  d$y <- as.integer(d$A) + rnorm(nrow(d))              # no block effect: the block variance is 0
  m <- suppressMessages(lme4::lmer(y ~ A + (1 | block), data = d))
  skip_if_not(lme4::isSingular(m))
  expect_message(pm <- predictmeans(m, "A", plot = FALSE, prtnum = FALSE), "singular")
  expect_true(pm$settings$singular)
  expect_message(contrastmeans(m, "A", ctrmatrix = rbind(c(-1, 1, 0))), "singular")
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  m2 <- lme4::lmer(yield ~ nitro + (1 | Block), data = Oats)
  expect_false(predictmeans(m2, "nitro", plot = FALSE, prtnum = FALSE)$settings$singular)
  expect_true(is.na(predictmeans(lm(yield ~ nitro, data = Oats), "nitro", plot = FALSE,
                                 prtnum = FALSE)$settings$singular))
})

test_that("Tukey and Scheffe adjustments use each pair's own df (SAS ADJDFE=ROW)", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  Ou <- Oats[-c(1, 2, 30), ]                          # unbalanced: the pairs' KR df differ
  m <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Ou)
  pr <- predictmeans(m, "Variety", adj = "tukey", plot = FALSE, prtnum = FALSE)$pairs
  expect_gt(diff(range(pr$df)), 0.01)
  expect_equal(pr$p.adj, ptukey(sqrt(2) * abs(pr$statistic), 3, pr$df, lower.tail = FALSE))
  ps <- predictmeans(m, "Variety", adj = "scheffe", plot = FALSE, prtnum = FALSE)$pairs
  expect_equal(ps$p.adj, pf(ps$statistic^2 / 2, 2, ps$df, lower.tail = FALSE))
})
