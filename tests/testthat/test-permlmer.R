## permlmer(): an LRT of 0 up to optimizer noise must not look extreme
test_that("permlmer() gives p = 1 when the tested variance is estimated at 0", {
  set.seed(22)
  d <- expand.grid(sub = 1:3, wp = 1:4, block = 1:6)
  d$A <- factor(d$wp); d$block <- factor(d$block); d$wpid <- interaction(d$block, d$wp)
  d$y <- rnorm(6)[d$block] + rnorm(nrow(d))
  ctl <- lmerControl(check.conv.singular = "ignore")
  m0 <- lmer(y ~ A + (1 | block), data = d, control = ctl)
  m1 <- lmer(y ~ A + (1 | block) + (1 | wpid), data = d, control = ctl)
  lrt <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  expect_lt(abs(lrt), 1e-4)
  res <- suppressMessages(permlmer(m0, m1, nperm = 19, ncore = 1, seed = 1))
  expect_equal(res$`Perm-p`[2], 1)
})
