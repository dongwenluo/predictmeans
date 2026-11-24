test_that("permmodels() says when it replaces the user's test_statistic", {
  skip_if_not_installed("nlme")
  dat <- data.frame(y = c(5.1, 6.3, 4.8, 7.2, 6.6, 5.9, 8.1, 7.4, 6.8, 7.9, 5.5, 6.1),
                    A = factor(rep(c("a", "b", "c"), 4)),
                    x = c(1, 3, 2, 5, 4, 3, 6, 5, 4, 6, 2, 3))
  m_gls <- nlme::gls(y ~ A + x, data = dat)
  expect_message(permmodels(m_gls, nperm = 9, type = "II", test_statistic = "F", method = "raw", seed = 3017),
                 "\"Chisq\" is used")
  m_lm <- lm(y ~ A + x, data = dat)
  expect_message(permmodels(m_lm, nperm = 9, type = "I", test_statistic = "Chisq", method = "raw", seed = 3017),
                 "Type I tests are F tests")
  expect_message(permmodels(m_lm, nperm = 9, type = "II", test_statistic = "Chisq", seed = 3017),
                 "only with method = \"raw\"")
  # default left alone: no message
  expect_no_message(permmodels(m_lm, nperm = 9, type = "I", method = "raw", seed = 3017))
  expect_no_message(permmodels(m_lm, nperm = 9, type = "II", seed = 3017))
  # asked statistic used: no message
  m_glm <- glm(y ~ A + x, data = dat)
  expect_no_message(permmodels(m_glm, nperm = 9, type = "II", test_statistic = "LR", method = "raw", seed = 3017))
})
