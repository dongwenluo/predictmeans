
test_that("curves along a covariate aliased with a factor stop clearly; differences between curves still work", {
  skip_on_cran()
  data(DataExam5.2, package = "predictmeans")
  m <- lm(ht ~ site + sitemean*seedlot, DataExam5.2)          # sitemean is constant within site
  expect_error(covariatemeans(m, "seedlot", "sitemean", plot = FALSE), "cannot be estimated.*contr")
  lv <- levels(DataExam5.2$seedlot)
  d <- suppressMessages(covariatemeans(m, "seedlot", "sitemean", contr = lv[1:2], plot = FALSE))
  dd <- as.data.frame(d)
  expect_false(anyNA(dd$estimate))
  b <- coef(m)                                                     # level 1 minus level 2 at each sitemean
  expect_equal(dd$estimate, unname(-b[paste0("seedlot", lv[2])] - b[paste0("sitemean:seedlot", lv[2])] * dd$sitemean),
               tolerance = 1e-8)
})
