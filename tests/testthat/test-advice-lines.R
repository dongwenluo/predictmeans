## Advice rules of covariatemeans() and pm_slopes() (R/advice.R)

ids_of <- function(x) pm_advice_get(x)$id

lines_data <- function() {
  set.seed(6142)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)), x = runif(90, 0, 10))
  d$y <- 1 + c(a = 0.5, b = 1, c = 1.5)[d$g] * d$x + rnorm(90)
  d
}

test_that("pm_slopes(): equal by model, at a value, extrapolation", {
  d <- lines_data()
  mi <- lm(y ~ g * x, d)
  expect_length(ids_of(pm_slopes(mi, "x", "g")), 0)
  expect_identical(ids_of(pm_slopes(lm(y ~ g + x, d), "x", "g")), "slopes_equal_by_model")
  expect_identical(ids_of(pm_slopes(lm(y ~ g * poly(x, 2), d), "x", "g", adj = "holm")), "slopes_at_value")
  out <- pm_slopes(mi, "x", "g", at = 15)
  expect_identical(ids_of(out), "extrapolation")
  expect_true(any(grepl("^Advice:", capture.output(print(out)))))
})

test_that("covariatemeans(): extrapolation and pointwise regions", {
  skip_on_cran()
  d <- lines_data()
  mi <- lm(y ~ g * x, d)
  expect_length(ids_of(covariatemeans(mi, "g", "x")), 0)
  out <- covariatemeans(mi, "g", "x", at = seq(0, 15, 1))
  expect_identical(ids_of(out), "extrapolation")
  expect_true(any(grepl("^Advice:", capture.output(print(out)))))
  expect_identical(ids_of(covariatemeans(mi, "g", "x", contr = c("a", "b"))), "pointwise_regions")
  expect_length(ids_of(covariatemeans(mi, "g", "x", contr = c("a", "b"), band = "simultaneous")), 0)
})
