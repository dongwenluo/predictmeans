## Advice rules of semireg() and semipred() (R/advice.R)

ids_of <- function(x) pm_advice_get(x)$id

sr_data <- function() {
  set.seed(3318)
  d <- data.frame(x = runif(120, 0, 10), g = factor(rep(1:6, 20)))
  d$y <- sin(d$x) + rnorm(120, 0, 0.3)
  d$z <- d$x + rnorm(120)
  d
}

test_that("semireg(): a smooth estimated as a straight line is advised, a curve is not", {
  d <- sr_data()
  curved <- suppressMessages(semireg(y ~ s(x) + (1 | g), data = d))
  expect_length(ids_of(curved), 0)
  straight <- suppressMessages(semireg(z ~ s(x), data = d))
  expect_identical(ids_of(straight), "smooth_linear")
  expect_true(any(grepl("^Advice:", capture.output(print(straight)))))
})

test_that("semipred(): extrapolation and pointwise regions", {
  skip_on_cran()
  d <- sr_data()
  fit <- suppressMessages(semireg(y ~ s(x) + (1 | g), data = d))
  expect_length(ids_of(semipred(fit, "x", "x")), 0)
  out <- semipred(fit, "x", "x", at = seq(-2, 12, length = 30))
  expect_identical(ids_of(out), "extrapolation")
  expect_true(any(grepl("^Advice:", capture.output(print(out)))))
  expect_identical(ids_of(semipred(fit, "x", "x", deriv = 1)), "pointwise_regions")
  expect_length(ids_of(semipred(fit, "x", "x", deriv = 1, band = "simultaneous")), 0)
})
