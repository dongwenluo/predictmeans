test_that("internal point-in-polygon matches HRW::pointsInPoly for closed polygons", {
  skip_if_not_installed("HRW")
  set.seed(5)
  pts <- cbind(runif(500, -1, 2), runif(500, -1, 2))
  poly_closed <- rbind(c(0, 0), c(1.5, 0.2), c(1.2, 1.4), c(0.3, 1.1), c(-0.2, 0.6), c(0, 0))
  expect_identical(predictmeans:::pm_points_in_poly(pts, poly_closed), HRW::pointsInPoly(pts, poly_closed))
  # An open polygon (last vertex not repeated) is closed automatically.
  expect_identical(predictmeans:::pm_points_in_poly(pts, poly_closed[-6, ]), predictmeans:::pm_points_in_poly(pts, poly_closed))
})

test_that("optional spline types and lme variance components work when installed", {
  x <- seq(0, 1, length = 50)
  skip_if_not_installed("splines2")
  Z <- smZ(x, k = 5, type = "bSpline")
  expect_equal(nrow(Z), 50)
  skip_if_not_installed("lmeInfo")
  fm <- lme(yield ~ Variety, random = ~1|Block, data = nlme::Oats)
  vc <- varcomp(fm)
  expect_true(all(c("vcov", "SE") %in% colnames(vc)))
})
