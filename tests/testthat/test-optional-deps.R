test_that("internal point-in-polygon matches HRW::pointsInPoly for closed polygons", {
  skip_if_not_installed("HRW")
  set.seed(5)
  pts <- cbind(runif(500, -1, 2), runif(500, -1, 2))
  poly_closed <- rbind(c(0, 0), c(1.5, 0.2), c(1.2, 1.4), c(0.3, 1.1), c(-0.2, 0.6), c(0, 0))
  expect_identical(predictmeans:::pm_points_in_poly(pts, poly_closed), HRW::pointsInPoly(pts, poly_closed))
  # An open polygon (last vertex not repeated) is closed automatically.
  expect_identical(predictmeans:::pm_points_in_poly(pts, poly_closed[-6, ]), predictmeans:::pm_points_in_poly(pts, poly_closed))
})

test_that("smZ() was removed in 2.0: it stops with the formula form", {
  expect_error(smZ(seq(0, 1, length = 50)), "s(x, bs = \"os\")", fixed = TRUE, class = "defunctError")
})
