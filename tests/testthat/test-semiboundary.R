## semiboundary(): the region covered by two-dimensional data, and
## pm_points_in_poly() with several rings (R/semiboundary.R, R/refgrid.R)

in_b <- function(P, b) pm_points_in_poly(as.matrix(P), b)

test_that("the scallop survey: one piece, every point inside, close to the boundary of the original analysis", {
  data(scallop, package = "predictmeans", envir = environment())
  data(scallopBdry, package = "predictmeans", envir = environment())
  b <- semiboundary(scallop[, c("longitude", "latitude")])
  expect_identical(names(b), c("longitude", "latitude"))
  expect_identical(attr(b, "scale"), "lonlat")
  expect_false(anyNA(b))                                   # one ring
  expect_true(all(in_b(scallop[, c("longitude", "latitude")], b)))
  fine <- expand.grid(longitude = seq(-74, -71, length.out = 120), latitude = seq(38.5, 41.2, length.out = 120))
  a <- in_b(fine, b)
  r <- in_b(fine, scallopBdry)
  expect_gt(sum(a & r) / sum(a | r), 0.82)                 # intersection over union (0.85)
  # a larger mult: a larger region, still one piece
  b3 <- semiboundary(scallop[, c("longitude", "latitude")], mult = 4)
  expect_gt(sum(in_b(fine, b3)), sum(a))
})

test_that("clusters apart give separate pieces and an empty middle gives a hole", {
  set.seed(1)
  two <- data.frame(x = c(runif(80, 0, 1), runif(80, 5, 6)), y = runif(160, 0, 1))
  b <- semiboundary(two)
  expect_identical(attr(b, "scale"), "sd")
  expect_equal(sum(is.na(b$x)), 1)                         # two rings
  expect_true(all(in_b(two, b)))
  expect_false(in_b(data.frame(3, 0.5), b))                 # between the clusters
  # points on an annulus: the centre is a hole
  th <- runif(400, 0, 2 * pi)
  rr <- runif(400, 2, 3)
  ann <- data.frame(x = rr * cos(th), y = rr * sin(th))
  ba <- semiboundary(ann)
  expect_equal(sum(is.na(ba$x)), 1)                        # the outer ring and the hole
  expect_true(all(in_b(ann, ba)))
  expect_false(in_b(data.frame(0, 0), ba))
  expect_true(in_b(data.frame(2.5, 0), ba))
})

test_that("pm_points_in_poly(): one polygon as before; rings by NA rows with the even-odd rule", {
  sq <- cbind(c(0, 4, 4, 0), c(0, 0, 4, 4))
  hole <- cbind(c(1, 3, 3, 1), c(1, 1, 3, 3))
  P <- rbind(c(0.5, 0.5), c(2, 2), c(5, 5))
  expect_identical(pm_points_in_poly(P, sq), c(TRUE, TRUE, FALSE))
  expect_identical(pm_points_in_poly(P, rbind(sq, NA, hole)), c(TRUE, FALSE, FALSE))
  far <- sq + 10
  expect_identical(pm_points_in_poly(rbind(P, c(12, 12)), rbind(sq, NA, far)), c(TRUE, TRUE, FALSE, TRUE))
})

test_that("scale: coordinates in the same unit are not stretched", {
  # points on a 4 x 1 strip in metres: with "same" the closing distance is in metres and
  # the same along both axes; "sd" divides each axis by its own SD
  set.seed(3)
  st <- data.frame(east = runif(300, 0, 400), north = runif(300, 0, 100))
  bs <- semiboundary(st, scale = "same")
  expect_identical(attr(bs, "scale"), "same")
  expect_true(all(in_b(st, bs)))
  expect_equal(range(bs$east, na.rm = TRUE), c(0, 400), tolerance = 0.05)
  expect_identical(attr(semiboundary(st), "scale"), "sd")
  expect_error(semiboundary(st, scale = "metres"), "should be one of")
})

test_that("small holes are filled, an area without data stays a hole", {
  g <- expand.grid(x = 1:30, y = 1:30)
  lake <- g[!(g$x %in% 10:15 & g$y %in% 10:15), ]               # 36 plots missing
  b <- semiboundary(lake, scale = "same")
  expect_false(in_b(data.frame(12.5, 12.5), b))
  gap <- g[!(g$x %in% 10:12 & g$y %in% 10:12), ]                # 9 plots missing: filled
  expect_true(in_b(data.frame(11, 11), semiboundary(gap, scale = "same")))
  set.seed(1)
  u <- data.frame(x = runif(10000), y = runif(10000))           # no spurious holes in dense random data
  bu <- semiboundary(u)
  expect_false(anyNA(bu))
  expect_true(all(in_b(u, bu)))
})

test_that("arguments are checked", {
  expect_error(semiboundary(data.frame(x = 1:5)), "two numeric columns")
  expect_error(semiboundary(1:2, 1:2), "three distinct points")
  expect_error(semiboundary(data.frame(x = rep(1, 5), y = 1:5)), "must vary")
  expect_error(semiboundary(data.frame(x = runif(9), y = runif(9)), d = -1), "positive")
  b <- semiboundary(runif(50), runif(50))
  expect_identical(names(b), c("x", "y"))
})

test_that("semipred(boundary = TRUE) keeps the grid points inside the data's region", {
  set.seed(2)
  d <- data.frame(x1 = c(runif(60, 0, 1), runif(60, 2, 3)), x2 = runif(120, 0, 1))
  d$y <- sin(2 * d$x1) + d$x2 + rnorm(120, 0, 0.2)
  fit <- semireg(y ~ s(x1, x2, bs = "ztps", k = 4), data = d, prt = FALSE)
  p_all <- semipred(fit, covariate = c("x1", "x2"), prt = FALSE)
  p_b <- semipred(fit, covariate = c("x1", "x2"), boundary = TRUE, prt = FALSE)
  p_s <- semipred(fit, covariate = c("x1", "x2"), boundary = semiboundary(d[, c("x1", "x2")]), prt = FALSE)
  expect_lt(nrow(p_b$pred_df), nrow(p_all$pred_df))
  expect_identical(p_b$pred_df, p_s$pred_df)
  expect_false(any(p_b$pred_df$x1 > 1.3 & p_b$pred_df$x1 < 1.7))   # the gap between the clusters
  expect_error(semipred(fit, covariate = "x1", boundary = TRUE, prt = FALSE), "two covariates")
})
