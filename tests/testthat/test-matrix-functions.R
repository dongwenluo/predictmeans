test_that("pm_place_blocks() and pm_sparse_components() rebuild crossed block structures", {
  set.seed(5512)
  n <- 12
  vg <- factor(rep(1:4, each = 3))
  zg <- factor(c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5))
  mk <- function(k) { A <- matrix(rnorm(k * k), k); crossprod(A) + diag(k) }
  Vl <- stats::setNames(lapply(table(vg), mk), levels(vg))
  Zl <- stats::setNames(lapply(table(zg), mk), levels(zg))
  pos <- function(g) split(seq_len(n), g)
  S <- pm_place_blocks(Vl, pos(vg), n) + pm_place_blocks(Zl, pos(zg), n)
  D <- matrix(0, n, n)
  for (k in levels(vg)) { i <- vg == k; D[i, i] <- D[i, i] + Vl[[k]] }
  for (k in levels(zg)) { i <- zg == k; D[i, i] <- D[i, i] + Zl[[k]] }
  expect_equal(as.matrix(S), D, tolerance = 0, ignore_attr = TRUE)
  # residual groups 1 and 2 are linked by random-effect group 2
  expect_identical(pm_sparse_components(S), rep(1:3, c(6, 3, 3)))
  expect_identical(pm_sparse_components(Matrix::Diagonal(3)), 1:3)
})

test_that("matrix_list() splits rows, columns, blocks and vectors", {
  x <- matrix(1:16, 4)
  f <- factor(c("a", "a", "b", "b"))
  expect_identical(matrix_list(x, f, "row")[[2]], x[3:4, ])
  expect_identical(matrix_list(x, f, "both")[[1]], x[1:2, 1:2])
  expect_identical(matrix_list(1:4, f, "both")$b, diag(3:4))
  expect_error(matrix_list(1:4, f, "row"), "must be a matrix")
})
