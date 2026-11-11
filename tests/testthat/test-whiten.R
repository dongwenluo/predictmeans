## Order-invariant square roots used by the residual-permutation tests

whiten_case <- function() {
  set.seed(1)
  d <- expand.grid(u = 1:6, block = factor(1:5))
  d$g <- factor(rep(1:10, 3)); d$x <- rnorm(30)
  d$y <- 2 * rnorm(5)[d$block] + rnorm(10)[d$g] + rnorm(30)
  m <- lmer(y ~ x + (1 | block) + (1 | g), data = d)
  list(d = d, A = getME(m, "Z") %*% getME(m, "Lambda"))
}

test_that("pm_whiten() gives the symmetric square root of I + A A'", {
  cs <- whiten_case()
  V <- as.matrix(Matrix::tcrossprod(cs$A)) + diag(30)
  W <- pm_whiten(cs$A)
  Wm <- W$mult(diag(30))
  expect_equal(Wm %*% t(Wm), V, tolerance = 1e-10)
  expect_equal(Wm, t(Wm), tolerance = 1e-10)
  e <- eigen(V, symmetric = TRUE)
  expect_equal(Wm, e$vectors %*% (sqrt(e$values) * t(e$vectors)), tolerance = 1e-10)
  expect_equal(W$solve(W$mult(cs$d$y)), cs$d$y, tolerance = 1e-10)
})

test_that("pm_whiten() with a diagonal D is a square root and does not depend on the order", {
  cs <- whiten_case()
  dd <- seq(0.5, 3, length.out = 30)
  W <- pm_whiten(cs$A, dd)$mult(diag(30))
  expect_equal(W %*% t(W), diag(dd) + as.matrix(Matrix::tcrossprod(cs$A)), tolerance = 1e-10)
  p <- c(30:16, 1:15)
  Wp <- pm_whiten(cs$A[p, ], dd[p])$mult(diag(30))
  expect_equal(Wp, W[p, p], tolerance = 1e-10)
  V <- as.matrix(Matrix::tcrossprod(cs$A)) + diag(30)
  # the Cholesky factor used before 2.0 does depend on the order
  expect_gt(max(abs(t(chol(V[p, p])) - t(chol(V))[p, p])), 0.1)
})

test_that("pm_whiten_blocks() gives symmetric roots of each block", {
  Vl <- list(matrix(c(2, 1, 1, 2), 2), matrix(c(3, 1, 0.5, 1, 3, 1, 0.5, 1, 3), 3))
  W <- pm_whiten_blocks(Vl)
  Wm <- W$mult(diag(5))
  expect_equal(Wm %*% t(Wm), as.matrix(Matrix::bdiag(Vl)), tolerance = 1e-10)
  expect_equal(Wm, t(Wm), tolerance = 1e-10)
  expect_equal(W$solve(W$mult(1:5)), 1:5, tolerance = 1e-10)
})

test_that("permlmer() and permmodels() run with the new square root", {
  set.seed(2)
  d <- expand.grid(sub = 1:3, wp = 1:4, block = 1:6)
  d$A <- factor(d$wp); d$block <- factor(d$block); d$wpid <- interaction(d$block, d$wp)
  d$y <- rnorm(6)[d$block] + 0.5 * rnorm(24)[d$wpid] + rnorm(nrow(d))
  m0 <- lmer(y ~ A + (1 | block), data = d)
  m1 <- lmer(y ~ A + (1 | block) + (1 | wpid), data = d)
  p1 <- suppressMessages(permlmer(m0, m1, nperm = 19, ncore = 1, seed = 1))$`Perm-p`[2]
  expect_true(p1 > 0 && p1 <= 1)
  fm <- lmer(y ~ A + (1 | block), data = d)
  pm <- permmodels(fm, nperm = 9, prt = FALSE, ncore = 1, seed = 1)
  expect_length(pm$permlist, 9)
})
