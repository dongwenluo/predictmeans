## Square roots of mixed-model covariance matrices for permutation tests.
##
## Residual-permutation tests (permmodels() for lme/gls/lmer, permlmer())
## decorrelate residuals with a square root W of their covariance V
## (W W' = V), permute them, and correlate them again with W. predictmeans
## 1.x used the Cholesky factor. That factor depends on the order of the data
## (within a block it weights the first unit differently from the others), so
## the results changed when the data were sorted differently, and tests of
## terms that vary within blocks rejected too often (6-8% at the 5% level in
## simulations; 5% with the root below). See dev/permutation-plan.md.
##
## pm_whiten() returns functions to multiply by W and by W^(-1), where W is
## unchanged in form when the units are reordered (W for P V P' is P W P').

## V = D + A A' with D diagonal (d = its diagonal, NULL for the identity) and
## A = Z Lambda (n x q, usually sparse). W = D^(1/2) (I + B K B'), with
## B = D^(-1/2) A and K = R diag((sqrt(1 + s) - 1) / s) R' from the eigen
## decomposition B'B = R diag(s) R'. For D = I, W is the symmetric square
## root of V. Only q x q matrices are formed.
pm_whiten <- function(A, d = NULL) {
  n <- nrow(A)
  sd <- if (is.null(d)) rep(1, n) else sqrt(d)
  B <- Matrix::Diagonal(x = 1 / sd) %*% A
  q <- ncol(B)
  if (q == 0) {
    K_up <- K_dn <- matrix(0, 0, 0)
  } else {
    e <- eigen(as.matrix(Matrix::crossprod(B)), symmetric = TRUE)
    keep <- e$values > max(e$values, 1) * 1e-10
    R <- e$vectors[, keep, drop = FALSE]
    s <- e$values[keep]
    K_up <- R %*% (((sqrt(1 + s) - 1) / s) * t(R))
    K_dn <- R %*% (((1 / sqrt(1 + s) - 1) / s) * t(R))
  }
  low <- function(x, K) {
    if (q == 0) return(x)
    x + as.matrix(B %*% (K %*% as.matrix(Matrix::crossprod(B, x))))
  }
  out <- function(x, v) if (is.null(dim(x))) as.vector(v) else v
  list(
    mult = function(x) out(x, sd * low(as.matrix(x), K_up)),
    solve = function(x) out(x, low(as.matrix(x) / sd, K_dn))
  )
}

## The same for a block-diagonal V given as a list of dense blocks (nlme
## models): the symmetric square root of each block.
pm_whiten_blocks <- function(V_list) {
  roots <- lapply(V_list, function(V) {
    e <- eigen(as.matrix(V), symmetric = TRUE)
    val <- pmax(e$values, 0)
    list(r = e$vectors %*% (sqrt(val) * t(e$vectors)),
         ri = e$vectors %*% (ifelse(val > 0, 1 / sqrt(val), 0) * t(e$vectors)))
  })
  W <- Matrix::bdiag(lapply(roots, `[[`, "r"))
  Wi <- Matrix::bdiag(lapply(roots, `[[`, "ri"))
  out <- function(x, v) if (is.null(dim(x))) as.vector(v) else as.matrix(v)
  list(mult = function(x) out(x, W %*% x), solve = function(x) out(x, Wi %*% x))
}
