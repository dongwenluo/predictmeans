## Finite differences: the package's own numerical derivatives, one rule.
##
## Central differences with the step h = c * max(|x|, scale): c = 1e-5 for
## first derivatives, 1e-4 for second, 1e-3 for third (about eps^(1/(d + 2)),
## which balances the truncation error, of order h^2, against rounding,
## eps / h^d), and 'scale' says on what scale x changes:
##   1      parameters on an unconstrained scale (nlme's modelStruct
##          coefficients and log sigma; the link scale of a back-
##          transformation);
##   0.01   lme4's theta and sigma (semireg() fits), whose size can be small;
##   spread a covariate, with x = 0 (no location: a covariate measured in
##          years from 0 or from 2000 gets the same step): the covariate's
##          SD (pm_slopes()) or the range of the grid (semipred(deriv)).
## The error is about 1e-10 (relative) for first derivatives and 1e-8 for
## second ones of functions computed to full precision. These are used
## where derivatives are taken many times (df of many points, the pieces of
## a covariance, back-transformations of many means). numDeriv (Richardson
## extrapolation, more evaluations) is kept where a derivative is computed
## once per model: the Satterthwaite df of lmer fits (equal to lmerTest's),
## the Hessian of the variance parameters of lme4 fits and varcomp()'s delta
## method.

pm_fd_h <- function(x, scale = 1, order = 1) 10^(order - 6) * pmax(abs(x), scale)

## The derivative of f (a number, vector or matrix, also sparse) with respect
## to x[i].
pm_fd_d1 <- function(f, x, i, h = pm_fd_h(x[i])) {
  xp <- xm <- x
  xp[i] <- x[i] + h
  xm[i] <- x[i] - h
  (f(xp) - f(xm)) / (2 * h)
}

## The Jacobian of a vector-valued f at x (a row per value of f).
pm_fd_jac <- function(f, x, scale = 1) {
  h <- pm_fd_h(x, scale)
  G <- matrix(0, length(f(x)), length(x))
  for (j in seq_along(x)) G[, j] <- pm_fd_d1(f, x, j, h[j])
  G
}

## The Hessian of a scalar f over the elements idx of x (steps h for all of x).
pm_fd_hess <- function(f, x, idx = seq_along(x), h = pm_fd_h(x, order = 2)) {
  at <- function(j, d) replace(x, j, x[j] + d)
  p <- length(idx)
  H <- matrix(0, p, p)
  f0 <- f(x)
  for (a in seq_len(p)) {
    i <- idx[a]
    H[a, a] <- (f(at(i, h[i])) - 2 * f0 + f(at(i, -h[i]))) / h[i]^2
    for (b in seq_len(a - 1)) {
      j <- idx[b]
      H[a, b] <- H[b, a] <- (f(at(c(i, j), c(h[i], h[j]))) - f(at(c(i, j), c(h[i], -h[j]))) -
                               f(at(c(i, j), c(-h[i], h[j]))) + f(at(c(i, j), c(-h[i], -h[j])))) / (4 * h[i] * h[j])
    }
  }
  H
}

## The derivative of a vectorised function of one variable, as a function
## ('order': of the derivative this gives when f is itself one, so that
## nested derivatives take the step of their order).
pm_fd_fun <- function(f, scale = 1, order = 1) {
  function(x) {
    h <- pm_fd_h(x, scale, order)
    (f(x + h) - f(x - h)) / (2 * h)
  }
}
