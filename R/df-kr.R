## Kenward-Roger denominator df for many one-row contrasts at once.
##
## Same result as calling pbkrtest::Lb_ddf(L[k, ], V0, Vadj) for every row k,
## but vectorised. For a single row l, pbkrtest's two sums are equal:
## A1 = A2 = a' W a, with a_i = (l V0 P_i V0 l') / (l V0 l'), where P_i and W
## are the "P" and "W" attributes of the adjusted covariance from
## pbkrtest::vcovAdj(). The rest is pbkrtest's formula with q = 1, applied
## element by element.
pm_kr_ddf <- function(L, V0, Vadj) {
  if (!is.matrix(L)) L <- matrix(L, nrow = 1)
  P <- attr(Vadj, "P")
  W <- as.matrix(attr(Vadj, "W"))
  V0 <- as.matrix(V0)
  # pbkrtest sums over the upper triangle of W; mirror it to match exactly.
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  s <- rowSums((L %*% V0) * L)
  a <- vapply(P, function(Pi) {
    M <- V0 %*% as.matrix(Pi) %*% V0
    rowSums((L %*% M) * L) / s
  }, numeric(nrow(L)))
  a <- matrix(a, nrow = nrow(L))
  A1 <- A2 <- rowSums((a %*% W) * a)
  q <- 1
  B <- (1/(2 * q)) * (A1 + 6 * A2)
  g <- ((q + 1) * A1 - (q + 4) * A2)/((q + 2) * A2)
  c1 <- g/(3 * q + 2 * (1 - g))
  c2 <- (q - g)/(3 * q + 2 * (1 - g))
  c3 <- (q + 2 - g)/(3 * q + 2 * (1 - g))
  V0s <- 1 + c1 * B
  V1 <- 1 - c2 * B
  V2 <- 1 - c3 * B
  V0s <- ifelse(abs(V0s) < 1e-10, 0, V0s)
  x <- 1 - A2/q
  ratio <- ifelse(abs(x) < 1e-14 & abs(V1) < 1e-14, 1, x/V1)
  rho <- 1/q * ratio^2 * V0s/V2
  df2 <- 4 + (q + 2)/(q * rho - 1)
  stats::setNames(df2, rownames(L))
}
