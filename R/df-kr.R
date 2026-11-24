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
  A <- rowSums((a %*% W) * a)
  ddf <- pm_kr_ddf_formula(A, A, q = 1)
  stats::setNames(ddf$ddf, rownames(L))
}

## Kenward-Roger denominator df from pbkrtest's two sums A1 and A2 (vectors
## allowed) for a hypothesis of rank q: the one place of this formula (used by
## pm_kr_ddf() and pm_kr_ftest()). Returns list(ddf, A2). A df that is not
## positive and finite means the approximation broke down (e.g. an indefinite
## information matrix of the variance parameters): it stops with a clear
## message rather than giving a test with a nonsense df. NA/NaN (rows that
## cannot be estimated) are passed through.
pm_kr_ddf_formula <- function(A1, A2, q) {
  B <- (1/(2 * q)) * (A1 + 6 * A2)
  g <- ((q + 1) * A1 - (q + 4) * A2)/((q + 2) * A2)
  den <- 3 * q + 2 * (1 - g)
  c1 <- g/den
  c2 <- (q - g)/den
  c3 <- (q + 2 - g)/den
  V0 <- 1 + c1 * B
  V1 <- 1 - c2 * B
  V2 <- 1 - c3 * B
  V0 <- ifelse(abs(V0) < 1e-10, 0, V0)
  x <- 1 - A2/q
  ratio <- ifelse(abs(x) < 1e-14 & abs(V1) < 1e-14, 1, x/V1)
  rho <- 1/q * ratio^2 * V0/V2
  ddf <- 4 + (q + 2)/(q * rho - 1)
  bad <- !is.na(ddf) & (!is.finite(ddf) | ddf <= 0)
  if (any(bad)) {
    stop("Kenward-Roger df could not be computed (", format(ddf[bad][1], digits = 4), "): the approximation ",
         "breaks down for this model, e.g. a variance parameter near its boundary or a poorly identified one. ",
         "Use df_method = \"Satterthwaite\" or \"KR-expected\".", call. = FALSE)
  }
  list(ddf = ddf, A2 = A2)
}
