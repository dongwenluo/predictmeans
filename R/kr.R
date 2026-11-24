## Kenward-Roger (1997) adjusted covariance and F test for lmer models.
##
## Ported from package pbkrtest (Halekoh & Hojsgaard, GPL >= 2): vcovAdj()
## and the F test of KRmodcomp(). Same results (see the tests), so pbkrtest
## is not needed; the adjusted covariance is computed without n x n matrices.

## Mirror the upper triangle onto the lower one (as Matrix::forceSymmetric()).
pm_upper_sym <- function(A) {
  A <- as.matrix(A)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A
}

## Kenward-Roger adjusted covariance of the fixed effects (pbkrtest::vcovAdj()),
## with attributes "P" and "W" used for the degrees of freedom.
##
## W, the covariance of the variance parameters, is the inverse of their
## observed information (the default, as SAS PROC MIXED DDFM=KR; valid also
## with missing values, see dev/sim-df/README.md), or with information =
## "expected" of their expected information (pbkrtest, Kenward and Roger
## 1997), or with information = "average" of the average information
## y'P G_r P G_s P y / 2 (Gilmour, Thompson and Cullis 1995: the matrix of
## Genstat's and ASReml's REML iterations, which gives Genstat's df). Observed:
## for V linear in the parameters, -d2 l_R / ds_r ds_s =
## y'P G_r P G_s P y - tr(P G_r P G_s)/2 (the trace term is IE2/2 below), so
## the average information is the mean of the observed and expected ones.
## Both agree in balanced designs; with one plot missing in the split plot of
## dev/sas the observed one reproduces SAS's df (9.79, 31.3) where the
## expected one gives 9.87, 31.46.
##
## pbkrtest inverts the n x n covariance of the data, Sigma (O(n^3) time,
## O(n^2) memory). Here Sigma = sigma^2 (I + U U'), U = Z Lambda (n x q), so by
## the Woodbury identity Sigma^-1 = (I - U M U') / sigma^2 with
## M = (I + U'U)^-1, and every term pbkrtest needs reduces to q x q or q x p
## matrices: with E_k the derivative of the random-effect covariance for
## parameter k (G_k = Z E_k Z'), A = Z' Sigma^-1 Z and B = Z' Sigma^-2 Z,
##   tr(G_r Sigma^-1 G_s Sigma^-1) = tr(E_r A E_s A),
##   tr(G_r Sigma^-2) = tr(E_r B), and so on. Same results as pbkrtest.
pm_vcov_kr <- function(model, information = c("observed", "expected", "average")) {
  information <- match.arg(information)
  if (!lme4::getME(model, "is_REML")) model <- pm_refit(model, REML = TRUE)
  Phi <- as.matrix(stats::vcov(model))
  X <- as.matrix(lme4::getME(model, "X"))
  Z <- lme4::getME(model, "Z")
  Lambda <- lme4::getME(model, "Lambda")
  s2 <- stats::sigma(model)^2
  n <- nrow(X)
  q <- ncol(Z)

  # prior weights w: Sigma = sigma^2 (W^-1 + U U'), W = diag(w), and the
  # residual variance parameter has derivative R = W^-1 (not I). Whitened with
  # D = W^1/2: Sigma^-1 = D (I + Ut Ut')^-1 D / sigma^2, Ut = D U. Below,
  # "Sigma^-2" stands for Sigma^-1 R Sigma^-1 (the two agree when w = 1).
  w <- stats::weights(model)
  sw <- sqrt(w)
  U <- Z %*% Lambda
  Ut <- sw * U
  K <- as.matrix(Matrix::crossprod(Ut))
  M <- solve(diag(q) + K)
  Sinv <- function(Y) {                       # Sigma^-1 %*% Y
    Y <- sw * as.matrix(Y)
    sw * (Y - as.matrix(Ut %*% (M %*% as.matrix(Matrix::crossprod(Ut, Y)))))/s2
  }
  SX <- Sinv(X)                               # Sigma^-1 X        (n x p)
  SZ <- Sinv(Z)                               # Sigma^-1 Z        (n x q)
  A <- as.matrix(Matrix::crossprod(Z, SZ))     # Z' Sigma^-1 Z     (q x q)
  B <- crossprod(SZ, SZ/w)                    # Z' Sigma^-1 R Sigma^-1 Z (q x q)
  ZSX <- as.matrix(Matrix::crossprod(Z, SX))  # Z' Sigma^-1 X     (q x p)
  ZS2X <- crossprod(SZ, SX/w)                 # Z' Sigma^-1 R Sigma^-1 X (q x p)
  S2X <- Sinv(SX/w)                           # Sigma^-1 R Sigma^-1 X (n x p)
  MK <- M %*% K
  trS2 <- (n - 2 * sum(diag(MK)) + sum(MK * t(MK)))/s2^2   # tr((Sigma^-1 R)^2)

  # derivative matrices E_k (q x q) for each random-effect (co)variance
  Gp <- lme4::getME(model, "Gp")
  n_comp <- vapply(lme4::getME(model, "cnms"), length, integer(1))
  upper_entry <- function(k, N) {
    aa <- cumsum(N:1)
    aa_low <- c(0, aa[-length(aa)])
    i <- which(aa_low < k & k <= aa)
    j <- k - N * i + N - i * (3 - i)/2 + i
    c(i, j)
  }
  E <- list()
  # a variance estimated as zero (a singular fit): its parameter, and the
  # covariances that involve it, are left out of the adjustment (treated as
  # known, as SAS PROC MIXED); their near-singular information otherwise gave
  # nonsense df (e.g. 0.04 or 485371)
  vc <- lme4::VarCorr(model)
  boundary <- logical(0)
  for (ss in seq_along(n_comp)) {
    S_ss <- as.matrix(vc[[ss]])
    tol_ss <- 1e-8 * max(1, max(abs(diag(S_ss))), s2)
    for (rr in seq_len((n_comp[ss] + 1) * n_comp[ss]/2)) {
      ij <- unique(upper_entry(rr, n_comp[ss]))
      boundary <- c(boundary, any(diag(S_ss)[ij] <= tol_ss))
    }
  }
  for (ss in seq_along(n_comp)) {
    nc <- n_comp[ss]
    n_lev <- (Gp[ss + 1] - Gp[ss])/nc
    for (rr in seq_len((nc + 1) * nc/2)) {
      ij <- unique(upper_entry(rr, nc))
      EE <- matrix(0, nc, nc)
      EE[ij[1], ij[length(ij)]] <- EE[ij[length(ij)], ij[1]] <- 1
      blk <- Matrix::kronecker(Matrix::Diagonal(n_lev), Matrix::Matrix(EE, sparse = TRUE))
      idx <- (Gp[ss] + 1):Gp[ss + 1]
      Ek <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(q, q))
      Ek[idx, idx] <- blk
      E <- c(E, list(Ek))
    }
  }
  n_re <- length(E)
  n_g <- n_re + 1                             # last parameter: residual variance

  # P_k = -X' Sigma^-1 G_k Sigma^-1 X
  PP <- vector("list", n_g)
  EZSX <- lapply(E, function(Ek) as.matrix(Ek %*% ZSX))
  for (k in seq_len(n_re)) PP[[k]] <- pm_upper_sym(-crossprod(ZSX, EZSX[[k]]))
  PP[[n_g]] <- pm_upper_sym(-crossprod(SX, SX/w))

  # Q_rs = X' Sigma^-1 G_r Sigma^-1 G_s Sigma^-1 X and trace terms
  index <- function(i, j) if (i <= j) (i - 1) * (n_g - i/2) + j else (j - 1) * (n_g - j/2) + i
  QQ <- vector("list", n_g * (n_g + 1)/2)
  Ktrace <- matrix(NA, n_g, n_g)
  AE <- lapply(E, function(Ek) as.matrix(A %*% Ek))
  for (r in seq_len(n_g)) for (s in r:n_g) {
    if (s <= n_re) {
      QQ[[index(r, s)]] <- crossprod(EZSX[[r]], A %*% EZSX[[s]])
      Ktrace[r, s] <- Ktrace[s, r] <- sum(AE[[r]] * t(AE[[s]]))
    } else if (r <= n_re) {
      QQ[[index(r, s)]] <- crossprod(EZSX[[r]], ZS2X)
      Ktrace[r, s] <- Ktrace[s, r] <- sum(as.matrix(E[[r]]) * B)
    } else {
      QQ[[index(r, s)]] <- crossprod(SX/w, S2X)
      Ktrace[r, s] <- trS2
    }
  }

  IE2 <- matrix(NA, n_g, n_g)
  for (ii in seq_len(n_g)) {
    PhiP <- Phi %*% PP[[ii]]
    for (jj in ii:n_g) {
      IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] - 2 * sum(Phi * QQ[[index(ii, jj)]]) +
        sum(PhiP * (PP[[jj]] %*% Phi))
    }
  }
  if (information != "expected") {
    # observed information: y'P G_r P G_s P y - IE2/2, with P y = Sigma^-1 (y - X b);
    # average information: y'P G_r P G_s P y / 2
    Pv <- function(v) {                        # P %*% v
      sv <- Sinv(v)
      sv - SX %*% (Phi %*% crossprod(X, sv))
    }
    Py <- Sinv(lme4::getME(model, "y") - lme4::getME(model, "offset") - X %*% lme4::fixef(model))
    ZPy <- as.matrix(Matrix::crossprod(Z, Py))
    a <- c(lapply(E, function(Ek) as.matrix(Z %*% (Ek %*% ZPy))), list(Py/w))   # G_r P y
    Pa <- lapply(a, Pv)
    Yq <- matrix(NA, n_g, n_g)
    for (r in seq_len(n_g)) for (s in r:n_g) Yq[r, s] <- Yq[s, r] <- sum(a[[r]] * Pa[[s]])
    IO <- if (information == "observed") Yq - IE2/2 else Yq/2
  } else {
    IO <- IE2/2
  }
  # the expected information is positive (semi)definite; the observed and
  # average ones need not be away from the maximum: then the expected one is
  # used, with a message (pm_kr_info_check()); before 2.0.0 an indefinite
  # matrix was inverted anyway (a "covariance" that is not one, nonsense df)
  keep <- c(!boundary, TRUE)                   # the residual variance is always kept
  information_used <- information
  if (information != "expected" && !pm_kr_info_check(IO[keep, keep, drop = FALSE], information)) {
    IO <- IE2/2
    information_used <- "expected"
  }
  IOk <- IO[keep, keep, drop = FALSE]
  condi <- min(abs(eigen(IOk, symmetric = TRUE, only.values = TRUE)$values))
  WW <- matrix(0, n_g, n_g)
  WW[keep, keep] <- pm_upper_sym(if (condi > 1e-10) solve(IOk) else MASS::ginv(IOk))
  UU <- matrix(0, ncol(X), ncol(X))
  if (n_g > 1) {
    for (ii in 1:(n_g - 1)) for (jj in (ii + 1):n_g) {
      Qij <- QQ[[index(ii, jj)]]
      UU <- UU + WW[ii, jj] * (Qij - PP[[ii]] %*% Phi %*% PP[[jj]])
    }
  }
  UU <- UU + t(UU)
  for (ii in seq_len(n_g)) UU <- UU + WW[ii, ii] * (QQ[[index(ii, ii)]] - PP[[ii]] %*% Phi %*% PP[[ii]])
  PhiA <- Phi + 2 * (Phi %*% UU %*% Phi)
  dimnames(PhiA) <- dimnames(Phi)
  attr(PhiA, "P") <- PP
  attr(PhiA, "W") <- WW
  attr(PhiA, "condi") <- condi
  attr(PhiA, "information") <- information_used
  attr(PhiA, "boundary") <- boundary
  PhiA
}

## Kenward-Roger F test of L beta = 0 (as pbkrtest::KRmodcomp(model, L)$test["Ftest", ]).
## L is first reduced to linearly independent rows (q = its rank), as
## pbkrtest does; a df that is not positive and finite stops
## (pm_kr_ddf_formula() in R/df-kr.R).
pm_kr_ftest <- function(model, L, PhiA = pm_vcov_kr(model)) {
  if (!is.matrix(L)) L <- matrix(L, nrow = 1)
  r <- qr(t(L))
  L <- L[r$pivot[seq_len(r$rank)], , drop = FALSE]
  Phi <- as.matrix(stats::vcov(model))
  beta <- if (inherits(model, "merMod")) lme4::fixef(model) else if (inherits(model, "lme")) nlme::fixef(model) else
    stats::coef(model)
  Theta <- crossprod(L, solve(L %*% tcrossprod(Phi, L), L))
  P <- attr(PhiA, "P")
  W <- attr(PhiA, "W")
  A1 <- A2 <- 0
  ThetaPhi <- Theta %*% Phi
  n_g <- length(P)
  for (ii in seq_len(n_g)) for (jj in ii:n_g) {
    e <- if (ii == jj) 1 else 2
    ui <- ThetaPhi %*% P[[ii]] %*% Phi
    uj <- ThetaPhi %*% P[[jj]] %*% Phi
    A1 <- A1 + e * W[ii, jj] * (sum(diag(ui)) * sum(diag(uj)))
    A2 <- A2 + e * W[ii, jj] * sum(ui * t(uj))
  }
  q <- nrow(L)
  ddf <- pm_kr_ddf_formula(A1, A2, q)$ddf
  F_scaling <- if (abs(ddf - 2) < 0.01) 1 else ddf * (1 - A2/q)/(ddf - 2)
  Wald <- pm_wald_quad(L, beta, PhiA)
  Fstat <- F_scaling * Wald/q
  c(stat = Fstat, ndf = q, ddf = ddf, F.scaling = F_scaling,
    p.value = stats::pf(Fstat, q, ddf, lower.tail = FALSE))
}

## Is the observed (or average) information of the variance parameters
## positive definite? The one check for lmer (pm_vcov_kr()) and nlme
## (pm_vcov_kr_nlme()) models. FALSE, with a message, when an eigenvalue is
## not positive (relative to the largest): the fit is not at a clear maximum
## in some direction of the variance parameters (a flat or poorly identified
## likelihood), and its inverse is no covariance. The caller then uses the
## expected information, which is positive definite.
pm_kr_info_check <- function(I, information, tol = 1e-8) {
  ev <- eigen((I + t(I))/2, symmetric = TRUE, only.values = TRUE)$values
  if (length(ev) && min(ev) > tol * max(abs(ev), 1e-300)) return(TRUE)
  message("Kenward-Roger: the ", information, " information of the variance parameters is not positive definite ",
          "(smallest eigenvalue ", format(min(ev), digits = 3), "); using the expected information (as \"KR-expected\").")
  FALSE
}
