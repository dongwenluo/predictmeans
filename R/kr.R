## Kenward-Roger (1997) adjusted covariance and F test for lmer models.
##
## Ported from package pbkrtest (Halekoh & Hojsgaard, GPL >= 2):
## vcovAdj(), its internal .get_SigmaG(), and the F test of KRmodcomp().
## Same algorithm and results (see the tests), so pbkrtest is not needed.

## Covariance matrix of the data, Sigma = sum_i ggamma_i G_i, and its
## derivatives G_i with respect to the variance parameters ggamma (each
## random-effect term's (co)variances, then the residual variance).
pm_kr_sigma_g <- function(model) {
  if (!inherits(model, "lmerMod")) stop("'model' is not a Gaussian linear mixed model (lmer).")
  GGamma <- lme4::VarCorr(model)
  Gp <- lme4::getME(model, "Gp")
  n_rt <- length(Gp) - 1
  n_comp <- vapply(lme4::getME(model, "cnms"), length, integer(1))
  n_parm <- (n_comp + 1) * n_comp / 2
  n_lev <- diff(Gp) / n_comp
  ggamma <- NULL
  for (ii in seq_len(n_rt)) {
    Lii <- GGamma[[ii]]
    ggamma <- c(ggamma, Lii[lower.tri(Lii, diag = TRUE)])
  }
  ggamma <- c(ggamma, stats::sigma(model)^2)
  # (i, j) of the k-th entry of an N x N upper triangle, row by row
  upper_entry <- function(k, N) {
    aa <- cumsum(N:1)
    aa_low <- c(0, aa[-length(aa)])
    i <- which(aa_low < k & k <= aa)
    j <- k - N * i + N - i * (3 - i)/2 + i
    c(i, j)
  }
  Zt <- lme4::getME(model, "Zt")
  G <- list()
  for (ss in seq_len(n_rt)) {
    ZZ <- Zt[(Gp[ss] + 1):Gp[ss + 1], , drop = FALSE]
    Ig <- Matrix::Diagonal(n_lev[ss])
    for (rr in seq_len(n_parm[ss])) {
      ij <- unique(upper_entry(rr, n_comp[ss]))
      EE <- if (length(ij) == 1) {
        Matrix::sparseMatrix(ij, ij, x = 1, dims = rep(n_comp[ss], 2))
      } else {
        Matrix::sparseMatrix(ij, ij[2:1], x = 1, dims = rep(n_comp[ss], 2))
      }
      G <- c(G, list(Matrix::t(ZZ) %*% (Ig %x% EE) %*% ZZ))
    }
  }
  n_obs <- nrow(lme4::getME(model, "X"))
  G <- c(G, list(Matrix::Diagonal(n_obs)))
  Sigma <- ggamma[1] * G[[1]]
  for (ii in seq_along(G)[-1]) Sigma <- Sigma + ggamma[ii] * G[[ii]]
  list(Sigma = Sigma, G = G, n.ggamma = length(ggamma))
}

## Mirror the upper triangle onto the lower one (as Matrix::forceSymmetric()).
pm_upper_sym <- function(A) {
  A <- as.matrix(A)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A
}

## Kenward-Roger adjusted covariance of the fixed effects (pbkrtest::vcovAdj()),
## with attributes "P" and "W" used for the degrees of freedom.
pm_vcov_kr <- function(model) {
  if (!lme4::getME(model, "is_REML")) model <- pm_refit(model, REML = TRUE)
  Phi <- as.matrix(stats::vcov(model))
  SG <- pm_kr_sigma_g(model)
  X <- as.matrix(lme4::getME(model, "X"))
  n_g <- SG$n.ggamma
  Sigma <- as.matrix(SG$Sigma)
  SigmaInv <- chol2inv(chol(pm_upper_sym(Sigma)))
  TT <- SigmaInv %*% X
  HH <- OO <- vector("list", n_g)
  for (ii in seq_len(n_g)) {
    HH[[ii]] <- as.matrix(SG$G[[ii]] %*% SigmaInv)
    OO[[ii]] <- HH[[ii]] %*% X
  }
  index <- function(i, j) if (i <= j) (i - 1) * (n_g - i/2) + j else (j - 1) * (n_g - j/2) + i
  PP <- vector("list", n_g)
  QQ <- list()
  for (rr in seq_len(n_g)) {
    Ot <- t(OO[[rr]])
    PP[[rr]] <- pm_upper_sym(-1 * Ot %*% TT)
    for (ss in rr:n_g) QQ <- c(QQ, list(Ot %*% SigmaInv %*% OO[[ss]]))
  }
  Ktrace <- matrix(NA, n_g, n_g)
  for (rr in seq_len(n_g)) {
    Ht <- t(HH[[rr]])
    for (ss in rr:n_g) Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(Ht * HH[[ss]])
  }
  IE2 <- matrix(NA, n_g, n_g)
  for (ii in seq_len(n_g)) {
    PhiP <- Phi %*% PP[[ii]]
    for (jj in ii:n_g) {
      IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] - 2 * sum(Phi * QQ[[index(ii, jj)]]) +
        sum(PhiP * (PP[[jj]] %*% Phi))
    }
  }
  condi <- min(abs(eigen(IE2, only.values = TRUE)$values))
  WW <- pm_upper_sym(if (condi > 1e-10) 2 * solve(IE2) else 2 * MASS::ginv(IE2))
  UU <- matrix(0, ncol(X), ncol(X))
  if (n_g > 1) {
    for (ii in 1:(n_g - 1)) for (jj in (ii + 1):n_g) {
      UU <- UU + WW[ii, jj] * (QQ[[index(ii, jj)]] - PP[[ii]] %*% Phi %*% PP[[jj]])
    }
  }
  UU <- UU + t(UU)
  for (ii in seq_len(n_g)) UU <- UU + WW[ii, ii] * (QQ[[index(ii, ii)]] - PP[[ii]] %*% Phi %*% PP[[ii]])
  PhiA <- Phi + 2 * (Phi %*% UU %*% Phi)
  dimnames(PhiA) <- dimnames(Phi)
  attr(PhiA, "P") <- PP
  attr(PhiA, "W") <- WW
  attr(PhiA, "condi") <- condi
  PhiA
}

## Kenward-Roger F test of L beta = 0 (as pbkrtest::KRmodcomp(model, L)$test["Ftest", ]).
pm_kr_ftest <- function(model, L, PhiA = pm_vcov_kr(model)) {
  if (!is.matrix(L)) L <- matrix(L, nrow = 1)
  Phi <- as.matrix(stats::vcov(model))
  beta <- lme4::fixef(model)
  Theta <- t(L) %*% solve(L %*% Phi %*% t(L), L)
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
  q <- as.numeric(Matrix::rankMatrix(L))
  B <- (1/(2 * q)) * (A1 + 6 * A2)
  g <- ((q + 1) * A1 - (q + 4) * A2)/((q + 2) * A2)
  c1 <- g/(3 * q + 2 * (1 - g))
  c2 <- (q - g)/(3 * q + 2 * (1 - g))
  c3 <- (q + 2 - g)/(3 * q + 2 * (1 - g))
  V0 <- 1 + c1 * B
  V1 <- 1 - c2 * B
  V2 <- 1 - c3 * B
  V0 <- ifelse(abs(V0) < 1e-10, 0, V0)
  x <- 1 - A2/q
  ratio <- if (abs(x) < 1e-14 && abs(V1) < 1e-14) 1 else x/V1
  rho <- 1/q * ratio^2 * V0/V2
  ddf <- 4 + (q + 2)/(q * rho - 1)
  F_scaling <- if (abs(ddf - 2) < 0.01) 1 else ddf * (1 - A2/q)/(ddf - 2)
  b <- L %*% beta
  Wald <- as.numeric(t(b) %*% solve(L %*% as.matrix(PhiA) %*% t(L), b))
  Fstat <- F_scaling * Wald/q
  c(stat = Fstat, ndf = q, ddf = ddf, F.scaling = F_scaling,
    p.value = stats::pf(Fstat, q, ddf, lower.tail = FALSE))
}
