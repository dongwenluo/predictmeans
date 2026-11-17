## Satterthwaite denominator df for gls and lme models (nlme), with any
## correlation (corStruct) and variance (varFunc) structure: the same
## ingredients as pm_satt_parts() for lmer, so pm_satt_ddf() and
## pm_satt_Fddf() work unchanged.
##
## The covariance of the data is rebuilt as a function of the model's
## variance parameters theta (coef(modelStruct), unconstrained, and
## log(sigma)): random effects Z D Z' from the reStruct, times the
## correlation matrices and variance weights, in the data order (nlme sorts
## the data by groups internally). Then
##   A = 2 H^-1 (H the Hessian of -2 REML log-likelihood: the observed
##       information, as SAS uses; from its formula, pm_nlme_info(), with
##       central differences of the covariance's pieces only),
##   G_k = d vcov(beta) / d theta_k = Phi X' S^-1 (dS/dtheta_k) S^-1 X Phi,
## and df = 2 (l V l')^2 / (g' A g) with g_k = l G_k l'. The rebuilt
## covariance is checked against vcov() of the fit. Checked against SAS
## PROC MIXED DDFM=SATTERTHWAITE (dev/sas): unstructured repeated measures
## (gls; 10, 10.1, 10.2) and the split plot with a plot missing (lme; as lmer).
## The data covariance is a sparse matrix, block diagonal by the outermost
## groups (subjects), so the work grows about linearly with the number of
## groups.

## Covariance of the data relative to sigma^2, as a function of the
## unconstrained modelStruct parameters, rows in the order of the data used:
## a sparse symmetric matrix (random effects Z (I x D) Z' by level, plus the
## correlation blocks scaled by the variance weights).
pm_nlme_sigma <- function(model) {
  dat <- pm_get_data(model)
  ms <- model$modelStruct
  n <- nrow(dat)
  gf <- tryCatch(nlme::getGroupsFormula(model), error = function(e) NULL)
  Zs <- NULL
  if (inherits(model, "lme")) {
    re_f <- stats::formula(ms$reStruct)
    grps <- nlme::getGroups(dat, gf)
    if (!is.data.frame(grps)) {
      grps <- data.frame(grps)
      names(grps) <- names(re_f)[1]
    }
    ord <- do.call(order, grps)
    # one sparse design per level: the columns of group g are (g - 1) q + 1:q
    Zs <- lapply(names(re_f), function(lv) {
      # the groups of a nested level: its factor within those of the levels above it
      gq <- interaction(grps[seq_len(match(lv, names(grps)))], drop = TRUE)
      Zq <- stats::model.matrix(re_f[[lv]], dat)
      q <- ncol(Zq)
      gi <- as.integer(gq)
      list(Z = Matrix::sparseMatrix(i = rep(seq_len(n), q), j = (rep(gi, q) - 1L) * q + rep(seq_len(q), each = n),
                                    x = as.vector(Zq), dims = c(n, nlevels(gq) * q)),
           I = Matrix::Diagonal(nlevels(gq)))
    })
    names(Zs) <- names(re_f)
  } else {
    grps <- if (is.null(gf)) NULL else nlme::getGroups(dat, gf)
    ord <- if (is.null(grps)) seq_len(n) else order(grps)
  }
  # residual part: the correlation blocks (nlme's sorted order) scaled by the
  # variance weights, put straight into data order: sorted position a is data
  # row ord[a]; the indices of the block entries are found once
  blk <- if (is.null(ms$corStruct)) NULL else {
    cm0 <- nlme::corMatrix(ms$corStruct)
    if (!is.list(cm0)) cm0 <- list(cm0)
    sz <- vapply(cm0, nrow, integer(1))
    off <- cumsum(c(0L, sz[-length(sz)]))
    ia <- unlist(lapply(seq_along(sz), function(b) rep(seq_len(sz[b]), sz[b]) + off[b]))
    jb <- unlist(lapply(seq_along(sz), function(b) rep(seq_len(sz[b]), each = sz[b]) + off[b]))
    list(ia = ia, jb = jb, up = ord[ia] <= ord[jb])
  }
  # the pieces at theta: variance weights w (sorted order), correlations r (the
  # entries of the blocks, as blk), relative random-effect covariances D
  pieces <- function(theta, what = c("w", "r", "D")) {
    if (length(theta)) ms <- nlme::`coef<-`(ms, value = theta)
    list(w = if ("w" %in% what) (if (is.null(ms$varStruct)) rep(1, n) else 1 / nlme::varWeights(ms$varStruct)),
         r = if ("r" %in% what && !is.null(blk)) {
           cm <- nlme::corMatrix(ms$corStruct)
           if (is.list(cm)) unlist(cm, use.names = FALSE) else as.vector(cm)
         },
         D = if ("D" %in% what && !is.null(Zs)) nlme::pdMatrix(ms$reStruct))
  }
  # the part each parameter belongs to (coef(modelStruct) is in the order of
  # the components, and of the levels within the reStruct)
  comp <- unlist(lapply(names(ms), function(nm) {
    if (nm == "reStruct") {
      unlist(lapply(names(ms$reStruct), function(lv) rep(paste0("re:", lv), length(stats::coef(ms$reStruct[[lv]])))))
    } else rep(nm, length(stats::coef(ms[[nm]])))
  }))
  if (length(comp) != length(stats::coef(ms))) comp <- rep("other", length(stats::coef(ms)))
  # attribute Zs: the random-effect designs by level (R2_model()); pieces and
  # layout: for the second derivatives of the covariance (pm_nlme_T())
  structure(function(theta) {
    pc <- pieces(theta)
    w <- pc$w
    S <- if (is.null(blk)) {
      Matrix::sparseMatrix(i = ord, j = ord, x = w^2, dims = c(n, n), symmetric = TRUE)
    } else {
      x <- pc$r * w[blk$ia] * w[blk$jb]
      u <- blk$up
      Matrix::sparseMatrix(i = ord[blk$ia][u], j = ord[blk$jb][u], x = x[u], dims = c(n, n), symmetric = TRUE)
    }
    if (!is.null(Zs)) {
      for (lv in names(Zs)) {
        Z <- Zs[[lv]]$Z
        S <- S + Z %*% Matrix::kronecker(Zs[[lv]]$I, Matrix::Matrix(pc$D[[lv]])) %*% Matrix::t(Z)
      }
    }
    Matrix::forceSymmetric(methods::as(S, "CsparseMatrix"), uplo = "U")
  }, Zs = Zs, pieces = pieces, layout = list(n = n, ord = ord, blk = blk, comp = comp))
}

## The likelihood pieces of a gls or lme fit as functions of theta =
## (unconstrained modelStruct coefficients, log sigma): the covariance of the
## data (sparse; vb: of the coefficients), -2 log-likelihood (reml2: REML for
## REML fits; ML for ML fits, as lmerTest's Satterthwaite for ML lmer fits),
## checked against vcov(). Kenward-Roger needs REML (reml = TRUE).
pm_nlme_reml <- function(model, reml = FALSE) {
  ml <- identical(model$method, "ML")
  if (reml && ml) {
    stop("Kenward-Roger df for gls and lme models need a REML fit.", call. = FALSE)
  }
  dat <- pm_get_data(model)
  b <- if (inherits(model, "lme")) nlme::fixef(model) else stats::coef(model)
  X <- stats::model.matrix(stats::formula(model), dat)[, names(b), drop = FALSE]
  y <- stats::model.response(stats::model.frame(stats::formula(model), dat))
  Srel <- pm_nlme_sigma(model)
  th0 <- c(stats::coef(model$modelStruct), lSigma = log(model$sigma))
  k <- length(th0)
  Sig <- function(p) exp(2 * p[k]) * Srel(p[-k])
  # nlme's vcov() of an ML gls fit is scaled by N / (N - p) (a constant, which
  # does not change Satterthwaite's df)
  vscale <- if (ml && inherits(model, "gls")) model$dims$N / (model$dims$N - model$dims$p) else 1
  vb <- function(p) vscale * solve(crossprod(X, as.matrix(Matrix::solve(Sig(p), X))))
  reml2 <- function(p) {                       # -2 REML log-likelihood (up to a constant)
    S <- Sig(p)
    C <- Matrix::Cholesky(S, LDL = FALSE)
    SiX <- as.matrix(Matrix::solve(C, X, system = "A"))
    Siy <- as.vector(Matrix::solve(C, y, system = "A"))
    XSX <- crossprod(X, SiX)
    bh <- solve(XSX, crossprod(X, Siy))
    r <- y - X %*% bh
    as.numeric(Matrix::determinant(C, sqrt = FALSE)$modulus) + (if (ml) 0 else as.numeric(determinant(XSX)$modulus)) +
      sum(r * (Siy - SiX %*% bh))
  }
  Vm <- as.matrix(stats::vcov(model))
  if (max(abs(vb(th0) - Vm)) > 1e-6 * max(abs(Vm))) {
    stop("Satterthwaite and Kenward-Roger df are not available for this ", class(model)[1], " model (its ",
         "covariance structure could not be rebuilt).", call. = FALSE)
  }
  list(X = X, y = y, th0 = th0, k = k, Sig = Sig, Srel = Srel, vb = vb, reml2 = reml2, Vm = Vm, ml = ml,
       vscale = vscale)
}

## The derivative pieces of a gls or lme fit at its estimates (R from
## pm_nlme_reml()): G_i = d Sigma / d theta_i (central differences of the
## covariance for the modelStruct parameters, exact for log sigma), Sigma^-1
## X, Sigma^-1 G_i, Phi = (X' Sigma^-1 X)^-1, A_i = X' S^-1 G_i S^-1 X,
## Q_ij = X' S^-1 G_i S^-1 G_j S^-1 X, P y and G_i P y with P G_i P y. All n x n
## products stay sparse (block diagonal); the rest is p x p.
pm_nlme_parts <- function(R) {
  X <- R$X
  th <- R$th0
  k <- R$k
  S0 <- R$Sig(th)
  G <- lapply(seq_len(k), function(i) {
    if (i == k) return(2 * S0)                 # d Sigma / d log sigma
    h <- 1e-4 * max(1, abs(th[i]))
    e <- replace(numeric(k), i, h)
    (R$Sig(th + e) - R$Sig(th - e)) / (2 * h)
  })
  C <- Matrix::Cholesky(S0, LDL = FALSE)
  SX <- as.matrix(Matrix::solve(C, X, system = "A"))
  SG <- lapply(G, function(Gi) Matrix::solve(C, methods::as(Gi, "CsparseMatrix"), system = "A"))   # Sigma^-1 G_i
  Phi <- R$Vm / R$vscale
  A <- lapply(G, function(Gi) as.matrix(crossprod(SX, Gi %*% SX)))
  Q <- lapply(seq_len(k), function(i) lapply(seq_len(k), function(j)
    as.matrix(crossprod(SX, G[[i]] %*% (SG[[j]] %*% SX)))))
  Pv <- function(v) {                          # P v = S^-1 v - S^-1 X Phi X' S^-1 v
    sv <- as.matrix(Matrix::solve(C, v, system = "A"))
    sv - SX %*% (Phi %*% crossprod(SX, v))
  }
  Py <- Pv(R$y)
  a <- lapply(G, function(Gi) as.matrix(Gi %*% Py))
  list(S0 = S0, G = G, C = C, SX = SX, SG = SG, Phi = Phi, A = A, Q = Q, Py = as.vector(Py), a = a,
       Pa = lapply(a, Pv))
}

## T_ij = sum_ab c_ab d2 Sigma_ab / d theta_i d theta_j, c = P - P y y' P
## (REML; Sigma^-1 in place of P for ML): the part of the Hessian of the -2
## log-likelihood that has the second derivatives of the covariance. The
## covariance is a sum of the random-effect levels, Z (I x D) Z', and the
## residual part, sigma^2 w_a w_b r_ab (weights from the varStruct,
## correlations from the corStruct), so parameters of different parts give 0
## and the residual part splits by the product rule; central differences of
## those pieces only (pdMatrix, corMatrix, varWeights), never of the
## likelihood. log sigma: d2 Sigma = 4 Sigma, and 2 G_i with the others.
pm_nlme_T <- function(R, pt, reml) {
  k <- R$k
  ctr <- function(M) {                         # sum_ab c_ab M_ab for a symmetric n x n M
    v <- sum(Matrix::diag(Matrix::solve(pt$C, methods::as(M, "CsparseMatrix"), system = "A"))) -
      sum(pt$Py * as.vector(M %*% pt$Py))
    if (reml) v - sum((pt$SX %*% pt$Phi) * as.matrix(M %*% pt$SX)) else v
  }
  Tm <- matrix(0, k, k)
  for (i in seq_len(k - 1)) Tm[i, k] <- Tm[k, i] <- 2 * ctr(pt$G[[i]])
  Tm[k, k] <- 4 * ctr(pt$S0)
  m <- k - 1
  if (!m) return(Tm)
  lay <- attr(R$Srel, "layout")
  pieces <- attr(R$Srel, "pieces")
  comp <- lay$comp
  th <- R$th0[-k]
  if (any(comp == "other")) {                  # parts not known here: second differences of the likelihood
    H <- numDeriv::hessian(R$reml2, R$th0)
    return(H - 2 * pm_nlme_info(R, pt, "expected", reml) + 4 * pm_nlme_info(R, pt, "average", reml))
  }
  n <- lay$n
  Sinv <- Matrix::solve(pt$C, Matrix::Diagonal(n), system = "A")
  XP <- pt$SX %*% pt$Phi
  cvec <- function(i, j) {                     # c at the entries (i, j)
    v <- Sinv[cbind(i, j)] - pt$Py[i] * pt$Py[j]
    if (reml) v - rowSums(XP[i, , drop = FALSE] * pt$SX[j, , drop = FALSE]) else v
  }
  h <- 1e-4 * pmax(1, abs(th))
  at <- function(idx, d) replace(th, idx, th[idx] + d)
  hess2 <- function(f, idx) {                  # central second differences of f over the parameters idx
    p <- length(idx)
    H <- matrix(0, p, p)
    f0 <- f(th)
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
  Ts <- matrix(0, m, m)
  # random-effect levels: sum_g Z_g' c Z_g (q x q) times d2 D
  for (lv in names(attr(R$Srel, "Zs"))) {
    idx <- which(comp == paste0("re:", lv))
    if (!length(idx)) next
    Z <- attr(R$Srel, "Zs")[[lv]]$Z
    ng <- nrow(attr(R$Srel, "Zs")[[lv]]$I)
    q <- ncol(Z) %/% ng
    K <- Matrix::crossprod(Z, Sinv %*% Z)
    U <- as.matrix(Matrix::crossprod(Z, pt$SX))
    u <- as.vector(Matrix::crossprod(Z, pt$Py))
    M <- matrix(0, q, q)
    for (r in seq_len(q)) for (s in seq_len(q)) {
      ir <- (seq_len(ng) - 1L) * q + r
      is <- (seq_len(ng) - 1L) * q + s
      M[r, s] <- sum(K[cbind(ir, is)]) - sum(u[ir] * u[is]) -
        (if (reml) sum((U[ir, , drop = FALSE] %*% pt$Phi) * U[is, , drop = FALSE]) else 0)
    }
    Ts[idx, idx] <- hess2(function(x) sum(M * pieces(x, "D")$D[[lv]]), idx)
  }
  # residual part
  ic <- which(comp == "corStruct")
  iv <- which(comp == "varStruct")
  if (length(ic) || length(iv)) {
    blk <- lay$blk
    ia <- if (is.null(blk)) seq_len(n) else blk$ia      # positions in nlme's sorted order
    jb <- if (is.null(blk)) seq_len(n) else blk$jb
    cp <- cvec(lay$ord[ia], lay$ord[jb])
    p0 <- pieces(th, c("w", "r"))
    r0 <- if (is.null(p0$r)) 1 else p0$r
    ww <- function(w) w[ia] * w[jb]
    if (length(ic)) {
      cw <- cp * ww(p0$w)
      Ts[ic, ic] <- hess2(function(x) sum(cw * pieces(x, "r")$r), ic)
    }
    if (length(iv)) {
      cr <- cp * r0
      Ts[iv, iv] <- hess2(function(x) sum(cr * ww(pieces(x, "w")$w)), iv)
    }
    if (length(ic) && length(iv)) {
      dr <- lapply(ic, function(i) (pieces(at(i, h[i]), "r")$r - pieces(at(i, -h[i]), "r")$r) / (2 * h[i]))
      dw <- lapply(iv, function(j) (ww(pieces(at(j, h[j]), "w")$w) - ww(pieces(at(j, -h[j]), "w")$w)) / (2 * h[j]))
      for (a in seq_along(ic)) for (b in seq_along(iv)) {
        Ts[ic[a], iv[b]] <- Ts[iv[b], ic[a]] <- sum(cp * dr[[a]] * dw[[b]])
      }
    }
  }
  Tm[seq_len(m), seq_len(m)] <- exp(2 * R$th0[k]) * Ts
  Tm
}

## Information about the variance parameters theta (pt from
## pm_nlme_parts()): "expected", 1/2 tr(P G_i P G_j); "average" (Genstat),
## 1/2 y'P G_i P G_j P y; "observed" (SAS), half the Hessian of the -2
## log-likelihood, 1/2 T - expected + 2 average (pm_nlme_T(); exact up to
## the differences of the covariance's pieces). ML fits (reml = FALSE): S^-1
## in place of P in the traces.
pm_nlme_info <- function(R, pt, type = c("observed", "expected", "average"), reml = TRUE) {
  type <- match.arg(type)
  k <- R$k
  if (type == "average") {
    AI <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in i:k) AI[i, j] <- AI[j, i] <- 0.5 * sum(pt$a[[i]] * pt$Pa[[j]])
    return(AI)
  }
  if (type == "expected") {
    IE <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in i:k) {
      tr <- sum(pt$SG[[i]] * Matrix::t(pt$SG[[j]]))
      if (reml) tr <- tr - 2 * sum(pt$Phi * t(pt$Q[[i]][[j]])) + sum((pt$Phi %*% pt$A[[i]]) * t(pt$Phi %*% pt$A[[j]]))
      IE[i, j] <- IE[j, i] <- 0.5 * tr
    }
    return(IE)
  }
  0.5 * pm_nlme_T(R, pt, reml) - pm_nlme_info(R, pt, "expected", reml) + 2 * pm_nlme_info(R, pt, "average", reml)
}

## The covariance of the variance parameters from the observed information
## (directions with a non-positive curvature left out).
pm_nlme_W_observed <- function(R, pt = pm_nlme_parts(R)) {
  I <- pm_nlme_info(R, pt, "observed", reml = !R$ml)
  I <- (I + t(I)) / 2
  eg <- eigen(I, symmetric = TRUE)
  pos <- eg$values > 5e-9
  eg$vectors[, pos, drop = FALSE] %*% diag(1 / eg$values[pos], nrow = sum(pos)) %*% t(eg$vectors[, pos, drop = FALSE])
}

## Satterthwaite: d vcov(beta) / d theta_i = Phi A_i Phi (exact, given G_i).
pm_satt_parts_nlme <- function(model) {
  R <- pm_nlme_reml(model)
  pt <- pm_nlme_parts(R)
  structure(list(vcov_beta = R$Vm, vcov_varpar = pm_nlme_W_observed(R, pt),
                 Jac_list = lapply(pt$A, function(a) R$vscale * pt$Phi %*% a %*% pt$Phi)),
            class = "pm_satt")
}

## Kenward-Roger adjusted covariance of the coefficients of a gls or lme
## model, with attributes "P" and "W" as pm_vcov_kr() (so pm_kr_ftest() and
## pm_kr_ddf() work unchanged). Derivatives of the data covariance with
## respect to theta are numerical; the second-derivative term of Kenward and
## Roger (1997) is left out: the rest does not depend on how the covariance
## is parameterised, and in the natural parameters of variance components
## and unstructured covariances (where SAS works) that term is zero, so this
## equals SAS DDFM=KR for those (and DDFM=KR(FIRSTORDER) for structures such
## as AR(1)). W: expected information, observed (SAS; information =
## "observed"), or average (Genstat; y'P G_i P G_j P y / 2, which transforms
## with the parameterisation as the expected information does). Checked: lme = lmer's KR (pbkrtest); unstructured repeated
## measures (gls) = SAS's KR F tests.
pm_vcov_kr_nlme <- function(model, information = c("observed", "expected", "average")) {
  information <- match.arg(information)
  R <- pm_nlme_reml(model, reml = TRUE)
  X <- R$X
  k <- R$k
  pt <- pm_nlme_parts(R)
  Phi <- pt$Phi
  Q <- pt$Q
  PP <- lapply(pt$A, function(a) -a)
  W <- if (information == "observed") pm_nlme_W_observed(R, pt) else solve(pm_nlme_info(R, pt, information))
  U <- matrix(0, ncol(X), ncol(X))
  for (i in seq_len(k)) for (j in seq_len(k)) {
    U <- U + W[i, j] * (Q[[i]][[j]] - PP[[i]] %*% Phi %*% PP[[j]])
  }
  PhiA <- Phi + 2 * Phi %*% U %*% Phi
  PhiA <- (PhiA + t(PhiA)) / 2
  dimnames(PhiA) <- dimnames(Phi)
  attr(PhiA, "P") <- PP
  attr(PhiA, "W") <- W
  attr(PhiA, "information") <- information
  PhiA
}

## Satterthwaite or Kenward-Roger df of the rows of L for a gls or lme model.
pm_nlme_df <- function(L, model, df_type, vcov_adj) {
  L <- matrix(as.numeric(as.matrix(L)), nrow = NROW(L))
  if (identical(df_type, "Kenward-Roger")) pm_kr_ddf(L, stats::vcov(model), vcov_adj) else pm_satt_ddf(L, vcov_adj)
}
