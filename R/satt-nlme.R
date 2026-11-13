## Satterthwaite denominator df for gls and lme models (nlme), with any
## correlation (corStruct) and variance (varFunc) structure: the same
## ingredients as pm_satt_parts() for lmer, so pm_satt_ddf() and
## pm_satt_Fddf() work unchanged.
##
## The covariance of the data is rebuilt as a function of the model's
## variance parameters theta (coef(modelStruct), unconstrained, and
## log(sigma)): random effects Z D Z' from the reStruct, times the
## correlation matrices and variance weights, in the data order (nlme sorts
## the data by groups internally). Then, numerically,
##   A = 2 H^-1 (H the Hessian of -2 REML log-likelihood: the observed
##       information, as SAS uses),
##   G_k = d vcov(beta) / d theta_k,
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
  function(theta) {
    if (length(theta)) ms <- nlme::`coef<-`(ms, value = theta)
    w <- if (is.null(ms$varStruct)) rep(1, n) else 1 / nlme::varWeights(ms$varStruct)
    S <- if (is.null(blk)) {
      Matrix::sparseMatrix(i = ord, j = ord, x = w^2, dims = c(n, n), symmetric = TRUE)
    } else {
      cm <- nlme::corMatrix(ms$corStruct)
      x <- (if (is.list(cm)) unlist(cm, use.names = FALSE) else as.vector(cm)) * w[blk$ia] * w[blk$jb]
      u <- blk$up
      Matrix::sparseMatrix(i = ord[blk$ia][u], j = ord[blk$jb][u], x = x[u], dims = c(n, n), symmetric = TRUE)
    }
    if (!is.null(Zs)) {
      D <- nlme::pdMatrix(ms$reStruct)
      for (lv in names(Zs)) {
        Z <- Zs[[lv]]$Z
        S <- S + Z %*% Matrix::kronecker(Zs[[lv]]$I, Matrix::Matrix(D[[lv]])) %*% Matrix::t(Z)
      }
    }
    Matrix::forceSymmetric(methods::as(S, "CsparseMatrix"), uplo = "U")
  }
}

## The REML pieces of a gls or lme fit as functions of theta = (unconstrained
## modelStruct coefficients, log sigma): the covariance of the data (sparse;
## vb: of the coefficients), -2 REML log-likelihood, checked against vcov().
pm_nlme_reml <- function(model) {
  if (!identical(model$method, "REML")) {
    stop("Satterthwaite and Kenward-Roger df for gls and lme models need a REML fit.", call. = FALSE)
  }
  dat <- pm_get_data(model)
  b <- if (inherits(model, "lme")) nlme::fixef(model) else stats::coef(model)
  X <- stats::model.matrix(stats::formula(model), dat)[, names(b), drop = FALSE]
  y <- stats::model.response(stats::model.frame(stats::formula(model), dat))
  Srel <- pm_nlme_sigma(model)
  th0 <- c(stats::coef(model$modelStruct), lSigma = log(model$sigma))
  k <- length(th0)
  Sig <- function(p) exp(2 * p[k]) * Srel(p[-k])
  vb <- function(p) solve(crossprod(X, as.matrix(Matrix::solve(Sig(p), X))))
  reml2 <- function(p) {                       # -2 REML log-likelihood (up to a constant)
    S <- Sig(p)
    C <- Matrix::Cholesky(S, LDL = FALSE)
    SiX <- as.matrix(Matrix::solve(C, X, system = "A"))
    Siy <- as.vector(Matrix::solve(C, y, system = "A"))
    XSX <- crossprod(X, SiX)
    bh <- solve(XSX, crossprod(X, Siy))
    r <- y - X %*% bh
    as.numeric(Matrix::determinant(C, sqrt = FALSE)$modulus) + as.numeric(determinant(XSX)$modulus) +
      sum(r * (Siy - SiX %*% bh))
  }
  Vm <- as.matrix(stats::vcov(model))
  if (max(abs(vb(th0) - Vm)) > 1e-6 * max(abs(Vm))) {
    stop("Satterthwaite and Kenward-Roger df are not available for this ", class(model)[1], " model (its ",
         "covariance structure could not be rebuilt).", call. = FALSE)
  }
  list(X = X, y = y, th0 = th0, k = k, Sig = Sig, vb = vb, reml2 = reml2, Vm = Vm)
}

## Inverse of half the Hessian of -2 REML log-likelihood: the covariance of
## the variance parameters from the observed information.
pm_nlme_W_observed <- function(R) {
  H <- numDeriv::hessian(R$reml2, R$th0)
  eg <- eigen(H, symmetric = TRUE)
  pos <- eg$values > 1e-8
  2 * eg$vectors[, pos, drop = FALSE] %*% diag(1 / eg$values[pos], nrow = sum(pos)) %*% t(eg$vectors[, pos, drop = FALSE])
}

pm_satt_parts_nlme <- function(model) {
  R <- pm_nlme_reml(model)
  J <- numDeriv::jacobian(function(p) as.vector(R$vb(p)), R$th0)
  p <- ncol(R$X)
  structure(list(vcov_beta = R$Vm, vcov_varpar = pm_nlme_W_observed(R),
                 Jac_list = lapply(seq_len(R$k), function(i) matrix(J[, i], p, p))),
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
  R <- pm_nlme_reml(model)
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
  # all n x n products stay sparse (block diagonal); the rest is p x p
  C <- Matrix::Cholesky(S0, LDL = FALSE)
  SX <- as.matrix(Matrix::solve(C, X, system = "A"))
  SG <- lapply(G, function(Gi) Matrix::solve(C, methods::as(Gi, "CsparseMatrix"), system = "A"))   # Sigma^-1 G_i
  Phi <- R$Vm
  A <- lapply(G, function(Gi) as.matrix(crossprod(SX, Gi %*% SX)))          # X' S^-1 G_i S^-1 X
  PP <- lapply(A, function(a) -a)
  Q <- lapply(seq_len(k), function(i) lapply(seq_len(k), function(j)
    as.matrix(crossprod(SX, G[[i]] %*% (SG[[j]] %*% SX)))))
  W <- if (information == "observed") pm_nlme_W_observed(R) else if (information == "average") {
    # average information: 1/2 (G_i P y)' P (G_j P y), with P v = S^-1 v - S^-1 X Phi X' S^-1 v
    Pv <- function(v) {
      sv <- as.matrix(Matrix::solve(C, v, system = "A"))
      sv - SX %*% (Phi %*% crossprod(SX, v))
    }
    Py <- Pv(R$y)
    a <- lapply(G, function(Gi) as.matrix(Gi %*% Py))
    Pa <- lapply(a, Pv)
    AI <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in i:k) AI[i, j] <- AI[j, i] <- 0.5 * sum(a[[i]] * Pa[[j]])
    solve(AI)
  } else {
    # expected information: 1/2 tr(P G_i P G_j), P = S^-1 - S^-1 X Phi X' S^-1, in pieces
    IE <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in i:k) {
      IE[i, j] <- IE[j, i] <- 0.5 * (sum(SG[[i]] * Matrix::t(SG[[j]])) - 2 * sum(Phi * t(Q[[i]][[j]])) +
                                       sum((Phi %*% A[[i]]) * t(Phi %*% A[[j]])))
    }
    solve(IE)
  }
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
