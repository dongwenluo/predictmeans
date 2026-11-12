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
## The data covariance is dense (n x n), so this suits up to a few thousand
## observations.

## Covariance of the data relative to sigma^2, as a function of the
## unconstrained modelStruct parameters, rows in the order of the data used.
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
    Zs <- lapply(names(re_f), function(lv) {
      # the groups of a nested level: its factor within those of the levels above it
      gq <- interaction(grps[seq_len(match(lv, names(grps)))], drop = TRUE)
      Zq <- stats::model.matrix(re_f[[lv]], dat)
      lapply(levels(gq), function(l) Zq * (gq == l))
    })
    names(Zs) <- names(re_f)
  } else {
    grps <- if (is.null(gf)) NULL else nlme::getGroups(dat, gf)
    ord <- if (is.null(grps)) seq_len(n) else order(grps)
  }
  back <- order(ord)                           # nlme's internal (sorted) order -> data order
  function(theta) {
    if (length(theta)) ms <- nlme::`coef<-`(ms, value = theta)
    S <- matrix(0, n, n)
    if (!is.null(Zs)) {
      D <- nlme::pdMatrix(ms$reStruct)
      for (lv in names(Zs)) for (Z in Zs[[lv]]) S <- S + Z %*% D[[lv]] %*% t(Z)
    }
    R <- diag(n)
    if (!is.null(ms$corStruct)) {
      cm <- nlme::corMatrix(ms$corStruct)
      R <- if (is.list(cm)) as.matrix(Matrix::bdiag(cm)) else cm
    }
    if (!is.null(ms$varStruct)) {
      w <- nlme::varWeights(ms$varStruct)
      R <- R / outer(w, w)
    }
    S + R[back, back]
  }
}

## The REML pieces of a gls or lme fit as functions of theta = (unconstrained
## modelStruct coefficients, log sigma): the covariance of the data (vb: of
## the coefficients), -2 REML log-likelihood, checked against vcov().
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
  vb <- function(p) solve(crossprod(X, solve(Sig(p), X)))
  reml2 <- function(p) {                       # -2 REML log-likelihood (up to a constant)
    R <- chol(Sig(p))
    Si <- chol2inv(R)
    XSX <- crossprod(X, Si %*% X)
    r <- y - X %*% solve(XSX, crossprod(X, Si %*% y))
    2 * sum(log(diag(R))) + as.numeric(determinant(XSX)$modulus) + sum(r * (Si %*% r))
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
## as AR(1)). W: expected information, or observed (SAS; information =
## "observed"). Checked: lme = lmer's KR (pbkrtest); unstructured repeated
## measures (gls) = SAS's KR F tests.
pm_vcov_kr_nlme <- function(model, information = c("expected", "observed")) {
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
  Si <- solve(S0)
  SX <- Si %*% X
  Phi <- R$Vm
  PP <- lapply(G, function(Gi) -crossprod(SX, Gi %*% SX))
  W <- if (information == "observed") pm_nlme_W_observed(R) else {
    Pm <- Si - SX %*% Phi %*% t(SX)
    PG <- lapply(G, function(Gi) Pm %*% Gi)
    IE <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in i:k) IE[i, j] <- IE[j, i] <- 0.5 * sum(PG[[i]] * t(PG[[j]]))
    solve(IE)
  }
  SG <- lapply(G, function(Gi) Si %*% Gi)
  U <- matrix(0, ncol(X), ncol(X))
  for (i in seq_len(k)) for (j in seq_len(k)) {
    Qij <- crossprod(SX, G[[i]] %*% (SG[[j]] %*% SX))
    U <- U + W[i, j] * (Qij - PP[[i]] %*% Phi %*% PP[[j]])
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
