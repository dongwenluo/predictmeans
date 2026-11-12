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

pm_satt_parts_nlme <- function(model) {
  if (!identical(model$method, "REML")) {
    stop("Satterthwaite df for gls and lme models need a REML fit.", call. = FALSE)
  }
  dat <- pm_get_data(model)
  b <- if (inherits(model, "lme")) nlme::fixef(model) else stats::coef(model)
  X <- stats::model.matrix(stats::formula(model), dat)[, names(b), drop = FALSE]
  y <- stats::model.response(stats::model.frame(stats::formula(model), dat))
  Sig <- pm_nlme_sigma(model)
  th0 <- c(stats::coef(model$modelStruct), lSigma = log(model$sigma))
  k <- length(th0)
  vb <- function(p) {
    Si <- solve(Sig(p[-k]))
    exp(2 * p[k]) * solve(crossprod(X, Si %*% X))
  }
  reml2 <- function(p) {                       # -2 REML log-likelihood (up to a constant)
    R <- chol(exp(2 * p[k]) * Sig(p[-k]))
    Si <- chol2inv(R)
    XSX <- crossprod(X, Si %*% X)
    r <- y - X %*% solve(XSX, crossprod(X, Si %*% y))
    2 * sum(log(diag(R))) + as.numeric(determinant(XSX)$modulus) + sum(r * (Si %*% r))
  }
  V0 <- vb(th0)
  Vm <- as.matrix(stats::vcov(model))
  if (max(abs(V0 - Vm)) > 1e-6 * max(abs(Vm))) {
    stop("Satterthwaite df are not available for this ", class(model)[1], " model (its covariance structure ",
         "could not be rebuilt).", call. = FALSE)
  }
  H <- numDeriv::hessian(reml2, th0)
  eg <- eigen(H, symmetric = TRUE)
  pos <- eg$values > 1e-8
  Hinv <- eg$vectors[, pos, drop = FALSE] %*% diag(1 / eg$values[pos], nrow = sum(pos)) %*%
    t(eg$vectors[, pos, drop = FALSE])
  J <- numDeriv::jacobian(function(p) as.vector(vb(p)), th0)
  p <- ncol(X)
  structure(list(vcov_beta = Vm, vcov_varpar = 2 * Hinv,
                 Jac_list = lapply(seq_len(k), function(i) matrix(J[, i], p, p))),
            class = "pm_satt")
}
