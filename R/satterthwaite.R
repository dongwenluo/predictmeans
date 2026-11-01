## Satterthwaite denominator df for lmer models, for many contrasts at once.
##
## Same computation as lmerTest (Kuznetsova et al., GPL >= 2): for a contrast
## l, df = 2 (l V l')^2 / (g' A g), where V = cov(beta), g = gradient of
## l V l' with respect to the variance parameters (theta, sigma) and A = their
## covariance (2 x inverse Hessian of the deviance). The ingredients are
## computed once per model by pm_satt_parts() and reused for every contrast.

pm_satt_parts <- function(model) {
  cl <- stats::getCall(model)
  cl[[1]] <- quote(lme4::lmer)
  cl$devFunOnly <- TRUE
  if (is.null(cl$control)) cl$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
  devfun <- eval(cl, pm_fit_env(model))
  is_reml <- lme4::getME(model, "is_REML")
  varpar <- unname(c(lme4::getME(model, "theta"), stats::sigma(model)))
  h <- numDeriv::hessian(func = devfun_vp, x = varpar, devfun = devfun, reml = is_reml)
  eig <- eigen(h, symmetric = TRUE)
  tol <- 1e-8
  pos <- eig$values > tol
  q <- sum(pos)
  h_inv <- eig$vectors[, pos, drop = FALSE] %*% diag(1/eig$values[pos], nrow = q) %*%
    t(eig$vectors[, pos, drop = FALSE])
  p <- length(lme4::fixef(model))
  Jac <- numDeriv::jacobian(func = get_covbeta, x = varpar, devfun = devfun)
  devfun(lme4::getME(model, "theta"))         # leave the deviance function at the optimum
  structure(list(vcov_beta = as.matrix(stats::vcov(model)),
                 vcov_varpar = 2 * h_inv,
                 Jac_list = lapply(seq_len(ncol(Jac)), function(i) array(Jac[, i], dim = c(p, p)))),
            class = "pm_satt")
}

pm_satt_ddf <- function(L, parts) {
  if (!is.matrix(L)) L <- matrix(L, nrow = 1)
  var_con <- rowSums((L %*% parts$vcov_beta) * L)
  grad <- vapply(parts$Jac_list, function(J) rowSums((L %*% J) * L), numeric(nrow(L)))
  grad <- matrix(grad, nrow = nrow(L))
  denom <- rowSums((grad %*% parts$vcov_varpar) * grad)
  stats::setNames(2 * var_con^2/denom, rownames(L))
}
