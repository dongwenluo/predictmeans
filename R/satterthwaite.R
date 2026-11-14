## Satterthwaite denominator df for lmer models, for many contrasts at once.
##
## Same computation as lmerTest (Kuznetsova et al., GPL >= 2): for a contrast
## l, df = 2 (l V l')^2 / (g' A g), where V = cov(beta), g = gradient of
## l V l' with respect to the variance parameters (theta, sigma) and A = their
## covariance (2 x inverse Hessian of the deviance). The ingredients are
## computed once per model by pm_satt_parts() and reused for every contrast.

pm_satt_parts <- function(model) {
  vp <- pm_lmer_varpar(model)                 # (theta, sigma) and their covariance (R/varcomp.R)
  p <- length(lme4::fixef(model))
  Jac <- numDeriv::jacobian(func = get_covbeta, x = vp$varpar, devfun = vp$devfun)
  vp$devfun(lme4::getME(model, "theta"))      # leave the deviance function at the optimum
  structure(list(vcov_beta = as.matrix(stats::vcov(model)),
                 vcov_varpar = vp$vcov_varpar,
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
