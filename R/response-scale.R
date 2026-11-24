## Back-transformed means: bias adjustment and averaging on the response
## scale (predictmeans(bias_adjust =, sigma =, average =)).
##
## h is the inverse link (or the user's back-transformation 'trans'). By
## default a mean averages over the other factors on the link scale and is
## back-transformed: h(eta), the response of a unit with random effects zero.
##   bias_adjust = TRUE: h(eta) + sigma^2/2 h''(eta), the second-order
##     approximation of E h(eta + u), u ~ N(0, sigma^2): the mean over the
##     population of random effects (as emmeans' bias.adjust); SE by the delta
##     method with h' + sigma^2/2 h'''.
##   average = "response": the cells of the reference grid are
##     back-transformed first (with the bias adjustment if asked) and then
##     averaged with the same weights (as emmeans' regrid()); SE by the delta
##     method, limits mean +- t SE.

## h and its derivatives (vectorised): h' from pm_back_deriv() (R/tidy.R:
## the family's mu.eta, exp, else central differences), h'' and h''' central
## differences of it (R/deriv.R).
pm_back_funs <- function(back, trans, model, link_scale) {
  d1 <- function(x) pm_back_deriv(trans, x, model, link_scale)
  d2 <- pm_fd_fun(d1, order = 2)
  d3 <- pm_fd_fun(d2, order = 3)
  list(h = back, d1 = d1, d2 = d2, d3 = d3)
}

## The bias adjustment asked for, checked, and its variance: the one place
## for predictmeans() (R/tidy.R) and covariatemeans() (pm_cov_back()).
## 'bias_adjust': TRUE, FALSE (or NULL) or "taylor"; 'sigma': the user's SD,
## else pm_bias_sigma2(), computed only when 'need' (an adjustment will be
## made). Returns list(on, method ("exact" or "taylor"), s2 (0 when off)).
pm_bias_setup <- function(bias_adjust, sigma, model, link_scale, need = TRUE, exclude = NULL) {
  on <- isTRUE(bias_adjust) || identical(bias_adjust, "taylor")
  if (!on && !isFALSE(bias_adjust) && !is.null(bias_adjust)) {
    stop("'bias_adjust' must be TRUE, FALSE or \"taylor\".", call. = FALSE)
  }
  s2 <- if (!on || !need) 0 else if (!is.null(sigma)) sigma^2 else pm_bias_sigma2(model, link_scale, exclude)
  list(on = on, method = if (identical(bias_adjust, "taylor")) "taylor" else "exact", s2 = s2)
}

## Default sigma^2 for the bias adjustment: the total variance of the random
## intercepts (GLMMs), plus the residual variance for a transformed response
## of a linear (mixed) model; aov() with Error(): the sum of the stratum
## variance components (negative estimates as zero; unrounded, from
## pm_vc_aov_est(), R/varcomp.R). Before 2.0.0 predictmeans() used the
## residual variance of the carrier lm here and covariatemeans() the strata
## sum rounded to 4 decimals; they differed by 3-4% on unbalanced data.
pm_bias_sigma2 <- function(model, link_scale, exclude = NULL) {
  if (inherits(model, "pm_aovlist")) {
    return(sum(pmax(suppressMessages(pm_vc_aov_est(attr(model, "aovlist")))$est, 0)))
  }
  s2 <- 0
  # 'exclude': random-effect terms that are not variation between clusters
  # (the smooths of a semireg fit: their random effects are the curve itself)
  if (inherits(model, "merMod")) {
    vc <- as.data.frame(lme4::VarCorr(model))
    re <- vc[vc$grp != "Residual" & !(vc$grp %in% exclude), , drop = FALSE]
    if (any(!is.na(re$var2) | re$var1 != "(Intercept)")) {
      stop("bias_adjust with random slopes or correlated random effects needs 'sigma' (the SD of the random effects ",
           "to average over).", call. = FALSE)
    }
    s2 <- sum(re$vcov)
  } else if (inherits(model, "glmmTMB")) {
    vc <- glmmTMB::VarCorr(model)$cond
    vc <- vc[setdiff(names(vc), exclude)]
    if (any(vapply(vc, function(m) nrow(as.matrix(m)) != 1, logical(1)))) {
      stop("bias_adjust with random slopes needs 'sigma'.", call. = FALSE)
    }
    s2 <- sum(vapply(vc, function(m) as.matrix(m)[1, 1], numeric(1)))
  } else if (inherits(model, "lme")) {
    D <- nlme::pdMatrix(model$modelStruct$reStruct)
    if (any(vapply(D, nrow, integer(1)) != 1)) stop("bias_adjust with random slopes needs 'sigma'.", call. = FALSE)
    s2 <- sum(unlist(D)) * model$sigma^2
  }
  if (!link_scale) s2 <- s2 + stats::sigma(model)^2      # transformed response: the residual variance too
  if (s2 == 0) warning("bias_adjust: the model has no random effects; give 'sigma' for an adjustment.", call. = FALSE)
  s2
}

## Response-scale means (and SEs) of the rows 'rows' of the reference grid
## average, cells back-transformed first. hb, db: the (bias-adjusted)
## back-transformation and its derivative.
pm_response_average <- function(rg, rows, b, V, hb, db, chunk = 20000L) {
  W <- rg$W[rows, , drop = FALSE]
  n <- nrow(rg$grid)
  mu <- numeric(length(rows))
  J <- matrix(0, length(rows), length(b))
  for (s in seq(1L, n, by = chunk)) {
    idx <- s:min(n, s + chunk - 1L)
    X <- as.matrix(rg$basis(rg$grid[idx, , drop = FALSE]))[, names(b), drop = FALSE]
    eta <- drop(X %*% b)
    Wi <- W[, idx, drop = FALSE]
    mu <- mu + as.vector(Wi %*% hb(eta))
    J <- J + as.matrix(Wi %*% (db(eta) * X))
  }
  list(mu = mu, se = sqrt(pmax(rowSums((J %*% as.matrix(V)) * J), 0)))
}
