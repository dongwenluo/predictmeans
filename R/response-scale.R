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

## h and its derivatives (vectorised): the family's mu.eta when there is one,
## else central differences.
pm_back_funs <- function(back, trans, model, link_scale) {
  fd <- function(f) function(x) {
    e <- 1e-4 * pmax(1, abs(x))
    (f(x + e) - f(x - e)) / (2 * e)
  }
  d1 <- NULL
  if (link_scale) {
    fam <- tryCatch(stats::family(model), error = function(e) NULL)
    if (is.list(fam) && is.function(fam$mu.eta)) d1 <- function(x) fam$mu.eta(x)
  }
  if (is.null(d1)) d1 <- if (identical(trans, exp)) exp else fd(back)
  d2 <- fd(d1)
  d3 <- fd(d2)
  list(h = back, d1 = d1, d2 = d2, d3 = d3)
}

## Default sigma^2 for the bias adjustment: the total variance of the random
## intercepts (GLMMs), plus the residual variance for a transformed response
## of a linear (mixed) model.
pm_bias_sigma2 <- function(model, link_scale) {
  s2 <- 0
  if (inherits(model, "merMod")) {
    vc <- as.data.frame(lme4::VarCorr(model))
    re <- vc[vc$grp != "Residual", , drop = FALSE]
    if (any(!is.na(re$var2) | re$var1 != "(Intercept)")) {
      stop("bias_adjust with random slopes or correlated random effects needs 'sigma' (the SD of the random effects ",
           "to average over).", call. = FALSE)
    }
    s2 <- sum(re$vcov)
  } else if (inherits(model, "glmmTMB")) {
    vc <- glmmTMB::VarCorr(model)$cond
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
