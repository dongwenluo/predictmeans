## The working (linearised) model of a GLM or GLMM, for permutation tests on
## the link scale (dev/permutation-plan.md, section 4.1).
##
## At the fit, with eta the linear predictor without the offset and mu the
## fitted (conditional) mean:
##   working response  z = eta + (y - mu) / (dmu/deta)
##   working variance  Var(z_i | b) = phi V(mu_i) / (a_i (dmu/deta)^2)
## (a: prior weights, number of trials for binomial; V and phi: the family's
## variance function and dispersion), so that Var(z) = A A' + diag(resvar)
## with A = Z L (L L' = G, the random-effect covariance on the absolute
## scale). This is the V of R2_glmm(). Family variances, checked against
## simulate() (tests/testthat/test-working.R):
##   poisson mu; binomial mu(1 - mu) (y a proportion, a = trials);
##   Gamma mu^2 sigma^2; gaussian sigma^2; nbinom2 mu(1 + mu/theta);
##   nbinom1 mu(1 + phi); beta mu(1 - mu)/(1 + phi)
##   (theta, phi and sigma as sigma() of the glmmTMB fit).
##
## pm_working(model) returns a list:
##   z, resvar   working response and variances (length n)
##   xb          X beta-hat (no offset, no random effects)
##   X, A        fixed-effect design (non-aliased columns) and Z L (n x q)
##   eta, mu, y  linear predictor without offset, mean, response (proportion
##               for binomial); size: trials (binomial) or 1
##   family, link, binary (TRUE for 0/1 responses with one trial each)
pm_working <- function(model) {
  if (inherits(model, "glmmTMB")) return(pm_working_tmb(model))
  if (inherits(model, "glmerMod")) return(pm_working_glmer(model))
  if (inherits(model, "glm")) return(pm_working_glm(model))
  stop("pm_working() handles glm, glmer and glmmTMB models.", call. = FALSE)
}

pm_working_out <- function(eta, mu, mu_eta, y, resvar, X, beta, A, size, family, link) {
  binary <- family == "binomial" && all(size == 1) && all(y %in% c(0, 1))
  list(z = as.vector(eta + (y - mu) / mu_eta), resvar = as.vector(resvar),
       xb = as.vector(X %*% beta), X = X, A = A, eta = as.vector(eta), mu = as.vector(mu),
       y = as.vector(y), size = size, family = family, link = link, binary = binary)
}

pm_working_glm <- function(model) {
  fam <- stats::family(model)
  eta_full <- model$linear.predictors
  off <- stats::model.offset(stats::model.frame(model))
  eta <- if (is.null(off)) eta_full else eta_full - off
  mu <- model$fitted.values
  mu_eta <- fam$mu.eta(eta_full)
  a <- model$prior.weights
  phi <- summary(model)$dispersion
  b <- stats::coef(model)
  X <- stats::model.matrix(model)[, !is.na(b), drop = FALSE]
  resvar <- phi * fam$variance(mu) / (a * mu_eta^2)
  size <- if (fam$family == "binomial") a else rep(1, length(mu))
  pm_working_out(eta, mu, mu_eta, model$y, resvar, X, b[!is.na(b)], matrix(0, length(mu), 0),
                 size, fam$family, fam$link)
}

pm_working_glmer <- function(model) {
  resp <- model@resp
  fam <- stats::family(model)
  eta_full <- resp$eta
  eta <- eta_full - resp$offset
  mu_eta <- fam$mu.eta(eta_full)
  # lme4: sqrtWrkWt()^2 = a (dmu/deta)^2 / V(mu); sigma is 1 for families
  # without a dispersion parameter, and G = sigma^2 Lambda Lambda'
  s <- stats::sigma(model)
  resvar <- s^2 / resp$sqrtWrkWt()^2
  A <- s * (lme4::getME(model, "Z") %*% lme4::getME(model, "Lambda"))
  size <- if (fam$family == "binomial") resp$weights else rep(1, length(eta))
  pm_working_out(eta, resp$mu, mu_eta, resp$y, resvar, as.matrix(lme4::getME(model, "X")), lme4::fixef(model),
                 A, size, fam$family, fam$link)
}

pm_working_tmb <- function(model) {
  info <- model$modelInfo
  if (!identical(deparse(info$allForm$ziformula), "~0") || !identical(deparse(info$allForm$dispformula), "~1")) {
    stop("Permutation tests on the working scale need a glmmTMB model without zero inflation and with the ",
         "default dispersion (dispformula = ~1).", call. = FALSE)
  }
  fam <- info$family
  fname <- fam$family
  ok <- c("poisson", "binomial", "nbinom1", "nbinom2", "Gamma", "gaussian", "beta")
  if (!fname %in% ok) {
    stop("Permutation tests on the working scale handle glmmTMB families ", paste(ok, collapse = ", "),
         "; not ", fname, ".", call. = FALSE)
  }
  dat <- model$obj$env$data
  X <- as.matrix(lme4::getME(model, "X"))
  beta <- glmmTMB::fixef(model)$cond[colnames(X)]
  Z <- lme4::getME(model, "Z")
  b <- lme4::getME(model, "b")[seq_len(ncol(Z))]
  eta <- as.vector(X %*% beta + Z %*% b)
  off <- if (length(dat$offset)) as.vector(dat$offset) else 0
  eta_full <- eta + off
  mu <- fam$linkinv(eta_full)
  mu_eta <- fam$mu.eta(eta_full)
  n <- length(eta)
  w <- if (length(dat$weights)) as.vector(dat$weights) else rep(1, n)
  size <- if (fname == "binomial" && length(dat$size)) as.vector(dat$size) else rep(1, n)
  y <- as.vector(dat$yobs)
  if (fname == "binomial") y <- y / size
  s <- stats::sigma(model)
  v <- switch(fname,
              poisson = mu, binomial = mu * (1 - mu), Gamma = mu^2 * s^2, gaussian = rep(s^2, n),
              nbinom2 = mu * (1 + mu / s), nbinom1 = mu * (1 + s), beta = mu * (1 - mu) / (1 + s))
  resvar <- v / (w * size * mu_eta^2)
  A <- if (ncol(Z)) Z %*% (pm_tmb_lambda(model) * s) else matrix(0, n, 0)
  pm_working_out(eta, mu, mu_eta, y, resvar, X, beta, A, size, fname, fam$link)
}

## A model's formula for the working linear mixed model: the response
## replaced by .pm_z and offset terms removed (the working response has no
## offset).
pm_working_formula <- function(f) {
  drop_off <- function(e) {
    if (is.call(e) && identical(e[[1]], as.name("+"))) {
      parts <- Filter(Negate(is.null), lapply(as.list(e)[-1], drop_off))
      if (!length(parts)) return(NULL)
      return(Reduce(function(a, b) call("+", a, b), parts))
    }
    if (is.call(e) && identical(e[[1]], as.name("offset"))) return(NULL)
    e
  }
  rhs <- drop_off(f[[3]])
  if (is.null(rhs)) rhs <- 1
  out <- stats::as.formula(call("~", as.name(".pm_z"), rhs), env = environment(f))
  out
}

## The rows of the data used to fit a glmer or glmmTMB model.
pm_model_rows <- function(model) {
  if (inherits(model, "glmmTMB")) {
    d <- eval(stats::getCall(model)$data, environment(stats::formula(model)))
    fr <- model$frame
  } else {
    d <- pm_mer_data(model)
    fr <- stats::model.frame(model)
  }
  if (!is.data.frame(d)) stop("Can't find the data used to fit the model; please refit it with a 'data' argument.",
                              call. = FALSE)
  d <- as.data.frame(d)
  if (nrow(d) != nrow(fr)) d <- d[rownames(fr), , drop = FALSE]
  d
}
