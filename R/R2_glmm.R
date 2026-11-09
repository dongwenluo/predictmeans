#' An adjusted coefficient of determination (R2) for generalized linear mixed models
#'
#' @description
#' This function produces adjusted R2 for generalized linear mixed models which was crafted following the guidance provided by Professor Hans-Peter Piepho.
#' @usage
#' R2_glmm(model, over_disp=FALSE)
#' @param model An object returned by \code{lmer}, \code{glmer} or \code{glmmTMB}.
#' @param over_disp A logical scalar to indicate whether \code{model} with over-dispersion or not. The default value is FALSE.
#' @return
#' Adjusted R2 in percentage for Total (fixed + random), Fiexd, Random and individual random term.
#' @references
#' Piepho HP. An adjusted coefficient of determination (R2 ) for generalized linear mixed models in one go. Biom J. 2023 Oct;65(7):e2200290. doi: 10.1002/bimj.202200290. Epub 2023 May 1. PMID: 37127864.
#' @examples
#'   library(predictmeans)
#'   Oats$nitro <- factor(Oats$nitro)
#'   (fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats))
#'   R2_glmm(fm)
#'   (gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#'               data = cbpp, family = binomial))
#'   R2_glmm(gm)
#' @name R2_glmm
NULL

R2_glmm <- function(model, over_disp=FALSE){
  cod_v <- pm_r2_piepho(model, over_disp)
  R2 <- paste(round(cod_v, 2), "%", sep="")
  names(R2) <- names(cod_v)
  cat("# Adjusted R2 for Mixed Models\n")
  return(R2)
}

## Piepho's R2 in percent: Total, Fixed, Random and each random component.
pm_r2_piepho <- function(model, over_disp=FALSE){

  stopifnot(any(
    inherits(model, "lmerMod"),
    inherits(model, "glmerMod"),
    inherits(model, "glmmTMB")))

  if (inherits(model, "glmmTMB")){
    # Rebuild the model structure (without fitting) where the model was fitted,
    # so the data are found even when they are not in the global environment.
    fit_env <- environment(formula(model))
    tmb_call <- model$call
    tmb_call$doFit <- FALSE
	tmbf <- try(eval(tmb_call, model$frame, fit_env), TRUE)
	if (inherits(tmbf, "try-error")) {
	  tmb_call$data <- base::merge(eval(tmb_call$data, fit_env), model$frame, sort=FALSE)
	  tmbf <- eval(tmb_call, fit_env)
	}
    weights <- tmbf$data.tmb$weights
    X <- as.matrix(getME(model, "X"))
    # conditional-model fixed effects (getME(, "beta") also holds the
    # dispersion and zero-inflation parameters); columns glmmTMB dropped as
    # aliased are not in X
    beta <- fixef(model)$cond[colnames(X)]
    covB <- vcov(model)$cond[colnames(X), colnames(X), drop = FALSE]
    family_n <- model$modelInfo$family$family
    if (family_n=="binomial") yobs <- unique(tmbf$data.tmb$yobs)
    link_n <- model$modelInfo$family$link
  }else{
    weights <- model@resp$weights
    Lambda <- getME(model, "Lambda")
    beta <- getME(model, "beta")
    if (inherits(model, "lmerMod")) family_n <- "gaussian" else family_n <- model@resp$family$family
    if (family_n=="binomial") yobs <- unique(model@resp$y)
    if (inherits(model, "glmerMod")) link_n <- model@resp$family$link
    covB <- vcov(model)
    X <- as.matrix(getME(model, "X"))
  }

  Z <- getME(model, "Z")
  re <- pm_re_terms(model)
  if (inherits(model, "glmmTMB")) Lambda <- pm_tmb_lambda(model, re)

  if (family_n != "gaussian"){
    if (family_n =="binomial" && all(yobs%in%c(0,1))){
      prior_weight <- switch(link_n,
                             "logit" = 3/(pi^2),
                             "probit" = 1,
                             "cloglog" = 6/(pi^2))
    }else{
      if (inherits(model, "glmmTMB")) {
        weight_fun <- function(eta, mod){
          family_fun <- mod$modelInfo$family
          family_fun$mu.eta(eta)^2/family_fun$variance(family_fun$linkinv(eta))
        }
        # getME(, "b") also holds the zero-inflation random effects (after
        # those of the conditional model)
        y_hat <- as.vector(cbind(X, Z)%*%c(beta, getME(model, "b")[seq_len(ncol(Z))]))
        prior_weight <- as.vector(weight_fun(y_hat, model))*tmbf$data.tmb$weights
        if (length(tmbf$data.tmb$size) > 0) prior_weight <- prior_weight*tmbf$data.tmb$size
      }else prior_weight <- model@resp$sqrtWrkWt()^2
    }
  }else prior_weight <- weights

  n <- nrow(X)
  s <- sigma(model)

  # The variance of the data is V = Z G Z' + diag(s^2/w), with G = s^2 Lambda Lambda',
  # and P = I - J/n centres it. Every quantity below is a trace of the form
  # tr(A P) = tr(A) - 1'A1/n, computed from Z Lambda (sparse, n x q) without
  # forming any n x n matrix.
  ZL <- Z %*% Lambda
  trP_ZGZ <- function(cols) {
    if (length(cols) == 0) return(0)
    B <- ZL[, cols, drop = FALSE]
    s^2 * (sum(B^2) - sum(Matrix::colSums(B)^2)/n)
  }
  resid_var <- s^2/prior_weight
  sum_resid <- if (length(resid_var) == 1) n*resid_var else sum(resid_var)
  all_cols <- seq_len(ncol(ZL))

  theta_V <- (trP_ZGZ(all_cols) + sum_resid*(1 - 1/n))/(n-1)
  XtPX <- as.matrix(crossprod(X) - tcrossprod(colSums(X))/n)
  theta_Xbeta <- (as.numeric(crossprod(beta, XtPX %*% beta)) - sum(XtPX*as.matrix(covB)))/(n-1)

  if (over_disp) {
    # Leave out an observation-level random effect (one level per observation)
    obs <- which(vapply(re, function(tm) length(tm$comps) == 1 && tm$nl == n, logical(1)))
    if (length(obs) == 0 || length(obs) == length(re)) {
      stop("'over_disp = TRUE' needs an observation-level random effect (one level per observation) in the model.")
    }
    theta_Zu <- trP_ZGZ(setdiff(all_cols, unlist(lapply(re[obs], `[[`, "cols"))))/(n-1)
  }else theta_Zu <- trP_ZGZ(all_cols)/(n-1)

  Omega_full <- theta_Xbeta+theta_V
  Omega_beta <- theta_Xbeta/Omega_full
  Omega_u <- theta_Zu/Omega_full
  Omega_beta_u <- Omega_beta+Omega_u

  Omega_beta_percent=Omega_beta*100
  Omega_u_percent=Omega_u*100
  Omega_beta_u_percent=Omega_beta_u*100

  if (!over_disp) {
    # Each random-effect component alone: sigma_k^2 tr(Z_k' P Z_k), with Z_k
    # every nc-th column of its term. Covariances between the components of
    # a term (e.g. intercept and slope) are only in the total.
    trP_cols <- function(cols) {
      B <- Z[, cols, drop = FALSE]
      sum(B^2) - sum(Matrix::colSums(B)^2)/n
    }
    Omega_u_list <- numeric(0)
    for (tm in re) {
      nc <- length(tm$comps)
      for (k in seq_len(nc)) {
        nm <- sub("\\.\\(Intercept\\)$", "", paste(tm$group, tm$comps[k], sep = "."))
        cols <- tm$cols[seq(k, length(tm$cols), by = nc)]
        Omega_u_list <- c(Omega_u_list, stats::setNames((tm$vcov[k, k]*trP_cols(cols)/(n-1))/Omega_full, nm))
      }
    }
  }else Omega_u_list <- NULL

  c(Total=Omega_beta_u_percent, Fixed=Omega_beta_percent, Random=Omega_u_percent,  unlist(Omega_u_list)*100)
}

#' R-squared for Mixed Models as Numbers
#'
#' @description
#' The R-squared of a mixed model as a numeric result. \code{method = "piepho"} (the default) gives the same values as
#' \code{\link{R2_glmm}} (Piepho, 2023), as proportions: total, fixed, random and each random-effect component.
#' \code{method = "nakagawa"} gives the marginal (fixed effects) and conditional (fixed and random effects) R-squared of
#' Nakagawa and Schielzeth (2013), with random slopes as Johnson (2014): the variance of each random-effect term is the
#' mean of \eqn{z_i' \Sigma z_i} over the observations.
#' @details
#' For \code{"nakagawa"} the observation-level variance is the residual variance for gaussian models and the variance
#' of the latent distribution for binomial models: \eqn{\pi^2/3} (logit), 1 (probit) or \eqn{\pi^2/6} (cloglog). Other
#' families have several competing approximations; use \code{method = "piepho"} for them.
#' @usage
#' R2_model(model, method = c("piepho", "nakagawa"), over_disp = FALSE)
#' @param model A model fitted by \code{lmer}, \code{glmer} or \code{glmmTMB}.
#' @param method \code{"piepho"} or \code{"nakagawa"}.
#' @param over_disp For \code{"piepho"}: as in \code{\link{R2_glmm}}.
#' @return
#' A named numeric vector of proportions (class \code{"pm_r2"}, printed as percentages). For \code{"nakagawa"},
#' \code{Marginal} and \code{Conditional}, with attribute \code{"variances"}: the fixed-effect, random-effect and
#' observation-level variances.
#' @references
#' Piepho, H.-P. (2023), An adjusted coefficient of determination (R2) for generalized linear mixed models in one go,
#' \emph{Biometrical Journal}, 65(7), 2200290.
#'
#' Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for obtaining R2 from generalized linear
#' mixed-effects models, \emph{Methods in Ecology and Evolution}, 4, 133-142.
#'
#' Johnson, P. C. D. (2014), Extension of Nakagawa & Schielzeth's R2GLMM to random slopes models, \emph{Methods in
#' Ecology and Evolution}, 5, 944-946.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#'   R2_model(fm)
#'   R2_model(fm, "nakagawa")
#' @name R2_model
NULL

R2_model <- function(model, method = c("piepho", "nakagawa"), over_disp = FALSE) {
  method <- match.arg(method)
  stopifnot("R2_model() needs a model fitted by lmer, glmer or glmmTMB" =
              inherits(model, c("lmerMod", "glmerMod", "glmmTMB")))
  if (method == "piepho") {
    return(structure(pm_r2_piepho(model, over_disp)/100, class = "pm_r2", method = "piepho"))
  }
  # Nakagawa & Schielzeth (2013), random slopes as Johnson (2014)
  X <- as.matrix(lme4::getME(model, "X"))
  beta <- if (inherits(model, "glmmTMB")) glmmTMB::fixef(model)$cond[colnames(X)] else lme4::fixef(model)[colnames(X)]
  var_fixed <- stats::var(drop(X %*% beta))
  Z <- lme4::getME(model, "Z")
  var_random <- sum(vapply(pm_re_terms(model), function(tm) {
    Zt <- as.matrix(Z[, tm$cols, drop = FALSE])
    G <- kronecker(diag(tm$nl), tm$vcov)
    mean(rowSums((Zt %*% G) * Zt))
  }, numeric(1)))
  fam <- if (inherits(model, "lmerMod")) list(family = "gaussian", link = "identity") else stats::family(model)
  var_obs <- if (identical(fam$family, "gaussian")) {
    stats::sigma(model)^2
  } else if (identical(fam$family, "binomial") && fam$link %in% c("logit", "probit", "cloglog")) {
    switch(fam$link, logit = pi^2/3, probit = 1, cloglog = pi^2/6)
  } else {
    stop("method = \"nakagawa\" is available for gaussian and binomial models; use method = \"piepho\".", call. = FALSE)
  }
  tot <- var_fixed + var_random + var_obs
  structure(c(Marginal = var_fixed/tot, Conditional = (var_fixed + var_random)/tot), class = "pm_r2",
            method = "nakagawa", variances = c(fixed = var_fixed, random = var_random, observation = var_obs))
}

print.pm_r2 <- function(x, digits = 2, ...) {
  m <- attr(x, "method")
  cat(if (identical(m, "nakagawa")) "# R2 for mixed models (Nakagawa and Schielzeth)\n" else
    "# Adjusted R2 for mixed models (Piepho)\n")
  print(stats::setNames(paste0(format(round(100*unclass(x), digits), nsmall = digits), "%"), names(x)), quote = FALSE)
  invisible(x)
}
