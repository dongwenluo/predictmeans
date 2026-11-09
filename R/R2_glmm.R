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

  cod_v <- c(Total=Omega_beta_u_percent, Fixed=Omega_beta_percent, Random=Omega_u_percent,  unlist(Omega_u_list)*100)
  R2 <- paste(round(cod_v, 2), "%", sep="")
  names(R2) <- names(cod_v)
  cat("# Adjusted R2 for Mixed Models\n")
  return(R2)
}
