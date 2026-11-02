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
    cond_reTrms <- tmbf$condList$reTrms
    if (is.null(cond_reTrms$Lambdat) || is.null(cond_reTrms$Lind)) {
      # glmmTMB >= 1.1.12 no longer stores Lambdat/Lind; rebuild them in
      # glmmTMB's own term order
      cond_reTrms <- lme4::mkReTrms(pm_findbars(stats::formula(model)), model$frame,
                                    reorder.terms = FALSE)
    }
    theta <- exp(getME(model, "theta"))/sigma(model)
    Lambda <- t(cond_reTrms$Lambdat)
    Lambda@x <- theta[cond_reTrms$Lind]
    beta <- getME(model, "beta")
    beta <- beta[setdiff(names(beta), c("betad", "betazi"))]
    family_n <- model$modelInfo$family$family
    if (family_n=="binomial") yobs <- unique(tmbf$data.tmb$yobs)
    link_n <- model$modelInfo$family$link
    covB <- vcov(model)$cond
  }else{
    weights <- model@resp$weights
    Lambda <- getME(model, "Lambda")
    beta <- getME(model, "beta")
    if (inherits(model, "lmerMod")) family_n <- "gaussian" else family_n <- model@resp$family$family
    if (family_n=="binomial") yobs <- unique(model@resp$y)
    if (inherits(model, "glmerMod")) link_n <- model@resp$family$link
    covB <- vcov(model)
  }

  X <- as.matrix(getME(model, "X"))
  Z <- getME(model, "Z")

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
        y_hat <- as.vector(cbind(X, Z)%*%c(beta, getME(model, "b")))
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

  # Term sizes, in the order of the columns of Z
  if (inherits(model, "glmmTMB")){
    u_ind <- c(ini=0, cond_reTrms$nl)
    u_ind <- cumsum(u_ind)
    u_indn <- names(u_ind)
    q_terms <- diff(c(0, cumsum(sapply(cond_reTrms$Ztlist, nrow))))
  }else{
    u_ind <- c(ini=0, sapply(getME(model, "Ztlist"), function(x) x@Dim[1]))
    u_ind <- cumsum(u_ind)
    u_indn <- sub("\\.\\(Intercept\\)$", "",names(u_ind))
    names(u_ind) <- u_indn
    q_terms <- sapply(getME(model, "Ztlist"), function(x) x@Dim[1])
  }

  if (over_disp) {
    # Leave out an observation-level random effect (one level per observation),
    # which comes first in Z.
    if (length(q_terms) == 0 || q_terms[1] != n || ncol(ZL) <= n) {
      stop("'over_disp = TRUE' needs an observation-level random effect (one level per observation) in the model.")
    }
    theta_Zu <- trP_ZGZ((n+1):ncol(ZL))/(n-1)
  }else theta_Zu <- trP_ZGZ(all_cols)/(n-1)

  Omega_full <- theta_Xbeta+theta_V
  Omega_beta <- theta_Xbeta/Omega_full
  Omega_u <- theta_Zu/Omega_full
  Omega_beta_u <- Omega_beta+Omega_u

  Omega_beta_percent=Omega_beta*100
  Omega_u_percent=Omega_u*100
  Omega_beta_u_percent=Omega_beta_u*100

  if (!over_disp) {
    Omega_u_list <- vector("list", length(u_indn)-1)
    names(Omega_u_list) <- u_indn[-1]

    for (i in 1:length(Omega_u_list)) {
	  contInds <- (u_ind[i] + 1):(u_ind[i + 1])
      Omega_u_list[i] <- (trP_ZGZ(contInds)/(n-1))/Omega_full
    }
  }else Omega_u_list <- NULL

  cod_v <- c(Total=Omega_beta_u_percent, Fixed=Omega_beta_percent, Random=Omega_u_percent,  unlist(Omega_u_list)*100)
  R2 <- paste(round(cod_v, 2), "%", sep="")
  names(R2) <- names(cod_v)
  cat("# Adjusted R2 for Mixed Models\n")
  return(R2)
}
