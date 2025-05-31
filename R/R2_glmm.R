#' An adjusted coefficient of determination (R2) for generalized linear mixed
#' models
#'
#' This function produces adjusted R2 for generalized linear mixed models which
#' was crafted following the guidance provided by Professor Hans-Peter Piepho.
#'
#'
#' @param model An object returned by \code{lmer}, \code{glmer} or
#' \code{glmmTMB}.
#' @param over_disp A logical scalar to indicate whether \code{model} with
#' over-dispersion or not. The default value is FALSE.
#' @return Adjusted R2 in percentage for Total (fixed + random), Fiexd, Random
#' and individual random term.
#' @references Piepho HP. An adjusted coefficient of determination (R2 ) for
#' generalized linear mixed models in one go. Biom J. 2023 Oct;65(7):e2200290.
#' doi: 10.1002/bimj.202200290. Epub 2023 May 1. PMID: 37127864.
#' @examples
#'
#'   library(predictmeans)
#'   Oats$nitro <- factor(Oats$nitro)
#'   (fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats))
#'   R2_glmm(fm)
#'   (gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#'               data = cbpp, family = binomial))
#'   R2_glmm(gm)
#' @importFrom Matrix t diag
#' @export
R2_glmm <- function(model, over_disp=FALSE) {
  stopifnot(
    any(
      inherits(model, "lmerMod"),
      inherits(model, "glmerMod"),
      inherits(model, "glmmTMB")
    )
  )

  tr <- function (m) {
    if ((dim(m)[1] != dim(m)[2])) {
      stop("m must be a square matrix")
    }
    return(sum(diag(m), na.rm = TRUE))
  }

  if (inherits(model, "glmmTMB")) {
    tmb_call <- model$call
    tmb_call$doFit <- FALSE
    tmbf <- try(with(model$frame, eval(tmb_call)), TRUE)
    if (inherits(tmbf, "try-error")) {
      tmb_call$data <- base::merge(eval(tmb_call$data), model$frame, sort=FALSE)
      tmbf <- eval(tmb_call)
    }
    # tmbf <- with(model$frame, eval(tmb_call))
    weights <- tmbf$data.tmb$weights
    cond_reTrms <- reTrms_tmb(model)
    Lambda <- t(cond_reTrms$Lambdat)

    beta <- getME(model, "beta")
    beta_L <- grep("betadisp|betazi", names(beta))
    if (length(beta_L > 0)) {
      beta <- beta[-beta_L]
    }
    b <- getME(model, "b")
    b_L <- grep("bzi", names(b))
    if (length(b_L > 0)) {
      b <- b[-b_L]
    }
    family_n <- model$modelInfo$family$family
    if (family_n=="binomial") {
      yobs <- unique(tmbf$data.tmb$yobs)
    }
    link_n <- model$modelInfo$family$link
    covB <- vcov(model)$cond
  } else {
    weights <- model@resp$weights
    Lambda <- getME(model, "Lambda")
    beta <- getME(model, "beta")
    if (inherits(model, "lmerMod")) {
      family_n <- "gaussian"
    } else {
      family_n <- model@resp$family$family
    }
    if (family_n=="binomial") {
      yobs <- unique(model@resp$y)
    }
    if (inherits(model, "glmerMod")) {
      link_n <- model@resp$family$link
    }
    covB <- vcov(model)
  }

  X <- getME(model, "X")
  Z <- getME(model, "Z")
  XZ <- cbind(X, Z)

  if (family_n != "gaussian") {
    if (family_n =="binomial" && all(yobs%in%c(0,1))) {
      prior_weight <- switch(link_n,
                             "logit" = 3/(pi^2),
                             "probit" = 1,
                             "cloglog" = 6/(pi^2))
    } else {
      if (inherits(model, "glmmTMB")) {
        weight_fun <- function(eta, mod) {
          family_fun <- mod$modelInfo$family
          family_fun$mu.eta(eta)^2/family_fun$variance(family_fun$linkinv(eta))
        }
        y_hat <- as.vector(XZ %*% c(beta, b))
        prior_weight <- as.vector(weight_fun(y_hat, model))*tmbf$data.tmb$weights
        if (length(tmbf$data.tmb$size) > 0) {
          prior_weight <- prior_weight*tmbf$data.tmb$size
        }
      } else {
        prior_weight <- model@resp$sqrtWrkWt()^2
      }
    }
  } else {
    prior_weight <- weights
  }

  n <- nrow(X)
  Zt <- t(Z)
  s <- sigma(model)
  G <- Matrix::tcrossprod(Lambda)*(s^2)
  V <- Z %*% G %*% Zt+diag(n)*(s^2/prior_weight)
  P <- diag(n) - matrix(1/n, n, n)
  theta_V <- tr(V %*% P)/(n-1)
  theta_Xbeta <- t(beta) %*% t(X) %*% P %*% X %*% (beta/(n-1)) - tr(t(X) %*% P %*% X %*% covB)/(n-1)
  if (over_disp) {
    theta_Zu <- tr(P %*% Z %*% G %*% diag(c(rep(0, n), rep(1, nrow(G) - n))) %*% Zt)/(n-1)
  } else {
    theta_Zu <- tr(P %*% Z %*% G %*% Zt)/(n-1)
  }

  Omega_full <- theta_Xbeta+theta_V
  Omega_beta <- theta_Xbeta/Omega_full
  Omega_u <- theta_Zu/Omega_full
  Omega_beta_u <- Omega_beta+Omega_u

  Omega_beta_percent=Omega_beta*100
  Omega_u_percent=Omega_u*100
  Omega_beta_u_percent=Omega_beta_u*100

  if (!over_disp) {
    if (inherits(model, "glmmTMB")) {
      u_ind <- c(ini=0, sapply(cond_reTrms$Ztlist, function(x) x@Dim[1]))
      u_ind <- cumsum(u_ind)
      u_indn <- names(u_ind)
    } else {
      u_ind <- c(ini=0, sapply(getME(model, "Ztlist"), function(x) x@Dim[1]))
      u_ind <- cumsum(u_ind)
      u_indn <- sub("\\.\\(Intercept\\)$", "",names(u_ind))
      names(u_ind) <- u_indn
    }

    Omega_u_list <- vector("list", length(u_indn)-1)
    names(Omega_u_list) <- u_indn[-1]

    for (i in 1:length(Omega_u_list)) {
      contInds <- (u_ind[i] + 1):(u_ind[i + 1])
      theta_Zu_i <- tr(P %*% Z[, contInds] %*% G[contInds, contInds] %*% Zt[contInds,])/(n-1)
      Omega_u_list[i] <- theta_Zu_i/Omega_full
    }
  } else {
    Omega_u_list <- NULL
  }

  cod_v <- c(Total=Omega_beta_u_percent, Fixed=Omega_beta_percent, Random=Omega_u_percent,  unlist(Omega_u_list)*100)
  R2 <- paste(round(cod_v, 2), "%", sep="")
  names(R2) <- names(cod_v)
  cat("# Adjusted R2 for Mixed Models\n")
  return(R2)
}

