mymodelparm <- function(model, coef., vcov., df, ...)
UseMethod("mymodelparm")

mymodelparm.default <- function(model, coef. = pm_coef, vcov. = pm_vcov, df = NULL, ...)
{
  ### extract coefficients and their covariance matrix
	if (inherits(model, "glmmTMB")) {
	  beta <- try(coef.(model)[["cond"]])
	} else {
	  beta <- try(coef.(model))
	}
  if (inherits(beta, "try-error")) {
    stop("no ", sQuote("coef"), " method for ",
         sQuote("model"), " found!")
  }

  if (inherits(model, "glmmTMB")) {
    sigma <- try(Matrix::as.matrix(vcov.(model)[["cond"]]))
  } else {
    sigma <- try(vcov.(model))
  }
  if (inherits(sigma, "try-error")) {
    stop("no ", sQuote("vcov"), " method for ",
         sQuote("model"), " found!")
  }
  sigma <- as.matrix(sigma)

  if (any(length(beta) != dim(sigma))) {
    beta = na.omit(beta)
    # stop("dimensions of coefficients and covariance matrix don't match")
  }

  ### determine degrees of freedom
  if (is.null(df)) {
    df <- 0
    ### check if a linear model was supplied
    if (inherits(model, "aov") || inherits(model, "lm") || inherits(model, "glm")) {
      class(model) <- "lm"
      df <- summary(model)$df[2]
    }
		if (inherits(model, "gls")) {
			dd <- model$dims
			df <- dd[["N"]] - dd[["p"]]
		}
		if (inherits(model, "lmerMod")) {
		    extra <- list(...)
		    df_type <- if (!is.null(extra$df_type)) extra$df_type else if (getME(model, "is_REML")) "Kenward-Roger" else "Satterthwaite"
		    L_list <- get_contrasts_type1(model)
            df <- if (identical(df_type, "containment")) pm_containment_df(model)[names(L_list)] else
              sapply(L_list, function(L) mean(df_term(model, ctrmatrix=L, type=df_type, vcov_adj=extra$vcov_adj)))
		}
		if (inherits(model, "glmerMod") || inherits(model, "glmmTMB")) {
			df <- summary(model)$AICtab["df.resid"]
		}
    if (inherits(model, "parm")) {
      df <- model$df
    }
    known <- c("lm", "gls", "lme", "lmerMod", "glmerMod", "glmmTMB", "parm")
    if (identical(df, 0) && !inherits(model, known)) {
      # other model classes: df.residual(), else normal-based (see ?pm_coef)
      df <- pm_residual_df(model)
    }
  } else {
    if (df < 0) {
      stop(sQuote("df"), " is not positive")
    }
  }

  ### try to identify non-estimable coefficients
  ### coef.aov removes NAs, thus touch coefficients
  ### directly
	if(inherits(model, "glmmTMB")) {
	  ocoef <- coef.(model)[["cond"]]
	} else {
	  ocoef <- coef.(model)
	}
  if (inherits(model, "aov")) {
    ocoef <- model$coefficients
  }
  estimable <- rep(TRUE, length(ocoef))
  if (any(is.na(ocoef))) {
    estimable[is.na(ocoef)] <- FALSE
    beta <- ocoef[estimable]
    if (dim(sigma)[1]==length(estimable)) {
      sigma <- sigma[estimable, estimable]
    }
  }

  ### just in case...
  if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable)) {
    stop("could not extract coefficients and covariance matrix from ",
         sQuote("model"))
  }

  RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
  class(RET) <- "mymodelparm"
  RET
}

mymodelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    stop("This function does not support objects of class ", sQuote("aovlist"))

mymodelparm.lme <- function(model, coef. = nlme::fixef, vcov. = vcov, df = NULL, ...)
    mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.lmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
   mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.glmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
   mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.gls <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.glmmTMB <- function(model, coef. = glmmTMB::fixef, vcov. = vcov, df = NULL, ...)
   mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

# mgcv gam/bam: all coefficients (parametric and smooth basis), the Bayesian
# posterior covariance matrix Vp, and the residual degrees of freedom.
mymodelparm.gam <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...) {
  beta <- coef.(model)
  sigma <- as.matrix(vcov.(model))
  if (is.null(df)) {
    df <- model$df.residual
  }
  RET <- list(coef = beta, vcov = sigma, df = df, estimable = rep(TRUE, length(beta)))
  class(RET) <- "mymodelparm"
  RET
}
