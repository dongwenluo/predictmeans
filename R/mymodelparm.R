mymodelparm <- function(model, coef., vcov., df, ...) {
  UseMethod("mymodelparm")
}

#' @export
#' @keywords internal
mymodelparm.default <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...) {
  ### extract coefficients and their covariance matrix
  if (inherits(model, "glmmTMB")) {
    beta <- try(coef.(model)[["cond"]])
  } else {
    beta <- try(coef.(model))
  }
  if (inherits(beta, "try-error")) {
    stop("no ", sQuote("coef"), " method for ", sQuote("model"), " found!")
  }

  if (inherits(model, "glmmTMB")) {
    sigma <- try(Matrix::as.matrix(vcov.(model)[["cond"]]))
  } else {
    sigma <- try(vcov.(model))
  }
  if (inherits(sigma, "try-error")) {
    stop("no ", sQuote("vcov"), " method for ", sQuote("model"), " found!")
  }
  sigma <- as.matrix(sigma)

  if (any(length(beta) != dim(sigma))) {
    beta <- na.omit(beta)
  }
  # stop("dimensions of coefficients and covariance matrix don't match")

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
    # if (inherits(model, "lmerMod")) {
    # L_list <- get_contrasts_type1(model)
    # df <- sapply(L_list, function(L) mean(df_term(model, ctrmatrix=L)))
    # }
    if (inherits(model, "glmerMod") || inherits(model, "glmmTMB")) {
      df <- summary(model)$AICtab["df.resid"]
    }
    if (inherits(model, "parm")) {
      df <- model$df
    }
  } else {
    if (df < 0) {
      stop(sQuote("df"), " is not positive")
    }
  }

  ### try to identify non-estimable coefficients
  ### coef.aov removes NAs, thus touch coefficients
  ### directly
  if (inherits(model, "glmmTMB")) {
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
    if (dim(sigma)[1] == length(estimable)) {
      sigma <- sigma[estimable, estimable]
    }
  }

  ### just in case...
  if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable)) {
    stop("could not extract coefficients and covariance matrix from ", sQuote("model"))
  }
  RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
  class(RET) <- "mymodelparm"
  RET
}

#' @export
#' @keywords internal
mymodelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...) {
  stop("This function does not support objects of class ", sQuote("aovlist"))
}

#' @export
#' @keywords internal
mymodelparm.lme <- function(model, coef. = nlme::fixef, vcov. = vcov, df = NULL, ...) {
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
}

#' @export
#' @keywords internal
mymodelparm.lmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...) {
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
}

#' @export
#' @keywords internal
mymodelparm.glmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...) {
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
}

#' @export
#' @keywords internal
mymodelparm.gls <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...) {
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
}

#' @export
#' @keywords internal
mymodelparm.glmmTMB <- function(model, coef. = glmmTMB::fixef, vcov. = vcov, df = NULL, ...) {
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
}
