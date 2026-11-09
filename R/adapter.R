## Model adapter: the few things predictmeans needs to know about a model.
##
## Kmatrix(), predictmeans(), contrastmeans() and covariatemeans() get these
## through the generics below. The default methods work for any model class
## with the usual methods (coef, vcov, terms, df.residual, family); to support
## a class that lacks them, or that needs something special (e.g. per-contrast
## degrees of freedom), write methods for these generics only. See ?pm_coef.

#' Supporting a New Model Class
#'
#' @description
#' \code{\link{Kmatrix}}, \code{\link{predictmeans}}, \code{\link{contrastmeans}} and \code{\link{covariatemeans}} get
#' everything they need from a model through these six generics. Their default methods work for any model class with
#' the usual methods (\code{coef}, \code{vcov}, \code{terms}, \code{df.residual}, \code{family}), so many classes need
#' no extra code. For a class that lacks one of these, or needs something special, write a method for the
#' corresponding generic only; nothing inside predictmeans needs to change.
#' @details
#' \describe{
#'   \item{\code{pm_coef(model)}}{Named fixed-effect estimates, with \code{NA} for aliased coefficients.
#'   Default: \code{coef(model)}; methods for \code{lme}, \code{merMod} and \code{glmmTMB} return the fixed effects.}
#'   \item{\code{pm_vcov(model)}}{Covariance matrix of the non-aliased estimates. Default: \code{vcov(model)}.}
#'   \item{\code{pm_basis(model, newdata, ...)}}{Design matrix of the fixed effects for \code{newdata}, with columns
#'   named as in \code{pm_coef(model)}. Default: \code{model.matrix()} of the model's terms (with their stored
#'   \code{predvars}, so \code{poly()}, \code{ns()}, \code{scale()} etc. use the fitted basis), contrasts and factor
#'   levels. The \code{gam} method uses the linear predictor matrix of package \pkg{mgcv}.}
#'   \item{\code{pm_data(model)}}{The data used to fit the model, or \code{NULL} (the default), in which case they are
#'   found from the \code{data} argument of the call, \code{getData()}, the model frame or the formula environment.}
#'   \item{\code{pm_df(model, L)}}{Denominator degrees of freedom for each row of \code{L}, or \code{NULL} (the
#'   default) to use the model's single residual df: \code{df.residual(model)}, or \code{Inf} (normal-based) when
#'   there is none. Kenward-Roger df for \code{lmer} and the rules of \pkg{nlme} for \code{lme} are built in.}
#'   \item{\code{pm_linkinv(model)}}{Inverse link function used to back-transform the means, or \code{NULL} for
#'   none. Default: the inverse link of \code{family(model)} unless it is the identity.}
#' }
#' The model also needs a \code{terms()} method for its fixed effects (used to find the variables in the model), and
#' its call (\code{getCall()}) or a \code{pm_data()} method so that the data can be found.
#' @usage
#' pm_coef(model, ...)
#' pm_vcov(model, ...)
#' pm_basis(model, newdata, ...)
#' pm_data(model, ...)
#' pm_df(model, L, ...)
#' pm_linkinv(model, ...)
#' @param model A fitted model.
#' @param newdata A data frame of predictor values (the reference grid).
#' @param L A matrix of linear functions of the coefficients, one row per contrast or mean.
#' @param \dots Further arguments. Kmatrix() passes \code{terms}, \code{contrasts} and \code{xlev} to
#'   \code{pm_basis()}, and \code{exclude} for \code{gam} models.
#' @return
#' As described in Details.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   # A made-up model class that wraps lm() and has only four methods.
#'   fit_mymodel <- function(formula, data) {
#'     structure(list(fit = lm(formula, data), call = match.call()), class = "mymodel")
#'   }
#'   coef.mymodel <- function(object, ...) coef(object$fit)
#'   vcov.mymodel <- function(object, ...) vcov(object$fit)
#'   terms.mymodel <- function(x, ...) terms(x$fit)
#'   df.residual.mymodel <- function(object, ...) object$fit$df.residual
#'
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   m <- fit_mymodel(yield ~ Block + nitro*Variety, data = Oats)
#'   predictmeans(m, "nitro", plot = FALSE)$means
#'
#'   # Degrees of freedom of its own, e.g. normal-based:
#'   pm_df.mymodel <- function(model, L, ...) rep(Inf, nrow(L))
#'   predictmeans(m, "nitro", plot = FALSE)$means
#' @name pm_coef
#' @aliases pm_vcov pm_basis pm_data pm_df pm_linkinv
NULL

## Fixed-effect estimates (named; NA for aliased coefficients).
pm_coef <- function(model, ...) UseMethod("pm_coef")
pm_coef.default <- function(model, ...) stats::coef(model)
pm_coef.lme <- function(model, ...) nlme::fixef(model)
pm_coef.merMod <- function(model, ...) lme4::fixef(model)
pm_coef.glmmTMB <- function(model, ...) glmmTMB::fixef(model)$cond

## Covariance matrix of the (non-aliased) fixed-effect estimates.
pm_vcov <- function(model, ...) UseMethod("pm_vcov")
pm_vcov.default <- function(model, ...) as.matrix(stats::vcov(model))
pm_vcov.glmmTMB <- function(model, ...) as.matrix(stats::vcov(model)$cond)

## Design matrix of the fixed effects for new data, with one column per
## coefficient (named as in pm_coef()). Kmatrix() passes the model's terms
## (with 'predvars', so poly(), ns() etc. use the fitted basis), contrasts and
## factor levels.
pm_basis <- function(model, newdata, ...) UseMethod("pm_basis")
pm_basis.default <- function(model, newdata, terms = NULL, contrasts = NULL, xlev = NULL, ...) {
  if (is.null(terms)) terms <- stats::delete.response(stats::terms(model))
  m <- stats::model.frame(terms, newdata, na.action = stats::na.pass, xlev = xlev)
  stats::model.matrix(terms, m, contrasts.arg = contrasts)
}
pm_basis.gam <- function(model, newdata, exclude = NULL, ...) {
  pm_need("mgcv", "'gam' models")
  mgcv::predict.gam(model, newdata = newdata, type = "lpmatrix", exclude = exclude)
}

## The data the model was fitted to, or NULL to let predictmeans find them
## (the 'data' argument of the call, getData(), the model frame, ...).
pm_data <- function(model, ...) UseMethod("pm_data")
pm_data.default <- function(model, ...) NULL

## Denominator degrees of freedom for each row of the contrast matrix L, or
## NULL to use the model's single residual df (the default). Kenward-Roger for
## lmer and the nlme rules for lme are handled inside predictmeans.
pm_df <- function(model, L, ...) UseMethod("pm_df")
pm_df.default <- function(model, L, ...) NULL

## Inverse link for back-transforming the means, or NULL for none.
pm_linkinv <- function(model, ...) UseMethod("pm_linkinv")
pm_linkinv.default <- function(model, ...) {
  fam <- tryCatch(stats::family(model), error = function(e) NULL)
  if (!inherits(fam, "family") || is.null(fam$linkinv) || identical(fam$link, "identity")) return(NULL)
  fam$linkinv
}
pm_linkinv.lm <- function(model, ...) NULL
pm_linkinv.glm <- function(model, ...) model$family$linkinv
pm_linkinv.glmerMod <- function(model, ...) slot(model, "resp")$family$linkinv
pm_linkinv.glmmTMB <- function(model, ...) model$modelInfo$family$linkinv

## The inverse link when it is not the identity (the means are then
## back-transformed with it), else NULL. With an identity link, e.g. a
## gaussian glm or glmmTMB of log(y), the user's 'trans' is used.
pm_link_trans <- function(model) {
  f <- pm_linkinv(model)
  if (is.null(f)) return(NULL)
  x <- c(-2, -0.5, 0, 0.5, 2)
  y <- tryCatch(suppressWarnings(as.numeric(f(x))), error = function(e) NULL)
  if (isTRUE(all.equal(y, x))) NULL else f
}

## Contrasts used to build the model matrix (internal).
pm_contrasts <- function(model) {
  if (inherits(model, "lme") || inherits(model, "gls")) return(model$contrasts)
  ctr <- tryCatch(attr(stats::model.matrix(model), "contrasts"), error = function(e) NULL)
  if (is.null(ctr)) ctr <- tryCatch(model$contrasts, error = function(e) NULL)
  ctr
}

## Residual degrees of freedom of a model class predictmeans has no rule for:
## df.residual() when it gives a positive number, else Inf (normal-based).
pm_residual_df <- function(model) {
  rd <- tryCatch(stats::df.residual(model), error = function(e) NULL)
  if (is.numeric(rd) && length(rd) == 1 && is.finite(rd) && rd > 0) rd else Inf
}
