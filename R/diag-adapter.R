## What the diagnostics (residplot(), pm_check(), CookD()) need to know about a
## fit, one method per class. Internal generics: a new model class adds
## methods here, the diagnostic functions stay as they are.

## Standardized residuals: deviance residuals over sqrt(scale) for gam,
## standardized deviance residuals for glm, rstandard() for lm, normalized for
## lme/gls, scaled (LMM) or Pearson (GLMM) for lme4, Pearson for glmmTMB,
## residuals over their SD otherwise.
pm_std_resid <- function(model) UseMethod("pm_std_resid")
#' @exportS3Method NULL
pm_std_resid.default <- function(model) {
  rr <- stats::residuals(model)
  rr / stats::sd(rr, na.rm = TRUE)
}
#' @exportS3Method NULL
pm_std_resid.gam <- function(model) stats::residuals(model, type = "deviance") / sqrt(model$scale)
#' @exportS3Method NULL
pm_std_resid.glm <- function(model) stats::rstandard(model, type = "deviance")
#' @exportS3Method NULL
pm_std_resid.lm <- function(model) stats::rstandard(model)
#' @exportS3Method NULL
pm_std_resid.lme <- function(model) {
  r <- if (!is.null(model$modelStruct$corStruct)) tryCatch(pm_nlme_normalized(model), error = function(e) NULL)
  if (is.null(r)) stats::residuals(model, type = "normalized") else r
}

## Normalized residuals of a gls or lme fit with a correlation structure,
## whitened group by group with the structure's own matrices (corMatrix()):
## nlme's residuals(type = "normalized") are misplaced when the data are not
## sorted by the correlation groups (checked: gls with corAR1(~ t | rep) on
## data sorted by field position differs by up to 2.2 from the hand
## computation, and is exact once sorted by rep). corAR1xAR1: the plots of a
## group in row, then column order, so a plot's residual does not depend on
## the order of the data. NULL when the residuals cannot be matched to the
## groups (nlme's are used).
pm_nlme_normalized <- function(model) {
  cs <- model$modelStruct$corStruct
  r <- stats::residuals(model, type = "pearson")
  keep <- !is.na(r)
  rk <- as.vector(r[keep])
  R <- nlme::corMatrix(cs)
  if (!is.list(R)) {
    R <- list(R)
    idx <- list(seq_along(rk))
  } else {
    g <- if (inherits(model, "lme")) model$groups[[ncol(model$groups)]] else model$groups
    if (length(g) != length(rk) || is.null(names(R))) return(NULL)
    idx <- lapply(names(R), function(nm) which(as.character(g) == nm))
  }
  if (!identical(unname(lengths(idx)), unname(vapply(R, NROW, integer(1))))) return(NULL)
  pos <- if (inherits(cs, "corAR1xAR1")) attr(cs, "covariate")
  if (!is.null(pos) && !is.list(pos)) pos <- list(pos)
  out <- rk
  for (k in seq_along(R)) {
    i <- idx[[k]]
    o <- if (is.null(pos)) seq_along(i) else order(pos[[k]][, 1], pos[[k]][, 2])
    U <- chol(as.matrix(R[[k]])[o, o, drop = FALSE])
    out[i[o]] <- backsolve(U, rk[i[o]], transpose = TRUE)
  }
  r[keep] <- out
  structure(as.vector(r), names = names(r), label = "Normalized residuals")
}
#' @exportS3Method NULL
pm_std_resid.gls <- pm_std_resid.lme
#' @exportS3Method NULL
pm_std_resid.merMod <- function(model) {
  if (lme4::isLMM(model)) stats::residuals(model, scaled = TRUE) else stats::residuals(model, type = "pearson")
}
#' @exportS3Method NULL
pm_std_resid.glmmTMB <- function(model) stats::residuals(model, type = "pearson")

## The predicted random effects as a list of data frames, one per grouping
## term (NULL for models without random effects).
pm_ranef_list <- function(model) UseMethod("pm_ranef_list")
#' @exportS3Method NULL
pm_ranef_list.default <- function(model) NULL
#' @exportS3Method NULL
pm_ranef_list.lme <- function(model) {
  r <- nlme::ranef(model)
  if (is.data.frame(r)) list(r) else r
}
#' @exportS3Method NULL
pm_ranef_list.merMod <- function(model) lme4::ranef(model)
#' @exportS3Method NULL
pm_ranef_list.glmmTMB <- function(model) lme4::ranef(model)[["cond"]]

## Number of levels of each grouping factor of the random effects (NULL if none).
pm_re_levels <- function(model) UseMethod("pm_re_levels")
#' @exportS3Method NULL
pm_re_levels.default <- function(model) NULL
#' @exportS3Method NULL
pm_re_levels.merMod <- function(model) lme4::ngrps(model)
#' @exportS3Method NULL
pm_re_levels.lme <- function(model) {
  g <- model$dims$ngrps
  g[setdiff(names(g), c("X", "y"))]
}
#' @exportS3Method NULL
pm_re_levels.glmmTMB <- function(model) vapply(model$modelInfo$reTrms$cond$flist, nlevels, 1L)

## The terms of the fixed part of the model.
pm_fixed_terms <- function(model) UseMethod("pm_fixed_terms")
#' @exportS3Method NULL
pm_fixed_terms.default <- function(model) stats::terms(model)
#' @exportS3Method NULL
pm_fixed_terms.merMod <- function(model) stats::terms(stats::formula(model, fixed.only = TRUE))

## Warnings about the fit (singular, not converged, Hessian), as short phrases.
pm_fit_warnings <- function(model) UseMethod("pm_fit_warnings")
#' @exportS3Method NULL
pm_fit_warnings.default <- function(model) character(0)
#' @exportS3Method NULL
pm_fit_warnings.merMod <- function(model) {
  out <- character(0)
  if (isTRUE(pm_singular(model))) out <- c(out, "singular (a variance component is zero)")
  msg <- model@optinfo$conv$lme4$messages
  if (length(msg)) out <- c(out, paste("convergence:", msg[1]))
  out
}
#' @exportS3Method NULL
pm_fit_warnings.glmmTMB <- function(model) if (isFALSE(model$sdr$pdHess)) "Hessian not positive definite" else character(0)
#' @exportS3Method NULL
pm_fit_warnings.lme <- function(model) {
  if (is.character(model$apVar)) "variance parameters: approximate covariance not positive definite" else character(0)
}
