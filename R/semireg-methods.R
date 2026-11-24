## A semireg fit as a model: a short print, a table of the smooth terms, and
## the usual extractors, all from the fitted mixed model ('semer'), so code
## written for lmer/glmmTMB fits works on a semireg fit too.

## The smooth terms of a fit, one row per random-effect block: the smooth as
## written (label), its random-effect term, basis, number of basis functions,
## dimension of the basis (columns of its random-effect design), smoothing
## variance (on the scale of the fit) and whether it is estimated as a straight
## line (variance 0).
pm_sr_smooth_table <- function(fit) {
  sm <- fit$smooths
  blocks <- names(sm$maps)
  if (!length(blocks)) return(NULL)
  v <- pm_sr_block_vars(fit)
  bas <- vapply(seq_along(blocks), function(i) {
    s <- sm$sms[[min(i, length(sm$sms))]]
    sub("\\.smooth.spec$|\\.smooth$", "", class(s)[1])
  }, "")
  dim <- vapply(seq_along(blocks), function(i) {
    s <- sm$sms[[min(i, length(sm$sms))]]
    if (is.null(s$bs.dim)) NA_real_ else as.numeric(s$bs.dim)
  }, numeric(1))
  var <- if (is.null(v)) rep(NA_real_, length(blocks)) else unname(v[blocks])
  flat <- !is.na(var) & var < 1e-8 * max(1, abs(var), na.rm = TRUE)
  # the covariates of each block: what semipred(covariate = ) takes, and what
  # sm_term = "<covariate>" selects (every smooth of that covariate)
  covs <- vapply(blocks, function(b) {
    s <- sm$sms[vapply(sm$maps[[b]], `[[`, 1L, "sm")]
    paste(unique(unlist(lapply(s, `[[`, "term"))), collapse = ", ")
  }, "")
  data.frame(smooth = if (length(sm$labels) == length(blocks)) sm$labels else blocks, term = blocks,
             covariate = unname(covs), basis = bas, basis_dim = dim, variance = var, linear = flat, row.names = NULL)
}

#' @exportS3Method print semireg
print.semireg <- function(x, ...) {
  m <- x$semer
  tmb <- inherits(m, "glmmTMB")
  fam <- tryCatch(stats::family(m), error = function(e) NULL)
  cat("Semiparametric ", if (tmb) "glmmTMB" else if (inherits(m, "glmerMod")) "glmer" else "lmer", " fit",
      if (!is.null(fam)) paste0(" (", fam$family, ", ", fam$link, ")"),
      if (isTRUE(tryCatch(lme4::isREML(m), error = function(e) tmb && isTRUE(m$modelInfo$REML)))) ", REML" else "",
      "\n", sep = "")
  cat("Formula: ", paste(trimws(deparse(stats::formula(x))), collapse = " "), "\n", sep = "")
  ll <- tryCatch(stats::logLik(m), error = function(e) NULL)
  cat("Observations: ", stats::nobs(m), if (!is.null(ll)) paste0("; logLik ", format(signif(as.numeric(ll), 6)),
      ", AIC ", format(signif(stats::AIC(m), 6))), "\n", sep = "")
  st <- pm_sr_smooth_table(x)
  if (!is.null(st)) {
    cat("\nSmooth terms (variance of each smooth's random effects):\n")
    pm_sr_print_smooths(st)
  }
  cat("\nFixed effects:\n")
  print(signif(pm_coef(m), 5))
  vc <- tryCatch(pm_sr_other_vc(x), error = function(e) NULL)
  if (!is.null(vc) && nrow(vc)) {
    cat("\nOther variance components:\n")
    print(vc, row.names = FALSE)
  }
  pm_print_advice(x)
  invisible(x)
}

## The smooth table as printed, with a line on how semipred() selects them.
pm_sr_print_smooths <- function(st) {
  p <- st[c("smooth", "term", "covariate", "basis", "basis_dim", "variance", "linear")]
  p$variance <- signif(p$variance, 4)
  p$linear <- ifelse(p$linear, "yes", "")
  print(p, row.names = FALSE)
  ex <- st$smooth[1]
  cat(strwrap(paste0("semipred(): covariate = a 'covariate' above; sm_term = a 'smooth' or 'term' (e.g. \"", ex,
                      "\"), or a covariate for all its smooths; default every smooth."), width = 0.95 * getOption("width")),
      sep = "\n")
}

## Variance components other than the smooths (random effects of the formula,
## residual), as a small table.
pm_sr_other_vc <- function(x) {
  m <- x$semer
  blocks <- names(x$smooths$maps)
  if (inherits(m, "glmmTMB")) {
    vc <- glmmTMB::VarCorr(m)$cond
    d <- data.frame(group = names(vc), variance = vapply(vc, function(v) as.matrix(v)[1, 1], numeric(1)))
    s <- stats::sigma(m)
    fam <- stats::family(m)$family
    if (fam %in% c("gaussian")) d <- rbind(d, data.frame(group = "Residual", variance = s^2))
  } else {
    v <- as.data.frame(lme4::VarCorr(m))
    v <- v[is.na(v$var2), ]
    d <- data.frame(group = ifelse(v$var1 %in% c(NA, "(Intercept)"), v$grp, paste(v$grp, v$var1)), variance = v$vcov)
    d <- d[!(v$grp %in% blocks), , drop = FALSE]
  }
  d <- d[!(d$group %in% blocks), , drop = FALSE]
  d$variance <- signif(d$variance, 4)
  d
}

## summary(): the summary of the fitted model with the smooth table.
#' @exportS3Method summary semireg
summary.semireg <- function(object, ...) {
  structure(list(model = summary(object$semer, ...), smooths = pm_sr_smooth_table(object),
                 advice = pm_advice_get(object)), class = "summary.semireg")
}

#' @exportS3Method print summary.semireg
print.summary.semireg <- function(x, ...) {
  print(x$model, ...)
  if (!is.null(x$smooths)) {
    cat("\nSmooth terms:\n")
    pm_sr_print_smooths(x$smooths)
  }
  if (!is.null(x$advice)) cat(pm_format_advice(x$advice))
  invisible(x)
}

## Extractors: those of the fitted mixed model.
#' @exportS3Method stats::coef semireg
coef.semireg <- function(object, ...) pm_coef(object$semer)
#' @exportS3Method stats::vcov semireg
vcov.semireg <- function(object, ...) pm_vcov(object$semer)
#' @exportS3Method stats::fitted semireg
fitted.semireg <- function(object, ...) stats::fitted(object$semer, ...)
#' @exportS3Method stats::residuals semireg
residuals.semireg <- function(object, ...) stats::residuals(object$semer, ...)
#' @exportS3Method stats::logLik semireg
logLik.semireg <- function(object, ...) stats::logLik(object$semer, ...)
#' @exportS3Method stats::nobs semireg
nobs.semireg <- function(object, ...) stats::nobs(object$semer, ...)
## The formula as written (with s() terms); the expanded formula that was
## fitted is x$refit$formula.
#' @exportS3Method stats::formula semireg
formula.semireg <- function(x, ...) {
  f <- x$call$formula
  if (is.null(f)) return(x$refit$formula)
  f <- eval(f, environment(x$refit$formula))
  if (inherits(f, "formula")) f else stats::as.formula(f)
}
#' @exportS3Method stats::family semireg
family.semireg <- function(object, ...) stats::family(object$semer, ...)
#' @exportS3Method nlme::fixef semireg
fixef.semireg <- function(object, ...) pm_coef(object$semer)
#' @exportS3Method nlme::ranef semireg
ranef.semireg <- function(object, ...) nlme::ranef(object$semer, ...)
#' @exportS3Method nlme::VarCorr semireg
VarCorr.semireg <- function(x, sigma = 1, ...) {
  if (inherits(x$semer, "glmmTMB")) glmmTMB::VarCorr(x$semer, ...) else lme4::VarCorr(x$semer, ...)
}

## tidy(): one table of the smooth terms (as print shows); tidy(fit, "fixed")
## the fixed effects with SEs.
#' @exportS3Method generics::tidy semireg
tidy.semireg <- function(x, what = c("smooths", "fixed"), ...) {
  what <- match.arg(what)
  if (what == "smooths") return(pm_sr_smooth_table(x))
  b <- pm_coef(x$semer)
  se <- sqrt(diag(pm_vcov(x$semer)))[names(b)]
  data.frame(term = names(b), estimate = unname(b), std.error = unname(se), statistic = unname(b / se), row.names = NULL)
}
