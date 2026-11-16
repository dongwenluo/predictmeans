#' Tests for a Semi-Parametric Model
#'
#' @description
#' Tests of the fixed effects and of the smooth terms of a \code{\link{semireg}} or \code{\link{semireg_tmb}} fit,
#' with the package's own code (no \pkg{lmerTest} needed).
#' @details
#' \strong{Fixed effects}: Type III (or II) tests from \code{pm_anova()}: for gaussian fits Wald F tests with
#' Kenward-Roger df and the adjusted covariance (the df equal those of the same model fitted by \code{nlme::lme});
#' for generalized fits Wald chi-square tests. The fixed parts of the smooths (e.g. \code{x} and \code{f:x} for
#' \code{s(x, by = f)}) are tested here: they are the straight-line parts of the curves.
#'
#' \strong{Smooth terms}: for each block of smooth random effects (one variance), the likelihood ratio test of the
#' variance being zero, i.e. of the curve being its straight-line part: the model is refitted without the block (with
#' the same fixed effects, so REML likelihoods can be compared for gaussian fits). As the variance is on the boundary
#' of its range under the null hypothesis, the p-value is half that of the chi-square with 1 df (Self and Liang,
#' 1987); \code{lmerTest::ranova()} used the full chi-square, which is conservative. For an \code{s(x, f)} smooth the
#' block of the deviation curves tests whether the curves of the levels of \code{f} have the same shape. For small
#' samples \code{\link{permlmer}} or \code{\link{bootmodels}} on \code{fit$semer} give resampling tests.
#' @param object A \code{\link{semireg}} or \code{\link{semireg_tmb}} fit.
#' @param type \code{"III"} (default) or \code{"II"}: the type of the tests of the fixed effects.
#' @param smooth Logical: test the smooth terms (refits the model once per smooth block)? The default is TRUE.
#' @param ... Not used.
#' @return A list of class \code{"pm_semi_anova"} with \code{fixed} (the table of the fixed effects) and
#'   \code{smooth} (a data frame: the variance of each smooth block, the log-likelihoods of the model with and without
#'   it, the likelihood ratio statistic and its p-value), printed together.
#' @references
#' Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum likelihood estimators and likelihood ratio
#' tests under nonstandard conditions. \emph{Journal of the American Statistical Association}, 82, 605-610.
#' @examples
#' library(predictmeans)
#' data(Dialyzer, package = "nlme")
#' fm <- semireg(rate ~ s(pressure, QB, k = 6) + (1 | Subject), data = Dialyzer, prt = FALSE)
#' anova(fm)
#' @name anova.semireg
NULL

anova.semireg <- function(object, type = c("III", "II"), smooth = TRUE, ...) {
  type <- match.arg(type)
  semer <- object$semer
  fixed <- if (inherits(semer, "lmerMod")) pm_anova(semer, type = type, test.statistic = "F") else
    pm_anova(semer, type = type, test.statistic = "Chisq")
  sm <- if (isTRUE(smooth)) pm_sr_smooth_tests(object)
  structure(list(fixed = fixed, smooth = sm), class = "pm_semi_anova")
}

print.pm_semi_anova <- function(x, digits = 4, ...) {
  print(x$fixed)
  if (!is.null(x$smooth)) {
    cat("\nSmooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);\n",
        "p-values halved for the boundary (chi-square 0:1)\n", sep = "")
    tab <- x$smooth
    tab$p.value <- format.pval(tab$p.value, digits = 3)
    print(format(tab, digits = digits), row.names = FALSE)
  }
  invisible(x)
}

## Likelihood ratio tests of the smooth blocks of a semireg fit (variance =
## 0), by refitting without each block (pm_sr_refit()).
pm_sr_smooth_tests <- function(object) {
  if (is.null(object$refit)) stop("This semireg fit cannot be refitted; fit it again with predictmeans >= 2.0.",
                                  call. = FALSE)
  semer <- object$semer
  blocks <- names(object$refit$smoothZ)
  ll_full <- as.numeric(stats::logLik(semer))
  vc <- if (inherits(semer, "glmmTMB")) {
    v <- glmmTMB::VarCorr(semer)$cond
    stats::setNames(vapply(blocks, function(b) as.matrix(v[[b]])[1, 1], numeric(1)), blocks)
  } else {
    v <- as.data.frame(lme4::VarCorr(semer))
    stats::setNames(v$vcov[match(blocks, v$grp)], blocks)
  }
  rows <- lapply(blocks, function(b) {
    red <- pm_sr_refit(object, drop = b)
    ll <- as.numeric(stats::logLik(red))
    stat <- max(0, 2*(ll_full - ll))
    data.frame(term = b, variance = unname(vc[b]), logLik = ll_full, logLik.without = ll, LRT = stat,
               p.value = 0.5*stats::pchisq(stat, 1, lower.tail = FALSE))
  })
  do.call(rbind, rows)
}

## Refit a semireg fit without the smooth blocks 'drop' (same fixed effects,
## the other smooth blocks and random terms kept): the fitted model (merMod,
## glmmTMB, or lm/glm when no random term is left).
pm_sr_refit <- function(object, drop) {
  r <- object$refit
  keep <- setdiff(names(r$smoothZ), drop)
  f <- r$formula
  specials <- pm_special_terms(f)        # glmmTMB ar1(), ou(), ...: kept as written
  bars <- pm_plain_bars(f, specials)
  grp <- vapply(bars, function(b) paste(all.vars(b[[3]]), collapse = ":"), "")
  other <- bars[!grp %in% names(r$smoothZ)]
  fixed <- pm_nobars(f)
  e <- new.env(parent = asNamespace("predictmeans"))
  rhs <- c(paste(deparse(fixed[[3]], width.cutoff = 500L), collapse = " "),
           vapply(other, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"), ""),
           vapply(specials, function(b) paste(deparse(b, width.cutoff = 500L), collapse = " "), ""))
  f_new <- stats::as.formula(paste(paste(deparse(f[[2]], width.cutoff = 500L), collapse = " "), "~",
                                   paste(rhs, collapse = " + ")), env = e)
  e$sr_data <- r$data
  e$sr_Z <- r$smoothZ[keep]
  e$sr_w <- r$weights
  e$sr_off <- r$offset
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  if (identical(r$engine, "glmmTMB")) {
    cl <- call("semireg_tmb", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z), family = r$family,
               ziformula = r$ziformula, dispformula = r$dispformula, REML = r$REML, prt = FALSE, predict_info = FALSE)
  } else if (!length(keep) && !length(other)) {
    # no random term left: lm / glm (REML log-likelihood for a gaussian REML fit)
    environment(fixed) <- e
    dd <- r$data
    dd$.pm_w <- if (is.null(r$weights)) rep(1, nrow(dd)) else r$weights
    dd$.pm_off <- if (is.null(r$offset)) rep(0, nrow(dd)) else r$offset
    # the weights and offset as columns of the data (found there by the fitting function)
    w <- as.name(".pm_w")
    off <- as.name(".pm_off")
    if (is.null(r$family) && !isTRUE(r$ngbinomial)) {
      fit <- eval(call("lm", fixed, data = quote(dd), weights = w, offset = off))
      return(structure(list(fit = fit), class = "pm_lm_reml", REML = isTRUE(r$REML)))
    }
    return(if (isTRUE(r$ngbinomial)) eval(call("glm.nb", fixed, data = quote(dd), weights = w),
                                          list(glm.nb = MASS::glm.nb, dd = dd)) else
      eval(call("glm", fixed, family = r$family, data = quote(dd), weights = w, offset = off)))
  } else {
    cl <- call("semireg", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z), REML = r$REML,
               ncenter = FALSE, nscale = FALSE, prt = FALSE, predict_info = FALSE, control = r$control)
    if (!is.null(r$family)) cl$family <- r$family
    if (isTRUE(r$ngbinomial)) cl$ngbinomial <- TRUE
    if (!is.null(r$weights) && any(r$weights != 1)) cl$weights <- quote(sr_w)
    if (!is.null(r$offset)) cl$offset <- quote(sr_off)
  }
  fit <- quiet(eval(cl, e))
  if (is.list(fit) && !is.null(fit$semer)) fit$semer else fit
}

## logLik of an lm on the REML scale when the full fit was REML.
logLik.pm_lm_reml <- function(object, ...) stats::logLik(object$fit, REML = isTRUE(attr(object, "REML")))
