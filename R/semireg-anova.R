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
#' the same fixed effects, so REML likelihoods can be compared for gaussian fits). The tests are those of
#' \code{varcomp(fit, test = TRUE)}: as the variance is on the boundary of its range under the null hypothesis, the
#' reference is the mixture \eqn{\frac12\chi^2_0 + \frac12\chi^2_1} (Self and Liang, 1987), and for the only random
#' effect of a gaussian fit the exact null distribution of the restricted likelihood ratio (Crainiceanu and Ruppert,
#' 2004), where the mixture is conservative; \code{lmerTest::ranova()} used the full chi-square, more conservative
#' still. Column \code{test} names the distribution. For an \code{s(x, f)} smooth the
#' block of the deviation curves tests whether the curves of the levels of \code{f} have the same shape. For small
#' samples \code{\link{permlmer}} or \code{\link{bootmodels}} on \code{fit$semer} give resampling tests.
#' @param object A \code{\link{semireg}} or \code{\link{semireg_tmb}} fit.
#' @param type \code{"III"} (default), \code{"II"} or \code{"I"} (also \code{3}, \code{2}, \code{1}): the type of
#'   the tests of the fixed effects. Type I is sequential (Satterthwaite F for gaussian fits, Wald chi-square
#'   otherwise); see \code{\link{anova.lmerMod}} for which type to use.
#' @param smooth Logical: test the smooth terms (refits the model once per smooth block)? The default is TRUE.
#' @param ... Not used.
#' @return A list of class \code{"pm_semi_anova"} with \code{fixed} (the table of the fixed effects) and
#'   \code{smooth} (a data frame: the variance of each smooth block, the log-likelihoods of the model with and without
#'   it, the likelihood ratio statistic, its p-value and the reference distribution), printed together. Advice on the
#'   type of tests is kept as for \code{\link{anova.lmerMod}} (see \code{\link{advice}}).
#' @seealso \code{\link{anova.lmerMod}}, \code{\link{varcomp}}, \code{\link{semipred}}
#' @references
#' Crainiceanu, C. M. and Ruppert, D. (2004). Likelihood ratio tests in linear mixed models with one variance
#' component. \emph{Journal of the Royal Statistical Society B}, 66, 165-185.
#'
#' Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of maximum likelihood estimators and likelihood ratio
#' tests under nonstandard conditions. \emph{Journal of the American Statistical Association}, 82, 605-610.
#' @examples
#' library(predictmeans)
#' data(Dialyzer, package = "nlme")
#' fm <- semireg(rate ~ s(pressure, QB, k = 6) + (1 | Subject), data = Dialyzer)
#' anova(fm)
#' @exportS3Method anova semireg
anova.semireg <- function(object, type = c("III", "II", "I"), smooth = TRUE, ...) {
  type <- pm_anova_type(type[1])
  semer <- object$semer
  fixed <- if (type == "I") pm_anova_type1(semer) else if (inherits(semer, "lmerMod"))
    pm_anova(semer, type = type, test_statistic = "F") else pm_anova(semer, type = type, test_statistic = "Chisq")
  fixed <- pm_anova_add_advice(fixed, semer, type)
  sm <- if (isTRUE(smooth)) pm_sr_smooth_tests(object)
  structure(list(fixed = fixed, smooth = sm), class = "pm_semi_anova")
}

#' @exportS3Method print pm_semi_anova
print.pm_semi_anova <- function(x, digits = 4, ...) {
  print(x$fixed)
  if (!is.null(x$smooth)) {
    cat("\nSmooth terms: likelihood ratio tests of variance = 0 (the curve is its straight-line part);\n",
        "reference: exact restricted LRT for the only random effect of a gaussian fit, otherwise 0.5 chi2(0) + ",
        "0.5 chi2(1)\n", sep = "")
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
  vc <- if (inherits(semer, "glmmTMB")) {
    v <- glmmTMB::VarCorr(semer)$cond
    stats::setNames(vapply(blocks, function(b) as.matrix(v[[b]])[1, 1], numeric(1)), blocks)
  } else {
    v <- as.data.frame(lme4::VarCorr(semer))
    stats::setNames(v$vcov[match(blocks, v$grp)], blocks)
  }
  # the tests of varcomp(test = TRUE) (R/vc-test.R): the same p-values
  tt <- pm_vc_test(object, rows = paste0(blocks, ".(Intercept)"))
  tt <- tt[match(paste0(blocks, ".(Intercept)"), tt$row), , drop = FALSE]
  data.frame(term = blocks, variance = unname(vc[blocks]), logLik = tt$logLik, logLik.without = tt$logLik.null,
             LRT = tt$LRT, p.value = tt$p.value, test = tt$test, stringsAsFactors = FALSE)
}

## Refit a semireg fit without the smooth blocks 'drop' (same fixed effects,
## the other smooth blocks and random terms kept; 'random': other random
## terms, as texts, instead of the fit's, for the tests of varcomp(),
## R/vc-test.R): the fitted model (merMod, glmmTMB, or lm/glm when no random
## term is left).
pm_sr_refit <- function(object, drop = character(0), random = NULL, whole = FALSE) {
  r <- object$refit
  keep <- setdiff(names(r$smoothZ), drop)
  f <- r$formula
  specials <- pm_special_terms(f)        # glmmTMB ar1(), ou(), ...: kept as written
  bars <- pm_plain_bars(f, specials)
  grp <- vapply(bars, function(b) paste(all.vars(b[[3]]), collapse = ":"), "")
  other <- bars[!grp %in% names(r$smoothZ)]
  fixed <- pm_nobars(f)
  e <- new.env(parent = asNamespace("predictmeans"))
  other_txt <- c(vapply(other, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"), ""),
                 vapply(specials, function(b) paste(deparse(b, width.cutoff = 500L), collapse = " "), ""))
  if (!is.null(random)) other_txt <- random
  rhs <- c(paste(deparse(fixed[[3]], width.cutoff = 500L), collapse = " "), other_txt)
  f_new <- stats::as.formula(paste(paste(deparse(f[[2]], width.cutoff = 500L), collapse = " "), "~",
                                   paste(rhs, collapse = " + ")), env = e)
  e$sr_data <- r$data
  e$sr_Z <- r$smoothZ[keep]
  e$sr_w <- r$weights
  e$sr_off <- r$offset
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  if (identical(r$engine, "glmmTMB")) {
    cl <- call("pm_semireg_tmb", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z), family = r$family,
               ziformula = r$ziformula, dispformula = r$dispformula, REML = r$REML, predict_info = FALSE)
    if (!is.null(r$weights)) cl$weights <- quote(sr_w)
    if (!is.null(r$offset)) cl$offset <- quote(sr_off)
  } else if (!length(keep) && !length(other_txt)) {
    # no random term left: lm / glm (REML log-likelihood for a gaussian REML fit)
    environment(fixed) <- e
    dd <- r$data
    dd$.pm_w <- if (is.null(r$weights)) rep(1, nrow(dd)) else r$weights
    dd$.pm_off <- if (is.null(r$offset)) rep(0, nrow(dd)) else r$offset
    # the weights and offset as columns of the data (found there by the fitting function)
    w <- as.name(".pm_w")
    off <- as.name(".pm_off")
    if (is.null(r$family)) {
      fit <- eval(call("lm", fixed, data = quote(dd), weights = w, offset = off))
      return(structure(list(fit = fit), class = "pm_lm_reml", REML = isTRUE(r$REML)))
    }
    return(eval(call("glm", fixed, family = r$family, data = quote(dd), weights = w, offset = off)))
  } else {
    cl <- call("pm_semireg_lme4", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z), REML = r$REML,
               predict_info = FALSE, control = r$control)
    if (!is.null(r$family)) cl$family <- r$family
    if (!is.null(r$weights) && any(r$weights != 1)) cl$weights <- quote(sr_w)
    if (!is.null(r$offset)) cl$offset <- quote(sr_off)
  }
  fit <- quiet(eval(cl, e))
  # whole = TRUE: the semireg result (its glmmTMB structure 'tmbf' is needed
  # by the parametric bootstrap, R/vc-test.R)
  if (whole) return(fit)
  if (is.list(fit) && !is.null(fit$semer)) fit$semer else fit
}

## logLik of an lm on the REML scale when the full fit was REML.
#' @exportS3Method logLik pm_lm_reml
logLik.pm_lm_reml <- function(object, ...) stats::logLik(object$fit, REML = isTRUE(attr(object, "REML")))
