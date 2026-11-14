#' Variance Components with Standard Errors and Confidence Intervals
#'
#' @description
#' The variance components of a fitted model, with standard errors and confidence intervals, for every model class
#' that \code{\link{predictmeans}} handles: the random-effect variances and covariances and the residual variance
#' of \code{lmer}, \code{lme} and \code{glmmTMB} models (also with correlation and variance structures and random
#' slopes), the random-effect (co)variances of \code{glmer} models, the correlation and variance-function parameters
#' of \code{gls} and \code{lme} models, the stratum variance components of \code{aov} models with \code{Error()},
#' the smoothing variance components of \code{mgcv::gam} models and \code{\link{semireg}} fits, and the error term
#' of \code{lm}, \code{aov} and \code{glm} models (residual variance, dispersion, or the negative binomial
#' \eqn{\theta}).
#' @usage
#' varcomp(model, ci = TRUE, level = 0.95, ci_method = c("Satterthwaite", "Wald", "profile"),
#'         information = c("observed", "average", "expected"))
#' @param model A fitted model (see Description).
#' @param ci \code{TRUE} (the default): add confidence intervals.
#' @param level Confidence level of the intervals; the default is 0.95.
#' @param ci_method How the intervals of variances are made: \code{"Satterthwaite"} (the default) chi-square
#'   intervals with \eqn{\nu = 2 (estimate/SE)^2} df, as SAS PROC MIXED (\code{COVTEST CL}); \code{"Wald"} Wald
#'   intervals on the log scale; \code{"profile"} (\code{lmer}, \code{glmer}) lme4's profile-likelihood intervals,
#'   which for a REML fit are those of the ML fit (lme4 profiles the ML deviance), so they do not go with the REML
#'   estimates. The other components always get the Wald intervals of the Details.
#' @param information For REML fits of \code{lmer}, \code{lme} and \code{gls} models: the information matrix of
#'   the variance parameters behind the SEs: \code{"observed"} (the default, as SAS PROC MIXED), \code{"average"}
#'   (the average information of Genstat's and ASReml's REML: Genstat's s.e.) or \code{"expected"}; the same
#'   choices as for the Kenward-Roger df (\code{df_method} of \code{\link{predictmeans}}). Other fits use the
#'   observed information.
#' @details
#' Standard errors are from the observed information of the variance parameters (the inverse Hessian of the
#' REML, or ML, deviance at the estimates, as SAS PROC MIXED and the Kenward-Roger df of predictmeans), carried to
#' the variance components by the delta method: \code{lmer} from the deviance in (theta, sigma), \code{glmer} from
#' the Hessian lme4 stores for (theta, beta) (\code{nAGQ = 0} fits: the deviance in theta), \code{lme} and
#' \code{gls} from the REML (ML) log-likelihood in their unconstrained \code{modelStruct} parameters (any
#' \code{pdMat}, \code{corStruct} and \code{varFunc}), \code{glmmTMB} from \code{vcov(full = TRUE)}.
#'
#' Intervals of variances: Satterthwaite's chi-square interval
#' \eqn{[\nu \hat\sigma^2/\chi^2_{\nu, 1-\alpha/2}, \nu \hat\sigma^2/\chi^2_{\nu, \alpha/2}]} with
#' \eqn{\nu = 2 (\hat\sigma^2/SE)^2} (positive and skewed as the sampling distribution of a variance; for the
#' residual variance of an \code{lm} it is the exact interval, \eqn{\nu} the residual df; for the stratum
#' components of an \code{aov} Satterthwaite's interval from the mean squares). Other components: Wald intervals
#' on the scale where the parameter is unbounded: atanh for correlations (a covariance: the interval of its
#' correlation times the two SDs), log for other positive parameters (\code{varIdent} ratios, ranges,
#' \eqn{\theta}), and the natural scale otherwise (e.g. the power of \code{varPower}). A variance estimated at its
#' bound, zero, has no SE (the interval is (0, NA)), as SAS prints it; an analysis of variance estimate can be
#' negative.
#'
#' \code{lm}, \code{aov}, gaussian \code{glm}: the residual variance \eqn{s^2} with SE \eqn{s^2 \sqrt{2/df}} and
#' the exact chi-square interval. Other \code{glm} families with an estimated dispersion (Gamma, inverse gaussian,
#' quasi): the Pearson dispersion, with the same (approximate) SE and interval. Poisson and binomial: the
#' dispersion is fixed at 1 (no SE). Negative binomial (\code{MASS::glm.nb}): \eqn{\theta} with its SE.
#'
#' \code{aov} with \code{Error()}: the analysis of variance (method of moments) estimates of the variance of each
#' error stratum, from the stratum mean squares and their expectations (as Genstat's stratum variances; they equal
#' the REML estimates in balanced designs when positive, and can be negative), with SEs from
#' \eqn{var(MS) = 2 MS^2/df} and Satterthwaite's chi-square intervals; with missing values, those of the design
#' completed as in \code{\link{predictmeans}}.
#'
#' \code{gam}: \code{mgcv::gam.vcomp()} (variances of the smooths, and the scale), its intervals squared, and the
#' SEs from their width on the log scale; it needs \code{method = "REML"} or \code{"ML"}.
#' @return A data frame with one row per component (random-effect variances and covariances by group, named
#'   \code{group.term} and \code{group.term1.term2}; correlation and variance-function parameters named
#'   \code{cor.name} and \code{var.name}; \code{residual}) and columns \code{vcov} (the estimate), \code{SE} and,
#'   with \code{ci = TRUE}, the limits. Values are rounded to 4 decimals; attribute \code{"ci_method"} says how
#'   the intervals were made.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' varcomp(fm)
#' fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' varcomp(fm1)
#' # the same components by the analysis of variance
#' varcomp(aov(yield ~ nitro*Variety + Error(Block/Variety), data=Oats))
#'
#' data(Orthodont, package="nlme")
#' mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
#' varcomp(mod1)
#' # correlation and variance structures
#' mod2 <- gls(distance ~ age*Sex, correlation=corAR1(form=~1|Subject),
#'             weights=varIdent(form=~1|Sex), data=Orthodont)
#' varcomp(mod2)
#' # the error term of a linear model
#' varcomp(lm(yield ~ Block + nitro*Variety, data=Oats))
#' @name varcomp
NULL

varcomp <- function(model, ci = TRUE, level = 0.95, ci_method = c("Satterthwaite", "Wald", "profile"),
                    information = c("observed", "average", "expected")) {
  ci_method <- match.arg(ci_method)
  information <- match.arg(information)
  if (inherits(model, "semireg")) model <- model$semer
  if (ci_method == "profile" && !inherits(model, c("lmerMod", "glmerMod"))) {
    message("Profile intervals are for lmer and glmer models; Satterthwaite intervals instead.")
    ci_method <- "Satterthwaite"
  }
  ci <- if (isTRUE(ci)) ci_method else "none"     # passed on: how to make the intervals
  if (inherits(model, "aovlist")) return(pm_vc_aov(model, ci, level))
  if (information != "observed") {
    reml <- if (inherits(model, "lmerMod")) lme4::getME(model, "is_REML") else
      inherits(model, c("lme", "gls")) && identical(model$method, "REML")
    if (!reml) {
      message("information = \"", information, "\" is for REML fits of lmer, lme and gls models; observed instead.")
      information <- "observed"
    }
  }
  if (inherits(model, "lmerMod")) return(pm_vc_lmer(model, ci, level, information))
  if (inherits(model, "glmerMod")) return(pm_vc_glmer(model, ci, level))
  if (inherits(model, c("lme", "gls"))) return(pm_vc_nlme(model, ci, level, information))
  if (inherits(model, "glmmTMB")) return(pm_tmb_varcomp(model, ci = ci, level = level))
  if (inherits(model, "gam")) return(pm_vc_gam(model, ci, level))
  if (inherits(model, "glm")) return(pm_vc_glm(model, ci, level))
  if (inherits(model, "lm") && !inherits(model, "mlm")) return(pm_vc_lm(model, ci, level))
  stop("varcomp() handles lm, aov (also with Error()), glm, gls, lme, lmer, glmer, glmmTMB, gam and semireg models.",
       call. = FALSE)
}

## ---- shared pieces ----

## The result table: estimates 'est' with their covariance 'V' (or SEs 'se'),
## and the intervals ('ci': "Satterthwaite", "Wald", "profile" or "none") by
## kind of component ('type'): "var" variances, Satterthwaite's chi-square
## with nu = 2 (est/SE)^2 or Wald on the log scale; "cov" covariances, the
## correlation's atanh interval times the two SDs ('pair': the names of the
## two variances); "cor" correlation parameters, atanh; "pos" other positive
## parameters, log; "par" natural scale. 'lower'/'upper' given: intervals
## made elsewhere (profile).
pm_vc_table <- function(est, V = NULL, se = NULL, type, pair = NULL, ci = "Satterthwaite", level = 0.95,
                        lower = NULL, upper = NULL) {
  if (is.null(se)) se <- sqrt(pmax(diag(as.matrix(V)), 0))
  se <- unname(se)
  # a variance estimated at its bound, zero (or below 1e-6 of the largest:
  # nlme only approaches it): no SE, the interval (0, NA), as SAS prints it
  vr <- type == "var"
  bound <- vr & is.finite(est) & est >= 0 & est <= 1e-6 * max(c(abs(est[vr]), 1e-300), na.rm = TRUE)
  est[bound] <- 0
  se[bound] <- NA_real_
  out <- data.frame(vcov = as.numeric(est), SE = se, row.names = names(est))
  if (ci != "none") {
    if (is.null(lower)) {
      z <- stats::qnorm(1 - (1 - level)/2)
      lower <- upper <- rep(NA_real_, length(est))
      for (i in seq_along(est)) {
        e <- est[[i]]; s <- se[i]
        if (bound[i]) {
          lower[i] <- 0
          next
        }
        if (!is.finite(e) || !is.finite(s)) next
        if (type[i] == "var" && ci != "Wald") {
          nu <- if (s > 0) 2 * (e/s)^2 else Inf
          lower[i] <- if (e > 0) nu * e / stats::qchisq(1 - (1 - level)/2, nu) else 0
          upper[i] <- if (e > 0) nu * e / stats::qchisq((1 - level)/2, nu) else NA_real_
          if (!is.finite(nu)) lower[i] <- upper[i] <- e
        } else if (type[i] %in% c("var", "pos")) {
          lower[i] <- if (e > 0) e * exp(-z * s/e) else 0
          upper[i] <- if (e > 0) e * exp(z * s/e) else NA_real_
        } else if (type[i] == "cor" && abs(e) < 1) {
          lower[i] <- tanh(atanh(e) - z * s/(1 - e^2))
          upper[i] <- tanh(atanh(e) + z * s/(1 - e^2))
        } else if (type[i] == "cov" && !is.null(V)) {
          ab <- match(pair[[names(est)[i]]], names(est))
          va <- est[[ab[1]]]; vb <- est[[ab[2]]]
          if (va > 0 && vb > 0) {
            sdp <- sqrt(va * vb)
            r <- e/sdp
            g <- numeric(length(est))                # gradient of r = c / sqrt(va vb)
            g[i] <- 1/sdp
            g[ab[1]] <- -r/(2 * va)
            g[ab[2]] <- -r/(2 * vb)
            sr <- sqrt(max(drop(t(g) %*% as.matrix(V) %*% g), 0))
            if (abs(r) < 1) {
              lower[i] <- tanh(atanh(r) - z * sr/(1 - r^2)) * sdp
              upper[i] <- tanh(atanh(r) + z * sr/(1 - r^2)) * sdp
            }
          }
        } else {
          lower[i] <- e - z * s
          upper[i] <- e + z * s
        }
      }
    }
    out$lower <- unname(lower)
    out$upper <- unname(upper)
    names(out)[3:4] <- paste(format(100 * c((1 - level)/2, 1 - (1 - level)/2), trim = TRUE, scientific = FALSE,
                                    digits = 3), "%")
  }
  out <- round(out, 4)
  attr(out, "ci_method") <- ci
  out
}

## The delta method: components fun(x) (named, with attributes "type" and
## "pair") of parameters x with covariance Vx; returns the table.
pm_vc_delta <- function(fun, x0, Vx, ci, level, lower = NULL, upper = NULL) {
  val <- fun(x0)
  J <- numDeriv::jacobian(function(x) as.numeric(fun(x)), x0)
  V <- J %*% as.matrix(Vx) %*% t(J)
  dimnames(V) <- list(names(val), names(val))
  pm_vc_table(val, V, type = attr(val, "type"), pair = attr(val, "pair"), ci = ci, level = level,
              lower = lower, upper = upper)
}

## The (co)variances of covariance matrices: list of matrices named by group;
## rows "group.term" (variances) and "group.term_a.term_b" (covariances, a > b),
## by column of the lower triangle; 'diag_only' leaves the covariances out.
pm_vc_entries <- function(mats, diag_only = NULL) {
  est <- numeric(0); type <- character(0); pair <- list()
  for (k in seq_along(mats)) {               # (by position: a group can have several terms)
    g <- names(mats)[k]
    M <- as.matrix(mats[[k]])
    cn <- colnames(M)
    if (is.null(cn)) cn <- paste0("V", seq_len(ncol(M)))
    for (b in seq_len(ncol(M))) for (a in b:nrow(M)) {
      if (a != b && isTRUE(diag_only[[k]])) next
      nm <- paste(c(g, cn[unique(c(a, b))]), collapse = ".")
      est[nm] <- M[a, b]
      type[nm] <- if (a == b) "var" else "cov"
      if (a != b) pair[[nm]] <- c(paste(g, cn[a], sep = "."), paste(g, cn[b], sep = "."))
    }
  }
  structure(est, type = unname(type), pair = pair)
}

## Join component vectors (each with "type" and "pair").
pm_vc_c <- function(...) {
  parts <- list(...)
  structure(unlist(lapply(parts, function(p) stats::setNames(as.numeric(p), names(p)))),
            type = unlist(lapply(parts, function(p) attr(p, "type"))),
            pair = do.call(c, lapply(parts, function(p) attr(p, "pair"))))
}

## Inverse of a symmetric information matrix, by its eigenvalues above 1e-8
## (as the Satterthwaite df always did).
pm_vc_inv <- function(H) {
  eg <- eigen((H + t(H))/2, symmetric = TRUE)
  pos <- eg$values > 1e-8
  eg$vectors[, pos, drop = FALSE] %*% (t(eg$vectors[, pos, drop = FALSE]) / eg$values[pos])
}

## ---- lmer and glmer ----

## Random-effect covariance matrices of a merMod from theta (and sigma):
## Lambda's lower-triangular blocks, by term, in lme4's theta order.
pm_vc_mer_fun <- function(model, residual) {
  cnms <- lme4::getME(model, "cnms")
  grp <- names(cnms)                           # (|| terms repeat a group, with other terms)
  nc <- lengths(cnms)
  nth <- nc * (nc + 1)/2
  diag_only <- as.list(rep(FALSE, length(grp)))
  function(x) {
    s <- if (residual) x[length(x)] else 1
    th <- if (residual) x[-length(x)] else x
    pos <- 0
    mats <- list()
    for (k in seq_along(cnms)) {
      Tk <- matrix(0, nc[k], nc[k])
      Tk[lower.tri(Tk, diag = TRUE)] <- th[pos + seq_len(nth[k])]
      pos <- pos + nth[k]
      M <- s^2 * tcrossprod(Tk)
      dimnames(M) <- list(cnms[[k]], cnms[[k]])
      mats[[k]] <- M
    }
    names(mats) <- grp
    re <- pm_vc_entries(mats, diag_only)
    if (!residual) return(re)
    pm_vc_c(re, structure(c(residual = s^2), type = "var", pair = list()))
  }
}

## Variance parameters (theta, sigma) of an lmer fit, with their covariance
## from the observed information of its REML (ML) deviance; also used for the
## Satterthwaite df (R/satterthwaite.R).
pm_lmer_varpar <- function(model) {
  cl <- stats::getCall(model)
  cl[[1]] <- quote(lme4::lmer)
  cl$devFunOnly <- TRUE
  if (is.null(cl$control)) cl$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
  devfun <- eval(cl, pm_fit_env(model))
  is_reml <- lme4::getME(model, "is_REML")
  varpar <- unname(c(lme4::getME(model, "theta"), stats::sigma(model)))
  h <- numDeriv::hessian(func = devfun_vp, x = varpar, devfun = devfun, reml = is_reml)
  devfun(lme4::getME(model, "theta"))         # leave the deviance function at the optimum
  list(devfun = devfun, varpar = varpar, vcov_varpar = 2 * pm_vc_inv(h))
}

pm_vc_lmer <- function(model, ci, level, information = "observed") {
  fun <- pm_vc_mer_fun(model, residual = TRUE)
  x0 <- unname(c(lme4::getME(model, "theta"), stats::sigma(model)))
  lims <- if (ci == "profile") pm_vc_mer_profile(model, fun(x0), level)
  if (ci == "profile" && is.null(lims)) ci <- "Satterthwaite"
  if (information != "observed") {
    # the covariance of the components themselves (their order: by term, the
    # lower triangle by column, then the residual), from R/kr.R
    val <- fun(x0)
    W <- attr(pm_vcov_kr(model, information), "W")
    return(pm_vc_table(val, W, type = attr(val, "type"), pair = attr(val, "pair"), ci = ci, level = level,
                       lower = lims$lower, upper = lims$upper))
  }
  vp <- pm_lmer_varpar(model)
  pm_vc_delta(fun, vp$varpar, vp$vcov_varpar, ci, level, lims$lower, lims$upper)
}

pm_vc_glmer <- function(model, ci, level) {
  th <- lme4::getME(model, "theta")
  H <- model@optinfo$derivs$Hessian
  Vth <- if (!is.null(H) && nrow(H) == length(th) + length(lme4::fixef(model))) {
    (2 * pm_vc_inv(H))[seq_along(th), seq_along(th), drop = FALSE]     # Hessian of the deviance in (theta, beta)
  } else {
    # (nAGQ = 0: no stored Hessian) the deviance in theta, beta profiled
    devfun <- tryCatch(stats::update(model, devFunOnly = TRUE), error = function(e) NULL)
    h <- if (is.function(devfun)) tryCatch(numDeriv::hessian(devfun, th), error = function(e) NULL)
    if (!is.null(devfun)) tryCatch(devfun(th), error = function(e) NULL)
    if (is.null(h) || nrow(h) != length(th)) matrix(NA_real_, length(th), length(th)) else 2 * pm_vc_inv(h)
  }
  fun <- pm_vc_mer_fun(model, residual = FALSE)
  lims <- if (ci == "profile") pm_vc_mer_profile(model, fun(th), level)
  if (ci == "profile" && is.null(lims)) ci <- "Satterthwaite"
  pm_vc_delta(fun, th, Vth, ci, level, lims$lower, lims$upper)
}

## Profile-likelihood intervals of a merMod (confint()): SD intervals squared
## for variances, the correlation's interval times the two SDs for
## covariances. NULL (Wald instead) when profiling fails.
pm_vc_mer_profile <- function(model, val, level) {
  pr <- tryCatch(suppressWarnings(suppressMessages(
    stats::confint(model, parm = "theta_", method = "profile", level = level, oldNames = TRUE))),
    error = function(e) NULL)
  if (is.null(pr) || nrow(pr) != length(val)) {
    message("Profile-likelihood intervals failed; Satterthwaite intervals instead.")
    return(NULL)
  }
  type <- attr(val, "type")
  pair <- attr(val, "pair")
  lower <- upper <- numeric(length(val))
  for (i in seq_along(val)) {
    if (type[i] == "var") {
      lower[i] <- max(pr[i, 1], 0)^2
      upper[i] <- pr[i, 2]^2
    } else {
      ab <- match(pair[[names(val)[i]]], names(val))
      sdp <- sqrt(val[[ab[1]]] * val[[ab[2]]])
      lower[i] <- pr[i, 1] * sdp
      upper[i] <- pr[i, 2] * sdp
    }
  }
  list(lower = lower, upper = upper)
}

## ---- lme and gls ----

## Components of an nlme fit as functions of its unconstrained modelStruct
## parameters and log(sigma) (R/satt-nlme.R): random-effect (co)variances by
## level, correlation parameters ("cor."), variance-function parameters
## ("var."), the residual variance (of the reference group with a varFunc).
pm_vc_nlme_fun <- function(model) {
  ms <- model$modelStruct
  lme <- inherits(model, "lme")
  function(x) {
    k <- length(x)
    if (k > 1) ms <- nlme::`coef<-`(ms, value = x[-k])
    s2 <- exp(2 * x[k])
    parts <- list()
    if (lme) {
      mats <- lapply(as.matrix(ms$reStruct), function(M) M * s2)
      parts$re <- pm_vc_entries(mats, lapply(ms$reStruct, function(p) inherits(p, c("pdDiag", "pdIdent"))))
    }
    if (!is.null(ms$corStruct) && length(stats::coef(ms$corStruct))) {
      p <- stats::coef(ms$corStruct, unconstrained = FALSE)
      if (is.null(names(p))) {
        names(p) <- if (inherits(ms$corStruct, "corSymm")) {   # by row of the upper triangle
          d <- (1 + sqrt(1 + 8 * length(p)))/2
          ij <- which(upper.tri(diag(d)), arr.ind = TRUE)
          ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
          paste(ij[, 1], ij[, 2], sep = ".")
        } else seq_along(p)
      }
      parts$cor <- structure(stats::setNames(p, paste0("cor.", names(p))),
                             type = ifelse(grepl("range", names(p), ignore.case = TRUE), "pos", "cor"),
                             pair = list())
    }
    if (!is.null(ms$varStruct) && length(stats::coef(ms$varStruct))) {
      p <- stats::coef(ms$varStruct, unconstrained = FALSE)
      if (is.null(names(p))) names(p) <- seq_along(p)
      parts$var <- structure(stats::setNames(p, paste0("var.", names(p))),
                             type = if (inherits(ms$varStruct, "varIdent")) rep("pos", length(p)) else
                               ifelse(names(p) == "const", "pos", "par"),
                             pair = list())
    }
    parts$res <- structure(c(residual = s2), type = "var", pair = list())
    do.call(pm_vc_c, unname(parts))
  }
}

pm_vc_nlme <- function(model, ci, level, information = "observed") {
  fun <- pm_vc_nlme_fun(model)
  x0 <- c(stats::coef(model$modelStruct), log(model$sigma))
  W <- if (information != "observed") {
    tryCatch(attr(pm_vcov_kr_nlme(model, information), "W"), error = function(e) NULL)   # R/satt-nlme.R
  } else {
    R <- tryCatch(pm_nlme_reml(model), error = function(e) NULL)
    if (!is.null(R)) tryCatch(pm_nlme_W_observed(R), error = function(e) NULL)
  }
  if (is.null(W)) {
    message("The covariance structure of this ", class(model)[1], " model could not be rebuilt: no SEs.")
    W <- matrix(NA_real_, length(x0), length(x0))
  }
  pm_vc_delta(fun, unname(x0), W, ci, level)
}

## ---- glmmTMB ----

## Variance components of a glmmTMB fit (conditional model, then the
## zero-inflation model with prefix "zi."): variances and covariances of each
## term ("us"), variances only ("diag"), and for structured terms the
## variance(s) and the correlation parameter(s) ("ar1", "hetar1", "cs",
## "homcs", "toep": rows "cor.group" or "cor.group.lag"); other structures:
## every (co)variance. Then the residual variance of a gaussian model. SEs
## by the delta method from vcov(full = TRUE): the components are recomputed
## by VarCorr() from perturbed parameters, so this works for every structure.
pm_tmb_varcomp <- function(model, ci = "Satterthwaite", level = 0.95) {
  pf <- model$fit$parfull
  keep <- !names(pf) %in% c("b", "bzi")   # vcov(full = TRUE) has these, in this order
  pn <- names(pf)[keep]
  V <- vcov(model, full = TRUE)
  if (nrow(V) != length(pn)) stop("cannot match the parameters of this glmmTMB model to vcov(full = TRUE)")
  gaus <- stats::family(model)$family == "gaussian"
  use <- pn %in% c("theta", "thetazi", if (gaus) c("betad", "betadisp"))
  pos <- which(keep)[use]
  on.exit(model$obj$env$report(model$fit$parfull))   # leave the TMB object as it was

  fun <- function(x) {
    m <- model
    m$fit$parfull[pos] <- x
    for (nm in unique(pn[use])) m$fit$par[names(m$fit$par) == nm] <- x[pn[use] == nm]
    vc <- glmmTMB::VarCorr(m)
    parts <- list()
    for (part in c("cond", "zi")) {
      for (j in seq_along(vc[[part]])) {
        v <- vc[[part]][[j]]
        code <- names(attr(v, "blockCode"))
        grp <- paste0(if (part == "zi") "zi.", names(vc[[part]])[j])
        M <- matrix(as.numeric(v), nrow(v), dimnames = list(rownames(v), colnames(v)))
        if (identical(code, "diag") || nrow(M) == 1) {
          parts[[grp]] <- pm_vc_entries(stats::setNames(list(M), grp), stats::setNames(list(TRUE), grp))
        } else if (code %in% c("ar1", "hetar1", "cs", "homcs", "toep")) {
          vars <- diag(M)
          if (code %in% c("ar1", "homcs")) vars <- vars[1]
          vnm <- if (length(vars) == 1) paste0(grp, ".var") else paste(grp, rownames(M), sep = ".")
          R <- stats::cov2cor(M)
          rho <- if (code == "toep") R[-1, 1] else R[2, 1]
          rnm <- if (code == "toep") paste0("cor.", grp, ".lag", seq_along(rho)) else paste0("cor.", grp)
          parts[[grp]] <- structure(c(stats::setNames(vars, vnm), stats::setNames(rho, rnm)),
                                    type = c(rep("var", length(vars)), rep("cor", length(rho))), pair = list())
        } else {
          parts[[grp]] <- pm_vc_entries(stats::setNames(list(M), grp))
        }
      }
    }
    if (gaus) parts$residual <- structure(c(residual = stats::sigma(m)^2), type = "var", pair = list())
    do.call(pm_vc_c, unname(parts))
  }
  pm_vc_delta(fun, pf[pos], V[use, use, drop = FALSE], ci, level)
}

## ---- the error term of lm, aov and glm ----

## A residual variance (dispersion) on df degrees of freedom: SE s2 sqrt(2/df);
## Satterthwaite's interval then has nu = df: the (exact) chi-square interval.
pm_vc_chisq <- function(s2, df, nm, ci, level) {
  se <- if (df > 0) s2 * sqrt(2/df) else NA_real_
  pm_vc_table(stats::setNames(s2, nm), se = se, type = "var", ci = ci, level = level)
}

pm_vc_lm <- function(model, ci, level) {
  df <- stats::df.residual(model)
  pm_vc_chisq(sum(stats::weighted.residuals(model)^2)/df, df, "residual", ci, level)
}

pm_vc_glm <- function(model, ci, level) {
  fam <- stats::family(model)$family
  if (inherits(model, "negbin")) {
    return(pm_vc_table(c(theta = model$theta), se = model$SE.theta, type = "pos", ci = ci, level = level))
  }
  if (pm_fixed_dispersion(model)) {
    out <- pm_vc_table(c(dispersion = 1), se = NA_real_, type = "par", ci = ci, level = level,
                       lower = NA_real_, upper = NA_real_)
    attr(out, "note") <- paste0("The dispersion of a ", fam, " model is fixed at 1.")
    return(out)
  }
  df <- stats::df.residual(model)
  phi <- sum(stats::residuals(model, type = "pearson")^2)/df
  pm_vc_chisq(phi, df, if (fam == "gaussian") "residual" else "dispersion", ci, level)
}

## ---- gam ----

pm_vc_gam <- function(model, ci, level) {
  vc <- NULL
  utils::capture.output(vc <- tryCatch(mgcv::gam.vcomp(model, conf.lev = level), error = function(e) NULL))
  if (is.null(vc)) stop("mgcv::gam.vcomp() failed for this model.", call. = FALSE)
  vc <- pm_gam_vcomp_matrix(vc)
  rn <- rownames(vc)
  rn[rn == "scale"] <- "residual"
  est <- stats::setNames(vc[, "std.dev"]^2, rn)
  # SE of a variance from the width of gam.vcomp()'s Wald interval of log(sd)
  z <- stats::qnorm(1 - (1 - level)/2)
  se_log_sd <- (log(vc[, "upper"]) - log(vc[, "lower"]))/(2 * z)
  pm_vc_table(est, se = est * 2 * se_log_sd, type = rep("var", length(est)), ci = ci, level = level)
}

## gam.vcomp() returns a matrix (std.dev, lower, upper) or, without REML/ML,
## a vector; since mgcv 1.9-3 either one sits in element 'vc' of a list
pm_gam_vcomp_matrix <- function(vc) {
  if (is.list(vc)) vc <- vc$vc
  if (is.matrix(vc)) return(unclass(vc)[, c("std.dev", "lower", "upper"), drop = FALSE])
  out <- cbind(std.dev = as.numeric(vc), lower = NA_real_, upper = NA_real_)
  rownames(out) <- names(vc)
  out
}

## ---- aov with Error() ----

## Analysis of variance estimates of the stratum variance components: the
## stratum mean squares (R/aovlist.R, with missing values estimated) and
## their expectations E(MS_s) = sum_j c_sj sigma2_j, c_sj = tr(R_s Z_j Z_j')/df_s
## (R_s the residual projector of stratum s, Z_j the indicators of error term
## j; Within: Z = I), solved for the components whose own stratum has
## residual df. SEs from var(MS_s) = 2 MS_s^2/df_s (independent strata);
## Satterthwaite's chi-square intervals for positive estimates.
pm_vc_aov <- function(model, ci, level) {
  carrier <- pm_aov_carrier(model)
  res <- suppressMessages(pm_aov_strata(carrier))
  strata <- res$strata
  ef <- attr(carrier, "error_formula")
  comp <- if (!is.null(attr(strata, "missing"))) tryCatch(pm_aov_complete(carrier), error = function(e) NULL)
  dat <- if (!is.null(comp)) comp$full else stats::model.frame(model)
  sq <- pm_aov_qr(ef, dat)
  X <- stats::model.matrix(stats::delete.response(stats::terms(carrier)), dat, contrasts.arg = carrier$contrasts)
  b <- stats::coef(carrier)
  X <- X[, intersect(colnames(X), names(b)[!is.na(b)]), drop = FALSE]
  S <- pm_aov_bases(sq$qe, sq$grp, X)
  st_names <- c("(Intercept)", attr(stats::terms(stats::as.formula(paste("~", paste(deparse(ef), collapse = "")))),
                                   "term.labels"), "Within")[sort(unique(sq$grp)) + 1L]
  comps <- setdiff(st_names, "(Intercept)")
  Z <- lapply(comps, function(t) {
    if (t == "Within") return(NULL)
    f <- interaction(dat[all.vars(str2lang(t))], drop = TRUE)
    m <- matrix(0, nrow(dat), nlevels(f))
    m[cbind(seq_len(nrow(dat)), as.integer(f))] <- 1
    qr.qty(sq$qe, m)
  })
  ms <- stats::setNames(strata$ms, strata$stratum)
  dfs <- stats::setNames(strata$df, strata$stratum)
  rows <- comps[comps %in% names(dfs) & dfs[comps] > 0]
  C <- matrix(0, length(rows), length(rows), dimnames = list(rows, rows))
  for (s in rows) {
    i <- match(s, st_names)
    for (j in rows) {
      C[s, j] <- if (j == "Within") 1 else {
        T <- Z[[match(j, comps)]][S[[i]]$sel, , drop = FALSE]
        sum((T - S[[i]]$U %*% crossprod(S[[i]]$U, T))^2) / S[[i]]$df
      }
    }
  }
  A <- solve(C)
  est <- drop(A %*% ms[rows])
  vms <- 2 * ms[rows]^2 / dfs[rows]
  V <- A %*% diag(vms, nrow = length(rows)) %*% t(A)
  # (Satterthwaite's interval from the mean squares, nu = (sum a_i MS_i)^2 /
  # sum (a_i MS_i)^2/df_i, is the common rule: nu = 2 (est/SE)^2)
  nm <- ifelse(rows == "Within", "residual", rows)
  ord <- c(setdiff(seq_along(rows), which(rows == "Within")), which(rows == "Within"))
  out <- pm_vc_table(stats::setNames(est, nm)[ord], V[ord, ord, drop = FALSE], type = rep("var", length(rows)),
                     ci = ci, level = level)
  if (any(est < 0)) attr(out, "note") <- "A negative estimate: that stratum varies less than the one below it."
  out
}
