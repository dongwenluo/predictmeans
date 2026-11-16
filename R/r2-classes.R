## Piepho's (2023) adjusted R2 for every model class of predictmeans().
##
## One calculation, pm_r2_shares(), for all classes. With V the covariance of
## the data (of the working response for a GLM), P = I - J/n and n
## observations:
##   theta_V     = tr(P V)/(n - 1)                    average variance
##   theta_Xbeta = (b'X'PXb - tr(X'PX cov(b)))/(n - 1)  variance of the fixed
##                 part, less the bias from estimating b
##   Fixed  = theta_Xbeta/(theta_Xbeta + theta_V)
##   Random = tr(P Z G Z')/(n - 1)/(theta_Xbeta + theta_V)
## For lm (V = s^2 I) Fixed is the usual adjusted R2. Each class only says
## what its X, b, cov(b), residual and random parts are (pm_r2_parts());
## lmer, glmer and glmmTMB keep pm_r2_piepho() of R/R2_glmm.R (the same
## formulas, checked equal for lme and lmer in tests/testthat/test-r2-model.R).

## tr(A P) for a diagonal A (its diagonal v), a product B B' and a matrix S
pm_trP_diag <- function(v, n) (if (length(v) == 1) n * v else sum(v)) * (1 - 1/n)
pm_trP_BB <- function(B) {
  if (!length(B)) return(0)
  sum(B^2) - sum(Matrix::colSums(B)^2)/nrow(B)
}
pm_trP_mat <- function(S) sum(Matrix::diag(S)) - sum(S)/nrow(S)

## parts: X, beta, covB; tr_resid = tr(P R); tr_random: named, tr(P Z_k G_kk
## Z_k') of each random component; tr_random_total: of all of them with
## their covariances (default: the sum).
pm_r2_shares <- function(p) {
  n <- nrow(p$X)
  tr_random <- if (is.null(p$tr_random)) numeric(0) else p$tr_random
  tr_rt <- if (is.null(p$tr_random_total)) sum(tr_random) else p$tr_random_total
  XtPX <- crossprod(p$X) - tcrossprod(colSums(p$X))/n
  theta_Xbeta <- (as.numeric(crossprod(p$beta, XtPX %*% p$beta)) - sum(XtPX * as.matrix(p$covB)))/(n - 1)
  full <- theta_Xbeta + (tr_rt + p$tr_resid)/(n - 1)
  fixed <- theta_Xbeta/full
  random <- tr_rt/(n - 1)/full
  c(Total = fixed + random, Fixed = fixed, Random = random, tr_random/(n - 1)/full)
}

pm_r2_parts <- function(model) {
  if (inherits(model, "aovlist")) return(pm_r2_parts_aovlist(model))
  if (inherits(model, "gam")) return(pm_r2_parts_gam(model))
  if (inherits(model, "glm")) return(pm_r2_parts_glm(model))
  if (inherits(model, "mlm")) stop("R2_model() does not handle models with several responses.", call. = FALSE)
  if (inherits(model, "lm")) return(pm_r2_parts_lm(model))
  if (inherits(model, c("lme", "gls"))) return(pm_r2_parts_nlme(model))
  stop("R2_model() handles lm, aov (also with Error()), glm, gls, lme, lmer, glmer, glmmTMB, gam and semireg models.",
       call. = FALSE)
}

pm_r2_parts_lm <- function(model) {
  b <- model$coefficients
  est <- !is.na(b)
  X <- stats::model.matrix(model)[, est, drop = FALSE]
  w <- model$weights
  keep <- if (is.null(w)) rep(TRUE, nrow(X)) else w > 0
  s2 <- stats::sigma(model)^2
  list(X = X[keep, , drop = FALSE], beta = b[est], covB = stats::vcov(model, complete = FALSE),
       tr_resid = pm_trP_diag(if (is.null(w)) s2 else s2/w[keep], sum(keep)))
}

## The latent-scale variance of a binary response (as pm_r2_piepho()); for
## other links, and other responses, the working variance.
pm_r2_resvar <- function(wk) {
  if (wk$binary && wk$link %in% c("logit", "probit", "cloglog")) {
    return(switch(wk$link, logit = pi^2/3, probit = 1, cloglog = pi^2/6))
  }
  wk$resvar
}

pm_r2_parts_glm <- function(model) {
  wk <- pm_working_glm(model)
  keep <- model$prior.weights > 0
  rv <- pm_r2_resvar(wk)
  list(X = wk$X[keep, , drop = FALSE], beta = stats::coef(model)[colnames(wk$X)],
       covB = stats::vcov(model, complete = FALSE),
       tr_resid = pm_trP_diag(if (length(rv) == 1) rv else rv[keep], sum(keep)))
}

## gam: the smooths are part of the mean (Fixed; cov(b) is mgcv's Bayesian
## Vp, so the bias term takes their effective df), except random-effect
## smooths (bs = "re"), which are Random with the variance of gam.vcomp().
pm_r2_parts_gam <- function(model) {
  wk <- tryCatch(pm_working_glm(model), error = function(e) NULL)
  if (is.null(wk) || !all(is.finite(wk$resvar[model$prior.weights > 0]))) {
    stop("R2_model() does not handle this gam family (", stats::family(model)$family, ").", call. = FALSE)
  }
  keep <- model$prior.weights > 0
  X <- wk$X[keep, , drop = FALSE]
  n <- nrow(X)
  b <- stats::coef(model)
  re <- Filter(function(sm) inherits(sm, "random.effect"), model$smooth)
  tr_random <- numeric(0)
  rcols <- integer(0)
  if (length(re)) {
    vc <- tryCatch({
      utils::capture.output(v <- mgcv::gam.vcomp(model))
      pm_gam_vcomp_matrix(v)
    }, error = function(e) NULL)
    for (sm in re) {
      if (is.null(vc) || !sm$label %in% rownames(vc)) next
      cols <- sm$first.para:sm$last.para
      tr_random <- c(tr_random, stats::setNames(vc[sm$label, "std.dev"]^2 * pm_trP_BB(X[, cols, drop = FALSE]), sm$label))
      rcols <- c(rcols, cols)
    }
  }
  fx <- setdiff(seq_along(b), rcols)
  rv <- pm_r2_resvar(wk)
  list(X = X[, fx, drop = FALSE], beta = b[fx], covB = model$Vp[fx, fx, drop = FALSE],
       tr_resid = pm_trP_diag(if (length(rv) == 1) rv else rv[keep], n), tr_random = tr_random)
}

## gls and lme: the covariance of the data from pm_nlme_sigma() (any
## correlation and variance structure); the random effects of lme by level
## and component.
pm_r2_parts_nlme <- function(model) {
  dat <- pm_get_data(model)
  b <- if (inherits(model, "lme")) nlme::fixef(model) else stats::coef(model)
  X <- stats::model.matrix(stats::formula(model), dat)[, names(b), drop = FALSE]
  n <- nrow(X)
  Srel <- pm_nlme_sigma(model)
  s2 <- model$sigma^2
  tr_all <- s2 * pm_trP_mat(Srel(stats::coef(model$modelStruct)))
  tr_random <- numeric(0)
  tr_rt <- 0
  Zs <- attr(Srel, "Zs")
  if (!is.null(Zs)) {
    D <- nlme::pdMatrix(model$modelStruct$reStruct)
    for (lv in names(Zs)) {
      Z <- Zs[[lv]]$Z
      Dl <- as.matrix(D[[lv]])
      q <- nrow(Dl)
      tr_rt <- tr_rt + s2 * pm_trP_BB(Z %*% Matrix::kronecker(Zs[[lv]]$I, Matrix::Matrix(pm_chol_lower(Dl))))
      for (k in seq_len(q)) {
        nm <- sub("\\.\\(Intercept\\)$", "", paste(lv, colnames(Dl)[k], sep = "."))
        tr_random <- c(tr_random, stats::setNames(
          s2 * Dl[k, k] * pm_trP_BB(Z[, seq(k, ncol(Z), by = q), drop = FALSE]), nm))
      }
    }
  }
  list(X = X, beta = b, covB = as.matrix(stats::vcov(model)), tr_resid = tr_all - tr_rt, tr_random = tr_random,
       tr_random_total = tr_rt)
}

## aov with Error(): the estimates combined over the strata and the stratum
## variance components of varcomp() (the analysis of variance estimates,
## which may be negative).
pm_r2_parts_aovlist <- function(model) {
  carrier <- pm_aov_carrier(model)
  res <- suppressMessages(pm_aov_strata(carrier))
  b <- res$mp$coef
  X <- stats::model.matrix(carrier)[, names(b), drop = FALSE]
  n <- nrow(X)
  vc <- suppressMessages(varcomp(model))
  dat <- stats::model.frame(model)
  comps <- setdiff(rownames(vc), "residual")
  tr_random <- vapply(comps, function(t) {
    f <- interaction(dat[all.vars(str2lang(t))], drop = TRUE)
    vc[t, "vcov"] * (sum(tabulate(f)) - sum(tabulate(f)^2)/n)
  }, numeric(1))
  list(X = X, beta = b, covB = res$mp$vcov, tr_resid = pm_trP_diag(vc["residual", "vcov"], n),
       tr_random = tr_random)
}
