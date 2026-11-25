## Cumulative link models: MASS::polr() and ordinal::clm() (decisions D1-D6 of
## dev/review-classes.md, 2026-11-25).
##
## Both fit P(Y <= j) = F(theta_j - eta), eta = x'b without an intercept, cut
## points theta = J alpha (J = I for polr and for clm's flexible thresholds;
## clm's tJac otherwise). The means are on the latent scale, x'b - mean(theta)
## (as emmeans' mode = "latent"): the coefficients are given as
## ((Intercept) = -mean(theta), b), so the model's own design (with its
## intercept column) gives the means, and the covariance is that of this
## linear map of (b, alpha). Differences of latent means are log cumulative
## odds ratios (logit link), the same for every cut point. Wald tests are z
## tests (pm_df_rule(), as fixed-dispersion GLMs).

## b, alpha, J, their joint covariance (order b, alpha), the link and the
## class labels of a polr or clm fit.
pm_ord_parts <- function(model) {
  if (inherits(model, "polr")) {
    b <- stats::coef(model)
    alpha <- model$zeta
    J <- diag(length(alpha))
    link <- model$method
    lev <- model$lev
    V <- stats::vcov(model)
  } else {
    if (length(model$nom.terms) || length(model$S.terms) || length(model$zeta)) {
      stop("predictmeans handles clm fits with location effects only; nominal and scale effects ",
           "(nominal = , scale = ) are not supported.", call. = FALSE)
    }
    b <- model$beta
    alpha <- model$alpha
    J <- model$tJac
    link <- model$link
    lev <- model$y.levels
    V <- stats::vcov(model)
  }
  if (is.null(V) || anyNA(V)) {
    stop("The ", class(model)[1], " fit has no covariance matrix", if (inherits(model, "polr"))
      ": refit with polr(..., Hess = TRUE)." else ".", call. = FALSE)
  }
  nm <- c(names(b), names(alpha))
  list(b = b, alpha = alpha, J = J, V = V[nm, nm, drop = FALSE], link = link, levels = lev)
}

## The coefficients on the latent scale and their covariance: ((Intercept) =
## -mean(theta), b), A (b, alpha) with A the linear map.
pm_ord_latent <- function(model) {
  p <- pm_ord_parts(model)
  nb <- length(p$b)
  A <- rbind(c(rep(0, nb), -colMeans(p$J)), cbind(diag(nb), matrix(0, nb, length(p$alpha))))
  coefs <- c("(Intercept)" = sum(A[1, ] * c(p$b, p$alpha)), p$b)
  V <- A %*% p$V %*% t(A)
  dimnames(V) <- list(names(coefs), names(coefs))
  list(coef = coefs, vcov = V)
}

#' @exportS3Method pm_coef polr
pm_coef.polr <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_coef clm
pm_coef.clm <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_vcov polr
pm_vcov.polr <- function(model, ...) pm_ord_latent(model)$vcov
#' @exportS3Method pm_vcov clm
pm_vcov.clm <- function(model, ...) pm_ord_latent(model)$vcov
## no back-transformation of the latent means (class probabilities: mean_type = "prob")
#' @exportS3Method pm_linkinv polr
pm_linkinv.polr <- function(model, ...) NULL
#' @exportS3Method pm_linkinv clm
pm_linkinv.clm <- function(model, ...) NULL

## z tests: the latent scale has no dispersion to estimate (D3; as emmeans and
## the fixed-dispersion GLM rule, dev/sim-df-glm/).
#' @exportS3Method NULL
pm_df_rule.polr <- function(model, df_method, kr_info) {
  if (df_method != "auto") {
    warning("'df_method = \"", df_method, "\"' does not apply to ", class(model)[1], " models; ignored.", call. = FALSE)
  }
  list(df_type = "asymptotic", vcov_adj = NULL)
}
#' @exportS3Method NULL
pm_df_rule.clm <- pm_df_rule.polr

## Class probabilities P(Y = j) (mean_type = "prob") or the mean class
## sum_j j P(Y = j) (mean_type = "mean_class") at each mean: rows of K (the
## latent means; averages over other factors are taken on the latent scale,
## D2), with delta-method SEs from the covariance of (b, alpha). Limits:
## probabilities on the logit scale (they stay in (0, 1)), the mean class on
## its own scale. emmeans' mode = "prob" / "mean.class" give the same values
## for the cells of the model (they average probabilities over other factors).
pm_ord_class <- function(model, K, what = c("prob", "mean_class"), conf_level = 0.95) {
  what <- match.arg(what)
  p <- pm_ord_parts(model)
  Fd <- switch(p$link, logistic = , logit = list(stats::plogis, stats::dlogis),
               probit = list(stats::pnorm, stats::dnorm),
               cloglog = list(function(x) 1 - exp(-exp(x)), function(x) exp(x - exp(x))),
               loglog = list(function(x) exp(-exp(-x)), function(x) exp(-x - exp(-x))),
               cauchit = list(stats::pcauchy, stats::dcauchy),
               stop("link '", p$link, "' is not supported.", call. = FALSE))
  X <- as.matrix(K)[, names(p$b), drop = FALSE]
  th <- as.vector(p$J %*% p$alpha)
  nj <- length(th)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  rows <- lapply(seq_len(nrow(X)), function(k) {
    x <- X[k, ]
    if (anyNA(x)) return(NULL)
    zz <- th - sum(x * p$b)
    f <- Fd[[2]](zz)
    cum <- c(0, Fd[[1]](zz), 1)
    # gradients of P(Y <= j) with respect to (b, alpha), j = 0, ..., J
    G <- rbind(0, cbind(-outer(f, x), f * p$J), 0)
    P <- diff(cum)
    GP <- G[-1, , drop = FALSE] - G[-(nj + 2), , drop = FALSE]
    if (what == "mean_class") {
      g <- colSums(seq_along(P) * GP)
      m <- sum(seq_along(P) * P)
      se <- sqrt(sum((g %*% p$V) * g))
      return(data.frame(row = k, estimate = m, std.error = se, conf.low = m - z * se, conf.high = m + z * se))
    }
    se <- sqrt(rowSums((GP %*% p$V) * GP))
    lg <- stats::qlogis(P)
    sl <- se / (P * (1 - P))                       # SE of logit(P)
    data.frame(row = k, class = factor(p$levels, levels = p$levels), estimate = P, std.error = se,
               conf.low = stats::plogis(lg - z * sl), conf.high = stats::plogis(lg + z * sl))
  })
  do.call(rbind, rows)
}

## Add the class table to a predictmeans() result: $probs, with the factor
## columns of the means in front.
pm_ord_attach <- function(res, model, K, what, conf_level) {
  tab <- pm_ord_class(model, K, what, conf_level)
  means <- res$means
  fac <- setdiff(names(means), c("estimate", "std.error", "df", "conf.low", "conf.high", "comparison.low",
                                 "comparison.high", "response", "response.conf.low", "response.conf.high", "group",
                                 "n", "estimable"))
  out <- cbind(means[tab$row, fac, drop = FALSE], tab[setdiff(names(tab), "row")])
  rownames(out) <- NULL
  attr(out, "what") <- what
  res$probs <- out
  res
}
