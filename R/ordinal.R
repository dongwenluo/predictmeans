## Cumulative link models: MASS::polr(), ordinal::clm() and ordinal::clmm()
## (decisions D1-D6 of dev/review-classes.md, 2026-11-25). clmm: the fixed part,
## random effects at 0 (a typical subject, as emmeans; D4); its covariance of
## (b, alpha) is the block of vcov() without the random-effect parameters.

## The classes handled here (one list for every check in the package).
pm_ord_classes <- c("polr", "clm", "clmm")
##
## Both fit P(Y <= j) = F(theta_j - eta), eta = x'b without an intercept, cut
## points theta = J alpha (J = I for polr and for clm's flexible thresholds;
## clm's tJac otherwise). The means are on the latent scale, x'b - mean(theta)
## (as emmeans' mode = "latent"): the coefficients are given as
## ((Intercept) = -mean(theta), b), so the model's own design (with its
## intercept column) gives the means, and the covariance is that of this
## linear map of (b, alpha). A difference of latent means (logit link) is the
## log odds ratio of a higher class, log[P(Y > j) / P(Y <= j)], the same for
## every cut point j. Wald tests are z tests (pm_df_rule(), as fixed-dispersion
## GLMs). Aliased coefficients are NA in pm_coef() and left out of pm_vcov()
## (the adapter contract, ?pm_coef).

## b (aliased coefficients NA), alpha, J, the joint covariance of the
## estimable b and alpha (in that order), the link and the class labels of a
## polr or clm fit; stops for what is not handled.
pm_ord_parts <- function(model) {
  if (inherits(model, "polr")) {
    if (is.null(model$Hessian)) {
      stop("The polr fit has no Hessian: refit with polr(..., Hess = TRUE).", call. = FALSE)
    }
    b <- stats::coef(model)
    alpha <- model$zeta
    J <- diag(length(alpha))
    link <- model$method
    lev <- model$lev
  } else {
    if (inherits(model, "clmm") && !is.null(model$call$nominal)) {
      stop("clmm() does not fit nominal effects: it ignored 'nominal ='; refit without it.", call. = FALSE)
    }
    if (length(model$nom.terms) || length(model$S.terms) || length(model$zeta)) {
      stop("predictmeans handles clm fits with location effects only; nominal and scale effects ",
           "(nominal = , scale = ) are not supported.", call. = FALSE)
    }
    if (!is.null(model$lambda) || model$link %in% c("Aranda-Ordaz", "log-gamma")) {
      stop("The ", model$link, " link has an estimated parameter (lambda), which the standard errors would treat as ",
           "known; use a logit, probit, cloglog, loglog or cauchit link.", call. = FALSE)
    }
    b <- model$beta
    alpha <- model$alpha
    J <- model$tJac
    link <- model$link
    lev <- model$y.levels
  }
  V <- stats::vcov(model)
  nm <- c(names(b)[!is.na(b)], names(alpha))
  if (is.null(V) || anyNA(V) || !all(nm %in% rownames(V))) {
    stop("The ", class(model)[1], " fit has no usable covariance matrix of its coefficients.", call. = FALSE)
  }
  list(b = b, alpha = alpha, J = J, V = V[nm, nm, drop = FALSE], link = link, levels = lev)
}

## The checks of an ordinal fit and of the options of a call, for every front
## end (predictmeans(), covariatemeans()): pm_ord_parts() stops for fits not
## handled; latent means have no back-transformation and no residual
## variance, so 'trans', 'bias_adjust', 'sigma' (and 'vcov', 'permlist', not
## built for the latent parametrisation) stop: the class scale is mean_type's.
pm_ord_check <- function(model, trans = NULL, bias_adjust = FALSE, sigma = NULL, vcov = NULL, permlist = NULL) {
  pm_ord_parts(model)
  bad <- c(trans = !is.null(trans) && !identical(trans, I), bias_adjust = !isFALSE(bias_adjust),
           sigma = !is.null(sigma), vcov = !is.null(vcov), permlist = !is.null(permlist))
  if (any(bad)) {
    stop("'", paste(names(bad)[bad], collapse = "', '"), "' cannot be used with ", class(model)[1],
         " models: the means are on the latent scale; for the class scale use mean_type = \"prob\" or ",
         "\"mean_class\" in predictmeans().", call. = FALSE)
  }
  invisible(TRUE)
}

## The coefficients on the latent scale and their covariance: ((Intercept) =
## -mean(theta), b), A (b, alpha) with A the linear map (estimable b only).
pm_ord_latent <- function(model) {
  p <- pm_ord_parts(model)
  be <- p$b[!is.na(p$b)]
  nb <- length(be)
  A <- rbind(c(rep(0, nb), -colMeans(p$J)), cbind(diag(nb), matrix(0, nb, length(p$alpha))))
  icpt <- sum(A[1, ] * c(be, p$alpha))
  V <- A %*% p$V %*% t(A)
  dimnames(V) <- list(c("(Intercept)", names(be)), c("(Intercept)", names(be)))
  list(coef = c("(Intercept)" = icpt, p$b), vcov = V)
}

#' @exportS3Method pm_coef polr
pm_coef.polr <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_coef clm
pm_coef.clm <- function(model, ...) pm_ord_latent(model)$coef
#' @exportS3Method pm_vcov polr
pm_vcov.polr <- function(model, ...) pm_ord_latent(model)$vcov
#' @exportS3Method pm_vcov clm
pm_vcov.clm <- function(model, ...) pm_ord_latent(model)$vcov
#' @exportS3Method pm_coef clmm
pm_coef.clmm <- pm_coef.clm
#' @exportS3Method pm_vcov clmm
pm_vcov.clmm <- pm_vcov.clm
## no back-transformation of the latent means (class probabilities: mean_type = "prob")
## the contrasts the fit used (clmm has no model.matrix() method: the default
## would rebuild the design with the current options)
#' @exportS3Method pm_contrasts clmm
pm_contrasts.clmm <- function(model, ...) model$contrasts
#' @exportS3Method pm_linkinv polr
pm_linkinv.polr <- function(model, ...) NULL
#' @exportS3Method pm_linkinv clm
pm_linkinv.clm <- function(model, ...) NULL
#' @exportS3Method pm_linkinv clmm
pm_linkinv.clmm <- function(model, ...) NULL

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
## clmm: inner-outer df by default (pm_df_method(), as GLMMs: z rejected 11% of
## true nulls at 5% with 6 clusters for a between-cluster term); this rule
## (z) is the fallback when they cannot be found.
#' @exportS3Method NULL
pm_df_rule.clmm <- pm_df_rule.polr

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
  be <- p$b[!is.na(p$b)]
  X <- as.matrix(K)[, names(be), drop = FALSE]
  th <- as.vector(p$J %*% p$alpha)
  nj <- length(th)
  z <- pm_ci_crit(1 - conf_level, Inf)                          # z: the df of these models (pm_df_rule)
  rows <- lapply(seq_len(nrow(X)), function(k) {
    x <- X[k, ]
    if (anyNA(x)) return(NULL)
    zz <- th - sum(x * be)
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
      # limits on the logit of (m - 1) / (J - 1), so they stay within the classes 1 to J
      u <- (m - 1) / nj
      lu <- stats::qlogis(u)
      su <- se / nj / (u * (1 - u))
      return(data.frame(row = k, estimate = m, std.error = se, conf.low = 1 + nj * stats::plogis(lu - z * su),
                        conf.high = 1 + nj * stats::plogis(lu + z * su)))
    }
    se <- sqrt(rowSums((GP %*% p$V) * GP))
    lg <- stats::qlogis(P)
    sl <- se / (P * (1 - P))                       # SE of logit(P)
    data.frame(row = k, class = factor(p$levels, levels = p$levels), estimate = P, std.error = se,
               conf.low = stats::plogis(lg - z * sl), conf.high = stats::plogis(lg + z * sl))
  })
  do.call(rbind, rows)
}

## Add the class table to a predictmeans() result: $probs, rows in the order
## of the means (matched by their labels, not position: atvar re-sorts them),
## with the factor columns of modelterm in front.
pm_ord_attach <- function(res, model, K, what, conf_level) {
  means <- res$means
  vars <- intersect(strsplit(res$settings$modelterm, ":")[[1]], names(means))
  key <- do.call(paste, c(lapply(means[vars], as.character), sep = ":"))
  rows <- match(key, rownames(K))
  if (anyNA(rows)) stop("internal: the means do not match the rows of K.", call. = FALSE)
  tab <- pm_ord_class(model, K[rows, , drop = FALSE], what, conf_level)
  out <- cbind(means[tab$row, vars, drop = FALSE], tab[setdiff(names(tab), "row")])
  rownames(out) <- NULL
  attr(out, "what") <- what
  res$probs <- out
  pm_add_advice(res, pm_ordinal_advice(model, res$settings$modelterm, what))
}
