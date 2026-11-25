## Back-transformation of means: the transformation of the response found in
## the model formula, Box-Cox and Yeo-Johnson transformations, and the mean
## of the back-transformed response (bias adjustment).

#' Box-Cox and Yeo-Johnson Transformations of a Response
#'
#' @description
#' Power transformations to use for the response in a model formula, e.g.
#' \code{lm(pm_boxcox(y, 0.25) ~ A)}. \code{\link{predictmeans}} and \code{\link{covariatemeans}} recognise them
#' (as they recognise \code{log(y)}, \code{log(y + c)}, \code{log10(y)}, \code{log2(y)}, \code{log1p(y)},
#' \code{sqrt(y)}, \code{y^p}, \code{1/y}, \code{qlogis(y)} and \code{asin(sqrt(y))}) and give the means
#' back-transformed to the original scale as well; \code{bias_adjust = TRUE} there gives the mean rather than the
#' median of the original scale.
#' @usage
#' pm_boxcox(y, lambda)
#' pm_yeojohnson(y, lambda)
#' @param y A numeric vector; positive for \code{pm_boxcox()}, any values for \code{pm_yeojohnson()}.
#' @param lambda The power: one number. Box-Cox: 1 leaves the shape of \code{y} as it is, 0.5 is like a square
#'   root, 0 is the logarithm, -1 like a reciprocal. \code{\link{pm_check}} reports the value the data suggest for
#'   an \code{lm} fit.
#' @details
#' Box-Cox: \eqn{(y^\lambda - 1)/\lambda}, and \eqn{\log y} for \eqn{\lambda = 0}. Yeo-Johnson (2000): the same for
#' \eqn{y + 1} when \eqn{y \ge 0}, and \eqn{-((1 - y)^{2 - \lambda} - 1)/(2 - \lambda)} (\eqn{-\log(1 - y)} for
#' \eqn{\lambda = 2}) when \eqn{y < 0}, so zero and negative values are allowed.
#'
#' Choose a transformation from the residuals of the model, not from the raw response, and prefer a generalized
#' linear (mixed) model for counts and proportions.
#' @return The transformed values.
#' @references
#' Box, G. E. P. and Cox, D. R. (1964), \emph{An analysis of transformations}, Journal of the Royal Statistical
#' Society B, 26(2), 211-252.
#'
#' Yeo, I.-K. and Johnson, R. A. (2000), \emph{A new family of power transformations to improve normality or
#' symmetry}, Biometrika, 87(4), 954-959.
#' @seealso \code{\link{predictmeans}}, \code{\link{pm_check}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lm(pm_boxcox(yield, 0.5) ~ Block + nitro, data = Oats)
#'   predictmeans(fm, "nitro", plot = FALSE)$means
#' @name pm_boxcox
#' @aliases pm_yeojohnson
NULL

#' @export
pm_boxcox <- function(y, lambda) {
  pm_check_lambda(lambda)
  if (any(y <= 0, na.rm = TRUE)) stop("pm_boxcox() needs positive values; pm_yeojohnson() takes any.", call. = FALSE)
  if (abs(lambda) < 1e-12) log(y) else (y^lambda - 1)/lambda
}
pm_boxcox_inv <- function(lambda) {
  if (abs(lambda) < 1e-12) return(exp)
  function(x) pmax(lambda*x + 1, 0)^(1/lambda)
}

#' @export
pm_yeojohnson <- function(y, lambda) {
  pm_check_lambda(lambda)
  out <- y
  pos <- !is.na(y) & y >= 0
  neg <- !is.na(y) & y < 0
  out[pos] <- if (abs(lambda) < 1e-12) log1p(y[pos]) else ((y[pos] + 1)^lambda - 1)/lambda
  out[neg] <- if (abs(lambda - 2) < 1e-12) -log1p(-y[neg]) else -((1 - y[neg])^(2 - lambda) - 1)/(2 - lambda)
  out
}
pm_yeojohnson_inv <- function(lambda) {
  function(x) {
    out <- x
    pos <- !is.na(x) & x >= 0
    neg <- !is.na(x) & x < 0
    out[pos] <- if (abs(lambda) < 1e-12) expm1(x[pos]) else pmax(lambda*x[pos] + 1, 0)^(1/lambda) - 1
    out[neg] <- if (abs(lambda - 2) < 1e-12) -expm1(-x[neg]) else 1 - pmax(1 - (2 - lambda)*x[neg], 0)^(1/(2 - lambda))
    out
  }
}
pm_check_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda)) stop("'lambda' must be one number.", call. = FALSE)
  invisible(lambda)
}

## The transformation of the response in the model formula and its inverse:
## NULL when the response is not transformed, or not in one of the forms
## known here. Returns list(fun (the inverse; exp itself for log(y), so that
## differences are also given as ratios), label (the response as written)).
pm_response_trans <- function(model) {
  tt <- tryCatch(stats::terms(model), error = function(e) NULL)
  if (is.null(tt) || !isTRUE(attr(tt, "response") > 0)) return(NULL)
  e <- attr(tt, "variables")[[attr(tt, "response") + 1]]
  env <- attr(tt, ".Environment")
  if (is.null(env)) env <- globalenv()
  inv <- tryCatch(pm_parse_trans(e, env), error = function(err) NULL)
  if (is.null(inv)) return(NULL)
  list(fun = inv, label = paste(deparse(e), collapse = ""))
}

## The inverse of a response expression: f(y) or f(y + c), f one of log (any
## base), log10, log2, log1p, sqrt, a power, a reciprocal, qlogis, asin(sqrt()),
## pm_boxcox, pm_yeojohnson; I() and brackets are looked through. NULL when
## the expression is not of these forms.
pm_parse_trans <- function(e, env = globalenv()) {
  num <- function(x) {
    v <- tryCatch(eval(x, env), error = function(err) NULL)
    if (is.numeric(v) && length(v) == 1 && is.finite(v)) v else NULL
  }
  # y, (y), y + c, c + y, y - c: the constant added to the response
  offset <- function(a) {
    while (is.call(a) && identical(a[[1]], as.name("("))) a <- a[[2]]
    if (is.name(a)) return(0)
    if (is.call(a) && length(a) == 3 && as.character(a[[1]]) %in% c("+", "-")) {
      sgn <- if (identical(as.character(a[[1]]), "-")) -1 else 1
      if (is.name(a[[2]]) && !is.null(num(a[[3]])) && !length(all.vars(a[[3]]))) return(sgn*num(a[[3]]))
      if (sgn == 1 && is.name(a[[3]]) && !is.null(num(a[[2]])) && !length(all.vars(a[[2]]))) return(num(a[[2]]))
    }
    NULL
  }
  while (is.call(e) && as.character(e[[1]])[1] %in% c("(", "I")) e <- e[[2]]
  if (!is.call(e)) return(NULL)
  fn <- as.character(e[[1]])
  fn <- fn[length(fn)]                                  # (pkg::fun)
  off <- offset(e[[2]])
  minus <- function(f) if (isTRUE(off == 0)) f else function(x) f(x) - off
  if (fn %in% c("pm_boxcox", "pm_yeojohnson") && length(e) == 3 && !is.null(off)) {
    lambda <- num(e[[3]])
    if (is.null(lambda)) return(NULL)
    return(minus(if (fn == "pm_boxcox") pm_boxcox_inv(lambda) else pm_yeojohnson_inv(lambda)))
  }
  if (fn == "log" && !is.null(off)) {
    if (length(e) == 2) return(if (off == 0) exp else function(x) exp(x) - off)
    base <- num(e[[3]])
    if (is.null(base)) return(NULL)
    return(minus(function(x) base^x))
  }
  if (fn == "asin" && length(e) == 2) {
    s <- e[[2]]
    if (is.call(s) && identical(s[[1]], as.name("sqrt")) && isTRUE(offset(s[[2]]) == 0)) return(function(x) sin(x)^2)
    return(NULL)
  }
  if (fn == "^" && length(e) == 3 && !is.null(off)) {
    p <- num(e[[3]])
    if (is.null(p) || p == 0 || length(all.vars(e[[3]]))) return(NULL)
    return(minus(function(x) pmax(x, 0)^(1/p)))
  }
  if (fn == "/" && length(e) == 3) {                    # k/y, k/(y + c)
    k <- num(e[[2]])
    off <- offset(e[[3]])
    if (is.null(k) || length(all.vars(e[[2]])) || is.null(off)) return(NULL)
    return(minus(function(x) k/x))
  }
  if (length(e) != 2 || is.null(off)) return(NULL)
  switch(fn,
         log10 = minus(function(x) 10^x),
         log2 = minus(function(x) 2^x),
         log1p = minus(expm1),
         sqrt = minus(function(x) x^2),
         qlogis = minus(stats::plogis),
         logit = minus(stats::plogis),
         NULL)
}

## Nodes and weights of Gauss-Hermite quadrature for a standard normal
## variable (Golub and Welsch): E f(Z) = sum(w f(z)).
pm_gauss_hermite <- function(n = 31) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}

## The mean of the back-transformed response, E h(eta + u) with
## u ~ N(0, s2), and its derivative in eta (for the delta-method SE).
##   "exact"   Gauss-Hermite quadrature: exp(eta + s2/2) for a log, eta^2 +
##             s2 for a square root, and accurate for any smooth h
##   "taylor"  h + s2/2 h'' (second order, as emmeans' bias.adjust): too low
##             for a log when s2 is large
## bf: pm_back_funs(). Returns list(h, d).
pm_bias_funs <- function(bf, s2, method = c("exact", "taylor")) {
  method <- match.arg(method)
  taylor <- list(h = function(x) bf$h(x) + s2/2 * bf$d2(x), d = function(x) bf$d1(x) + s2/2 * bf$d3(x))
  if (method == "taylor" || !isTRUE(s2 > 0)) return(taylor)
  gh <- pm_gauss_hermite(31)
  gh2 <- pm_gauss_hermite(21)
  s <- sqrt(s2)
  avg <- function(f, q) function(x) {
    vapply(x, function(xi) if (is.na(xi)) NA_real_ else sum(q$w * f(xi + s*q$z)), numeric(1))
  }
  # the integral exists and is found when two rules agree; not for a
  # back-transformation with a pole within the random variation (a
  # reciprocal near zero), where the mean does not exist
  settled <- function(a, b, x) {
    ok <- is.finite(x)
    all(is.finite(a[ok]) & is.finite(b[ok])) && all(abs(a[ok] - b[ok]) <= 1e-4 * pmax(1, abs(a[ok])))
  }
  list(h = function(x) {
         v <- avg(bf$h, gh)(x)
         if (!settled(v, avg(bf$h, gh2)(x), x)) {
           warning("bias_adjust: the mean of the back-transformed response is not finite over the range of the ",
                   "random variation; the second-order formula is used.", call. = FALSE)
           return(taylor$h(x))
         }
         v
       },
       d = function(x) {
         v <- avg(bf$d1, gh)(x)
         if (!settled(v, avg(bf$d1, gh2)(x), x) || !settled(avg(bf$h, gh)(x), avg(bf$h, gh2)(x), x)) return(taylor$d(x))
         v
       })
}

## What the back-transformed means are (settings$back_transform): NULL when
## there are none, else list(by, estimate):
##   a GLM: the mean; a GLMM: the response of a typical cluster (random
##   effects zero); a transformed response: the median of the original scale
##   (the geometric mean for a log); with bias_adjust: the mean.
pm_back_label <- function(means, by, link_scale, model, bias_adjusted, average = "link") {
  if (is.null(means$response)) return(NULL)
  random <- inherits(model, c("merMod", "glmmTMB", "lme"))
  user_on_link <- !link_scale && !is.null(pm_link_trans(model))
  estimate <- if (inherits(model, "pm_catfit")) {
    "probability of the class"
  } else if (bias_adjusted) {
    "mean (bias adjusted)"
  } else if (user_on_link) {
    "the user's 'trans' of the link-scale mean (not the model's inverse link)"
  } else if (link_scale && !random) {
    "mean"
  } else if (link_scale) {
    "typical cluster (random effects at zero; bias_adjust = TRUE gives the mean over the clusters)"
  } else {
    "median of the original scale (bias_adjust = TRUE gives the mean)"
  }
  list(by = if (is.na(by)) "trans" else by, estimate = estimate, average = average)
}

## The Box-Cox power the residuals of a linear model suggest: the profile
## likelihood estimate of lambda and its 95% interval (as MASS::boxcox()), for an
## lm fit with a positive, untransformed response. NULL otherwise.
pm_boxcox_profile <- function(model) {
  if (!inherits(model, "lm") || inherits(model, c("glm", "mlm", "aovlist", "gam"))) return(NULL)
  if (!is.null(pm_response_trans(model))) return(NULL)
  y <- tryCatch(stats::model.response(stats::model.frame(model)), error = function(e) NULL)
  if (!is.numeric(y) || is.matrix(y) || any(y <= 0) || length(unique(y)) < 5) return(NULL)
  if (!is.null(model$weights) || !is.null(model$offset)) return(NULL)
  # profile log-likelihood of lambda: -n/2 log(RSS(lambda)/n) + (lambda - 1) sum(log y)
  qrx <- model$qr
  if (is.null(qrx) || length(y) != nrow(qrx$qr)) return(NULL)
  logy <- log(y)
  n <- length(y)
  lam <- seq(-3, 3, by = 0.01)
  ll <- vapply(lam, function(l) {
    z <- if (abs(l) < 1e-8) logy else (y^l - 1)/l
    -n/2 * log(sum(qr.resid(qrx, z)^2)/n) + (l - 1) * sum(logy)
  }, numeric(1))
  if (any(!is.finite(ll))) return(NULL)
  bc <- list(x = lam, y = ll)
  best <- which.max(bc$y)
  inside <- bc$x[bc$y >= bc$y[best] - stats::qchisq(0.95, 1)/2]
  lo <- min(inside)
  hi <- max(inside)
  # the usual powers, when the interval holds one: the closest to the estimate
  usual <- c(-1, -0.5, 0, 0.5, 1, 2)
  ok <- usual[usual >= lo & usual <= hi]
  pick <- if (length(ok)) ok[which.min(abs(ok - bc$x[best]))] else NA_real_
  name <- c("1/y", "1/sqrt(y)", "log(y)", "sqrt(y)", "none", "y^2")[match(pick, usual)]
  list(lambda = bc$x[best], lower = lo, upper = hi, pick = pick, name = name,
       edge = lo <= -3 + 1e-8 || hi >= 3 - 1e-8)
}
