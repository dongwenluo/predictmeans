## Further uses of the permutations of permmodels() (plan step 14):
## standardised effect sizes and a check of the parametric test against the
## permutation null (permstats()), and nonparametric combination of tests
## made with the same permutations (permcombine()).

## Satterthwaite denominator df of the Wald F test of L b = 0 (lmerTest's
## rule for several rows, as in pm_joint()).
pm_satt_Fddf <- function(L, parts) {
  q <- nrow(L)
  if (q == 1) return(unname(pm_satt_ddf(L, parts)))
  eg <- eigen(L %*% parts$vcov_beta %*% t(L), symmetric = TRUE)
  nu <- pm_satt_ddf(crossprod(eg$vectors, L), parts)
  E <- sum((nu / (nu - 2))[nu > 2])
  if (E > q) 2 * E / (E - q) else NA_real_
}

## Denominator df of the usual parametric test of each permutation test:
## the residual df (lm, aov, glm, gls), Satterthwaite (lmer), nlme's df
## (lme), Inf (Wald chi-square; glmer, glmmTMB).
pm_perm_df2 <- function(model, nms, Ls, df_res) {
  if (inherits(model, c("glmerMod", "glmmTMB"))) return(rep(Inf, length(nms)))
  if (inherits(model, "lmerMod")) {
    parts <- tryCatch(suppressWarnings(pm_satt_parts(model)), error = function(e) NULL)
    if (is.null(parts)) return(rep(NA_real_, length(nms)))
    return(vapply(Ls, function(L) tryCatch(pm_satt_Fddf(L, parts), error = function(e) NA_real_), numeric(1)))
  }
  if (inherits(model, "lme")) {
    aT <- stats::anova(model)
    tT <- summary(model)$tTable
    return(vapply(nms, function(k) {
      out <- if (startsWith(k, "coef:")) tT[sub("^coef:", "", k), "DF"] else aT[k, "denDF"]
      if (length(out) == 1 && !is.na(out)) out else NA_real_
    }, numeric(1)))
  }
  rep(df_res, length(nms))
}

## What permstats() and permcombine() need: observed statistics (F for
## terms, |t| for coefficients), the permuted ones (tests x permutations, NA
## for failed refits), the permutation p-values of the tables, numerator and
## denominator df of the parametric test, and the number of observations.
pm_perm_statset <- function(model, obs, T, p, df1, kind, Ls, df_res, n) {
  nms <- names(obs)
  o <- order(startsWith(nms, "coef:"))            # terms first, then coefficients
  nms <- nms[o]
  T <- T[o, , drop = FALSE]
  rownames(T) <- nms
  df2 <- tryCatch(pm_perm_df2(model, nms, Ls[nms], df_res), error = function(e) rep(NA_real_, length(nms)))
  list(obs = unname(obs[o]), T = T, p = unname(p[o]), df1 = unname(df1[o]), df2 = unname(df2),
       kind = unname(kind[o]), n = n, names = nms)
}

pm_perm_get_stats <- function(x, what = "x") {
  st <- if (inherits(x, "pm_perm")) x$stats else if (is.list(x)) x[["perm"]]$stats
  if (is.null(st)) {
    stop(sprintf(paste0("'%s' must be a result of permmodels() with Freedman-Lane permutations ",
                        "(method = \"auto\" or \"freedman-lane\")."), what), call. = FALSE)
  }
  st
}

#' Effect Sizes and Model Checks from Permutations
#'
#' @description
#' Uses the permuted statistics kept by \code{\link{permmodels}} to give, for each term and coefficient, a
#' standardised effect size and a check of the usual parametric test against the permutation distribution.
#'
#' @param x A result of \code{\link{permmodels}} with Freedman-Lane permutations (\code{method = "auto"}, the
#'   default, or \code{"freedman-lane"}).
#' @param level The significance level of the check of the parametric test (default 0.05).
#' @return A data frame with one row per term and coefficient (rows \code{"coef:<name>"}):
#' \item{statistic}{the observed statistic: F for terms, |t| for coefficients (robust ones with
#'   \code{robust = TRUE}).}
#' \item{df1, df2}{numerator and denominator df of the usual parametric test: the residual df for \code{lm},
#'   \code{aov}, \code{glm} and \code{gls}, Satterthwaite df for \code{lmer}, nlme's df for \code{lme}, and
#'   \code{Inf} (Wald chi-square) for \code{glmer} and \code{glmmTMB}.}
#' \item{p.param}{the parametric p-value of the observed statistic with these df.}
#' \item{p.perm}{the permutation p-value (as in the tables of \code{permmodels()}).}
#' \item{ses}{the standardised effect size, (observed - mean of the permuted statistics) / their standard
#'   deviation (Gotelli and McCabe 2002): how far the observed statistic lies from what the permutations give
#'   when the effect is absent, in standard deviations.}
#' \item{param.size}{the proportion of permuted data sets in which the parametric test rejects at
#'   \code{level}. The permuted data sets have no effect for the tested term, so this is the type I error rate of
#'   the parametric test for data like these: close to \code{level} when the parametric reference
#'   distribution suits the data, clearly above it when the parametric test is too liberal here (for example
#'   with unequal variances, overdispersion or few groups), below it when it is conservative. With few
#'   permutations it is rough (its standard error is about \code{sqrt(level (1 - level) / nperm)}).}
#' @details
#' Every permuted data set comes from the model without the tested term (Freedman-Lane), so the permuted
#' statistics show how the statistic behaves when that term has no effect but everything else in the data is
#' kept. Comparing the parametric p-value with the permutation p-value, and \code{param.size} with
#' \code{level}, shows whether the usual test can be trusted for these data. For mixed models the
#' permuted statistics are Wald statistics (unadjusted covariance), so \code{param.size} checks the Wald F
#' with the given df, not the Kenward-Roger test.
#' @references
#' Gotelli, N. J. and McCabe, D. J. (2002). Species co-occurrence: a meta-analysis of J. M. Diamond's
#' assembly rules model. \emph{Ecology}, 83, 2091-2096.
#' @seealso \code{\link{permmodels}}, \code{\link{permcombine}}
#' @examples
#' mtcars$cyl <- factor(mtcars$cyl)
#' fm <- lm(mpg ~ wt + cyl, data = mtcars)
#' pfl <- permmodels(fm, nperm = 999, type = "II", seed = 1, prt = FALSE)
#' permstats(pfl)
#' @usage permstats(x, level = 0.05)
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @name permstats
NULL

permstats <- function(x, level = 0.05) {
  st <- pm_perm_get_stats(x)
  T <- st$T
  ok <- !is.na(T)
  nok <- rowSums(ok)
  mu <- rowSums(ifelse(ok, T, 0)) / nok
  s <- sqrt(rowSums(ifelse(ok, (T - mu)^2, 0)) / (nok - 1))
  pref <- function(stat, df1, df2) stats::pf(stat^ifelse(st$kind == "t", 2, 1), df1, df2, lower.tail = FALSE)
  p_param <- pref(st$obs, st$df1, st$df2)
  pT <- pref(T, st$df1, st$df2)                  # vectors recycle down the rows of T
  size <- rowSums(ifelse(ok, pT <= level, FALSE)) / nok
  out <- data.frame(statistic = st$obs, df1 = st$df1, df2 = st$df2, p.param = p_param, p.perm = st$p,
                    ses = (st$obs - mu) / s, param.size = size, row.names = st$names)
  attr(out, "nperm") <- ncol(T)
  attr(out, "level") <- level
  out
}

#' Combine Permutation Tests of Several Responses
#'
#' @description
#' Nonparametric combination (NPC) of permutation tests made with the same permutations, for example the same
#' model fitted to several responses (traits, time points) measured on the same units: for each term or
#' coefficient, one p-value for the hypothesis that it has no effect on any of the responses, and p-values
#' adjusted for testing it on every response.
#'
#' @param ... Results of \code{\link{permmodels}} with Freedman-Lane permutations (the default), all made with
#'   the same permutations: the same observations (rows of the data) and the same \code{seed} (or the same
#'   \code{perms}) and \code{nperm}. Names given to the arguments name the responses. A single list of such
#'   results can be given instead.
#' @param method How the p-values of the responses are combined: \code{"fisher"} (the default),
#'   \eqn{-2 \sum \log p}{-2 sum log p}, good when several responses have effects; \code{"tippett"}, the smallest
#'   p-value, good when one response has a large effect.
#' @param tests Names of the terms and coefficients (\code{"coef:<name>"}) to combine; the default is every
#'   one that all the results have.
#' @return A data frame with one row per test: the permutation p-value of each response, the combined
#'   p-value (\code{combined}), and the p-values adjusted for the number of responses (columns
#'   \code{adj.<response>}): single-step min-p adjusted p-values (Westfall and Young 1993), which control the
#'   familywise error rate over the responses and allow for their correlation.
#' @details
#' Each response's statistic is turned into a p-value within each permuted data set (its rank among all the
#' permutations, the observed data counted as one of them); the p-values of the responses are combined in
#' every permuted data set, and the combined p-value is the proportion of permuted data sets whose combined
#' statistic is at least the observed one (Pesarin and Salmaso 2010). Because the same permutation of the
#' units is applied to every response, the correlation between the responses is kept and needs no model.
#' Permuted data sets in which a refit failed for any response are left out.
#'
#' Several tests within one family over time or over levels of a factor (a max-statistic test over time
#' points) are available from \code{\link{contrastmeans}} with \code{adj = "westfall"} and a \code{permlist}.
#' @references
#' Pesarin, F. and Salmaso, L. (2010). \emph{Permutation Tests for Complex Data: Theory, Applications and
#' Software}. Wiley.
#'
#' Westfall, P. H. and Young, S. S. (1993). \emph{Resampling-Based Multiple Testing}. Wiley.
#'
#' Winkler, A. M., Webster, M. A., Brooks, J. C., Tracey, I., Smith, S. M. and Nichols, T. E. (2016).
#' Non-parametric combination and related permutation tests for neuroimaging. \emph{Human Brain Mapping},
#' 37, 1486-1511.
#' @seealso \code{\link{permmodels}}, \code{\link{permstats}}
#' @examples
#' mtcars$cyl <- factor(mtcars$cyl)
#' p1 <- permmodels(lm(mpg ~ wt + cyl, data = mtcars), nperm = 999, type = "II", seed = 1, prt = FALSE)
#' p2 <- permmodels(lm(qsec ~ wt + cyl, data = mtcars), nperm = 999, type = "II", seed = 1, prt = FALSE)
#' permcombine(mpg = p1, qsec = p2)
#' @usage permcombine(..., method = c("fisher", "tippett"), tests = NULL)
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @name permcombine
NULL

permcombine <- function(..., method = c("fisher", "tippett"), tests = NULL) {
  method <- match.arg(method)
  res <- list(...)
  if (length(res) == 1 && is.list(res[[1]]) && is.null(res[[1]][["perm"]])) res <- res[[1]]
  if (length(res) < 2) stop("give at least two results of permmodels().", call. = FALSE)
  nm <- names(res)
  if (is.null(nm)) nm <- rep("", length(res))
  nm[nm == ""] <- paste0("r", which(nm == ""))
  nm <- make.unique(nm)
  sts <- lapply(seq_along(res), function(i) pm_perm_get_stats(res[[i]], nm[i]))
  perms <- lapply(res, function(r) r[["perm"]])
  same <- function(f) all(vapply(perms[-1], function(p) identical(f(p), f(perms[[1]])), logical(1)))
  if (!same(function(p) p$nperm) || !same(function(p) p$exact) || !same(function(p) p$start) ||
      !same(function(p) isTRUE(p$flip)) ||
      !same(function(p) p$idx) || !all(vapply(sts[-1], function(s) identical(s$n, sts[[1]]$n), logical(1)))) {
    stop("the results must use the same permutations: the same observations and the same 'seed' ",
         "(or 'perms') and 'nperm'.", call. = FALSE)
  }
  common <- Reduce(intersect, lapply(sts, `[[`, "names"))
  if (is.null(tests)) tests <- common
  miss <- setdiff(tests, common)
  if (length(miss)) stop("not in every result: ", paste(miss, collapse = ", "), call. = FALSE)
  exact <- isTRUE(perms[[1]]$exact)
  one <- function(tst) {
    # rows: the observed data, then the permutations; columns: the responses
    S <- vapply(sts, function(s) {
      i <- match(tst, s$names)
      round(c(s$obs[i], s$T[i, ]), 6)
    }, numeric(ncol(sts[[1]]$T) + 1))
    S <- S[stats::complete.cases(S), , drop = FALSE]
    B1 <- nrow(S)
    # p-value of each response in each data set: its rank among all of them
    P <- apply(S, 2, function(v) (B1 + 1 - rank(v, ties.method = "min")) / B1)
    C <- if (method == "fisher") -2 * rowSums(log(P)) else -apply(P, 1, min)
    comb <- (sum(C[-1] >= C[1]) + !exact) / (B1 - 1 + !exact)
    minp <- apply(P[-1, , drop = FALSE], 1, min)
    adj <- vapply(P[1, ], function(p) (sum(minp <= p) + !exact) / (B1 - 1 + !exact), numeric(1))
    c(vapply(sts, function(s) s$p[match(tst, s$names)], numeric(1)), comb, adj)
  }
  out <- as.data.frame(do.call(rbind, lapply(tests, one)))
  names(out) <- c(nm, "combined", paste0("adj.", nm))
  rownames(out) <- tests
  attr(out, "method") <- method
  out
}
