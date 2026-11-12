#' Slopes of a Covariate by Group
#'
#' @description
#' The slope (rate of change) of the fitted response in a numeric covariate within each level of a model term, with
#' standard errors, degrees of freedom, confidence limits and tests, and pairwise comparisons of the slopes, as
#' \code{emmeans::emtrends()} and SAS \code{ESTIMATE} statements on slope terms.
#' @details
#' The slopes are on the scale of the linear predictor. They come from the reference grid of \code{\link{Kmatrix}}:
#' the change in the K matrix between \code{at - h} and \code{at + h}, divided by \code{2h} with \code{h} = 0.0001 times
#' the standard deviation of the covariate. This is exact when the covariate enters the model linearly (also through
#' interactions) and a central difference for \code{poly()}, splines and smooths, whose slope depends on \code{at}.
#' Standard errors and degrees of freedom follow \code{\link{predictmeans}}: for \code{lmer} models fitted by REML
#' Kenward-Roger (with the adjusted covariance when \code{vcov_adjust = TRUE}), for \code{lme} the smallest \pkg{nlme} df
#' of the terms with the covariate, otherwise the residual df.
#' @usage
#' pm_slopes(model, covariate, modelterm, at = NULL, level = 0.05, adj = "none",
#'           df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-SAS"),
#'           vcov_adjust = TRUE,
#'           weights = c("equal", "proportional", "outer", "cells", "flat"), vcov = NULL)
#' @param model A fitted model, as for \code{\link{predictmeans}}.
#' @param covariate Name (in "quotes") of the numeric covariate.
#' @param modelterm Name (in "quotes") of the factor term whose levels get a slope each, e.g. \code{"Variety"}.
#' @param at The value of the covariate where the slopes are taken (matters only when the covariate enters the model
#'   non-linearly); the default is its mean.
#' @param level Significance level for the confidence limits, the default is 0.05.
#' @param adj Adjustment of the p-values of the pairwise comparisons of slopes, as in \code{\link{predictmeans}}.
#' @param df_method,vcov_adjust,weights,vcov As in \code{\link{predictmeans}}.
#' @return
#' A list of two data frames:
#' \item{slopes}{one row per level of \code{modelterm}: the factor columns, \code{slope}, \code{std.error}, \code{df},
#'   \code{conf.low}, \code{conf.high}, \code{statistic} (t) and \code{p.value} (of slope = 0).}
#' \item{pairs}{one row per pair of levels: \code{level1}, \code{level2}, \code{estimate} (slope of \code{level1} minus
#'   slope of \code{level2}), \code{std.error}, \code{df}, \code{statistic}, \code{p.value} and \code{p.adj}.}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm_slopes(fm, "nitro", "Variety")
#' @name pm_slopes
NULL

pm_slopes <- function(model, covariate, modelterm, at = NULL, level = 0.05, adj = "none",
                      df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-SAS"),
                      vcov_adjust = TRUE, weights = c("equal", "proportional", "outer", "cells", "flat"),
                      vcov = NULL) {
  vcov_adjust_default <- missing(vcov_adjust)
  df_method <- match.arg(df_method)
  # as in predictmeans(): containment and between-within df with the unadjusted covariance
  if (df_method %in% c("containment", "between-within") && vcov_adjust_default) vcov_adjust <- FALSE
  weights <- match.arg(weights)
  adj <- pm_match_adj(adj)
  if (inherits(model, "aovlist")) model <- aovlist_lmer(model)
  pm_check_modelterm(model, modelterm)
  dat <- pm_get_data(model)
  if (!covariate %in% names(dat) || !is.numeric(dat[[covariate]])) {
    stop("'", covariate, "' must be a numeric variable of the model's data.", call. = FALSE)
  }
  x <- dat[[covariate]]
  if (is.null(at)) at <- mean(x, na.rm = TRUE)
  h <- 1e-4 * (if (isTRUE(stats::sd(x, na.rm = TRUE) > 0)) stats::sd(x, na.rm = TRUE) else 1)
  kp <- Kmatrix(model, modelterm, covariate, covariateV = at + h, weights = weights)
  km <- Kmatrix(model, modelterm, covariate, covariateV = at - h, weights = weights)
  vars <- unlist(strsplit(modelterm, ":"))
  lab <- kp$fctnames[vars]
  rownames(lab) <- NULL
  levs <- do.call(paste, c(lapply(lab, as.character), sep = ":"))

  parms <- pm_model_parms(model, df_method, vcov_adjust, vcov_adjust_default, vcov)
  mp <- parms$mp
  L <- (kp$K - km$K)/(2*h)
  L <- L[, mp$estimable, drop = FALSE]
  rownames(L) <- levs
  V <- as.matrix(mp$vcov)
  df_of <- function(LL) pm_slope_df(model, LL, modelterm, covariate, mp, parms)

  est <- drop(L %*% mp$coef)
  se <- sqrt(rowSums((L %*% V) * L))
  dfs <- df_of(L)
  crit <- stats::qt(1 - level/2, dfs)
  slopes <- cbind(lab, slope = est, std.error = se, df = dfs, conf.low = est - crit*se, conf.high = est + crit*se,
                  statistic = est/se, p.value = 2*stats::pt(-abs(est/se), dfs))

  pairs <- NULL
  n <- length(est)
  if (n > 1) {
    pp <- pm_pair_index(n)
    rL <- pm_pair_matrix(n) %*% L
    pe <- drop(rL %*% mp$coef)
    pV <- rL %*% V %*% t(rL)
    pse <- sqrt(diag(pV))
    pdf <- df_of(rL)
    stat <- pe/pse
    p <- 2*stats::pt(-abs(stat), pdf)
    p_adj <- if (adj == "tukey") {
      stats::ptukey(sqrt(2)*abs(stat), n, pdf, lower.tail = FALSE)
    } else if (adj %in% pm_adjust_extra_methods) {
      pm_adjust_extra(adj, stat, p, pdf, rank = n - 1, R = if (adj == "mvt") stats::cov2cor(pV))
    } else stats::p.adjust(p, adj)
    pairs <- data.frame(level1 = levs[pp$j], level2 = levs[pp$i], estimate = pe, std.error = pse, df = pdf,
                        statistic = stat, p.value = p, p.adj = p_adj, stringsAsFactors = FALSE)
  }
  rownames(slopes) <- NULL
  structure(list(slopes = slopes, pairs = pairs), class = "pm_slopes",
            covariate = covariate, at = at, adjust = adj)
}

## Degrees of freedom of rows L (coefficient space) of slope contrasts.
pm_slope_df <- function(model, L, modelterm, covariate, mp, parms) {
  if (identical(parms$df_type, "between-within")) {
    # the smallest between-within df of the terms with the covariate
    bw <- pm_betwithin_df(model)
    with_cov <- vapply(strsplit(names(bw), ":"), function(v) covariate %in% v, logical(1))
    return(rep(min(bw[with_cov]), nrow(L)))
  }
  if (inherits(model, c("lme", "gls")) && any(parms$df_type %in% c("Satterthwaite", "Kenward-Roger"))) {
    return(as.numeric(pm_nlme_df(L, model, parms$df_type, parms$vcov_adj)))
  }
  if (inherits(model, "lmerMod")) {
    return(as.numeric(df_term(model, modelterm, ctrmatrix = L, type = parms$df_type, vcov_adj = parms$vcov_adj)))
  }
  if (inherits(model, "lme")) {
    tdf <- stats::terms(model$fixDF)
    with_cov <- vapply(strsplit(names(tdf), ":"), function(v) covariate %in% v, logical(1))
    return(rep(min(tdf[with_cov]), nrow(L)))
  }
  custom <- pm_df(model, L)
  if (!is.null(custom)) return(rep_len(as.numeric(custom), nrow(L)))
  d <- mp$df
  rep(if (is.numeric(d) && length(d) == 1 && is.finite(d) && d > 0) d else Inf, nrow(L))
}

print.pm_slopes <- function(x, digits = 4, ...) {
  cat("Slopes of", sQuote(attr(x, "covariate")), "at", format(attr(x, "at"), digits = digits), "\n\n")
  print(x$slopes, digits = digits, row.names = FALSE)
  if (!is.null(x$pairs)) {
    cat("\nPairwise differences of slopes (p.adj: ", attr(x, "adjust"), ")\n\n", sep = "")
    print(x$pairs, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
