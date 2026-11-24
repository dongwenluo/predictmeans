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
#' Standard errors and degrees of freedom follow \code{\link{predictmeans}}: for models fitted by REML (\code{lmer},
#' \code{lme}, \code{gls}) Kenward-Roger (with the adjusted covariance when \code{vcov_adjust = TRUE}); with one df per
#' term (containment, between-within, \pkg{nlme}'s df) the smallest df of the terms with the covariate; otherwise the
#' residual df.
#' @usage
#' pm_slopes(model, covariate, modelterm, at = NULL, conf_level = 0.95, adj = "none",
#'           df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment",
#'             "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
#'           vcov_adjust = TRUE,
#'           weights = c("equal", "proportional", "outer", "cells", "flat"),
#'           vcov = NULL, ...)
#' @param model A fitted model, as for \code{\link{predictmeans}}.
#' @param covariate Name (in "quotes") of the numeric covariate.
#' @param modelterm Name (in "quotes") of the factor term whose levels get a slope each, e.g. \code{"Variety"}.
#' @param at The value of the covariate where the slopes are taken (matters only when the covariate enters the model
#'   non-linearly); the default is its mean.
#' @param conf_level Confidence level of the limits, the default is 0.95 (\code{level = 0.05} before 2.0).
#' @param adj Adjustment of the p-values of the pairwise comparisons of slopes, as in \code{\link{predictmeans}}.
#' @param df_method,vcov_adjust,weights,vcov As in \code{\link{predictmeans}}.
#' @param ... Not used: \code{level} (before 2.0) stops with its 2.0 spelling.
#' @return
#' An object of class \code{"pm_slopes"}, a list in the layout of \code{\link{predictmeans}} (\code{tidy()} and
#' \code{as.data.frame()} return its tables, \code{what = "slopes"} or \code{"pairs"}):
#' \item{slopes}{one row per level of \code{modelterm}: the factor columns, \code{estimate} (the slope),
#'   \code{std.error}, \code{df}, \code{statistic} (t), \code{p.value} (of slope = 0), \code{conf.low} and
#'   \code{conf.high}.}
#' \item{pairs}{one row per pair of levels: \code{level1}, \code{level2}, \code{estimate} (slope of \code{level1} minus
#'   slope of \code{level2}), \code{std.error}, \code{df}, \code{statistic}, \code{p.value}, \code{p.adj},
#'   \code{conf.low}, \code{conf.high} and, with the adjustments that have them (tukey, bonferroni, sidak, scheffe,
#'   mvt), \code{adj.conf.low}, \code{adj.conf.high}: simultaneous limits that leave out 0 exactly when
#'   \code{p.adj} is below \code{1 - conf_level}.}
#' \item{settings}{\code{modelterm}, \code{covariate}, \code{at}, \code{level}, \code{adjust}, \code{df_method},
#'   \code{vcov_adjust}, \code{vcov}, \code{model_class}, \code{version}.}
#' @seealso \code{\link{covariatemeans}}, \code{\link{predictmeans}}, \code{\link{advice}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   sl <- pm_slopes(fm, "nitro", "Variety", adj = "tukey")
#'   sl
#'   as.data.frame(sl, what = "pairs")
#' @name pm_slopes
NULL

#' @export
pm_slopes <- function(model, covariate, modelterm, at = NULL, conf_level = 0.95, adj = "none",
                      df_method = c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
                      vcov_adjust = TRUE, weights = c("equal", "proportional", "outer", "cells", "flat"),
                      vcov = NULL, ...) {
  pm_defunct("pm_slopes", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  level <- pm_conf_alpha(conf_level)
  vcov_adjust_default <- missing(vcov_adjust)
  df_method <- match.arg(df_method)
  weights <- match.arg(weights)
  adj <- pm_match_adj(adj)
  if (inherits(model, "aovlist")) model <- pm_aov_carrier(model)
  pm_check_modelterm(model, modelterm)
  dat <- pm_get_data(model)
  if (!covariate %in% names(dat) || !is.numeric(dat[[covariate]])) {
    stop("'", covariate, "' must be a numeric variable of the model's data.", call. = FALSE)
  }
  x <- dat[[covariate]]
  if (is.null(at)) at <- mean(x, na.rm = TRUE)
  h <- pm_fd_h(0, if (isTRUE(stats::sd(x, na.rm = TRUE) > 0)) stats::sd(x, na.rm = TRUE) else 1)   # R/deriv.R
  kp <- Kmatrix(model, modelterm, covariate, at = at + h, weights = weights)
  km <- suppressWarnings(Kmatrix(model, modelterm, covariate, at = at - h, weights = weights))   # (warned once)
  vars <- unlist(strsplit(modelterm, ":"))
  lab <- kp$fctnames[vars]
  rownames(lab) <- NULL
  levs <- do.call(paste, c(lapply(lab, as.character), sep = ":"))

  inf <- pm_inference(model, df_method, vcov_adjust, vcov_adjust_default, vcov)     # R/inference.R
  mp <- inf$mp
  L <- (kp$K - km$K)/(2*h)
  L <- pm_K_cols(L, mp)
  rownames(L) <- levs
  V <- as.matrix(mp$vcov)
  # one df per term: the smallest of the terms with the covariate
  cov_terms <- pm_terms_with(model, covariate)
  df_of <- function(LL) inf$df(LL, cov_terms)

  est <- drop(L %*% mp$coef)
  se <- sqrt(rowSums((L %*% V) * L))
  dfs <- df_of(L)
  crit <- stats::qt(1 - level/2, dfs)
  # the layout of predictmeans(): estimate, std.error, df, statistic, p.value, limits
  slopes <- cbind(lab, estimate = est, std.error = se, df = dfs, statistic = est/se,
                  p.value = 2*stats::pt(-abs(est/se), dfs), conf.low = est - crit*se, conf.high = est + crit*se)

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
    fams <- pm_families(rep("all", length(stat)), est = !is.na(stat), ij = cbind(pp$j, pp$i), L = rL)
    R <- if (adj == "mvt") stats::cov2cor(pV)
    p_adj <- pm_adjust_p(adj, stat, p, pdf, fams, R = R)   # R/p-adjust.R
    ct <- stats::qt(1 - level/2, pdf)
    pairs <- data.frame(level1 = levs[pp$j], level2 = levs[pp$i], estimate = pe, std.error = pse, df = pdf,
                        statistic = stat, p.value = p, p.adj = p_adj, conf.low = pe - ct*pse, conf.high = pe + ct*pse,
                        stringsAsFactors = FALSE)
    # simultaneous limits that agree with p.adj (the adjustments that have them)
    crit_adj <- pm_adjust_crit(adj, level, pdf, fams, R = R)
    if (!is.null(crit_adj)) {
      pairs$adj.conf.low <- pe - crit_adj*pse
      pairs$adj.conf.high <- pe + crit_adj*pse
    }
  }
  rownames(slopes) <- NULL
  settings <- list(modelterm = modelterm, covariate = covariate, at = at, level = level, adjust = adj,
                   df_method = inf$label, vcov_adjust = isTRUE(inf$vcov_adjust),
                   vcov = if (is.null(inf$vcov_label)) "model" else inf$vcov_label,
                   model_class = class(model)[1], version = as.character(utils::packageVersion("predictmeans")))
  out <- structure(list(slopes = slopes, pairs = pairs, settings = settings), class = "pm_slopes",
                   covariate = covariate, at = at, adjust = adj)
  pm_add_advice(out, tryCatch(pm_slopes_advice(out, model), error = function(e) NULL))
}

#' @exportS3Method generics::tidy pm_slopes
tidy.pm_slopes <- function(x, what = c("slopes", "pairs"), ...) {
  pm_result_table(x, if (!missing(what)) what, c("slopes", "pairs"))       # R/tidy.R
}

#' @exportS3Method as.data.frame pm_slopes
as.data.frame.pm_slopes <- function(x, row.names = NULL, optional = FALSE, what = c("slopes", "pairs"), ...) {
  pm_result_table(x, if (!missing(what)) what, c("slopes", "pairs"), row.names)
}

#' @exportS3Method print pm_slopes
print.pm_slopes <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("Slopes of", sQuote(attr(x, "covariate")), "at", format(attr(x, "at"), digits = digits),
      if (!is.null(s)) paste0("(", format(100*(1 - s$level)), "% limits; ", s$df_method, " df)"), "\n\n")
  print(x$slopes, digits = digits, row.names = FALSE)
  if (!is.null(x$pairs)) {
    cat("\nPairwise differences of slopes (p.adj: ", attr(x, "adjust"), ")\n\n", sep = "")
    print(x$pairs, digits = digits, row.names = FALSE)
  }
  pm_print_advice(x)
  invisible(x)
}
