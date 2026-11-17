## Results of the curve functions, covariatemeans() and semipred(), in the
## layout of predictmeans(): a classed list (class c(<function>, "pm_curves"))
## with
##   means     one row per point of a curve (or surface): the factor columns,
##             the covariate(s) under their own names, estimate, std.error,
##             df, conf.low, conf.high (model scale) and, with a
##             back-transformation, response, response.std.error,
##             response.conf.low, response.conf.high
##   slopes    semipred(deriv = 1, 2): the slopes (curvatures) in the same
##             columns, with statistic and p.value of slope = 0
##   pairs     with contr: one row per covariate value, the difference of
##             the two curves (level1 - level2) in the columns of
##             predictmeans()$pairs
##   settings  the options that produced them
## and the 1.x fields (plt, pltdf / pred_df, test, significant).
##
## Limits: conf.low, conf.high are pointwise (t with each row's df), as in
## predictmeans(); a simultaneous band (band = "simultaneous") adds
## adj.conf.low, adj.conf.high (and p.adj where a test goes with it), as the
## adjusted limits of predictmeans()$pairs. The prefix "response." marks the
## same quantity on the response scale.

## The table of a curve. lab: the factor and covariate columns; eta, se, df:
## the fitted values (or slopes) on the model scale with their SEs and df;
## crit: the pointwise critical values; crit_s: those of a simultaneous band
## (NULL: none); h, d: the back-transformation and its derivative (NULL:
## none), off: transOff. The response limits are h of the model-scale limits,
## the response SE by the delta method. test = TRUE (slopes): statistic and
## p.value of estimate = 0, and p_adj (with crit_s) the adjusted p-values.
pm_curve_means <- function(lab, eta, se, df, crit, h = NULL, d = NULL, off = 0, crit_s = NULL, test = FALSE,
                           p_adj = NULL) {
  out <- lab
  rownames(out) <- NULL
  est <- as.numeric(eta)
  se <- as.numeric(se)
  df <- rep_len(as.numeric(df), length(est))
  out$estimate <- est
  out$std.error <- se
  out$df <- df
  if (test) {
    out$statistic <- est/se
    out$p.value <- 2*stats::pt(abs(est/se), df, lower.tail = FALSE)
    if (!is.null(crit_s)) out$p.adj <- if (is.null(p_adj)) NA_real_ else p_adj
  }
  out$conf.low <- est - crit*se
  out$conf.high <- est + crit*se
  if (!is.null(crit_s)) {
    out$adj.conf.low <- est - crit_s*se
    out$adj.conf.high <- est + crit_s*se
  }
  if (!is.null(h)) {
    hl <- function(lo, hi) {
      lo <- as.numeric(h(lo)) - off
      hi <- as.numeric(h(hi)) - off
      list(pmin(lo, hi), pmax(lo, hi))
    }
    out$response <- as.numeric(h(est)) - off
    out$response.std.error <- abs(as.numeric(d(est)))*se
    l <- hl(out$conf.low, out$conf.high)
    out$response.conf.low <- l[[1]]
    out$response.conf.high <- l[[2]]
    if (!is.null(crit_s)) {
      l <- hl(out$adj.conf.low, out$adj.conf.high)
      out$response.adj.conf.low <- l[[1]]
      out$response.adj.conf.high <- l[[2]]
    }
  }
  out
}

## The table of the difference of two curves (pm_diff_table()$table, 'tab')
## in the columns of predictmeans()$pairs: level1, level2, the covariate(s)
## ('x', a data frame; default tab$xvar), estimate, std.error, df,
## statistic, p.value, conf.low, conf.high (pointwise t limits); with a
## simultaneous band (critical values crit_s) also p.adj (p_adj, which
## agrees with the band) and adj.conf.low, adj.conf.high; ratio or
## odds.ratio with their (pointwise) limits for a log or logit; the
## response-scale difference (response, ...) when there is a
## back-transformation.
pm_diff_pairs <- function(tab, covariate, levs, level, kind = "none", back = FALSE, crit_s = NULL, p_adj = NULL,
                          x = NULL) {
  n <- nrow(tab)
  if (is.null(x)) x <- stats::setNames(data.frame(tab$xvar), covariate[1])
  out <- cbind(data.frame(level1 = rep(levs[1], n), level2 = rep(levs[2], n), stringsAsFactors = FALSE), x)
  rownames(out) <- NULL
  est <- tab$diff
  se <- tab$ses
  df <- rep_len(as.numeric(tab$Df), n)
  t <- est/se
  out$estimate <- est
  out$std.error <- se
  out$df <- df
  out$statistic <- t
  out$p.value <- 2*stats::pt(abs(t), df, lower.tail = FALSE)
  simult <- !is.null(crit_s)
  if (simult) out$p.adj <- if (is.null(p_adj)) NA_real_ else p_adj
  ct <- stats::qt(1 - level/2, df)
  out$conf.low <- est - ct*se
  out$conf.high <- est + ct*se
  if (simult) {
    out$adj.conf.low <- est - crit_s*se
    out$adj.conf.high <- est + crit_s*se
  }
  if (back && kind %in% c("ratio", "odds.ratio")) {
    out[[kind]] <- exp(est)
    out[[paste0(kind, ".conf.low")]] <- exp(out$conf.low)
    out[[paste0(kind, ".conf.high")]] <- exp(out$conf.high)
  }
  if (back) {
    out$response <- tab$response.diff
    out$response.std.error <- tab$response.ses
    out$response.conf.low <- tab$response.diff - ct*tab$response.ses
    out$response.conf.high <- tab$response.diff + ct*tab$response.ses
    if (simult) {
      out$response.adj.conf.low <- tab$response.diff - crit_s*tab$response.ses
      out$response.adj.conf.high <- tab$response.diff + crit_s*tab$response.ses
    }
  }
  out
}

## Adjusted p-values that go with a simulated simultaneous band
## (pm_sim_crit(draws = TRUE)): the share of simulated maxima at least as
## large as each row's |t| (on the t probability scale when the rows have t
## df), so p.adj < level where the band leaves out zero (up to the
## resolution of the draws).
pm_sim_padj <- function(sim, t, df) {
  if (is.null(sim$M)) return(rep(NA_real_, length(t)))
  m <- if (isTRUE(sim$tdist)) stats::pt(abs(t), rep_len(df, length(t))) else abs(t)
  M <- sort(sim$M)
  # the share of draws; below 1/(number of draws) it is reported as that
  p <- pmax(1 - findInterval(m, M, left.open = TRUE)/length(M), 1/length(M))
  p[!is.finite(t)] <- NA_real_
  p
}

## What the back-transformation is, for settings$back_transform (as
## predictmeans(): pm_back_label(), R/back-transform.R). bt: pm_cov_back().
pm_curve_back_label <- function(model, trans_arg, bt) {
  if (is.null(bt) || is.null(bt$h)) return(NULL)
  by <- if (!is.null(trans_arg)) "trans" else if (isTRUE(bt$link_scale)) "the inverse link" else {
    rt <- pm_response_trans(model)
    if (is.null(rt)) "trans" else paste0("the inverse of ", rt$label, " (the response in the model formula)")
  }
  pm_back_label(list(response = TRUE), by, isTRUE(bt$link_scale), model, isTRUE(bt$adjusted))
}

## The settings of a curve result.
pm_curve_settings <- function(model, modelterm, covariate, level, df_method, mean_type, back, contr = NULL,
                              band = NULL, test = NULL, vcov_adjust = NA, vcov = "model", ...) {
  c(list(modelterm = modelterm, covariate = covariate, level = level, contrast = contr, band = band,
         df_method = df_method, vcov_adjust = vcov_adjust, vcov = vcov, back_transform = back, test = test),
    list(...),
    list(component = if (is.null(pm_zi_component(model))) NA_character_ else mean_type,
         model_class = class(model)[1],
         version = as.character(utils::packageVersion("predictmeans"))))
}

## The classed result: the stable tables first, then the 1.x fields.
pm_curve_result <- function(means = NULL, pairs = NULL, settings, legacy, fun = "covariatemeans", slopes = NULL) {
  out <- c(list(means = means, slopes = slopes, pairs = pairs, settings = settings), legacy)
  out <- out[!vapply(out, is.null, logical(1)) | names(out) %in% c("means", "pairs")]
  class(out) <- c(fun, "pm_curves")
  out
}

#' Curves and Their Differences as Data Frames
#'
#' @description
#' The results of \code{\link{covariatemeans}} and \code{\link{semipred}} hold their tables in the layout of
#' \code{\link{predictmeans}}: \code{means} (one row per point of a curve or surface), \code{slopes}
#' (\code{semipred(deriv = )}) and, with \code{contr}, \code{pairs} (the difference of two curves at each
#' covariate value), with the same column names (\code{estimate}, \code{std.error}, \code{df}, \code{conf.low},
#' \code{conf.high}, ...) and the covariates under their own names. \code{tidy()} and \code{as.data.frame()} return
#' them, \code{autoplot()} returns the plot and \code{print()} shows a short summary.
#' @usage
#' \method{tidy}{pm_curves}(x, what = c("means", "pairs", "slopes"), ...)
#'
#' \method{as.data.frame}{pm_curves}(x, row.names = NULL, optional = FALSE,
#'               what = c("means", "pairs", "slopes"), ...)
#'
#' \method{autoplot}{pm_curves}(object, ...)
#'
#' \method{print}{pm_curves}(x, n = 10, digits = 4, ...)
#' @param x,object A result of \code{\link{covariatemeans}} or \code{\link{semipred}}.
#' @param what Which table: \code{"means"}, \code{"pairs"} or \code{"slopes"}. The default is the one the result
#'   has (\code{"pairs"} for a result with \code{contr}, \code{"slopes"} for \code{semipred(deriv = 1)}).
#' @param row.names Optional row names for the result.
#' @param optional Not used.
#' @param n The number of rows printed.
#' @param digits Significant digits printed.
#' @param \dots Not used.
#' @return
#' \code{tidy()} and \code{as.data.frame()}: a data frame with unrounded values. In every table
#' \code{conf.low} and \code{conf.high} are pointwise limits (t with the row's \code{df}); with
#' \code{band = "simultaneous"} the simultaneous band is in \code{adj.conf.low} and \code{adj.conf.high}, with
#' \code{p.adj} where there is a test (Scheffe over the covariate values for \code{covariatemeans()}, the simulated
#' maximum for \code{semipred()}; the adjusted p-value is below \code{level} exactly where the band leaves out zero,
#' up to the resolution of the 10000 draws for \code{semipred()}). Columns beginning \code{response} are the same
#' quantities on the response scale (back-transformed limits for means, delta-method limits for differences).
#' \describe{
#' \item{\code{means}}{The factor columns of \code{modelterm}, the covariate(s), \code{estimate}, \code{std.error},
#'   \code{df}, \code{conf.low}, \code{conf.high} (the model scale: the link, or the transformed response) and,
#'   when the curves are back-transformed, \code{response}, \code{response.std.error} (delta method),
#'   \code{response.conf.low}, \code{response.conf.high}.}
#' \item{\code{slopes}}{(\code{semipred(deriv = 1, 2)}) The same columns for the slopes (curvatures) on the model
#'   scale, with \code{statistic} (t) and \code{p.value} of slope = 0.}
#' \item{\code{pairs}}{\code{level1}, \code{level2}, the covariate(s), \code{estimate} (\code{level1} minus
#'   \code{level2} on the model scale), \code{std.error}, \code{df}, \code{statistic} (t), \code{p.value},
#'   \code{conf.low}, \code{conf.high}; \code{ratio} or \code{odds.ratio} with their (pointwise) limits for a log or
#'   logit scale; and \code{response}, \code{response.std.error}, \code{response.conf.low},
#'   \code{response.conf.high}, the difference of the back-transformed curves.}
#' }
#' \code{settings} (a list in the result) records \code{modelterm}, \code{covariate}, \code{level},
#' \code{contrast}, \code{band}, \code{df_method}, \code{vcov_adjust}, \code{vcov}, \code{back_transform},
#' \code{test} (\code{covariatemeans()}: the joint F test that the two curves are the same, with the columns of
#' \code{\link{pm_joint}}: \code{df1}, \code{df2}, \code{F.ratio}, \code{p.value}), for
#' \code{semipred()} also \code{deriv}, \code{unconditional} and \code{boundary}, then \code{component}
#' (zero-inflated \code{glmmTMB} fits), \code{model_class} and \code{version}. \code{significant} (in the result)
#' holds the ranges of the covariate where the band asked for leaves out zero.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   data(Oats, package="nlme")
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'   cm <- covariatemeans(fm, "Variety", covariate="nitro", covariateV=c(0, 0.3, 0.6))
#'   cm
#'   as.data.frame(cm)
#'   cd <- covariatemeans(fm, "Variety", covariate="nitro", contr=c("Golden Rain", "Victory"),
#'                        band="simultaneous")
#'   head(as.data.frame(cd))
#'   if (requireNamespace("ggplot2", quietly = TRUE)) ggplot2::autoplot(cd)
#' @name tidy.pm_curves
#' @aliases as.data.frame.pm_curves autoplot.pm_curves print.pm_curves
NULL

## The table a result holds (pairs, then slopes, then means).
pm_curve_tables <- c("pairs", "slopes", "means")
pm_curve_what <- function(x) pm_curve_tables[!vapply(pm_curve_tables, function(t) is.null(x[[t]]), logical(1))][1]

tidy.pm_curves <- function(x, what = c("means", "pairs", "slopes"), ...) {
  pm_result_table(x, if (!missing(what)) what, pm_curve_tables)            # R/tidy.R
}

as.data.frame.pm_curves <- function(x, row.names = NULL, optional = FALSE, what = c("means", "pairs", "slopes"), ...) {
  pm_result_table(x, if (!missing(what)) what, pm_curve_tables, row.names)
}

autoplot.pm_curves <- function(object, ...) {
  if (is.null(object$plt)) stop("This result has no plot.", call. = FALSE)
  object$plt
}

print.pm_curves <- function(x, n = 10, digits = 4, ...) {
  s <- x$settings
  what <- pm_curve_what(x)
  tab <- x[[what]]
  pct <- paste0(format(100*(1 - s$level)), "%")
  cv <- paste(s$covariate, collapse = ", ")
  d <- if (is.null(s$deriv)) 0 else s$deriv
  qty <- c("curves", "slopes", "curvatures")[d + 1]
  grp <- if (!setequal(s$modelterm, s$covariate)) paste0(" of '", s$modelterm, "'") else ""
  along <- if (length(s$covariate) > 1) " over " else " along "
  simult <- identical(s$band, "simultaneous")
  cat(switch(what,
    pairs = paste0("Difference ", s$contrast[1], " - ", s$contrast[2], if (d) paste(" of the", qty), grp),
    slopes = paste0(if (d == 1) "Slopes" else "Curvatures", grp, " (model scale)"),
    means = paste0(if (length(s$covariate) > 1) "Surface" else "Curve", if (nzchar(grp)) "s", grp)),
    along, cv, "\n", sep = "")
  cat(pct, " limits", if (simult) paste0("; adj.conf.*: ", pct, " simultaneous band"), "; ",
      if (identical(s$df_method, "user")) "df given" else paste(s$df_method, "df"), "\n", sep = "")
  tst <- s$test
  if (!is.null(tst) && is.finite(tst$p.value)) {
    cat("Curves equal: F(", tst$df1, ", ", format(round(tst$df2, 1)), ") = ", format(round(tst$F.ratio, 3)),
        ", p = ", format.pval(tst$p.value, digits = 3), "\n", sep = "")
  }
  sig <- x$significant
  if (is.data.frame(sig) && length(s$covariate) == 1) {
    bw <- if (simult) "simultaneous band" else "pointwise limits"
    ref <- if (what == "pairs") "no difference" else "0"
    lev <- if (!is.null(sig$level) && any(nzchar(sig$level))) paste0(sig$level, ": ") else ""
    cat(if (!nrow(sig)) paste0("The ", bw, " include ", ref, " at every value of ", cv, " shown.\n") else
      paste0(if (what == "pairs") "They differ" else paste(if (d == 1) "Slope" else "Curvature", "differs from 0"),
             " (", bw, ") for ", cv, ": ", paste0(lev, format(sig$from, digits = 4), " to ",
                                                  format(sig$to, digits = 4), " (", sig$direction, ")", collapse = "; "),
             "\n"))
  }
  b <- s$back_transform
  if (!is.null(b)) cat("response: back-transformed by ", b$by, "; ", b$estimate, "\n", sep = "")
  if (!is.na(s$component)) cat("Zero inflation: ", s$component, " means\n", sep = "")
  if (is.null(tab) || !nrow(tab)) return(invisible(x))
  if (what == "pairs") tab$level1 <- tab$level2 <- NULL
  print(utils::head(tab, n), digits = digits, row.names = FALSE)
  if (nrow(tab) > n) cat("... ", nrow(tab) - n, " more rows: as.data.frame(x) gives the table.\n", sep = "")
  cat("autoplot(x) draws the plot.\n")
  invisible(x)
}
