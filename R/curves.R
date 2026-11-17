## Results of the curve functions (covariatemeans(); semipred() next) in the
## layout of predictmeans(): a classed list with
##   means     one row per point of a curve: the factor columns, the
##             covariate under its own name, estimate, std.error, df,
##             conf.low, conf.high (model scale) and, with a
##             back-transformation, response, response.std.error,
##             response.conf.low, response.conf.high
##   pairs     with contr: one row per covariate value, the difference of
##             the two curves (level1 - level2) in the columns of
##             predictmeans()$pairs
##   settings  the options that produced them
## and the 1.x fields (plt, pltdf, test, significant) unchanged.

## The table of curves. lab: the factor and covariate columns; eta, se, df:
## the fitted values on the model scale with their SEs and df; crit: the
## critical values of the limits; h, d: the back-transformation and its
## derivative (NULL: none), off: transOff. The response limits are h of the
## model-scale limits, the response SE by the delta method.
pm_curve_means <- function(lab, eta, se, df, crit, h = NULL, d = NULL, off = 0) {
  out <- lab
  rownames(out) <- NULL
  out$estimate <- as.numeric(eta)
  out$std.error <- as.numeric(se)
  out$df <- rep_len(as.numeric(df), length(eta))
  out$conf.low <- out$estimate - crit*out$std.error
  out$conf.high <- out$estimate + crit*out$std.error
  if (!is.null(h)) {
    lo <- as.numeric(h(out$conf.low)) - off
    hi <- as.numeric(h(out$conf.high)) - off
    out$response <- as.numeric(h(out$estimate)) - off
    out$response.std.error <- abs(as.numeric(d(out$estimate)))*out$std.error
    out$response.conf.low <- pmin(lo, hi)
    out$response.conf.high <- pmax(lo, hi)
  }
  out
}

## The table of the difference of two curves (pm_diff_table()$table, 'tab')
## in the columns of predictmeans()$pairs: level1, level2, the covariate,
## estimate, std.error, df, statistic, p.value, conf.low, conf.high
## (pointwise t limits); with a simultaneous band (r: the rank of the
## contrasts along the covariate, crit_s its Scheffe critical values) also
## p.adj and adj.conf.low, adj.conf.high, which agree (pm_adjust_p(),
## pm_adjust_crit(), R/p-adjust.R); ratio or odds.ratio with their limits for
## a log or logit; the response-scale difference (response, ...) when there
## is a back-transformation.
pm_diff_pairs <- function(tab, covariate, levs, level, kind = "none", back = FALSE, r = NULL, crit_s = NULL) {
  n <- nrow(tab)
  out <- data.frame(level1 = rep(levs[1], n), level2 = rep(levs[2], n), x = tab$xvar, stringsAsFactors = FALSE)
  names(out)[3] <- covariate[1]
  est <- tab$diff
  se <- tab$ses
  df <- rep_len(as.numeric(tab$Df), n)
  t <- est/se
  p <- 2*stats::pt(abs(t), df, lower.tail = FALSE)
  out$estimate <- est
  out$std.error <- se
  out$df <- df
  out$statistic <- t
  out$p.value <- p
  simult <- !is.null(crit_s)
  if (simult) {
    fam <- list(list(idx = seq_len(n), ok = which(is.finite(t)), m = n, k = 2, r = r))
    out$p.adj <- pm_adjust_p("scheffe", t, p, df, fam)
  }
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
  }
  out
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
pm_curve_settings <- function(model, modelterm, covariate, level, inf, Df, mean_type, back, contr = NULL,
                              band = NULL, test = NULL) {
  list(modelterm = modelterm, covariate = covariate, level = level,
       contrast = contr, band = band,
       df_method = if (!is.null(Df) && !all(Df %in% c("NULL", ""))) "user" else inf$label,
       vcov_adjust = isTRUE(inf$vcov_adjust),
       vcov = if (is.null(inf$vcov_label)) "model" else inf$vcov_label,
       back_transform = back,
       test = test,
       component = if (is.null(pm_zi_component(model))) NA_character_ else mean_type,
       model_class = class(model)[1],
       version = as.character(utils::packageVersion("predictmeans")))
}

## The classed result: the stable tables first, then the 1.x fields.
pm_curve_result <- function(means = NULL, pairs = NULL, settings, legacy, fun = "covariatemeans") {
  out <- c(list(means = means, pairs = pairs, settings = settings), legacy)
  class(out) <- c(fun, "pm_curves")
  out
}

#' Curves and Their Differences as Data Frames
#'
#' @description
#' The result of \code{\link{covariatemeans}} holds its tables in the layout of \code{\link{predictmeans}}:
#' \code{means} (one row per point of a curve) and, with \code{contr}, \code{pairs} (the difference of the two
#' curves at each covariate value), with the same column names (\code{estimate}, \code{std.error}, \code{df},
#' \code{conf.low}, \code{conf.high}, ...) and the covariate under its own name. \code{tidy()} and
#' \code{as.data.frame()} return them, \code{autoplot()} returns the plot and \code{print()} shows a short summary.
#' @usage
#' \method{tidy}{pm_curves}(x, what = c("means", "pairs"), ...)
#'
#' \method{as.data.frame}{pm_curves}(x, row.names = NULL, optional = FALSE,
#'               what = c("means", "pairs"), ...)
#'
#' \method{autoplot}{pm_curves}(object, ...)
#'
#' \method{print}{pm_curves}(x, n = 10, digits = 4, ...)
#' @param x,object A result of \code{\link{covariatemeans}}.
#' @param what Which table: \code{"means"} or \code{"pairs"}. The default is the one the result has
#'   (\code{"pairs"} for a result with \code{contr}).
#' @param row.names Optional row names for the result.
#' @param optional Not used.
#' @param n The number of rows printed.
#' @param digits Significant digits printed.
#' @param \dots Not used.
#' @return
#' \code{tidy()} and \code{as.data.frame()}: a data frame with unrounded values.
#' \describe{
#' \item{\code{means}}{The factor columns of \code{modelterm}, the covariate, \code{estimate}, \code{std.error},
#'   \code{df}, \code{conf.low}, \code{conf.high} (the model scale: the link, or the transformed response) and,
#'   when the curves are back-transformed, \code{response}, \code{response.std.error} (delta method),
#'   \code{response.conf.low}, \code{response.conf.high} (the back-transformed limits).}
#' \item{\code{pairs}}{\code{level1}, \code{level2}, the covariate, \code{estimate} (\code{level1} minus
#'   \code{level2} on the model scale), \code{std.error}, \code{df}, \code{statistic} (t), \code{p.value},
#'   \code{conf.low}, \code{conf.high} (pointwise limits); with \code{band = "simultaneous"} also \code{p.adj}
#'   and \code{adj.conf.low}, \code{adj.conf.high} (Scheffe over the covariate values; the adjusted p-value is
#'   below \code{level} exactly when these limits leave out zero); \code{ratio} or \code{odds.ratio} with
#'   their (pointwise) limits for a log or logit scale; and \code{response}, \code{response.std.error},
#'   \code{response.conf.low}, \code{response.conf.high}, the difference of the back-transformed curves with
#'   delta-method limits.}
#' }
#' \code{settings} (a list in the result) records \code{modelterm}, \code{covariate}, \code{level},
#' \code{contrast}, \code{band}, \code{df_method}, \code{vcov_adjust}, \code{vcov}, \code{back_transform},
#' \code{test} (the joint F test that the two curves are the same), \code{component} (zero-inflated
#' \code{glmmTMB} fits), \code{model_class} and \code{version}.
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

tidy.pm_curves <- function(x, what = c("means", "pairs"), ...) {
  what <- if (missing(what) && is.null(x$means)) "pairs" else match.arg(what)
  out <- x[[what]]
  if (is.null(out)) stop("There are no ", what, " in this result.", call. = FALSE)
  out
}

as.data.frame.pm_curves <- function(x, row.names = NULL, optional = FALSE, what = c("means", "pairs"), ...) {
  out <- if (missing(what)) tidy.pm_curves(x) else tidy.pm_curves(x, what = what)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

autoplot.pm_curves <- function(object, ...) {
  if (is.null(object$plt)) stop("This result has no plot.", call. = FALSE)
  object$plt
}

print.pm_curves <- function(x, n = 10, digits = 4, ...) {
  s <- x$settings
  pct <- paste0(format(100*(1 - s$level)), "%")
  cv <- paste(s$covariate, collapse = ", ")
  dfm <- if (identical(s$df_method, "user")) "df given" else paste(s$df_method, "df")
  if (!is.null(x$pairs)) {
    tab <- x$pairs
    cat("Difference ", s$contrast[1], " - ", s$contrast[2], " of '", s$modelterm, "' along ", cv, " (", pct, " ",
        s$band, " limits, ", dfm, ")\n", sep = "")
    tst <- s$test
    if (!is.null(tst) && is.finite(tst[["p.value"]])) {
      cat("Curves equal: F(", tst[["numDF"]], ", ", format(round(tst[["denDF"]], 1)), ") = ", format(round(tst[["F"]], 3)),
          ", p = ", format.pval(tst[["p.value"]], digits = 3), "\n", sep = "")
    }
    sig <- x$significant
    if (is.data.frame(sig)) {
      cat(if (!nrow(sig)) paste0("The ", s$band, " band includes no difference at every value of ", cv, " shown.\n") else
        paste0("They differ for ", cv, ": ", paste0(format(sig$from, digits = 4), " to ", format(sig$to, digits = 4),
                                                   " (", sig$direction, ")", collapse = "; "), "\n"))
    }
    tab$level1 <- tab$level2 <- NULL
  } else {
    tab <- x$means
    cat("Curves of '", s$modelterm, "' along ", cv, " (", pct, " limits, ", dfm, ")\n", sep = "")
  }
  b <- s$back_transform
  if (!is.null(b)) cat("response: back-transformed by ", b$by, "; ", b$estimate, "\n", sep = "")
  if (!is.na(s$component)) cat("Zero inflation: ", s$component, " means\n", sep = "")
  if (is.null(tab) || !nrow(tab)) return(invisible(x))
  print(utils::head(tab, n), digits = digits, row.names = FALSE)
  if (nrow(tab) > n) cat("... ", nrow(tab) - n, " more rows: as.data.frame(x) gives the table.\n", sep = "")
  cat("autoplot(x) draws the plot.\n")
  invisible(x)
}
