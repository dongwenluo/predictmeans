## Plots and summaries built from the stable result tables of predictmeans().
## They never compute statistics and never open devices or print, so they can
## be used in Shiny (renderPlot), reports or plotly::ggplotly().

#' Plots and Summary of Predicted Means
#'
#' @description
#' \code{autoplot()} draws the predicted means, their pairwise p-values or the back-transformed means of a
#' \code{\link{predictmeans}} result as a \pkg{ggplot2} object. \code{summary()} prints the means with a footer that
#' states the degrees of freedom, the SED and LSD ranges, the p-value adjustment and the covariate values used.
#' Both are built from the \code{means} and \code{pairs} tables of the result: they do not compute statistics, open
#' graphics devices or print on their own, so they can be used in reports and Shiny apps
#' (e.g. \code{renderPlot(autoplot(pm))}) or made interactive with \code{plotly::ggplotly()}.
#' @usage
#' \method{autoplot}{predictmeans}(object, type = c("means", "pvalue", "response"),
#'          interval = c("ci", "comparison", "lsd", "none"), values = FALSE, ...)
#'
#' \method{summary}{predictmeans}(object, digits = 4, ...)
#'
#' \method{print}{summary.predictmeans}(x, ...)
#' @param object A result of \code{\link{predictmeans}}.
#' @param type \code{"means"} (the default) plots the predicted means; \code{"pvalue"} a square level plot of the
#'   adjusted pairwise p-values, like a correlation plot (the plot of \code{\link{PMplot}}; with \code{atvar} the
#'   comparisons within each of its levels are blocks down the diagonal); \code{"response"} the back-transformed
#'   means. \code{"pvalue"} needs the \code{pairs} table (with more than 100 means: \code{pairwise = TRUE} in
#'   \code{predictmeans()}).
#' @param interval Bars around the means: confidence limits (\code{"ci"}, the default); comparison intervals
#'   (\code{"comparison"}); half the average LSD (\code{"lsd"}); or none. Confidence limits are for each mean and
#'   must not be read as tests: two means can differ significantly while their 95\% limits overlap, and with
#'   correlated means (blocks, mixed models) by a wide margin. Comparison intervals are made for that reading: mean
#'   \eqn{\pm c\,SE} with \eqn{c = crit \times SED/(SE_1 + SE_2)} averaged over the pairs, so that two means differ
#'   at the level (and adjustment) of the result when their bars do not overlap; exactly so when the SEDs are
#'   equal, and the subtitle says how many pairs are judged differently by their own tests otherwise (columns
#'   \code{comparison.low} and \code{comparison.high} of \code{means}; see \code{\link{predictmeans}}). Half-LSD
#'   bars are the same idea with the average LSD (adjusted as the p-values when \code{adj} has a critical value).
#' @param values For \code{type = "pvalue"}: print each p-value on its tile (the default is \code{FALSE}).
#' @param digits Number of significant digits printed (plus two).
#' @param x A result of \code{summary()}.
#' @param \dots Not used.
#' @return
#' \code{autoplot()} returns a \code{ggplot} object; \code{summary()} returns an object whose print method shows the
#' table and footer.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
#'   summary(pm)
#'   ggplot2::autoplot(pm)
#'   ggplot2::autoplot(pm, interval = "lsd")
#'   ggplot2::autoplot(pm, type = "pvalue")
#' @name autoplot.predictmeans
#' @aliases summary.predictmeans print.summary.predictmeans
NULL

autoplot.predictmeans <- function(object, type = c("means", "pvalue", "response"),
                                  interval = c("ci", "comparison", "lsd", "none"), values = FALSE, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  means <- object$means
  vars <- unlist(strsplit(object$settings$modelterm, ":"))
  level <- object$settings$level

  if (type == "pvalue") return(pm_pvalue_plot(object, values = values))

  if (type == "response") {
    if (is.null(means$response)) {
      stop("These means are not back-transformed; use type = \"means\".", call. = FALSE)
    }
    means$y <- means$response
    means$lo <- means$response.conf.low
    means$hi <- means$response.conf.high
    ylab <- "Back-transformed mean"
    if (interval %in% c("lsd", "comparison")) interval <- "ci"
  } else {
    means$y <- means$estimate
    means$lo <- means$conf.low
    means$hi <- means$conf.high
    ylab <- "Predicted mean"
    if (interval == "lsd") {
      # bars of half an LSD: two means whose bars do not overlap differ at
      # the LSD's level (exactly so when all LSDs are equal)
      half <- pm_average_lsd(object, means, vars) / 2
      means$lo <- means$y - half
      means$hi <- means$y + half
    }
    if (interval == "comparison") {
      if (is.null(means$comparison.low)) {
        stop("There are no comparison intervals: they need the pairs (one mean, or more than 100 means without ",
             "pairwise = TRUE, or permutation tests without limits).", call. = FALSE)
      }
      means$lo <- means$comparison.low
      means$hi <- means$comparison.high
    }
  }
  cmp <- object$settings$comparison
  ci_label <- switch(interval,
                     ci = paste0(" with ", (1 - level) * 100, "% confidence limits"),
                     lsd = paste0(" with half-LSD (", level * 100, "%) bars"),
                     comparison = " with comparison intervals",
                     none = "")
  subtitle <- if (interval == "comparison") {
    pm_wrap_title(paste0("Means differ at the ", level * 100, "% level",
                         if (isTRUE(cmp$adjusted)) paste0(" (", cmp$adjust, ")"),
                         " when their bars do not overlap", pm_comparison_caveat(cmp)), 11, rel = 1)
  }

  x <- vars[1]
  grp <- if (length(vars) >= 2) vars[2] else NULL
  dodge <- ggplot2::position_dodge(width = if (is.null(grp)) 0 else 0.3)
  p <- if (is.null(grp)) {
    ggplot2::ggplot(means, ggplot2::aes(x = .data[[x]], y = .data$y, group = 1))
  } else {
    ggplot2::ggplot(means, ggplot2::aes(x = .data[[x]], y = .data$y,
                                        colour = .data[[grp]], group = .data[[grp]]))
  }
  p <- p + ggplot2::geom_line(position = dodge) + ggplot2::geom_point(size = 2, position = dodge)
  if (interval != "none") {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                                    width = 0.15, position = dodge)
  }
  if (length(vars) >= 3) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", paste(sprintf("`%s`", vars[-(1:2)]), collapse = "+"))))
  }
  p + ggplot2::labs(x = x, y = ylab, subtitle = subtitle,
                    title = pm_wrap_title(paste0("Predicted means for '", object$settings$modelterm, "'", ci_label), 11)) +
    ggplot2::theme_bw()
}

## Average LSD used for the half-LSD bars: within atvar levels when given.
## Without a pairs table (no comparisons asked for) the average LSD of the
## summary.
pm_average_lsd <- function(object, means, vars) {
  pairs <- object$pairs
  if (is.null(pairs)) {
    lsd <- object$LSD
    return(rep(if (is.numeric(lsd) && "Aveg.LSD" %in% names(lsd)) unname(lsd["Aveg.LSD"]) else NA_real_, nrow(means)))
  }
  if (!is.null(pairs$by)) pairs <- pairs[!is.na(pairs$by), , drop = FALSE]
  # the LSDs that go with the p-values: adjusted when 'adj' has a critical value
  lsd <- if (!is.null(pairs$adj.conf.high)) (pairs$adj.conf.high - pairs$adj.conf.low)/2 else pairs$lsd
  rep(mean(lsd, na.rm = TRUE), nrow(means))
}

## Level plot of the adjusted pairwise p-values: the plot of PMplot()
## (pm_ptile_plot(), R/PMplot.R), from the pairs table. With 'atvar' the
## comparisons within each of its levels are blocks down the diagonal.
pm_pvalue_plot <- function(object, values = FALSE) {
  pairs <- object$pairs
  if (is.null(pairs) || nrow(pairs) == 0) {
    stop("There are no pairwise comparisons to plot (one mean, or more than 100 means without pairwise = TRUE).",
         call. = FALSE)
  }
  s <- object$settings
  vars <- unlist(strsplit(s$modelterm, ":"))
  all_levs <- do.call(paste, c(lapply(object$means[vars], as.character), sep = ":"))
  by <- if (is.null(pairs$by)) rep("all", nrow(pairs)) else pairs$by
  keep <- !is.na(by)
  groups <- unique(by[keep])
  mats <- lapply(groups, function(g) {
    pg <- pairs[keep & by == g, , drop = FALSE]
    levs <- intersect(all_levs, c(pg$level1, pg$level2))
    P <- matrix(NA_real_, length(levs), length(levs), dimnames = list(levs, levs))
    i <- match(pg$level1, levs)
    j <- match(pg$level2, levs)
    P[cbind(pmax(i, j), pmin(i, j))] <- pg$p.adj
    P
  })
  P <- do.call(adiag, c(unname(mats), pad = NA))
  title <- paste0("Pairwise p-values (adjusted by '", s$adjust, "') for '", s$modelterm, "'",
                  if (!is.null(pairs$by)) paste0("\nat each level of '", paste(s$atvar, collapse = ":"), "'"))
  pm_ptile_plot(P, unlist(lapply(mats, rownames)), level = s$level, title = title,
                legend_title = if (identical(s$adjust, "none")) "p-value" else "Adjusted p-value",
                blocks = vapply(mats, nrow, integer(1)), values = values)
}

summary.predictmeans <- function(object, digits = 4, ...) {
  s <- object$settings
  means <- object$means
  vars <- unlist(strsplit(s$modelterm, ":"))
  show <- c(setdiff(names(means), c("std.error", "df", "conf.low", "conf.high", "comparison.low", "comparison.high", "estimate", "response",
                                    "response.conf.low", "response.conf.high", "group", "n", "estimable")),
            "estimate", "std.error", "df", "conf.low", "conf.high",
            intersect(c("response", "response.conf.low", "response.conf.high", "group"), names(means)))
  tab <- means[show]
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], function(z) signif(z, digits + 2))
  pairs <- object$pairs
  if (!is.null(pairs) && !is.null(pairs$by)) pairs <- pairs[!is.na(pairs$by), , drop = FALSE]
  # the SED and LSD summaries of the result (also there without a pairs table)
  structure(list(table = tab, settings = s, pairs = pairs, nonestimable = sum(!means$estimable),
                 sed = object[["Standard Error of Differences"]], lsd = object$LSD, digits = digits),
            class = "summary.predictmeans")
}

print.summary.predictmeans <- function(x, ...) {
  s <- x$settings
  fmt <- function(v) format(signif(v, x$digits))
  conf <- paste0((1 - s$level) * 100, "%")
  cat("Predicted means for '", s$modelterm, "' (", s$model_class, " model; df: ", s$df_method,
      if (isTRUE(s$vcov_adjust)) ", Kenward-Roger adjusted SEs" else "",
      if (!is.null(s$vcov) && !identical(s$vcov, "model")) paste0(", covariance: ", s$vcov) else "", ")\n", sep = "")
  cat("Confidence limits: ", if (identical(s$df_method, "permutation")) "mean +/- 2 SE" else conf, "\n\n", sep = "")
  print(x$table, row.names = FALSE, digits = x$digits)
  cat("\n")
  p <- x$pairs
  if (is.null(p) && is.numeric(x$sed) && all(c("Max.SED", "Min.SED", "Aveg.SED") %in% names(x$sed))) {
    # no comparisons asked for: the summaries over all pairs of means
    n <- nrow(x$table)
    cat(n*(n - 1)/2, " pairs of means, not tabulated (pairwise = TRUE compares them).\n", sep = "")
    cat("SED: min ", fmt(x$sed[["Min.SED"]]), ", mean ", fmt(x$sed[["Aveg.SED"]]), ", max ", fmt(x$sed[["Max.SED"]]), "\n", sep = "")
    if (is.numeric(x$lsd) && "Aveg.LSD" %in% names(x$lsd)) {
      ddf <- attr(x$lsd, "Degree of freedom")
      cat("LSD (", s$level * 100, "%", if (!is.null(ddf)) paste0(", ", fmt(ddf), " df"), "): min ", fmt(x$lsd[["Min.LSD"]]),
          ", mean ", fmt(x$lsd[["Aveg.LSD"]]), ", max ", fmt(x$lsd[["Max.LSD"]]), "\n", sep = "")
    }
  }
  if (!is.null(p) && nrow(p) > 0) {
    within <- if (!is.null(s$atvar)) paste0(" within each level of ", paste(s$atvar, collapse = " and ")) else ""
    cat(nrow(p), " pairwise comparisons", within, ".\n", sep = "")
    cat("SED: min ", fmt(min(p$std.error)), ", mean ", fmt(mean(p$std.error)), ", max ", fmt(max(p$std.error)), "\n", sep = "")
    # the LSDs that go with the p-values: adjusted when 'adj' has a critical value
    adj_lsd <- !is.null(p$adj.conf.high)
    lsd <- if (adj_lsd) (p$adj.conf.high - p$adj.conf.low)/2 else p$lsd
    cat("LSD (", s$level * 100, "%", if (adj_lsd) paste0(", ", s$adjust), "): min ", fmt(min(lsd, na.rm = TRUE)), ", mean ",
        fmt(mean(lsd, na.rm = TRUE)), ", max ", fmt(max(lsd, na.rm = TRUE)), "\n", sep = "")
    method <- if (identical(s$df_method, "permutation")) paste0("permutation (", s$nperm, " permutations)") else "t tests"
    cat("p-values: ", method, ", adjusted by '", s$adjust, "'.\n", sep = "")
    if ("group" %in% names(x$table)) {
      cat("Means sharing a letter are not significantly different at the ", s$level * 100, "% level",
          within, ".\n", sep = "")
    }
    cmp <- s$comparison
    if (!is.null(cmp)) {
      cat("Comparison intervals: mean +/- ", fmt(cmp$mult), " SE (", round(100 * cmp$conf.level, 1),
          "% intervals); means differ at the ", s$level * 100, "% level when these do not overlap",
          pm_comparison_caveat(cmp), ". The ", (1 - s$level) * 100, "% confidence limits are not for comparing means.\n",
          sep = "")
    }
  }
  if (!is.null(s$covariates)) {
    cv <- s$covariates
    if (nrow(cv) == 1) {
      cat("Covariates fixed at: ", paste(names(cv), fmt(unlist(cv)), sep = " = ", collapse = ", "), "\n", sep = "")
    }
  }
  if (x$nonestimable > 0) cat(x$nonestimable, "mean(s) cannot be estimated.\n")
  invisible(x)
}
