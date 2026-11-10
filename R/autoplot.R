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
#'          interval = c("ci", "lsd", "none"), ...)
#'
#' \method{summary}{predictmeans}(object, digits = 4, ...)
#'
#' \method{print}{summary.predictmeans}(x, ...)
#' @param object A result of \code{\link{predictmeans}}.
#' @param type \code{"means"} (the default) plots the predicted means; \code{"pvalue"} a heatmap of the adjusted
#'   pairwise p-values (within each level of \code{atvar} when it was given); \code{"response"} the back-transformed
#'   means.
#' @param interval Bars around the means: confidence limits (\code{"ci"}, the default), half the average LSD
#'   (\code{"lsd"}; means whose bars do not overlap differ significantly when the LSDs are equal), or none.
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
                                  interval = c("ci", "lsd", "none"), ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  means <- object$means
  vars <- unlist(strsplit(object$settings$modelterm, ":"))
  level <- object$settings$level

  if (type == "pvalue") return(pm_pvalue_plot(object))

  if (type == "response") {
    if (is.null(means$response)) {
      stop("These means are not back-transformed; use type = \"means\".", call. = FALSE)
    }
    means$y <- means$response
    means$lo <- means$response.conf.low
    means$hi <- means$response.conf.high
    ylab <- "Back-transformed mean"
    if (interval == "lsd") interval <- "ci"
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
  }
  ci_label <- switch(interval,
                     ci = paste0(" with ", (1 - level) * 100, "% confidence limits"),
                     lsd = paste0(" with half-LSD (", level * 100, "%) bars"),
                     none = "")

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
  p + ggplot2::labs(x = x, y = ylab,
                    title = paste0("Predicted means for '", object$settings$modelterm, "'", ci_label)) +
    ggplot2::theme_bw()
}

## Average LSD used for the half-LSD bars: within atvar levels when given.
pm_average_lsd <- function(object, means, vars) {
  pairs <- object$pairs
  if (is.null(pairs)) return(rep(NA_real_, nrow(means)))
  if (!is.null(pairs$by)) pairs <- pairs[!is.na(pairs$by), , drop = FALSE]
  rep(mean(pairs$lsd, na.rm = TRUE), nrow(means))
}

## Heatmap of adjusted pairwise p-values in the significance bands of PMplot().
pm_pvalue_plot <- function(object) {
  pairs <- object$pairs
  if (is.null(pairs) || nrow(pairs) == 0) stop("There are no pairwise comparisons to plot.", call. = FALSE)
  if (!is.null(pairs$by)) pairs <- pairs[!is.na(pairs$by), , drop = FALSE]
  levs <- unique(c(rbind(pairs$level1, pairs$level2)))
  levs <- levs[order(match(levs, do.call(paste, c(lapply(object$means[unlist(strsplit(object$settings$modelterm, ":"))],
                                                          as.character), sep = ":"))))]
  pairs$row <- factor(pairs$level2, levels = rev(levs))
  pairs$col <- factor(pairs$level1, levels = levs)
  bands <- c("p <= 0.01", "0.01 < p <= 0.05", "0.05 < p <= 0.1", "p > 0.1")
  pairs$band <- factor(bands[findInterval(pairs$p.adj, c(0.01, 0.05, 0.1), left.open = TRUE) + 1], levels = bands)
  adj <- object$settings$adjust
  p <- ggplot2::ggplot(pairs, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$band)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_manual(values = stats::setNames(c("#0D0DFF", "#5D5DFF", "#A1A1FF", "#E4E4FF"), bands),
                               drop = FALSE, name = "Adjusted p-value") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("Pairwise p-values (adjusted by '", adj, "') for '",
                                 object$settings$modelterm, "'")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
  if (!is.null(pairs$by)) {
    p <- p + ggplot2::facet_wrap(~by, scales = "free")
  }
  p
}

summary.predictmeans <- function(object, digits = 4, ...) {
  s <- object$settings
  means <- object$means
  vars <- unlist(strsplit(s$modelterm, ":"))
  show <- c(setdiff(names(means), c("std.error", "df", "conf.low", "conf.high", "estimate", "response",
                                    "response.conf.low", "response.conf.high", "group", "n", "estimable")),
            "estimate", "std.error", "df", "conf.low", "conf.high",
            intersect(c("response", "response.conf.low", "response.conf.high", "group"), names(means)))
  tab <- means[show]
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], function(z) signif(z, digits + 2))
  pairs <- object$pairs
  if (!is.null(pairs) && !is.null(pairs$by)) pairs <- pairs[!is.na(pairs$by), , drop = FALSE]
  structure(list(table = tab, settings = s, pairs = pairs, nonestimable = sum(!means$estimable),
                 digits = digits),
            class = "summary.predictmeans")
}

print.summary.predictmeans <- function(x, ...) {
  s <- x$settings
  fmt <- function(v) format(signif(v, x$digits))
  conf <- paste0((1 - s$level) * 100, "%")
  cat("Predicted means for '", s$modelterm, "' (", s$model_class, " model; df: ", s$df_method,
      if (isTRUE(s$vcov_adjust)) ", Kenward-Roger adjusted SEs" else "", ")\n", sep = "")
  cat("Confidence limits: ", if (identical(s$df_method, "permutation")) "mean +/- 2 SE" else conf, "\n\n", sep = "")
  print(x$table, row.names = FALSE, digits = x$digits)
  cat("\n")
  p <- x$pairs
  if (!is.null(p) && nrow(p) > 0) {
    within <- if (!is.null(s$atvar)) paste0(" within each level of ", paste(s$atvar, collapse = " and ")) else ""
    cat(nrow(p), " pairwise comparisons", within, ".\n", sep = "")
    cat("SED: min ", fmt(min(p$std.error)), ", mean ", fmt(mean(p$std.error)), ", max ", fmt(max(p$std.error)), "\n", sep = "")
    cat("LSD (", s$level * 100, "%): min ", fmt(min(p$lsd)), ", mean ", fmt(mean(p$lsd)), ", max ", fmt(max(p$lsd)), "\n", sep = "")
    method <- if (identical(s$df_method, "permutation")) paste0("permutation (", s$nperm, " permutations)") else "t tests"
    cat("p-values: ", method, ", adjusted by '", s$adjust, "'.\n", sep = "")
    if ("group" %in% names(x$table)) {
      cat("Means sharing a letter are not significantly different at the ", s$level * 100, "% level",
          within, ".\n", sep = "")
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
