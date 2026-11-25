## Plots and summaries built from the stable result tables of predictmeans().
## They never compute statistics and never open devices or print, so they can
## be used in Shiny (renderPlot), reports or plotly::ggplotly().

#' Plots and Summary of Predicted Means
#'
#' @description
#' \code{autoplot()} draws a \code{\link{predictmeans}} result as a \pkg{ggplot2} object: the means with one LSD
#' bar (as predictmeans 1.x), as bars, as horizontal confidence intervals, or the pairwise p-values.
#' \code{plot()} draws the same and prints it. \code{summary()} prints the means with a footer that states the
#' degrees of freedom, the SED and LSD ranges, the p-value adjustment and the covariate values used. All are built
#' from the \code{means} and \code{pairs} tables of the result: they do not compute statistics, so they can be used
#' in reports and Shiny apps (e.g. \code{renderPlot(autoplot(pm))}) or made interactive with
#' \code{plotly::ggplotly()}. Titles, labels and font sizes are changed with \pkg{ggplot2}, e.g.
#' \code{autoplot(pm) + labs(title = "Yield") + theme_bw(14)}.
#' @usage
#' \method{autoplot}{predictmeans}(object, type = c("means", "bar", "ci", "pvalue"),
#'          interval = c("lsd", "sed", "ci", "comparison", "none"),
#'          scale = c("auto", "model", "response"),
#'          data = FALSE, x = NULL, group = NULL, panel = NULL,
#'          scales = c("fixed", "free_y"), lines = TRUE,
#'          jitter = 0.1, values = FALSE, ...)
#'
#' \method{plot}{predictmeans}(x, ...)
#'
#' \method{summary}{predictmeans}(object, digits = 4, ...)
#'
#' \method{print}{summary.predictmeans}(x, ...)
#' @param object A result of \code{\link{predictmeans}}.
#' @param type \code{"means"} (the default): points joined by lines for each level of \code{group};
#'   \code{"bar"}: a bar chart of the means; \code{"ci"}: each mean with its confidence limits, drawn horizontally
#'   (one row per mean; on the response scale when the means are back-transformed); \code{"pvalue"}: a square level
#'   plot of the adjusted pairwise p-values (the plot of \code{\link{PMplot}}; with \code{atvar} the comparisons
#'   within each of its levels are blocks down the diagonal), which needs the \code{pairs} table.
#'   (\code{"response"} of earlier 2.0 versions is \code{type = "ci"} or \code{scale = "response"}.)
#' @param interval For \code{"means"} and \code{"bar"}: \code{"lsd"} (the default) draws one bar of the average LSD
#'   beside the means of each panel, as predictmeans 1.x: two means differ at the level of the result when they are
#'   further apart than the bar (exactly so when all the LSDs are equal; the caption gives their range when they are
#'   not). It is the LSD that goes with the p-values: adjusted when \code{adj} has a critical value. With
#'   \code{atvar} each panel (level of \code{atvar}) gets the LSD of the pairs compared within it. \code{"sed"} draws
#'   the average SED instead; \code{"ci"} confidence limits on each mean (for each mean, not for comparing means: two
#'   means can differ while their limits overlap); \code{"comparison"} comparison intervals (mean
#'   \eqn{\pm c\,SE}, \eqn{c = crit \times SED/(SE_1 + SE_2)} averaged over the pairs: means differ when their bars
#'   do not overlap; see \code{\link{predictmeans}}); \code{"none"}. One LSD bar replaces the half-LSD bars on each
#'   mean of earlier 2.0 versions: those said the same when the SEDs are equal, and looked like confidence limits.
#' @param scale \code{"model"}: the scale of the model (e.g. \code{log(y + 1)} or the link); \code{"response"}: the
#'   back-transformed means (needs them); \code{"auto"} (the default): \code{"response"} for \code{type = "ci"}
#'   when the means are back-transformed, else \code{"model"}. LSD and SED bars exist only on the model scale.
#' @param scales For \code{"means"} and \code{"bar"} with panels: \code{"fixed"} (the default), one y axis for all
#'   panels, with the LSD (SED) bar at the same place in each, from the lowest point of the plot, so equal LSDs look
#'   equal; \code{"free_y"}, each panel its own y range, with its bar from its own lowest point (the bars are
#'   then read within their panel only: the caption says so).
#' @param data Add the observations (\code{TRUE}; jittered points), from \code{\link{pm_data}} of the model, on the
#'   scale shown: the transformed response on the model scale, the original response on the response scale. For a
#'   link function the data are shown only on the response scale. A data frame can be given instead (with the
#'   factors of \code{modelterm} and the response).
#' @param x,group,panel Names of factors of \code{modelterm} for the x axis, the lines (colours) and the panels.
#'   The defaults: with \code{atvar}, the panels are its levels and x the first other factor; else x is the first
#'   factor, group the second and the panels the rest. These replace \code{plotord} of predictmeans 1.x.
#' @param lines For \code{type = "means"}: join the means of a group by lines (\code{TRUE}) or not.
#' @param jitter Horizontal jitter of the data points, in units of the x categories.
#' @param values For \code{type = "pvalue"}: print each p-value on its tile (the default is \code{FALSE}).
#' @param digits Number of significant digits printed (plus two).
#' @param \dots For \code{plot()}: arguments of \code{autoplot()}, with \code{xvar} for its \code{x} (in \code{plot()} \code{x} is the
#'   result). Otherwise not used.
#' @return
#' \code{autoplot()} returns a \code{ggplot} object; \code{plot()} prints it and returns it invisibly;
#' \code{summary()} returns an object whose print method shows the table and footer.
#' @seealso \code{\link{predictmeans}}, \code{\link{tidy.predictmeans}}, \code{\link{pm_data}}, \code{\link{advice}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
#'   summary(pm)
#'   ggplot2::autoplot(pm)
#'   ggplot2::autoplot(pm, data = TRUE)
#'   ggplot2::autoplot(pm, type = "bar")
#'   ggplot2::autoplot(pm, type = "ci")
#'   ggplot2::autoplot(pm, type = "pvalue")
#'   plot(predictmeans(fm, "nitro:Variety"), xvar = "Variety", group = "nitro")
#'   fl <- lm(log(yield) ~ nitro * Variety, data = Oats)
#'   ggplot2::autoplot(predictmeans(fl, "nitro:Variety"), type = "ci", data = TRUE)
#' @name autoplot.predictmeans
#' @aliases plot.predictmeans summary.predictmeans print.summary.predictmeans
NULL

#' @exportS3Method ggplot2::autoplot predictmeans
autoplot.predictmeans <- function(object, type = c("means", "bar", "ci", "pvalue"),
                                  interval = c("lsd", "sed", "ci", "comparison", "none"),
                                  scale = c("auto", "model", "response"), data = FALSE, x = NULL, group = NULL,
                                  panel = NULL, scales = c("fixed", "free_y"), lines = TRUE, jitter = 0.1,
                                  values = FALSE, ...) {
  if (identical(type, "response")) {
    stop("type = \"response\" is now type = \"ci\" (horizontal intervals) or scale = \"response\".", call. = FALSE)
  }
  type <- match.arg(type)
  interval <- match.arg(interval)
  scale <- match.arg(scale)
  scales <- match.arg(scales)
  if (type == "pvalue") return(pm_pvalue_plot(object, values = values))
  s <- object$settings
  means <- object$means
  vars <- unlist(strsplit(s$modelterm, ":"))
  level <- s$level
  has_resp <- !is.null(means$response)
  if (scale == "auto") scale <- if (type == "ci" && has_resp) "response" else "model"
  if (scale == "response" && !has_resp) {
    stop("These means are not back-transformed; use scale = \"model\".", call. = FALSE)
  }
  if (scale == "response") {
    means$y <- means$response
    means$lo <- means$response.conf.low
    means$hi <- means$response.conf.high
  } else {
    means$y <- means$estimate
    means$lo <- means$conf.low
    means$hi <- means$conf.high
  }
  if (type == "ci") interval <- "ci"
  if (scale == "response" && interval %in% c("lsd", "sed", "comparison")) {
    message("LSD, SED and comparison bars are on the model scale; confidence limits are drawn on the response scale.")
    interval <- "ci"
  }
  if (interval == "comparison") {
    if (is.null(means$comparison.low)) {
      stop("There are no comparison intervals: they need the pairs (one mean, or more than 100 means without ",
           "pairwise = TRUE, or permutation tests without limits).", call. = FALSE)
    }
    means$lo <- means$comparison.low
    means$hi <- means$comparison.high
  }
  lay <- pm_plot_layout(vars, s$atvar, x, group, panel)
  # means over empty cells (column 'over'): open points; no comparison bars
  # for those whose pairs are not the differences of the shown means
  cmp <- s$comparison
  aff <- if (is.null(means$over)) rep(FALSE, nrow(means)) else !is.na(means$over)
  misread <- do.call(paste, c(lapply(means[vars], as.character), sep = ":")) %in% cmp$left_out
  if (any(misread) && interval == "comparison") means$lo[misread] <- means$hi[misread] <- NA_real_
  means$cells <- ifelse(aff, "observed cells", "all cells")
  obs <- pm_plot_data(object, data, vars, scale)

  ylab <- pm_scale_label(s, scale)
  bar <- if (interval %in% c("lsd", "sed")) pm_plot_bars(object, means, lay, interval)
  int_label <- switch(interval,
                      ci = paste0(" with ", (1 - level) * 100, "% confidence limits"),
                      lsd = paste0(" with the average LSD (", level * 100, "%",
                                   if (isTRUE(bar$adjusted)) paste0(", ", s$adjust), ")"),
                      sed = " with the average SED",
                      comparison = " with comparison intervals",
                      none = "")
  title <- pm_wrap_title(paste0(pm_scale_label(s, scale, "title"), int_label), 14)
  captions <- c(if (any(aff)) "Open points: means over their observed cells only (column 'over')",
                if (any(misread) && interval == "comparison") "no bar where a mean's pairs are compared on common cells",
                bar$caption, if (!is.null(obs)) attr(obs, "caption"),
                if (!is.null(bar) && scales == "free_y" && type != "ci" && length(lay$panel))
                  "free y scales: read each LSD bar within its own panel only")
  caption <- if (length(captions)) pm_wrap_title(paste(captions, collapse = "; "), 11, rel = 1)
  subtitle <- if (interval == "comparison") {
    pm_wrap_title(paste0("Means differ at the ", level * 100, "% level",
                         if (isTRUE(cmp$adjusted)) paste0(" (", cmp$adjust, ")"),
                         " when their bars do not overlap", pm_comparison_caveat(cmp)), 11, rel = 1)
  }

  p <- if (type == "ci") {
    pm_plot_ci(means, lay, obs, jitter)
  } else {
    pm_plot_means(means, lay, obs, type, interval, lines, jitter, bar, free = scales == "free_y" && length(lay$panel) > 0)
  }
  if (length(lay$panel)) {
    fs <- if (type == "ci") "free_y" else scales
    p <- p + if (length(lay$panel) == 2 && type != "ci") {
      ggplot2::facet_grid(stats::as.formula(paste(sprintf("`%s`", lay$panel), collapse = " ~ ")),
                          labeller = ggplot2::label_both, scales = fs)
    } else {
      ggplot2::facet_wrap(stats::as.formula(paste("~", paste(sprintf("`%s`", lay$panel), collapse = "+"))),
                          labeller = ggplot2::label_both, scales = fs)
    }
  }
  p + ggplot2::labs(x = if (type == "ci") NULL else lay$x, y = if (type == "ci") NULL else ylab,
                    title = title, subtitle = subtitle, caption = caption) +
    (if (type == "ci") ggplot2::labs(x = ylab)) +
    ggplot2::theme_bw() +
    (if (type == "bar") ggplot2::theme(panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_line(),
                                       panel.grid = ggplot2::element_blank(), legend.position = "top"))
}

#' @exportS3Method plot predictmeans
plot.predictmeans <- function(x, ...) {
  # 'x' of autoplot() (the factor on the x axis) cannot be given by name to
  # plot(): it is the object; use xvar =
  args <- list(...)
  if (!is.null(args$xvar)) { args$x <- args$xvar; args$xvar <- NULL }
  p <- do.call(autoplot.predictmeans, c(list(x), args))
  print(p)
  invisible(p)
}

## Which factor goes where: x, group (lines/colours) and panels. With atvar
## the panels are its levels (the comparisons are within them).
pm_plot_layout <- function(vars, atvar, x = NULL, group = NULL, panel = NULL) {
  chk <- function(v, what) {
    bad <- setdiff(v, vars)
    if (length(bad)) stop("'", what, "' must name factors of modelterm; not: '", paste(bad, collapse = "', '"), "'.",
                          call. = FALSE)
    v
  }
  chk(c(x, group, panel), "x, group and panel")
  if (is.null(panel)) {
    panel <- if (length(atvar)) setdiff(atvar, c(x, group)) else character(0)
  }
  rest <- setdiff(vars, c(panel, x, group))
  if (is.null(x)) {
    x <- rest[1]
    if (is.na(x)) x <- setdiff(vars, c(group, panel))[1]
    if (is.na(x)) x <- vars[1]
    rest <- setdiff(rest, x)
  }
  if (is.null(group) && length(rest)) {
    group <- rest[1]
    rest <- rest[-1]
  }
  if (!length(atvar) || length(panel)) panel <- c(panel, rest)
  list(x = x, group = group, panel = unique(panel))
}

## The average LSD (or SED) of each panel: of the pairs compared within it
## (with atvar: within its levels), the LSDs that go with the p-values
## (adjusted when 'adj' has a critical value). Without a pairs table, the
## average of the result's summary. Returns list(data (one row per panel:
## the panel factors, value), adjusted, caption).
pm_plot_bars <- function(object, means, lay, interval) {
  pairs <- object$pairs
  s <- object$settings
  key_of <- function(d) if (length(lay$panel)) do.call(paste, c(lapply(d[lay$panel], as.character), sep = "\r")) else
    rep("all", nrow(d))
  panels <- unique(means[lay$panel])
  if (!length(lay$panel)) panels <- data.frame(row.names = 1)
  if (is.null(pairs)) {
    v <- if (interval == "lsd") object$LSD else object[["Standard Error of Differences"]]
    nm <- if (interval == "lsd") "Aveg.LSD" else "Aveg.SED"
    val <- if (is.numeric(v) && nm %in% names(v)) unname(v[[nm]]) else NA_real_
    panels$value <- val
    return(list(data = panels, adjusted = FALSE,
                caption = paste0("one ", toupper(interval), " bar: the average over all pairs of means")))
  }
  vars <- unlist(strsplit(s$modelterm, ":"))
  lab <- do.call(paste, c(lapply(means[vars], as.character), sep = ":"))
  m1 <- match(pairs$level1, lab)
  m2 <- match(pairs$level2, lab)
  k1 <- key_of(means[m1, , drop = FALSE])
  k2 <- key_of(means[m2, , drop = FALSE])
  keep <- !is.na(pairs$estimate) & k1 == k2
  if (!is.null(pairs$by)) keep <- keep & !is.na(pairs$by)
  adjusted <- !is.null(pairs$adj.conf.high)
  val <- if (interval == "sed") pairs$std.error else if (adjusted) (pairs$adj.conf.high - pairs$adj.conf.low)/2 else
    pairs$lsd
  keep <- keep & is.finite(val)
  if (!any(keep)) {
    # no pairs within a panel (each panel one mean): all the pairs
    keep <- is.finite(val)
    k1[] <- "all"
  }
  agg <- tapply(val[keep], k1[keep], mean)
  panels$value <- unname(agg[if (length(lay$panel)) key_of(panels) else "all"])
  r <- if (any(keep)) range(val[keep]) else c(0, 0)
  caption <- if (diff(r) > 1e-8 * max(abs(r), 1)) {
    paste0(if (interval == "sed") "SEDs" else "LSDs", " range from ", format(signif(r[1], 4)), " to ",
           format(signif(r[2], 4)), "; the bar is their average: use the pairs for each comparison")
  }
  list(data = panels, adjusted = adjusted && interval == "lsd", caption = caption)
}

## The observations for data = TRUE (the rows kept in the result by
## predictmeans(), pm_plot_rows()) or a data frame given: the factors of
## modelterm and .yplot on the scale shown. NULL when not asked for.
pm_plot_data <- function(object, data, vars, scale) {
  if (isFALSE(data) || is.null(data)) return(NULL)
  r <- object$settings$response
  if (is.data.frame(data)) {
    d <- data
    if (is.null(d$.y)) {
      nm <- if (!is.null(r)) r$name else NA
      if (is.na(nm) || !nm %in% names(d)) stop("The data must have the response '", nm, "'.", call. = FALSE)
      d$.y_orig <- d$.y <- d[[nm]]
    }
  } else {
    d <- attr(object, "plot_data")
    if (is.null(d)) stop("This result does not keep its observations; give the data: autoplot(x, data = df).",
                         call. = FALSE)
  }
  if (is.null(d$.y_orig)) d$.y_orig <- d$.y
  miss <- setdiff(vars, names(d))
  if (length(miss)) stop("The data have no column '", paste(miss, collapse = "', '"), "'.", call. = FALSE)
  if (scale == "model" && identical(r$scale, "latent")) {
    message("The data are not shown with latent means (classes of an ordinal response); mean_type = \"prob\" gives ",
            "the class probabilities.")
    return(NULL)
  }
  link <- identical(r$scale, "link")
  if (scale == "model" && link) {
    out <- NULL
    attr(out, "caption") <- NULL
    message("The data are not shown on the link scale; use scale = \"response\" (type = \"ci\") to see them.")
    return(NULL)
  }
  d$.yplot <- if (scale == "response") d$.y_orig else d$.y
  for (v in vars) d[[v]] <- factor(d[[v]], levels = levels(factor(object$means[[v]])))
  d <- d[stats::complete.cases(d[c(vars, ".yplot")]), , drop = FALSE]
  attr(d, "caption") <- "open points: the observations"
  d
}

## Means plot ("means": points and lines; "bar": bars) with one LSD/SED bar
## per panel beside the means, or bars on each mean.
pm_plot_means <- function(means, lay, obs, type, interval, lines, jitter, bar, free = FALSE) {
  x <- lay$x
  grp <- lay$group
  one <- is.null(grp)
  ng <- if (one) 1L else length(unique(means[[grp]]))
  has_bar <- !is.null(bar) && any(!is.na(bar$data$value))
  # the LSD/SED bar at the lower left of each panel, half a category step
  # left of the first level, from the lowest point shown (0 for bars)
  bar_lab <- if (has_bar) toupper(interval) else NULL
  xlev <- levels(factor(means[[x]]))
  means[[x]] <- factor(as.character(means[[x]]), levels = xlev)
  if (!is.null(obs)) obs[[x]] <- factor(as.character(obs[[x]]), levels = xlev)
  dodge_w <- if (one) 0 else if (type == "bar") 0.9 else 0.25
  dodge <- ggplot2::position_dodge(width = dodge_w)
  p <- if (one) ggplot2::ggplot(means, ggplot2::aes(x = .data[[x]], y = .data$y, group = 1)) else
    if (type == "bar") ggplot2::ggplot(means, ggplot2::aes(x = .data[[x]], y = .data$y, fill = .data[[grp]],
                                                         group = .data[[grp]])) else
      ggplot2::ggplot(means, ggplot2::aes(x = .data[[x]], y = .data$y, colour = .data[[grp]], group = .data[[grp]]))
  p <- p + ggplot2::scale_x_discrete(limits = xlev, drop = FALSE,
                                     expand = ggplot2::expansion(add = c(if (has_bar) 0.85 else 0.6, 0.6)))
  if (type == "bar") {
    p <- p + if (one) ggplot2::geom_col(width = 0.7, fill = "#9ECAE1") else
      ggplot2::geom_col(position = dodge, width = 0.85)
    if (!one) p <- p + ggplot2::scale_fill_manual(values = grDevices::colorRampPalette(c("#DEEBF7", "#3182BD"))(ng))
  }
  if (!is.null(obs)) {
    pos <- if (one) ggplot2::position_jitter(width = jitter, height = 0, seed = 1) else
      ggplot2::position_jitterdodge(jitter.width = jitter, jitter.height = 0, dodge.width = dodge_w, seed = 1)
    # position_jitterdodge() dodges by fill or colour (ggplot2 < 3.5 ignores
    # 'group'); the open grey circles draw no fill, and it adds no legend
    aes_o <- if (one) ggplot2::aes(x = .data[[x]], y = .data$.yplot) else
      ggplot2::aes(x = .data[[x]], y = .data$.yplot, group = .data[[grp]], fill = .data[[grp]])
    p <- p + ggplot2::geom_point(aes_o, data = obs, inherit.aes = FALSE, position = pos, colour = "grey50",
                                 shape = 1, size = 1.4, alpha = 0.7, show.legend = FALSE)
  }
  if (type == "means") {
    if (lines) p <- p + if (one) ggplot2::geom_line(linewidth = 0.5) else
      ggplot2::geom_line(ggplot2::aes(linetype = .data[[grp]]), position = dodge, linewidth = 0.9)
    p <- p + (if (one) ggplot2::geom_point(ggplot2::aes(shape = .data$cells), size = 2, colour = "red") else
      ggplot2::geom_point(ggplot2::aes(shape = .data$cells), size = 2, position = dodge)) +
      ggplot2::scale_shape_manual(values = c("all cells" = 16, "observed cells" = 1), guide = "none")
  }
  if (interval %in% c("ci", "comparison")) {
    p <- p + if (one) ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0.2,
                                             colour = "blue", na.rm = TRUE) else
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0.2,
                             position = dodge, na.rm = TRUE)
  }
  if (has_bar) {
    b <- bar$data
    key <- if (length(lay$panel)) do.call(paste, c(lapply(means[lay$panel], as.character), sep = "\r")) else
      rep("all", nrow(means))
    bkey <- if (length(lay$panel)) do.call(paste, c(lapply(b[lay$panel], as.character), sep = "\r")) else "all"
    ylo <- if (interval %in% c("ci", "comparison")) means$lo else means$y
    # fixed scales: the same place in every panel, from the lowest point of the
    # plot; free_y: from each panel's own lowest point
    low <- tapply(ylo, key, min, na.rm = TRUE)
    if (!free) low[] <- min(ylo, na.rm = TRUE)
    if (!is.null(obs)) {
      okey <- if (length(lay$panel)) do.call(paste, c(lapply(obs[lay$panel], as.character), sep = "\r")) else
        rep("all", nrow(obs))
      olow <- tapply(obs$.yplot, okey, min, na.rm = TRUE)
      low[] <- if (free) pmin(low, olow[names(low)], na.rm = TRUE) else min(low, olow, na.rm = TRUE)
    }
    if (type == "bar") low[] <- 0
    b$xb <- 0.5
    b$label <- bar_lab
    b$ymin <- unname(low[bkey])
    b$ymax <- b$ymin + b$value
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(x = .data$xb, ymin = .data$ymin, ymax = .data$ymax), data = b,
                                    inherit.aes = FALSE, width = 0.15, linewidth = 0.8, colour = "blue",
                                    na.rm = TRUE) +
      ggplot2::geom_text(ggplot2::aes(x = .data$xb, y = .data$ymax, label = .data$label), data = b,
                         inherit.aes = FALSE, vjust = -0.5, hjust = 0.3, size = 3, colour = "blue", na.rm = TRUE)
    # room for the label above the bar (each panel's range with free scales)
    yr <- if (free) tapply(c(means$y, b$ymax), c(key, bkey), function(v) diff(range(v, na.rm = TRUE)))[bkey] else
      diff(range(c(means$y, b$ymin, b$ymax), na.rm = TRUE))
    b$ytop <- b$ymax + 0.12 * unname(yr)
    p <- p + ggplot2::geom_blank(ggplot2::aes(x = .data$xb, y = .data$ytop), data = b, inherit.aes = FALSE)
  }
  p
}

## Horizontal confidence intervals, one row per mean (the "back" plot of
## predictmeans 1.x): red means and limits, the observations as open blue
## points.
pm_plot_ci <- function(means, lay, obs, jitter) {
  rowvars <- setdiff(c(lay$x, lay$group), lay$panel)
  lab <- function(d) do.call(paste, c(lapply(d[rowvars], as.character), sep = ":"))
  levs <- unique(lab(means))
  means$row <- factor(lab(means), levels = rev(levs))
  p <- ggplot2::ggplot(means, ggplot2::aes(y = .data$row, x = .data$y))
  if (!is.null(obs)) {
    obs$row <- factor(lab(obs), levels = rev(levs))
    p <- p + ggplot2::geom_point(ggplot2::aes(y = .data$row, x = .data$.yplot), data = obs, inherit.aes = FALSE,
                                 position = ggplot2::position_jitter(width = 0, height = jitter, seed = 1),
                                 colour = "blue", shape = 1, alpha = 0.6)
  }
  p + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), width = 0.2, orientation = "y",
                             colour = "red", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$cells), colour = "red", size = 2) +
    ggplot2::scale_shape_manual(values = c("all cells" = 16, "observed cells" = 1), guide = "none")
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
  labs <- unlist(lapply(mats, rownames))
  if (!is.null(pairs$by)) {
    # the panel names the level of 'atvar': the axes show the other factors only
    other <- setdiff(vars, s$atvar)
    if (length(other) && length(other) < length(vars)) {
      key <- do.call(paste, c(lapply(object$means[other], as.character), sep = ":"))
      labs <- key[match(labs, all_levs)]
    }
  }
  title <- paste0("Pairwise p-values (adjusted by '", s$adjust, "') for '", s$modelterm, "'",
                  if (!is.null(pairs$by)) paste0("\nat each level of '", paste(s$atvar, collapse = ":"), "'"))
  pm_ptile_plot(P, labs, level = s$level, title = title,
                legend_title = if (identical(s$adjust, "none")) "p-value" else "Adjusted p-value",
                blocks = vapply(mats, nrow, integer(1)), values = values,
                block_names = if (!is.null(pairs$by)) as.character(groups))
}

#' @exportS3Method summary predictmeans
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
  structure(list(table = tab, settings = s, pairs = pairs, nonestimable = sum(!means$estimable), probs = object$probs,
                 sed = object[["Standard Error of Differences"]], lsd = object$LSD, digits = digits),
            class = "summary.predictmeans", advice = attr(object, "advice"))
}

#' @exportS3Method print summary.predictmeans
print.summary.predictmeans <- function(x, ...) {
  s <- x$settings
  fmt <- function(v) format(signif(v, x$digits))
  conf <- paste0((1 - s$level) * 100, "%")
  cat("Predicted means for '", s$modelterm, "' (", s$model_class, " model; df: ", s$df_method,
      if (isTRUE(s$vcov_adjust)) ", Kenward-Roger adjusted SEs" else "",
      if (!is.null(s$vcov) && !identical(s$vcov, "model")) paste0(", covariance: ", s$vcov) else "", ")\n", sep = "")
  cat("Confidence limits: ", conf, if (identical(s$df_method, "permutation"))
    " (the model's t limits of each mean; the permutation tests are of the comparisons)", "\n", sep = "")
  if (identical(s$response$scale, "latent")) {
    cat("Scale: latent mean of ", s$response$name, " (", s$response$link, "); differences are ",
        if (s$response$link == "logit") "log odds ratios of a higher class, the same at every cut point" else
          paste0("on the ", s$response$link, " scale (a positive difference: higher classes more likely)"),
        "; class probabilities: mean_type = \"prob\"\n", sep = "")
  }
  cat("\n")
  print(x$table, row.names = FALSE, digits = x$digits)
  cat("\n")
  if (!is.null(x$probs)) {
    cat(if (identical(attr(x$probs, "what"), "mean_class")) "Mean class" else "Class probabilities",
        " at each mean (delta-method SEs):\n", sep = "")
    print(x$probs, row.names = FALSE, digits = x$digits)
    cat("\n")
  }
  p <- x$pairs
  if (is.null(p) && is.numeric(x$sed) && all(c("Max.SED", "Min.SED", "Aveg.SED") %in% names(x$sed))) {
    # no comparisons asked for: the summaries over all pairs of means
    n <- nrow(x$table)
    cat(n*(n - 1)/2, " pairs of means, not tabulated (pairwise = TRUE compares them).\n", sep = "")
    cat("SED: min ", fmt(x$sed[["Min.SED"]]), ", mean ", fmt(x$sed[["Aveg.SED"]]), ", max ", fmt(x$sed[["Max.SED"]]), "\n", sep = "")
    if (is.numeric(x$lsd) && "Aveg.LSD" %in% names(x$lsd)) {
      ddf <- attr(x$lsd, "Degree of freedom")
      r <- unname(x$lsd[c("Min.LSD", "Max.LSD")])
      lab <- paste0("LSD (", s$level * 100, "%", if (!is.null(ddf)) paste0(", ", fmt(ddf), " df"), ")")
      if (isTRUE(attr(x$lsd, "equal"))) cat(lab, ": ", fmt(r[1]), " (equal for all pairs)\n", sep = "") else
        cat(lab, ": min ", fmt(r[1]), ", mean ", fmt(x$lsd[["Aveg.LSD"]]), ", max ", fmt(r[2]), "\n", sep = "")

    }
  }
  if (!is.null(p) && nrow(p) > 0) {
    within <- if (!is.null(s$atvar)) paste0(" within each level of ", paste(s$atvar, collapse = " and ")) else ""
    cat(nrow(p), " pairwise comparisons", within, ".\n", sep = "")
    cat("SED: min ", fmt(min(p$std.error)), ", mean ", fmt(mean(p$std.error)), ", max ", fmt(max(p$std.error)), "\n", sep = "")
    # the LSDs that go with the p-values: adjusted when 'adj' has a critical value
    adj_lsd <- !is.null(p$adj.conf.high)
    lsd <- if (adj_lsd) (p$adj.conf.high - p$adj.conf.low)/2 else p$lsd
    r <- suppressWarnings(range(lsd, na.rm = TRUE))
    lab <- paste0("LSD (", s$level * 100, "%", if (adj_lsd) paste0(", ", s$adjust), ")")
    if (!all(is.finite(r))) NULL else if (diff(r) <= 1e-8 * max(abs(r), 1)) cat(lab, ": ", fmt(r[1]), " (equal for all pairs)\n", sep = "") else
      cat(lab, ": min ", fmt(r[1]), ", mean ", fmt(mean(lsd, na.rm = TRUE)), ", max ", fmt(r[2]),
          " (they differ: use each pair's p-value or LSD)\n", sep = "")

    method <- if (identical(s$df_method, "permutation")) paste0("permutation (", s$nperm, " permutations)") else "t tests"
    cat("p-values: ", method, ", adjusted by '", s$adjust, "'.\n", sep = "")
    if ("group" %in% names(x$table)) {
      cat("Means sharing a letter are not significantly different at the ", s$level * 100, "% level",
          within, ".\n", sep = "")
    }
    cmp <- s$comparison
    if (!is.null(cmp)) {
      cat("Comparison intervals: mean +/- ", fmt(cmp$mult), " SE (", round(100 * cmp$conf_level, 1),
          "% intervals); means differ at the ", s$level * 100, "% level when these do not overlap",
          pm_comparison_caveat(cmp), ". The ", (1 - s$level) * 100, "% confidence limits are not for comparing means.\n",
          sep = "")
    }
  }
  bt <- s$back_transform
  if (!is.null(bt) && "response" %in% names(x$table)) {
    cat("Back-transformed by ", bt$by, ": 'response' is the ", bt$estimate,
        if (identical(bt$average, "response")) ", averaged over the other factors on the response scale" else "",
        ".\n", sep = "")
  }
  if (!is.null(s$covariates)) {
    cv <- s$covariates
    if (nrow(cv) == 1) {
      cat("Covariates fixed at: ", paste(names(cv), fmt(unlist(cv)), sep = " = ", collapse = ", "), "\n", sep = "")
    }
  }
  if (x$nonestimable > 0) cat(x$nonestimable, "mean(s) cannot be estimated.\n")
  pm_print_advice(x)
  invisible(x)
}
