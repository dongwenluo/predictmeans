## Plots of the curve results of covariatemeans() and semipred(), built from
## their tables (means, slopes, pairs), settings and the data rows kept in
## attr(x, "plot_data"). Nothing is computed here; the 1.x look is kept
## (theme_bw, bands under their lines, one panel per level by default).

#' Plots of Curves, Surfaces and Their Differences
#'
#' @description
#' \code{autoplot()} draws a result of \code{\link{covariatemeans}} or \code{\link{semipred}} as a \pkg{ggplot2}
#' object (a \pkg{plotly} widget for \code{type = "3d"}): the curves with their bands and the data, surfaces as
#' filled contours, slopes (\code{semipred(deriv = )}) with the line at zero, or the difference of two curves
#' (\code{contr}). It reads the tables of the result only, so the layout can be changed without computing again.
#' \code{plot = TRUE} in \code{covariatemeans()} and \code{semipred()} prints it, and \code{$plt} of the result
#' holds it; their \code{plot_args} (superseded) give the same arguments under their 1.x names. Those settings are
#' kept in the result (\code{attr(x, "plot_args")}) and are the defaults of \code{autoplot(x)}; arguments given to
#' \code{autoplot()} win (\code{autoplot(x, interval = "ci")} draws the bands of a result made with \code{ci = FALSE}).
#'
#' The curves and their bands share one legend. When each panel is one group (one factor, the default layout) the
#' strips name the groups and there is no legend. With more than 8 groups the colours are evenly spaced hues, and
#' with more than 12 the legend is set in columns with small keys; \code{panel = character(0)} puts all curves in
#' one panel.
#' @usage
#' \method{autoplot}{pm_curves}(object, type = NULL, group = NULL, panel = NULL,
#'          scales = c("fixed", "free", "free_x", "free_y"), scale = NULL,
#'          data = TRUE, jitter = 0, interval = c("ci", "none"), title = NULL,
#'          style = list(), ...)
#' @param object A result of \code{\link{covariatemeans}} or \code{\link{semipred}}.
#' @param type \code{"lines"} (one covariate: curves or slopes), \code{"surface"} (two covariates: filled contours,
#'   one panel per level), \code{"3d"} (two covariates: a \pkg{plotly} 3-D surface, one widget with a scene per level;
#'   needs \pkg{plotly}) or \code{"difference"} (a result with \code{contr}). The default is \code{"difference"} with
#'   \code{contr}, else \code{"lines"} or \code{"surface"} by the number of covariates.
#' @param group,panel Names of factors of \code{modelterm}: \code{group} colours the curves, \code{panel} makes the
#'   panels (one factor: \code{facet_wrap()}; two: \code{facet_grid()} rows by columns; more: the first gives the rows,
#'   the others together the columns). Factors in neither are joined to \code{group}, so each combination has its own
#'   curve. The defaults: with one factor, colour and one panel per level (as 1.x); with more, the first factor
#'   colours and the others make the panels. \code{panel = character(0)} draws one panel. These replace
#'   \code{plotord} and \code{trellis} of 1.x.
#' @param scales The axes of the panels: \code{"fixed"} (default), \code{"free"}, \code{"free_x"}, \code{"free_y"}.
#' @param scale \code{"model"} (the link or the transformed response, column \code{estimate}) or \code{"response"}
#'   (the back-transformed curves, columns \code{response...}). Default: \code{"response"} when the curves are
#'   back-transformed. For \code{contr}: \code{"model"} shows the difference (as a ratio or odds ratio for a log or
#'   logit), \code{"response"} the difference of the back-transformed curves; the default is \code{contr_scale} of
#'   the call. Slopes are on the model scale.
#' @param data Draw the data (\code{TRUE}, jittered by \code{jitter}): the original response on the response scale,
#'   the modelled one on the model scale (none for a link on the model scale). Not drawn for slopes and differences.
#' @param jitter Jitter of the data points, in x and y (units of the axes).
#' @param interval \code{"ci"}: the band the result has (simultaneous with \code{band = "simultaneous"});
#'   \code{"none"}.
#' @param title The title; \code{NULL} gives the default (\code{covariatemeans()}: that of 1.x), \code{FALSE} none.
#' @param style Colours and sizes, any of: \code{line} and \code{band} (one curve), \code{diff_line} and
#'   \code{diff_band} (a difference), \code{ref} (the reference line and the marks where a band leaves it out),
#'   \code{alpha}, \code{linewidth}, \code{points}, \code{point_alpha}, \code{point_size}, \code{palette} (curves by
#'   group; ggplot2's own beyond its length), \code{surface_points}.
#' @param \dots Not used.
#' @return A \code{ggplot} object, or for \code{type = "3d"} a \pkg{plotly} widget. Add to it with \pkg{ggplot2}
#'   (\code{+ labs(y = "Rate")}).
#' @seealso \code{\link{covariatemeans}}, \code{\link{semipred}}, \code{\link{tidy.pm_curves}}
#' @examples
#'   library(predictmeans)
#'   data(Oats, package = "nlme")
#'   Oats$B <- factor(Oats$Block %in% c("I", "II", "III"), labels = c("late", "early"))
#'   fm <- lm(yield ~ nitro * Variety * B, data = Oats)
#'   cm <- covariatemeans(fm, "Variety:B", covariate = "nitro")
#'   ggplot2::autoplot(cm)                                   # colour Variety, panels B
#'   ggplot2::autoplot(cm, group = "B", panel = "Variety")
#'   ggplot2::autoplot(cm, panel = character(0), data = FALSE)
#' @name autoplot.pm_curves
#' @exportS3Method ggplot2::autoplot pm_curves
autoplot.pm_curves <- function(object, type = NULL, group = NULL, panel = NULL,
                               scales = c("fixed", "free", "free_x", "free_y"), scale = NULL, data = TRUE, jitter = 0,
                               interval = c("ci", "none"), title = NULL, style = list(), ...) {
  # the plot settings of the call (plot_args, recorded in the result) are the
  # defaults; arguments given here win
  pa <- attr(object, "plot_args")
  if (length(pa)) {
    given <- names(match.call(expand.dots = FALSE))[-1]
    for (nm in setdiff(intersect(names(pa), c("type", "group", "panel", "scales", "scale", "data", "jitter",
                                              "interval", "title", "style")), given)) assign(nm, pa[[nm]])
  }
  s <- object$settings
  sty <- pm_curve_style(style)
  scales <- match.arg(scales)
  interval <- match.arg(interval)
  cv <- s$covariate
  has_pairs <- !is.null(object$pairs)
  if (is.null(type)) type <- if (has_pairs) "difference" else if (length(cv) > 1) "surface" else "lines"
  type <- match.arg(type, c("lines", "surface", "3d", "difference"))
  if (has_pairs && type != "difference") {
    stop("This result holds the difference of two curves (contr): type = \"difference\".", call. = FALSE)
  }
  if (!has_pairs && type == "difference") stop("type = \"difference\" needs a result with 'contr'.", call. = FALSE)
  if (type %in% c("surface", "3d") && length(cv) != 2) stop("type = \"", type, "\" needs two covariates.", call. = FALSE)
  if (type == "lines" && length(cv) != 1) stop("type = \"lines\" needs one covariate; use type = \"surface\".",
                                               call. = FALSE)
  if (type == "difference") return(pm_curve_diff_plot(object, scale, interval, title, sty))
  slopes <- !is.null(object$slopes)
  tab <- if (slopes) object$slopes else object$means
  vars <- if (setequal(s$modelterm, cv)) character(0) else strsplit(s$modelterm, ":", fixed = TRUE)[[1]]
  back <- "response" %in% names(tab)
  if (is.null(scale)) scale <- if (back && !slopes) "response" else "model"
  scale <- match.arg(scale, c("model", "response"))
  if (scale == "response" && !back) stop("scale = \"response\": the curves are not back-transformed.", call. = FALSE)
  simult <- "adj.conf.low" %in% names(tab)
  pre <- if (scale == "response") "response." else ""
  pd <- tab[c(vars, cv)]
  pd$.y <- tab[[if (scale == "response") "response" else "estimate"]]
  pd$.lo <- tab[[paste0(pre, if (simult) "adj.conf.low" else "conf.low")]]
  pd$.hi <- tab[[paste0(pre, if (simult) "adj.conf.high" else "conf.high")]]
  lay <- pm_curve_layout(vars, group, panel)
  pd <- pm_curve_layout_cols(pd, lay)
  obs <- if (isTRUE(data) && !slopes) pm_curve_obs(object, scale, c(vars, cv), lay)
  if (!is.null(obs)) for (k in intersect(c(".group", ".row", ".col"), names(obs))) {
    obs[[k]] <- factor(as.character(obs[[k]]), levels = levels(pd[[k]]))
  }
  ylab <- if (slopes) {
    paste0(if (s$deriv == 1) "Slope" else "Curvature", " (d", if (s$deriv == 2) "2", "/d", cv, if (s$deriv == 2) "2",
           if (!is.null(s$response$link)) ", link scale", ")")
  } else pm_scale_label(s, scale, "axis")
  if (is.null(title) && inherits(object, "covariatemeans")) {
    title <- paste0("Fitted", if (isTRUE(data)) " and observed", " relationship",
                    if (interval == "ci") paste0(" with ", (1 - s$level)*100, "% ", if (simult) "simultaneous band" else "CI"))
  }
  if (isFALSE(title)) title <- NULL
  if (type == "3d") return(pm_curve_3d(pd, obs, cv, ylab, sty))
  if (type == "surface") {
    plt <- ggplot2::ggplot(pd, ggplot2::aes(x = .data[[cv[1]]], y = .data[[cv[2]]])) +
      ggplot2::geom_contour_filled(ggplot2::aes(z = .data$.y)) +
      ggplot2::labs(x = cv[1], y = cv[2], fill = ylab) +
      ggplot2::theme_bw()
    if (!is.null(obs)) plt <- plt + ggplot2::geom_point(data = obs, ggplot2::aes(x = .data[[cv[1]]], y = .data[[cv[2]]]),
                                                        colour = sty$surface_points, size = 1.5, inherit.aes = FALSE)
  } else {
    grp <- !is.null(pd$.group)
    plt <- ggplot2::ggplot(pd, ggplot2::aes(x = .data[[cv]])) +
      ggplot2::labs(x = cv, y = ylab) + ggplot2::theme_bw()
    if (!is.null(obs)) {
      plt <- plt + if (grp && !is.null(obs$.group)) {
        ggplot2::geom_point(data = obs, ggplot2::aes(x = .data[[cv]], y = .data$.y, colour = .data$.group),
                            alpha = sty$point_alpha, size = sty$point_size, inherit.aes = FALSE,
                            position = ggplot2::position_jitter(width = jitter, height = jitter))
      } else {
        ggplot2::geom_point(data = obs, ggplot2::aes(x = .data[[cv]], y = .data$.y), colour = sty$points,
                            alpha = sty$point_alpha, size = sty$point_size, inherit.aes = FALSE,
                            position = ggplot2::position_jitter(width = jitter, height = jitter))
      }
    }
    if (interval == "ci") {
      plt <- plt + if (grp) {
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$.lo, ymax = .data$.hi, fill = .data$.group, group = .data$.group),
                             alpha = sty$alpha, colour = NA)
      } else ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$.lo, ymax = .data$.hi), fill = sty$band, alpha = sty$alpha,
                                  colour = NA)
    }
    plt <- plt + if (grp) {
      ggplot2::geom_line(ggplot2::aes(y = .data$.y, colour = .data$.group, group = .data$.group), linewidth = sty$linewidth)
    } else ggplot2::geom_line(ggplot2::aes(y = .data$.y), colour = sty$line, linewidth = sty$linewidth)
    if (grp) {
      ng <- nlevels(pd$.group)
      plt <- plt + pm_group_scales(sty, ng, lay$group_label) +
        ggplot2::labs(colour = lay$group_label)
      if (interval == "ci") plt <- plt + ggplot2::labs(fill = lay$group_label)
      if (setequal(lay$group, lay$panel)) {
        plt <- plt + ggplot2::theme(legend.position = "none")   # each panel is one group: its strip names it
      } else if (ng > 12) {
        # many groups: a compact legend at the side, in columns of up to 15 keys
        nc <- ceiling(ng / 15)
        plt <- plt + ggplot2::guides(colour = ggplot2::guide_legend(ncol = nc),
                                     fill = ggplot2::guide_legend(ncol = nc)) +
          ggplot2::theme(legend.text = ggplot2::element_text(size = 7),
                         legend.key.size = ggplot2::unit(0.75, "lines"),
                         legend.key.spacing.y = ggplot2::unit(0, "pt"))
      }
    }
    if (slopes) {
      plt <- plt + ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = sty$ref)
      if (!length(vars)) plt <- plt + pm_sig_bars(object$significant, pd[[cv]], sty)
    }
  }
  plt <- plt + pm_curve_facets(lay, scales)
  if (type == "lines") {
    # small panels: fewer breaks, and labels that would overlap are dropped
    nb <- if (length(lay$panel)) 4 else 6
    plt <- plt + ggplot2::scale_x_continuous(n.breaks = nb, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(n.breaks = nb, guide = ggplot2::guide_axis(check.overlap = TRUE))
  }
  if (!is.null(title)) {
    plt <- plt + ggplot2::labs(title = paste0(pm_wrap_title(title, 11), "\n"))
  }
  plt
}

## Which factors colour the curves and which make the panels. vars: the
## factors of modelterm; group, panel: the user's (NULL: default).
pm_curve_layout <- function(vars, group = NULL, panel = NULL) {
  if (!length(vars)) return(list(group = character(0), panel = character(0)))
  chk <- function(v, what) {
    bad <- setdiff(v, vars)
    if (length(bad)) stop("'", what, "' must be among the factors of modelterm: ", paste(vars, collapse = ", "),
                          "; not ", paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  chk(group, "group")
  chk(panel, "panel")
  if (is.null(group)) group <- setdiff(vars, panel)[1]
  if (is.na(group[1])) group <- character(0)
  if (is.null(panel)) panel <- if (length(vars) == 1) vars else setdiff(vars, group)
  # factors in neither: joined to the colours
  group <- c(group, setdiff(vars, c(group, panel)))
  list(group = group, panel = panel, group_label = paste(group, collapse = ":"))
}

## The columns .group, .row, .col of a table from the layout.
pm_curve_layout_cols <- function(d, lay) {
  lab <- function(v) {
    if (length(v) == 1) return(factor(d[[v]], levels = unique(d[[v]])))
    x <- do.call(paste, c(lapply(d[v], as.character), sep = ":"))
    factor(x, levels = unique(x))
  }
  if (length(lay$group)) d$.group <- lab(lay$group)
  p <- lay$panel
  if (length(p) == 1) d$.col <- lab(p)
  if (length(p) >= 2) {
    d$.row <- lab(p[1])
    d$.col <- lab(p[-1])
  }
  d
}

pm_curve_facets <- function(lay, scales) {
  p <- lay$panel
  if (!length(p)) return(NULL)
  if (length(p) == 1) return(ggplot2::facet_wrap(ggplot2::vars(.data$.col), scales = scales))
  ggplot2::facet_grid(rows = ggplot2::vars(.data$.row), cols = ggplot2::vars(.data$.col), scales = scales)
}

## The data under the curves: attr(x, "plot_data") on the scale shown, with
## the layout columns; NULL when there are none to show.
pm_curve_obs <- function(x, scale, cols, lay) {
  d <- attr(x, "plot_data")
  if (!is.data.frame(d) || !all(cols %in% names(d))) return(NULL)
  r <- x$settings$response
  link <- !is.null(r) && identical(r$scale, "link")
  if (scale == "model" && link) return(NULL)
  y <- if (scale == "response") d$.y_orig else d$.y
  if (is.factor(y) || is.logical(y)) y <- as.numeric(as.factor(y)) - 1
  if (!is.numeric(y) || all(is.na(y))) return(NULL)
  o <- d[cols]
  o$.y <- as.numeric(y)
  for (v in intersect(c(lay$group, lay$panel), names(o))) o[[v]] <- as.character(o[[v]])
  o <- pm_curve_layout_cols(o, lay)
  o
}

## A 3-D surface with plotly: one scene per level (subplot).
pm_curve_3d <- function(pd, obs, cv, zlab, sty) {
  pm_need("plotly", "3-D plots (type = \"3d\")")
  grp <- if (!is.null(pd$.group) || !is.null(pd$.col)) {
    interaction(pd[intersect(c(".group", ".row", ".col"), names(pd))], sep = ":", drop = TRUE, lex.order = TRUE)
  } else factor(rep("", nrow(pd)))
  ogrp <- if (!is.null(obs) && (!is.null(obs$.group) || !is.null(obs$.col))) {
    interaction(obs[intersect(c(".group", ".row", ".col"), names(obs))], sep = ":", drop = TRUE, lex.order = TRUE)
  }
  fml <- function(v) stats::as.formula(paste0("~`", v, "`"))
  levs <- levels(grp)
  k <- length(levs)
  sc <- function(i) paste0("scene", if (i > 1) i)
  # one widget, one scene per level side by side in a grid (plotly::subplot()
  # does not place 3-D scenes)
  nc <- ceiling(sqrt(k)); nr <- ceiling(k / nc)
  p <- plotly::plot_ly()
  lay <- list()
  for (i in seq_len(k)) {
    g <- levs[i]
    p <- plotly::add_trace(p, data = pd[grp == g, , drop = FALSE], x = fml(cv[1]), y = fml(cv[2]), z = ~.y,
                           intensity = ~.y, type = "mesh3d", opacity = 0.7, name = g, scene = sc(i),
                           showscale = i == 1, inherit = FALSE)
    if (!is.null(obs)) {
      o <- if (is.null(ogrp)) obs else obs[ogrp == g, , drop = FALSE]
      p <- plotly::add_trace(p, x = o[[cv[1]]], y = o[[cv[2]]], z = o$.y, mode = "markers", type = "scatter3d",
                             marker = list(size = 2.5, color = sty$surface_points), name = paste(g, "data"),
                             scene = sc(i), inherit = FALSE)
    }
    r <- (i - 1) %/% nc; cc <- (i - 1) %% nc
    lay[[sc(i)]] <- list(xaxis = list(title = cv[1]), yaxis = list(title = cv[2]), zaxis = list(title = zlab),
                         domain = list(x = c(cc, cc + 1) / nc, y = 1 - c(r + 1, r) / nr))
    if (k > 1) lay$annotations <- c(lay$annotations, list(list(text = g, showarrow = FALSE, xref = "paper",
      yref = "paper", x = (cc + 0.5) / nc, y = 1 - r / nr, xanchor = "center", yanchor = "bottom")))
  }
  do.call(plotly::layout, c(list(p), lay))
}

## The plot of the difference of two curves from the pairs table (one
## covariate: pm_diff_plot(); two: filled contours with marks where the band
## leaves out zero).
pm_curve_diff_plot <- function(x, scale, interval, title, sty) {
  s <- x$settings
  p <- x$pairs
  cv <- s$covariate
  simult <- "adj.conf.low" %in% names(p)
  kind <- intersect(c("ratio", "odds.ratio"), names(p))
  kind <- if (length(kind)) kind else "none"
  back <- "response" %in% names(p)
  if (is.null(scale)) scale <- if (is.null(s$contr_scale)) "model" else s$contr_scale
  scale <- match.arg(scale, c("model", "response"))
  lim <- function(pre) {
    a <- if (simult) "adj.conf." else "conf."
    list(p[[paste0(pre, a, "low")]], p[[paste0(pre, a, "high")]])
  }
  m <- lim("")
  tab <- data.frame(xvar = if (length(cv) == 1) p[[cv]] else seq_len(nrow(p)), diff = p$estimate, ses = p$std.error,
                    Df = p$df, LL = m[[1]], UL = m[[2]])
  tab$includes0 <- tab$LL <= 0 & tab$UL >= 0
  if (kind != "none") {
    tab[[kind]] <- p[[kind]]
    tab[[paste0(kind, ".LL")]] <- exp(tab$LL)
    tab[[paste0(kind, ".UL")]] <- exp(tab$UL)
  }
  if (back) {
    r <- lim("response.")
    tab$response.diff <- p$response
    tab$response.LL <- r[[1]]
    tab$response.UL <- r[[2]]
  }
  tab$significant <- !tab$includes0
  if (interval == "none") tab$LL <- tab$UL <- tab$response.LL <- tab$response.UL <- NULL
  if (length(cv) == 2) {
    resp <- back && scale == "response"
    pd <- p[cv]
    pd$.y <- if (resp) p$response else p$estimate
    sig <- pd[tab$significant %in% TRUE, , drop = FALSE]
    lev <- s$contrast
    ttl <- if (is.null(title)) paste("Factor", s$modelterm, "level", lev[1], "vs level", lev[2],
                                     "with Points at the Area out of", (1 - s$level) * 100,
                                     if (simult) "% simultaneous band" else "% CI") else if (!isFALSE(title)) title
    plt <- ggplot2::ggplot(pd, ggplot2::aes(x = .data[[cv[1]]], y = .data[[cv[2]]])) +
      ggplot2::geom_contour_filled(ggplot2::aes(z = .data$.y)) +
      ggplot2::geom_point(data = sig, ggplot2::aes(x = .data[[cv[1]]], y = .data[[cv[2]]]), colour = sty$surface_points,
                          size = 2.5, shape = 22, inherit.aes = FALSE) +
      ggplot2::labs(x = cv[1], y = cv[2], fill = paste0("Difference (", lev[1], " - ", lev[2], ")",
                                                        if (resp) ", response scale")) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if (!is.null(ttl)) plt <- plt + ggplot2::ggtitle(pm_wrap_title(ttl, 11))
    return(plt)
  }
  if (interval == "none") {
    tab$LL <- tab$UL <- tab$diff
    if (back) tab$response.LL <- tab$response.UL <- tab$response.diff
    if (kind != "none") tab[[paste0(kind, ".LL")]] <- tab[[paste0(kind, ".UL")]] <- tab[[kind]]
  }
  res <- list(table = tab, significant = x$significant, kind = if (back) kind else "none", back = back)
  test <- s$test
  test <- if (!is.null(test$F.ratio)) c(numDF = test$df1, denDF = test$df2, F = test$F.ratio, p.value = test$p.value) else
    if (!is.null(test$statistic)) c(maxt = test$statistic, p.value = test$p.value)
  bt <- s$back_transform
  adjusted <- !is.null(bt) && grepl("bias adjusted", bt$estimate)
  typical <- !is.null(bt) && grepl("typical cluster", bt$estimate)
  plt <- pm_diff_plot(res, cv, s$contrast, s$modelterm, s$level, if (simult) "simultaneous" else "pointwise", scale,
                      test, mtitle = if (!isFALSE(title)) title, typical = typical, adjusted = adjusted, style = sty)
  if (isFALSE(title)) plt <- plt + ggplot2::labs(title = NULL)
  d <- s$deriv
  if (!is.null(d) && d > 0) {
    q <- if (d == 1) "slopes" else "curvatures"
    plt <- plt + ggplot2::labs(y = paste0("The difference of the ", q, " (d", if (d == 2) "2", "/d", cv,
                                          if (d == 2) "2", ")"))
    if (!is.null(plt$labels$title) && is.null(title)) {
      plt <- plt + ggplot2::ggtitle(sub("^Difference between", paste("Difference between the", q, "of"),
                                        plt$labels$title))
    }
  }
  plt
}

## plot_args of covariatemeans() and semipred() (1.x names, superseded) as
## arguments of autoplot.pm_curves(): trellis, plotord, mtitle, ci, point,
## jitterv, threeD, scales and the colours. 'pa': pm_plot_args() output;
## 'vars': the factors of modelterm. Returns the list for do.call().
pm_curve_plot_args <- function(pa, vars, contr_scale = NULL) {
  out <- list()
  po <- pa$plotord
  if (!is.null(po) && !all(po %in% c("NULL", ""))) {
    if (length(po) != length(vars)) stop("plotord must be a vector of length ", length(vars), ".", call. = FALSE)
    if (is.character(po)) {
      i <- match(po, vars)
      if (anyNA(i)) stop("The entries ", paste(po[is.na(i)], collapse = ", "), " do not match any factor of ",
                         paste(vars, collapse = ":"), ".", call. = FALSE)
      po <- i
    }
    out$group <- vars[po][1]
    out$panel <- vars[po][-1]
  }
  if (isFALSE(pa$trellis)) out$panel <- character(0)
  if (!is.null(pa$mtitle) && !all(pa$mtitle %in% c("NULL", ""))) out$title <- pa$mtitle
  if (isFALSE(pa$ci)) out$interval <- "none"
  if (isFALSE(pa$point)) out$data <- FALSE
  if (!is.null(pa$jitterv)) out$jitter <- pa$jitterv
  if (isTRUE(pa$threeD)) out$type <- "3d"
  if (!is.null(pa$scales)) out$scales <- pa$scales
  sty <- pa[intersect(names(pa), names(pm_curve_style_defaults))]
  if (length(sty)) out$style <- sty
  if (!is.null(contr_scale)) out$scale <- contr_scale
  out
}

## The plot of a new curve result ($plt; printed when 'draw'): autoplot()
## with the mapped plot_args. A plot that cannot be made leaves $plt NULL with
## a warning; the tables are complete.
pm_curve_add_plot <- function(out, args, draw = FALSE) {
  plt <- tryCatch(do.call(autoplot.pm_curves, c(list(out), args)), error = function(e) {
    if (draw || identical(args$type, "3d")) stop(conditionMessage(e), call. = FALSE)
    warning("The plot could not be made (", conditionMessage(e), "); the tables are complete.", call. = FALSE)
    NULL
  })
  if (draw && !is.null(plt)) suppressWarnings(print(plt))
  if (length(args)) attr(out, "plot_args") <- args
  i <- match("plt", names(out))
  if (!is.na(i)) {
    if (is.null(plt)) out <- out[-i] else out[i] <- list(plt)
  }
  out
}

## The data rows kept for the plots: the factors and covariates with .y and
## .y_orig (pm_data()); NULL when they cannot be found.
pm_curve_plot_data <- function(model, cols, data = NULL, trans = NULL) {
  d <- tryCatch(suppressWarnings(pm_data(model, data = data, trans = trans)), error = function(e) NULL)
  if (!is.data.frame(d)) return(NULL)
  keep <- intersect(c(cols, ".y", ".y_orig"), names(d))
  out <- d[keep]
  rownames(out) <- NULL
  out
}

## 1.x scripts: covariatemeans()$pltdf is $pred_df since 2.0.
pm_once <- new.env(parent = emptyenv())
pm_pltdf_name <- function(name) {
  if (identical(name, "pltdf")) {
    if (is.null(pm_once$pltdf)) {
      pm_once$pltdf <- TRUE
      message("'pltdf' is called 'pred_df' since predictmeans 2.0 (as in semipred()); x$pltdf still gives it.")
    }
    return("pred_df")
  }
  name
}

#' @export
`$.pm_curves` <- function(x, name) {
  unclass(x)[[pm_pltdf_name(name), exact = FALSE]]
}

#' @export
`[[.pm_curves` <- function(x, i, ...) {
  if (is.character(i) && length(i) == 1) i <- pm_pltdf_name(i)
  unclass(x)[[i, ...]]
}
