#' Diagnostic Plots for a Fitted Model
#'
#' @description
#' Four diagnostic plots, drawn with \pkg{ggplot2}, for any model that \code{\link{predictmeans}} handles: \code{lm},
#' \code{aov} (also with \code{Error()} strata), \code{glm}, \code{gls}, \code{lme}, \code{lmer}, \code{glmer},
#' \code{glmmTMB}, \code{mgcv::gam}, \code{\link{semireg}} fits, and other classes with \code{fitted()} and
#' \code{residuals()} methods. With \code{interactive = TRUE} they are shown with \pkg{plotly}: hovering a point gives
#' its row of the data.
#' @usage
#' residplot(model, group = "none", re_level = 1, slope = FALSE, id = 3,
#'     interactive = FALSE, plot = TRUE, dense = c("auto", "points", "bins"), ...)
#' @param model A fitted model (see Description).
#' @param group Optional variable that colours the points: the name (in "quotes") of a variable in the data of the
#'   model, or a vector with one value per observation. The default \code{"none"} uses one colour. The plots have no
#'   legend; \code{interactive = TRUE} shows each point's row and group on hover, and the plots can be given one with
#'   \code{+ ggplot2::theme(legend.position = "right")}.
#' @param re_level For mixed models (\code{level} before 2.0): which random-effect term the random-effect plot shows, 1, 2, ... (in the order of
#'   \code{ranef()}: for \code{lme} the outermost grouping first, for \code{lmer} the innermost first).
#' @param slope For mixed models: \code{TRUE} plots the random slopes (the second column of the random effects) instead
#'   of the random intercepts.
#' @param id Number of the most extreme points labelled with their row names in the data (rows left out for missing
#'   values keep their numbers, so a label is the row to look at; largest absolute residuals; for the
#'   Cook's distance plot the largest distances); \code{FALSE} or 0 labels none. In 1.x \code{id = TRUE} let you
#'   click on points; now \code{TRUE} labels 3 points, and \code{interactive = TRUE} shows every point's row name.
#' @param ... Not used: \code{level}, \code{newwd}, \code{ask} (before 2.0) stop with their 2.0 spelling.
#' @param interactive \code{TRUE}: show the plots with \pkg{plotly} (\code{plotly::ggplotly()}), which must be
#'   installed.
#' @param plot \code{FALSE}: do not draw; only return the plots.
#' @param dense How overlapping points are drawn. \code{"points"} (also \code{"auto"}, the default): every
#'   point; beyond 200 points as open circles, more transparent as their number grows, so dense regions show darker
#'   and no point is hidden. \code{"bins"}: the scatter plots as shaded
#'   counts in a grid of bins (darker: more points; no colour bar), with the points that stand out (the \code{id}
#'   largest and those beyond 3 on the residual scale) drawn on top. The Q-Q plot of more than
#'   2000 residuals draws every point in the tails (outside the central 95\%) and an even selection of the others.
#' @details
#' The plots:
#' \enumerate{
#' \item Residuals against fitted values, with a loess smooth (with fewer than 7 distinct fitted values, the mean
#'   residual at each, joined): look for trends (a missing term, a wrong link) and
#'   for spread that changes with the fitted values.
#' \item Normal Q-Q plot of the residuals: the points should lie near the line.
#' \item Mixed models: normal Q-Q plot of the random effects of term \code{re_level}. \code{gls}: autocorrelation of the
#'   normalized residuals (with approximate 95\% limits). Other models: square root of the absolute residuals
#'   against fitted values (scale-location), for spread that changes with the fitted values.
#' \item \code{lm}, \code{aov} and \code{glm}: Cook's distances. Other models (also \code{gam}: Cook's distances of a
#'   penalized fit are only approximate): fitted against observed values of the
#'   response as modelled (from \code{\link{pm_data}}: the rows used in the fit; the axis is labelled with the response
#'   as written, e.g. \code{log(yield + 1)}, or the observed proportion for \code{cbind(successes, failures)}); for a
#'   0/1 response (also a two-level factor), the observed proportion against the mean fitted probability in tenths of
#'   the fitted values (calibration).
#' }
#' Residuals: standardized for \code{lm} (\code{rstandard()}), normalized for \code{gls} and \code{lme} (the
#' correlation structure removed), scaled by the residual standard deviation for \code{lmer}. For Poisson, binomial and
#' negative binomial models (\code{glm}, \code{glmer}, \code{glmmTMB}, \code{gam}) randomized quantile residuals
#' (Dunn and Smyth, 1996): exactly normal when the model is right, so plots 1 and 2 are read as for a normal model
#' (they use a fixed random seed, so the plots do not change between calls); fitted values are then on the link
#' scale. Other families: standardized deviance residuals. \code{aov} with \code{Error()}: the residuals of the lowest
#' stratum. The plots inform; they do not test. For a table of checks see \code{\link{pm_check}}; for Cook's
#' distances of mixed models \code{\link{CookD}}.
#' @return
#' Invisibly, a list of the four \code{ggplot} objects (class \code{"pm_residplot"}), with attribute \code{"data"}:
#' the data frame of the points (row label, fitted, residual, observed values, group). With \code{interactive = TRUE},
#' the \pkg{plotly} widget.
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996), Randomized quantile residuals, \emph{Journal of Computational and Graphical
#' Statistics}, 5, 236-244.
#' @seealso \code{\link{pm_check}}, \code{\link{CookD}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' residplot(fm, re_level = 2)    # lme: re_level = 2 for random effect "Block:Variety"
#'
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' residplot(fm, group = "Variety") # lmer: re_level = 1 for random effect "Block:Variety"
#'
#' # a Poisson GLM: randomized quantile residuals
#' gm <- glm(breaks ~ wool * tension, family = poisson, data = warpbreaks)
#' p <- residplot(gm, plot = FALSE)
#' p[[2]]                               # one of the plots
#' # interactive plots (hover for each point's row)
#' if (interactive() && requireNamespace("plotly", quietly = TRUE)) residplot(fm, interactive = TRUE)
#' @name residplot
NULL

#' @export
residplot <- function(model, group = "none", re_level = 1, slope = FALSE, id = 3, interactive = FALSE, plot = TRUE,
                      dense = c("auto", "points", "bins"), ...) {
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("residplot", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  model <- pm_fit(model)
  dense <- match.arg(dense)
  nlab <- if (isTRUE(id)) 3L else if (isFALSE(id) || is.null(id)) 0L else as.integer(id)
  d <- pm_resid_data(model)
  d$group <- pm_resid_group(model, group, d)
  re <- pm_resid_ranef(model, re_level, slope)
  plots <- pm_resid_plots(model, d, re, nlab, tooltip = isTRUE(interactive), dense = dense)
  if (isTRUE(interactive)) {
    pm_need("plotly", "interactive = TRUE")
    w <- plotly::subplot(lapply(plots, function(p) plotly::ggplotly(p, tooltip = c("text", "x", "y"))),
                         nrows = 2, titleX = TRUE, titleY = TRUE, margin = 0.07)
    w <- plotly::layout(w, showlegend = FALSE)
    if (plot) print(w)
    return(invisible(w))
  }
  if (plot) pm_print_grid(plots)
  invisible(structure(plots, class = c("pm_residplot", "list"), data = d))
}

## print() of residplot()'s value: the four plots on one page.
#' @exportS3Method print pm_residplot
print.pm_residplot <- function(x, ...) {
  pm_print_grid(unclass(x))
  invisible(x)
}

## Draw ggplots two by two on one page.
pm_print_grid <- function(plots) {
  grid::grid.newpage()
  nr <- ceiling(length(plots) / 2)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, 2)))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = (i - 1) %/% 2 + 1, layout.pos.col = (i - 1) %% 2 + 1))
  }
  grid::upViewport()
  invisible()
}

## The points of the plots: row label, fitted values (x axis of plot 1),
## residuals (standardized, normalized or randomized quantile), fitted and
## observed values on the response scale, with attributes saying which.
pm_resid_data <- function(model) {
  qfam <- pm_qresid_family(model)
  if (inherits(model, "aovlist")) {
    pr <- stats::proj(model)
    low <- pr[[length(pr)]]
    res <- low[, "Residuals"]
    y <- stats::model.response(stats::model.frame(model))
    dfr <- model[[length(model)]]$df.residual
    if (is.null(dfr) || !is.finite(dfr) || dfr <= 0) dfr <- length(res) - 1
    s <- sqrt(sum(res^2) / dfr)
    return(structure(data.frame(label = pm_resid_labels(res, rownames(pm_resid_aov_rows(model, length(res)))), x = y - res, resid = res / s, fitted = y - res,
                                observed = y),
                     xlab = "Fitted values", olab = pm_resid_ylab(model), rlab = paste0("Standardized residuals (", names(pr)[length(pr)],
                                                           " stratum)")))
  }
  fitted_r <- as.vector(stats::fitted(model))
  observed <- pm_resid_observed(model, length(fitted_r))
  if (!is.null(qfam)) {
    r <- pm_qresid(model, qfam, fitted_r)
    x <- as.vector(pm_link_predict(model))
    return(structure(data.frame(label = pm_resid_labels(stats::fitted(model), attr(observed, "rows")), x = x, resid = r, fitted = fitted_r,
                                observed = as.vector(observed)),
                     xlab = "Linear predictor", rlab = "Randomized quantile residuals", olab = attr(observed, "label")))
  }
  r <- if (inherits(model, "gam")) {                   # mgcv: before glm, which gam objects also are
    stats::residuals(model, type = "deviance") / sqrt(model$scale)
  } else if (inherits(model, "glm")) {
    stats::rstandard(model, type = "deviance")
  } else if (inherits(model, "lm")) {
    stats::rstandard(model)
  } else if (inherits(model, c("lme", "gls"))) {
    stats::residuals(model, type = "normalized")
  } else if (inherits(model, "merMod")) {
    if (lme4::isLMM(model)) stats::residuals(model, scaled = TRUE) else stats::residuals(model, type = "pearson")
  } else if (inherits(model, "glmmTMB")) {
    stats::residuals(model, type = "pearson")
  } else {
    rr <- stats::residuals(model)
    rr / stats::sd(rr, na.rm = TRUE)
  }
  glmlike <- inherits(model, c("glm", "glmerMod")) || (inherits(model, "glmmTMB") &&
                                                         !identical(stats::family(model)$family, "gaussian"))
  rlab <- if (inherits(model, "glm")) "Standardized deviance residuals" else
    if (inherits(model, c("lme", "gls"))) "Normalized residuals" else
      if (glmlike) "Pearson residuals" else "Standardized residuals"
  x <- if (glmlike) as.vector(pm_link_predict(model)) else fitted_r
  structure(data.frame(label = pm_resid_labels(r, attr(observed, "rows")), x = x, resid = as.vector(r), fitted = fitted_r,
                       observed = as.vector(observed)),
            xlab = if (glmlike) "Linear predictor" else "Fitted values", rlab = rlab, olab = attr(observed, "label"))
}

## Point labels: the row names of the data the model was fitted to (pm_data(),
## rows with missing values dropped, so the labels are the original row
## numbers), else the names of the residuals, else 1, 2, ...
pm_resid_labels <- function(v, rows = NULL) {
  if (length(rows) == length(v) && !anyDuplicated(rows)) return(as.character(rows))
  nm <- names(v)
  if (is.null(nm) || anyDuplicated(nm)) as.character(seq_along(v)) else nm
}

## The response as modelled (pm_data(): rows used in the fit, proportions for
## cbind(successes, failures), a two-level factor or logical as 0/1), with
## attribute "label": the response as written (e.g. "log(yield + 1)").
pm_resid_observed <- function(model, n) {
  pd <- tryCatch(suppressWarnings(pm_data(model)), error = function(e) NULL)
  r <- if (!is.null(pd)) attr(pd, "response")
  y <- if (!is.null(pd)) pm_resid_numeric(pd$.y)
  if (is.null(y) || length(y) != n || all(is.na(y))) y <- pm_resid_numeric(pm_resid_observed_fit(model, n))
  if (length(y) != n) y <- rep(NA_real_, n)
  lab <- if (is.null(r) || identical(r$expr, "response")) "Observed values" else
    if (grepl("^cbind\\(", r$expr)) paste0("Observed proportion of ", r$name) else r$expr
  rows <- if (!is.null(pd) && nrow(pd) == n) rownames(pd)
  structure(y, label = lab, rows = rows)
}

## The data rows of an aov with Error() (NULL when they do not match).
pm_resid_aov_rows <- function(model, n) {
  pd <- tryCatch(suppressWarnings(pm_data(model)), error = function(e) NULL)
  if (!is.null(pd) && nrow(pd) == n) pd
}

## The response as written, for axis labels.
pm_resid_ylab <- function(model) {
  e <- tryCatch(pm_response_info(model)$expr, error = function(e) "response")
  if (identical(e, "response")) "Observed values" else e
}

## A response as numbers: a two-level factor or logical as 0/1 (second level = 1).
pm_resid_numeric <- function(y) {
  if (is.null(y)) return(NULL)
  if (is.factor(y) || is.character(y)) {
    f <- factor(y)
    return(if (nlevels(f) == 2) as.numeric(f == levels(f)[2]) else rep(NA_real_, length(y)))
  }
  as.numeric(as.vector(y))
}

## The response from the fit itself, when pm_data() cannot find the data.
pm_resid_observed_fit <- function(model, n) {
  y <- tryCatch({
    if (inherits(model, "glm")) model$y
    else if (inherits(model, "merMod")) lme4::getME(model, "y")
    else if (inherits(model, c("lme", "gls"))) stats::fitted(model) + stats::residuals(model, type = "response")
    else {
      yy <- stats::model.response(stats::model.frame(model))
      if (is.matrix(yy) && ncol(yy) == 2) yy[, 1] / rowSums(yy) else yy
    }
  }, error = function(e) NULL)
  if (is.matrix(y) && ncol(y) == 2) y <- y[, 1] / rowSums(y)
  if (!is.factor(y)) y <- as.vector(y)
  if (length(y) != n) rep(NA_real_, n) else y
}

## The linear predictor (link scale, including the random effects).
pm_link_predict <- function(model) {
  if (inherits(model, "glmmTMB")) return(stats::predict(model, type = "link"))
  if (inherits(model, "gam")) return(stats::predict(model, type = "link"))
  stats::predict(model)
}

## The family for randomized quantile residuals, or NULL.
pm_qresid_family <- function(model) {
  if (!inherits(model, c("glm", "glmerMod", "glmmTMB", "gam"))) return(NULL)
  fam <- tryCatch(stats::family(model)$family, error = function(e) "")
  if (fam %in% c("poisson", "binomial")) return(fam)
  if (grepl("^Negative Binomial", fam) || identical(fam, "nbinom2")) return("negbin")
  NULL
}

## Randomized quantile residuals (Dunn and Smyth 1996), with a fixed seed.
pm_qresid <- function(model, fam, mu) {
  if (exists(".Random.seed", envir = globalenv())) {
    seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", seed, envir = globalenv()))
  } else on.exit(rm(".Random.seed", envir = globalenv()))
  set.seed(20261114)
  if (inherits(model, "glmmTMB")) {
    yy <- stats::model.response(stats::model.frame(model))
    if (is.matrix(yy)) { n <- rowSums(yy); y <- yy[, 1] } else {
      n <- rep(1, length(yy))
      y <- if (is.factor(yy)) as.numeric(yy != levels(yy)[1]) else yy
    }
  } else {
    y <- if (inherits(model, "merMod")) lme4::getME(model, "y") else model$y
    n <- stats::weights(model)
    if (is.null(n) || fam != "binomial") n <- rep(1, length(y))
    if (fam == "binomial") y <- round(y * n)
  }
  if (fam == "poisson") {
    a <- stats::ppois(y - 1, mu); b <- stats::ppois(y, mu)
  } else if (fam == "binomial") {
    a <- stats::pbinom(y - 1, n, mu); b <- stats::pbinom(y, n, mu)
  } else {
    theta <- if (inherits(model, "glmmTMB")) stats::sigma(model) else if (inherits(model, "merMod"))
      lme4::getME(model, "glmer.nb.theta") else model$theta
    a <- stats::pnbinom(y - 1, size = theta, mu = mu); b <- stats::pnbinom(y, size = theta, mu = mu)
  }
  u <- a + (b - a) * stats::runif(length(y))
  stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
}

## The grouping variable for the colours: NULL, a column of the model's data, or a vector.
pm_resid_group <- function(model, group, d) {
  if (is.null(group) || identical(group, "none")) return(NULL)
  if (length(group) == nrow(d)) return(factor(group))
  dat <- tryCatch(if (inherits(model, c("lme", "gls"))) pm_get_data(model) else stats::model.frame(model),
                  error = function(e) NULL)
  if (is.character(group) && length(group) == 1 && !is.null(dat) && group %in% names(dat)) {
    g <- dat[[group]]
    if (length(g) != nrow(d) && !is.null(rownames(dat)) && all(d$label %in% rownames(dat))) g <- dat[d$label, group]
    if (length(g) == nrow(d)) return(factor(g))
  }
  stop("'group' must be the name of a variable in the data of the model, or a vector with one value per ",
       "observation (", nrow(d), ").", call. = FALSE)
}

## The random effects of term number 're_level' (column 2 with slope = TRUE), or NULL.
pm_resid_ranef <- function(model, level, slope) {
  re <- tryCatch({
    if (inherits(model, "lme")) {
      r <- nlme::ranef(model)
      if (is.data.frame(r)) r <- list(r)
      r
    } else if (inherits(model, "merMod")) lme4::ranef(model)
    else if (inherits(model, "glmmTMB")) lme4::ranef(model)[["cond"]]
  }, error = function(e) NULL)
  if (is.null(re) || !length(re)) return(NULL)
  if (level > length(re)) stop("'re_level' is larger than the number of random-effect terms (", length(re), ").",
                               call. = FALSE)
  r <- as.data.frame(re[[level]])
  col <- if (slope) 2 else 1
  if (ncol(r) < col) stop("there is no random slope in random-effect term ", level, ".", call. = FALSE)
  list(value = r[[col]], label = rownames(r), term = names(re)[level], what = names(r)[col])
}

## The four plots.
pm_resid_plots <- function(model, d, re, nlab, tooltip = FALSE, dense = "auto") {
  th <- ggplot2::theme_bw(base_size = 10)
  n_all <- sum(is.finite(d$resid) & is.finite(d$x))
  # many points: smaller and more transparent, or shaded bins
  bins <- dense == "bins" && !tooltip
  rp_title <- function(x) pm_wrap_title(x, 10, panels = 2)      # (the panels are drawn two side by side)
  col <- "steelblue"
  grp <- !is.null(d$group)
  top <- function(v, k) if (k > 0) utils::head(order(abs(v), decreasing = TRUE), k) else integer(0)
  # points; with tooltip = TRUE the row label is the plotly hover text
  # ('text' is not a ggplot2 aesthetic: only plotly uses it)
  pts <- function(data, x, y, colour = grp) {
    # transparency and size for the number of points in this panel
    n <- nrow(data)
    p_alpha <- if (n <= 200) 0.9 else max(0.15, min(0.9, 20 / sqrt(n)))
    p_size <- if (n <= 500) 1.8 else if (n <= 5000) 1.4 else 1.1
    p_shape <- if (n <= 200) 16 else 1             # many points: open circles, overlaps show as darker rings
    m <- if (colour && tooltip) ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data$group, text = .data$label)
      else if (colour) ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data$group)
      else if (tooltip) ggplot2::aes(x = .data[[x]], y = .data[[y]], text = .data$label)
      else ggplot2::aes(x = .data[[x]], y = .data[[y]])
    suppressWarnings(if (colour) ggplot2::geom_point(m, data = data, alpha = p_alpha, size = p_size, shape = p_shape) else
      ggplot2::geom_point(m, data = data, colour = col, alpha = p_alpha, size = p_size, shape = p_shape))
  }
  # a scatter: points, or (bins) shaded counts with the points that stand out on top
  scatter <- function(data, x, y, out = integer(0)) {
    if (!bins) return(pts(data, x, y))
    list(ggplot2::geom_bin_2d(ggplot2::aes(x = .data[[x]], y = .data[[y]]), data = data, bins = 70),
         ggplot2::scale_fill_gradient(low = "#dbe7f3", high = "#08306b", trans = "log10", guide = "none"),
         if (length(out)) ggplot2::geom_point(ggplot2::aes(x = .data[[x]], y = .data[[y]]), data = data[out, , drop = FALSE],
                                              colour = "#08306b", size = 0.9))
  }
  lab <- function(data, x, y, rows) {
    if (!length(rows)) return(NULL)
    ggplot2::geom_text(ggplot2::aes(x = .data[[x]], y = .data[[y]], label = .data$label), data = data[rows, ],
                       size = 3, vjust = -0.6, colour = "grey20", check_overlap = TRUE)
  }
  xlab <- attr(d, "xlab")
  rlab <- attr(d, "rlab")
  ok <- is.finite(d$resid) & is.finite(d$x)
  dd <- d[ok, ]
  if (grp) dd$group <- d$group[ok]

  # a loess smooth of y against the fitted values; with few distinct fitted
  # values (a model of factors only) the mean at each of them, joined
  trend <- function(dd, y) {
    if (length(unique(dd$x)) >= 7) {
      big <- nrow(dd) > 2000
      ggplot2::geom_smooth(ggplot2::aes(x = .data$x, y = .data[[y]]), data = dd, method = if (big) "gam" else "loess",
                           formula = if (big) y ~ s(x, bs = "cs") else y ~ x, se = FALSE, colour = "firebrick",
                           linewidth = 0.7)
    } else {
      ggplot2::stat_summary(ggplot2::aes(x = .data$x, y = .data[[y]]), data = dd, fun = mean, geom = "line",
                            colour = "firebrick", linewidth = 0.6)
    }
  }

  # 1. residuals vs fitted
  out1 <- union(top(dd$resid, nlab), which(abs(dd$resid) > 3))
  p1 <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    scatter(dd, "x", "resid", out1) +
    trend(dd, "resid") +
    lab(dd, "x", "resid", top(dd$resid, nlab)) +
    ggplot2::labs(title = rp_title("Residuals vs fitted"), x = xlab, y = rlab) + th
  if (grp) p1 <- p1 + ggplot2::labs(colour = NULL)

  # 2. normal Q-Q of the residuals
  qq <- dd
  qq$theoretical <- stats::qnorm(stats::ppoints(nrow(qq)))[rank(qq$resid, ties.method = "first")]
  # the reference line through the quartiles of all residuals (before any thinning)
  qline <- stats::quantile(qq$resid, c(0.25, 0.75)); qn <- stats::qnorm(c(0.25, 0.75))
  qlab <- top(qq$resid, nlab)
  if (nrow(qq) > 2000) {
    # every point in the tails, an even selection of the central 95%
    r <- rank(qq$resid, ties.method = "first")
    tail <- r <= 0.025 * nrow(qq) | r > 0.975 * nrow(qq)
    keep <- tail | (r %% ceiling(nrow(qq) / 2000) == 0)
    keep[qlab] <- TRUE
    qlab <- match(qlab, which(keep))
    qq <- qq[keep, , drop = FALSE]
  }
  slope <- diff(qline) / diff(qn)
  p2 <- ggplot2::ggplot() + ggplot2::geom_abline(intercept = qline[1] - slope * qn[1], slope = slope, colour = "grey40") +
    pts(qq, "theoretical", "resid") + lab(qq, "theoretical", "resid", qlab) +
    ggplot2::labs(title = rp_title("Normal Q-Q of residuals"), x = "Normal quantiles", y = rlab) + th +
    ggplot2::theme(legend.position = "none")

  # 3. random effects / ACF (gls) / scale-location
  if (!is.null(re)) {
    rd <- data.frame(value = re$value, label = re$label)
    rd$theoretical <- stats::qnorm(stats::ppoints(nrow(rd)))[rank(rd$value, ties.method = "first")]
    ql <- stats::quantile(rd$value, c(0.25, 0.75)); s3 <- diff(ql) / diff(qn)
    p3 <- ggplot2::ggplot() + ggplot2::geom_abline(intercept = ql[1] - s3 * qn[1], slope = s3, colour = "grey40") +
      pts(rd, "theoretical", "value", colour = FALSE) + lab(rd, "theoretical", "value", top(rd$value, nlab)) +
      ggplot2::labs(title = rp_title(paste0("Normal Q-Q of random effects: ", re$what, " | ", re$term)),
                    x = "Normal quantiles", y = "Predicted random effects") + th
  } else if (inherits(model, "gls")) {
    form <- attr(model$modelStruct$corStruct, "formula")
    if (is.null(form)) form <- ~1
    ac <- tryCatch(nlme::ACF(model, form = form, resType = "normalized"), error = function(e) NULL)
    if (is.null(ac) || anyNA(ac$ACF)) {
      p3 <- ggplot2::ggplot() + ggplot2::labs(title = rp_title("Autocorrelation: not available")) + th
    } else {
      ac <- as.data.frame(ac)
      nu <- attr(nlme::ACF(model, form = form), "n.used")
      lim <- if (length(nu) == nrow(ac)) stats::qnorm(0.975) / sqrt(nu) else rep(NA_real_, nrow(ac))
      ac$lim <- lim
      p3 <- ggplot2::ggplot(ac, ggplot2::aes(x = .data$lag)) +
        ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, y = 0, yend = .data$ACF), colour = col, linewidth = 0.8) +
        ggplot2::geom_line(ggplot2::aes(y = .data$lim), linetype = 2, colour = "grey40", na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(y = -.data$lim), linetype = 2, colour = "grey40", na.rm = TRUE) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::labs(title = rp_title("Autocorrelation of normalized residuals"), x = "Lag", y = "Autocorrelation") + th
    }
  } else {
    dd$root <- sqrt(abs(dd$resid))
    p3 <- ggplot2::ggplot() + scatter(dd, "x", "root", out1) + trend(dd, "root") +
      ggplot2::labs(title = rp_title("Scale-location"), x = xlab, y = paste0("sqrt(|", rlab, "|)")) + th +
      ggplot2::theme(legend.position = "none")
  }

  # 4. Cook's distances (lm, aov, glm) / calibration (0/1 response) / fitted vs observed
  if (inherits(model, "lm") && !inherits(model, c("aovlist", "gam"))) {
    cd <- data.frame(case = seq_len(nrow(d)), cook = as.vector(stats::cooks.distance(model)), label = d$label)
    p4 <- ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(x = .data$case, xend = .data$case, y = 0, yend = .data$cook), data = cd,
                            colour = col) +
      pts(cd, "case", "cook", colour = FALSE) + lab(cd, "case", "cook", top(cd$cook, nlab)) +
      ggplot2::labs(title = rp_title("Cook's distance"), x = "Observation", y = "Cook's distance") + th
  } else {
    od <- d[is.finite(d$observed) & is.finite(d$fitted), ]
    binary <- nrow(od) > 0 && all(od$observed %in% c(0, 1))
    if (binary) {
      br <- unique(stats::quantile(od$fitted, seq(0, 1, 0.1)))
      bin <- cut(od$fitted, br, include.lowest = TRUE)
      cal <- data.frame(fitted = tapply(od$fitted, bin, mean), observed = tapply(od$observed, bin, mean),
                        n = as.vector(table(bin)))
      p4 <- ggplot2::ggplot(cal, ggplot2::aes(x = .data$fitted, y = .data$observed)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey40") +
        ggplot2::geom_point(ggplot2::aes(size = .data$n), colour = col) +
        ggplot2::scale_size_area(max_size = 4, guide = "none") + ggplot2::coord_equal(xlim = 0:1, ylim = 0:1) +
        ggplot2::labs(title = rp_title("Calibration (tenths of fitted values)"), x = "Mean fitted probability",
                      y = "Observed proportion") + th
    } else if (nrow(od)) {
      p4 <- ggplot2::ggplot() + scatter(od, "observed", "fitted") +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "grey40") +
        ggplot2::labs(title = rp_title("Fitted vs observed"), x = attr(d, "olab"), y = "Fitted values") + th +
        ggplot2::theme(legend.position = "none")
    } else {
      p4 <- ggplot2::ggplot() + ggplot2::labs(title = rp_title("Fitted vs observed: response not available")) + th
    }
  }
  # labels at the edge of a panel are not cut off
  out <- list(resid_fitted = p1, qq = p2, third = p3, fourth = p4)
  lapply(out, function(p) {
    p$coordinates$clip <- "off"          # (keeps the coordinates of each panel, e.g. the calibration's 0-1)
    p + ggplot2::theme(legend.position = "none")      # no legends: a clean page of four plots
  })
}
