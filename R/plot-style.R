## Colours and sizes of the curve plots of semipred() and of the difference
## plots of semipred(contr = ) and covariatemeans(contr = ), in one place.
## A band is drawn first, its line on top in a darker shade of the same hue
## (a band drawn over its line hides it); the reference line (0, or 1 for
## ratios) is dashed red, as in 1.x, and the stretches where the band leaves
## out the reference are marked under the curve. Groups take a palette that
## colour-blind readers can tell apart (Okabe and Ito).
pm_curve_style_defaults <- list(
  line = "#08519C", band = "#6BAED6",             # one curve: dark and light blue
  diff_line = "#00441B", diff_band = "#74C476",   # a difference of two curves: dark and light green
  ref = "#CB181D",                                # the reference line and the marks: red
  alpha = 0.35, linewidth = 1,
  points = "grey30", point_alpha = 0.6, point_size = 1.3,
  palette = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"),
  surface_points = "red"
)

## The plot settings of semipred() ('plot_args', with the colours above).
pm_semi_plot_defaults <- list(trellis = TRUE, scales = "fixed", plotord = NULL, ci = TRUE, point = TRUE, jitterv = 0,
                              threeD = FALSE)

pm_curve_style <- function(plot_args = list()) {
  if (!is.list(plot_args) || (length(plot_args) && (is.null(names(plot_args)) || !all(nzchar(names(plot_args)))))) {
    stop("'plot_args' must be a named list, e.g. plot_args = list(line = \"black\", band = \"grey70\").", call. = FALSE)
  }
  bad <- setdiff(names(plot_args), names(pm_curve_style_defaults))
  if (length(bad)) {
    stop("'plot_args' can hold ", paste(names(pm_curve_style_defaults), collapse = ", "), "; not ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  utils::modifyList(pm_curve_style_defaults, plot_args)
}

## The plot settings of a call ('plot_args' of covariatemeans() and
## semipred(); names may be abbreviated), with the 'defaults' of the calling
## function for those not given.
pm_plot_args <- function(plot_args = list(), defaults) {
  nms <- names(defaults)
  if (!is.list(plot_args) || (length(plot_args) && (is.null(names(plot_args)) || !all(nzchar(names(plot_args)))))) {
    stop("'plot_args' must be a named list, e.g. plot_args = list(mtitle = \"Yield\").", call. = FALSE)
  }
  k <- pmatch(names(plot_args), nms, duplicates.ok = TRUE)
  if (anyNA(k)) stop("'plot_args' can hold ", paste(nms, collapse = ", "), "; not ",
                     paste(names(plot_args)[is.na(k)], collapse = ", "), call. = FALSE)
  out <- defaults
  for (i in seq_along(k)) out[nms[k[i]]] <- list(plot_args[[i]])
  out
}

## Colour and fill scales for curves by group: the palette when it has enough
## colours, else ggplot2's own.
pm_group_scales <- function(style, n, name = ggplot2::waiver()) {
  # the style's palette while it has enough colours, else evenly spaced hues;
  # the same values and name for colour and fill, so the line and band keys
  # merge into one legend
  pal <- if (n <= length(style$palette)) style$palette[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
  list(ggplot2::scale_colour_manual(values = pal, name = name, drop = FALSE),
       ggplot2::scale_fill_manual(values = pal, name = name, drop = FALSE))
}

## Bars along the bottom of the panel for the stretches of x where a band
## leaves out the reference ('sig': data frame with from, to); each widened by
## half a grid step so a single point shows too. Drawn as a rug, which sits at
## the panel edge without a y value: y = -Inf became NaN on the log axis of
## ratio and odds-ratio plots, and the bars were dropped.
pm_sig_bars <- function(sig, x, style) {
  if (is.null(sig) || !nrow(sig)) return(NULL)
  ux <- sort(unique(x))
  step <- if (length(ux) > 1) stats::median(diff(ux)) / 2 else 0
  # the points of x in each stretch, widened by half a step, as a dense rug
  pts <- unlist(lapply(seq_len(nrow(sig)), function(i) {
    seq(sig$from[i] - step, sig$to[i] + step, length.out = max(2, 200 * (sig$to[i] - sig$from[i] + 2 * step) /
                                                                  max(diff(range(ux)), .Machine$double.eps)))
  }))
  ggplot2::geom_rug(data = data.frame(x = pts), ggplot2::aes(x = .data$x), sides = "b", colour = style$ref,
                    linewidth = 1.2, length = ggplot2::unit(0.02, "npc"), inherit.aes = FALSE)
}
