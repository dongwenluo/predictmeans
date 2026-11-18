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

## Colour and fill scales for curves by group: the palette when it has enough
## colours, else ggplot2's own.
pm_group_scales <- function(style, n) {
  if (n > length(style$palette)) return(NULL)
  list(ggplot2::scale_colour_manual(values = style$palette[seq_len(n)]),
       ggplot2::scale_fill_manual(values = style$palette[seq_len(n)]))
}

## Bars along the bottom of the panel for the stretches of x where a band
## leaves out the reference ('sig': data frame with from, to); each widened by
## half a grid step so a single point shows too.
pm_sig_bars <- function(sig, x, style) {
  if (is.null(sig) || !nrow(sig)) return(NULL)
  step <- if (length(unique(x)) > 1) stats::median(diff(sort(unique(x)))) / 2 else 0
  seg <- data.frame(from = sig$from - step, to = sig$to + step)
  ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$from, xend = .data$to, y = -Inf, yend = -Inf),
                        colour = style$ref, linewidth = 2.5, inherit.aes = FALSE)
}
