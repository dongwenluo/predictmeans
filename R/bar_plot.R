#' Produce a bar plot for predicted means
#'
#'
#' @param plot_mt Data frame of the mean table which is the output of
#'   \code{mean_table} or \code{mean_table$Table} from function
#'   \code{predictmeans}.
#'
#' @param x_var Name (in "quotes") for the x axis variable on the plot.
#' @param y_var Name (in "quotes") for the y axis variable on the plot.
#' @param se_var Name (in "quotes") for the variable with SE values.
#' @param col_var Name (in "quotes") for the color group variable on the plot.
#' @param panel_var Name (in "quotes") for the panel variable on the plot.
#' @param title The text for the title.
#' @param xlab The title of the respective axis (for xlab() or ylab()) or of the
#'   plot (for ggtitle()).
#' @param ylab The title of the respective axis (for xlab() or ylab()) or of the
#'   plot (for ggtitle()).
#' @param scales Should scales be fixed ("fixed", the default), free ("free"),
#'   or free in one dimension ("free_x", "free_y") in a trellis graph?
#' @param basesz Base font size, given in pts.
#'
#' @author Dongwen Luo
#'
#' @examples
#'
#' library(predictmeans)
#' data(ATP, package="predictmeans")
#' str(ATP)
#'
#' mod <- lmer(ATP^1.5 ~ A*B*time+(1| heart), ATP)
#' # residplot(mod)
#' #--------------------------------------------------------
#' predm_out <- predictmeansN(mod, "A:B:time", atvar=c("time"), adj="BH",
#' trans=function(x) x^(1/1.5))
#'
#' bar_plot(predm_out$mean_table$Table, x_var="time", y_var="Mean",
#' se_var="SE", col_var="A", panel_var ="B")
#' # Use plot method
#' plot(predm_out, mod_df = ATP, resp_name = "ATP")
#'
#' @importFrom ggplot2 element_line facet_grid geom_errorbarh
#' @importFrom ggplot2 geom_bar geom_errorbar guides guide_legend
#' @importFrom ggplot2 position_dodge scale_fill_brewer xlim ylim
#' @importFrom rlang .data
#' @importFrom stats aov dist
#' @export

bar_plot <- function(plot_mt,
                     x_var,
                     y_var,
                     se_var,
                     col_var = NULL,
                     panel_var = NULL,
                     title = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     scales = c("fixed", "free", "free_x", "free_y"),
                     basesz = 12) {
  # custom_colors <- c(
  #   "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  #   "#66a61e", "#e6ab02", "#a6761d", "#666666",
  #   "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"
  # )

  plot_mt <- na.omit(plot_mt)
  xlab <- ifelse(is.null(xlab) ||
                   any(c("NULL", "") == xlab), x_var, xlab)
  ylab <- ifelse(is.null(ylab) ||
                   any(c("NULL", "") == ylab), y_var, ylab)
  yMin <- min(0, (min(plot_mt[[y_var]]) - max(plot_mt[[se_var]])) * 1.2)
  yMax <- max(0, (max(plot_mt[[y_var]]) + max(plot_mt[[se_var]])) * 1.2)
  barMin <- (plot_mt[[y_var]] - plot_mt[[se_var]]) * (plot_mt[[y_var]] < 0) +
    pmax(plot_mt[[y_var]] - plot_mt[[se_var]], 0) * (plot_mt[[y_var]] >= 0)
  barMax <- (plot_mt[[y_var]] + plot_mt[[se_var]]) * (plot_mt[[y_var]] > 0) +
    pmin(plot_mt[[y_var]] + plot_mt[[se_var]], 0) * (plot_mt[[y_var]] <= 0)
  dodge <- position_dodge(width = 0.9)
  scales <- match.arg(scales)

  if (!is.null(x_var) && is.null(col_var) && is.null(panel_var)) {
    if (is.null(title) || any(c("NULL", "") == title)) {
      title <- paste("Predicted means for '",
                     ylab,
                     "' vs '",
                     xlab,
                     "' with SE Bars",
                     sep = "")
    }
    barPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]])) +
      labs(
        title = paste(title, "\n", sep = ""),
        x = paste("\n", xlab, sep = ""),
        y = paste(ylab, "\n", sep = "")
      ) +
      ylim(c(yMin, yMax)) +
      geom_point() +
      geom_bar(
        stat = "identity",
        fill = "papayawhip",
        alpha = 0.8,
        colour = "darkgreen"
      ) +
      geom_errorbar(
        aes(ymin = barMin, ymax = barMax),
        position = dodge,
        width = 0.25,
        colour = "blue"
      ) +
      theme_bw(basesz) +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid = element_blank()
      )
  }

  if (!is.null(x_var) && !is.null(col_var) && is.null(panel_var)) {
    if (is.null(title) || any(c("NULL", "") == title)) {
      title <- paste("Predicted means for '",
                     ylab,
                     "' vs '",
                     xlab,
                     "' by '",
                     col_var,
                     "' with SE Bars",
                     sep = "")
    }

    barPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]], group = .data[[col_var]], fill = .data[[col_var]])) +
      geom_point(position = dodge) +
      geom_bar(stat = "identity", position = dodge) +
      labs(
        title = paste(title, "\n", sep = ""),
        x = paste("\n", xlab, sep = ""),
        y = paste(ylab, "\n", sep = "")
      ) +
      ylim(c(yMin, yMax)) +
      geom_errorbar(
        aes(ymin = barMin, ymax = barMax),
        position = dodge,
        width = 0.25,
        colour = "blue"
      ) +
      theme_bw(basesz) +
      # scale_fill_brewer()+
      scale_fill_brewer(palette = "Set3") +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid = element_blank()
      ) +
      theme(legend.position = "top") + guides(fill = guide_legend(title = col_var))
  }

  if (!is.null(x_var) && !is.null(col_var) && !is.null(panel_var)) {
    if (is.null(title) || any(c("NULL", "") == title)) {
      title <- paste(
        "Predicted means for '",
        ylab,
        "' vs '",
        xlab,
        "' by '",
        col_var,
        "' for each '",
        panel_var,
        "' with SE Bars",
        sep = ""
      )
    }

    barPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]], group = .data[[col_var]], fill = .data[[col_var]])) +
      geom_point(position = dodge) +
      geom_bar(stat = "identity", position = dodge) +
      labs(
        title = paste(title, "\n", sep = ""),
        x = paste("\n", xlab, sep = ""),
        y = paste(ylab, "\n", sep = "")
      ) +
      ylim(c(yMin, yMax)) +
      geom_errorbar(
        aes(ymin = barMin, ymax = barMax),
        position = dodge,
        width = 0.25,
        colour = "blue"
      ) +
      facet_grid(as.formula(paste("~", panel_var)), scales = scales) +
      theme_bw(basesz) +
      # scale_fill_brewer()+
      scale_fill_brewer(palette = "Set3") +
      #  scale_fill_manual(values = custom_colors) +
      theme(axis.line = element_line(), panel.grid = element_blank()) +
      theme(legend.position = "right") + guides(fill = guide_legend(title = col_var))
  }
  return(barPlot)
}
