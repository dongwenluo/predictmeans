#'Produce a mean value plot for predicted means
#'
#'
#'@param plot_mt Data frame of the mean table which is the output of \code{mean_table} or \code{mean_table$Table} from
#'  function \code{predictmeans}.
#'
#'@param x_var Name (in "quotes") for the x axis variable on the plot.
#'@param y_var Name (in "quotes") for the y axis variable on the plot.
#'@param col_var Name (in "quotes") for the color group variable on the plot.
#'@param panel_var Name (in "quotes") for the panel variable on the plot.
#'@param title The text for the title.
#'@param xlab The title of the respective axis (for xlab() or ylab()) or of the plot (for ggtitle()).
#'@param ylab The title of the respective axis (for xlab() or ylab()) or of the plot (for ggtitle()).
#'@param scales Should scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x",
#'  "free_y") in a trellis graph?
#'@param bar_value The length of the bar with default value of zero.
#'@param bar_label Should bar label be "Aveg.LSD" (the default) or "Aveg.SED"?
#'@param level A significant level for calculating LSD, CI etc. The default value is 0.05.
#'@param basesz Base font size, given in pts.
#'@param line An option for drawing a line chart, or dot chart. The default is TRUE.
#'
#'@author Dongwen Luo
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
#' predm_out <- predictmeansN(mod, "A:B:time", atvar=c("time"), adj="BH", trans=function(x) x^(1/1.5))
#'
#' mean_plot(predm_out$mean_table$Table, x_var="time", y_var="Mean", col_var="A", panel_var ="B",
#' bar_value = 118.06)
#' # Use plot method
#' plot(predm_out, mod_df = ATP, resp_name = "ATP")
#'
#'@importFrom ggplot2 lims
#'@importFrom rlang .data
#'@export

mean_plot <- function(plot_mt,
                      x_var,
                      y_var,
                      col_var = NULL,
                      panel_var = NULL,
                      title = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      scales = c("fixed", "free", "free_x", "free_y"),
                      bar_value = 0,
                      bar_label = c("Aveg.LSD", "Aveg.SED"),
                      level = 0.05,
                      basesz = 12,
                      line = TRUE) {
  plot_mt <- na.omit(plot_mt)
  yMin <- min(plot_mt[, y_var])
  yMax <- max(plot_mt[, y_var])
  offSet <- 0.25 * (yMax - yMin)
  xlab <- ifelse(is.null(xlab) ||
                   any(c("NULL", "") == xlab), x_var, xlab)
  ylab <- ifelse(is.null(ylab) ||
                   any(c("NULL", "") == ylab), y_var, ylab)
  bar_label <- match.arg(bar_label)
  barMin <- yMin
  barMax <- barMin + bar_value
  scales <- match.arg(scales)

  if (bar_value > 0) {
    barMin <- yMin
    barMax <- barMin + bar_value
    plot_mt[, x_var] <- factor(plot_mt[, x_var], levels = c(bar_label, levels(plot_mt[, x_var])))
  }

  if (!is.null(x_var) && is.null(col_var) && is.null(panel_var)) {
    if (is.null(title) || any(c("NULL", "") == title)) {
      if (bar_value > 0) {
        title <- paste(
          "Predicted means for '",
          ylab,
          "' vs '",
          xlab,
          "' with ",
          bar_label,
          " Bar at ",
          level,
          " Significant Level\n",
          sep = ""
        )
      } else {
        title <- paste("Predicted means for '", ylab, "' vs '", xlab, "'\n", sep =
                         "")
      }
    }
    meanPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]], group =
                                      1)) +
      labs(
        title = paste(title, "\n", sep = ""),
        x = paste("\n", xlab, sep = ""),
        y = paste(ylab, "\n", sep = "")
      ) +
      lims(x = levels(plot_mt[, x_var]), y = c(yMin - offSet, max(yMax + offSet, yMin + bar_value + offSet))) +
      geom_point(colour = "red", size = 2) +
      geom_errorbar(
        aes(ymax = barMax, ymin = barMin, x = bar_label),
        width = 0.15,
        linewidth = 0.8,
        colour = "blue"
      ) +
      theme_bw(basesz)
    if (line) {
      meanPlot <- meanPlot + geom_line(linewidth = 0.5)
    }
  }

  if (!is.null(x_var) && !is.null(col_var) && is.null(panel_var)) {
    if (is.null(title) || any(c("NULL", "") == title)) {
      if (bar_value > 0) {
        title <- paste(
          "Predicted means for '",
          ylab,
          "' vs '",
          xlab,
          "' by '",
          col_var,
          "' with ",
          bar_label,
          " Bar at ",
          level,
          " Significant Level\n",
          sep = ""
        )
      } else {
        title <- paste("Predicted means for '",
                       ylab,
                       "' vs '",
                       xlab,
                       "' by '",
                       col_var,
                       "'",
                       sep = "")
      }
    }

    meanPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]], group =
                                      .data[[col_var]], col = .data[[col_var]])) +
      labs(
        title = paste(title, "\n", sep = ""),
        x = paste("\n", xlab, sep = ""),
        y = paste(ylab, "\n", sep = "")
      ) +
      lims(x = levels(plot_mt[, x_var]), y = c(yMin - offSet, max(yMax + offSet, yMin + bar_value + offSet))) +
      geom_point(size = 2) +
      geom_errorbar(
        aes(ymax = barMax, ymin = barMin, x = bar_label),
        width = 0.15,
        linewidth = 0.8,
        colour = "blue"
      ) +
      guides(col = guide_legend(title = col_var)) +
      theme_bw(basesz)
    if (line) {
      meanPlot <- meanPlot + geom_line(aes(linetype = .data[[col_var]], col =
                                             .data[[col_var]]), linewidth = 0.96) +
        guides(linetype = guide_legend(title = col_var))
    }
  }

  if (!is.null(x_var) && !is.null(col_var) && !is.null(panel_var)) {
    if (is.null(title) || any(c("NULL", "") == title)) {
      if (bar_value > 0) {
        title <- paste(
          "Predicted means for '",
          ylab,
          "' vs '",
          xlab,
          "' by '",
          col_var,
          "' for each '",
          panel_var,
          "'\n with ",
          bar_label,
          " Bar at ",
          level,
          " Significant Level\n",
          sep = ""
        )
      } else {
        title <- paste(
          "Predicted means for '",
          ylab,
          "' vs '",
          xlab,
          "' by '",
          col_var,
          "' for each '",
          panel_var,
          "'\n",
          sep = ""
        )
      }
    }

    meanPlot <- ggplot(plot_mt, aes(
      .data[[x_var]],
      .data[[y_var]],
      group = factor(.data[[col_var]]),
      col = factor(.data[[col_var]])
    )) +
      labs(
        title = paste(title, "\n", sep = ""),
        x = paste("\n", xlab, sep = ""),
        y = paste(ylab, "\n", sep = "")
      ) +
      lims(x = levels(plot_mt[, x_var]),
           y = c(
             yMin - 0.5 * offSet,
             max(yMax + 0.5 * offSet, yMin + bar_value + 0.5 * offSet)
           )) +
      geom_errorbar(
        aes(ymax = barMax, ymin = barMin, x = bar_label),
        width = 0.15,
        linewidth = 0.8,
        colour = "blue"
      ) +
      geom_point(size = 2) +
      facet_grid(as.formula(paste("~", panel_var)), scales = scales) +
      guides(group = guide_legend(col_var)) +
      guides(col = guide_legend(col_var)) +
      theme_bw(basesz)
    if (line) {
      meanPlot <- meanPlot + geom_line(aes(linetype = .data[[col_var]], col =
                                             .data[[col_var]]), linewidth = 0.8) +
        guides(linetype = guide_legend(title = col_var))
    }

  }
  return(meanPlot)
}
