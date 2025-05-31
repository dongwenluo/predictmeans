#' Level Plot of a Matrix of p-values.
#'
#' Creates a plot of p-values of pairwise comparisons.
#'
#'
#' @param pmatrix A matrix with p-values from pairwise comparisons. (This is a
#' lower triangle matrix.)
#' @param level The level of p-value to be highlighted. Default is 0.05.
#' @param mtitle The main title in the graph.
#' @param xylabel The x and y labels in the graph.
#' @param margin A value for specifying x and y margins in the graph. The
#' default value is 5.
#' @param legendx A value for specifying x coordinate of legend. The default
#' value is 0.73.
#' @param newwd A logical variable to indicate whether to print graph in a new
#' window. The default is FALSE.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
#'   library(predictmeans)
#'   set.seed(2013)
#'   pvalues <- runif(28)
#'   pmatrix <- matrix(0,8,8)
#'   pmatrix[lower.tri(pmatrix)] <- pvalues
#'   round(pmatrix, 4)
#'   PMplot(pmatrix)
#'
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", barplot=TRUE)
#'   PMplot(predictout$p_valueMatrix)
#'
#' @importFrom ggplot2 aes coord_fixed element_blank element_rect element_text
#' @importFrom ggplot2 ggplot geom_tile labs margin
#' @importFrom ggplot2 scale_fill_manual scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 theme theme_minimal
#'
#' @export
PMplot <- function(pmatrix, level=0.05, mtitle=NULL, xylabel=NULL, margin=5, legendx=0.73, newwd=FALSE) {

  # Simple block diagonal function to replace adiag
  # block_diag <- function(mat_list) {
  # if (length(mat_list) == 1) return(mat_list[[1]])

  # # Calculate dimensions
  # nrows <- sapply(mat_list, nrow)
  # ncols <- sapply(mat_list, ncol)
  # total_rows <- sum(nrows)
  # total_cols <- sum(ncols)

  # # Create result matrix filled with NA
  # result <- matrix(NA, nrow = total_rows, ncol = total_cols)

  # # Fill in blocks
  # row_start <- 1
  # col_start <- 1
  # for (i in seq_along(mat_list)) {
  # row_end <- row_start + nrows[i] - 1
  # col_end <- col_start + ncols[i] - 1
  # result[row_start:row_end, col_start:col_end] <- mat_list[[i]]
  # row_start <- row_end + 1
  # col_start <- col_end + 1
  # }

  # return(result)
  # }

  # Handle matrix input
  if (is.matrix(pmatrix)) {
    nr <- nrow(pmatrix)
    pmatrix[upper.tri(pmatrix, diag = TRUE)] <- NA
    if (is.null(rownames(pmatrix))) {
      rnpltm <- as.character(1:nrow(pmatrix))
    } else {
      rnpltm <- rownames(pmatrix)
    }
  }

  # Handle list input using base R block diagonal
  if (is.list(pmatrix)) {
    for (i in 1:length(pmatrix)) {
      pmatrix[[i]][upper.tri(pmatrix[[i]], diag = TRUE)] <- NA
    }
    # pmatrix <- block_diag(pmatrix)
    pmatrix <- do.call(adiag, c(pmatrix, pad=NA))
    nr <- nrow(pmatrix)
    if (is.null(xylabel)) {
      rnpltm <- as.character(1:nrow(pmatrix))
    } else {
      rnpltm <- xylabel
    }
  }

  # Check minimum size
  if (nr <= 3) {
    cat("\nThere is no plot for p-values matrix less than six values!\n")
    return(invisible(NULL))
  }

  # Set default title
  if (is.null(mtitle)) {
    mtitle <- "Level Plot of p-value Matrix"
  }

  # Prepare matrix for plotting (flip vertically to match base R image)
  pltmm <- pmatrix[nr:1, ]

  # Create significance categories
  if (level == 0.05) {
    # Create cut categories
    pltm_cut <- cut(as.numeric(pltmm),
                    breaks = c(-0.1, 0.01, 0.05, 0.1, 1),
                    labels = c("<= 0.01", "0.01 < p <= 0.05", "0.05 < p <= 0.1", "> 0.1"),
                    include.lowest = TRUE)
    fact_level <- rev(c("<= 0.01", "0.01 < p <= 0.05", "0.05 < p <= 0.1", "> 0.1"))

    # Define colors
    pcolr <- c("<= 0.01" = "#0D0DFF",
               "0.01 < p <= 0.05" = "#5D5DFF",
               "0.05 < p <= 0.1" = "#A1A1FF",
               "> 0.1" = "#E4E4FF")
  } else {
    # Simple significant/non-significant
    pltm_cut <- cut(as.numeric(pltmm),
                    breaks = c(-0.1, level, 1),
                    labels = c("significant", "insignificant"),
                    include.lowest = TRUE)

    fact_level <- rev(c("significant", "insignificant"))

    # Define colors
    pcolr <- c("significant" = "#0D0DFF",
               "insignificant" = "#A1A1FF")
  }

  # Convert to long format for ggplot using base R
  # # x_coords <- rep(1:ncol(pltmm), each = nrow(pltmm))
  # # y_coords <- rep(1:nrow(pltmm), times = ncol(pltmm))
  # # values <- as.numeric(pltmm)
  # # categories <- factor(as.character(pltm_cut), levels=fact_level)
  x <- rep(1:ncol(pltmm), each = nrow(pltmm))
  y <- rep(1:nrow(pltmm), times = ncol(pltmm))
  values <- as.numeric(pltmm)
  category <- factor(as.character(pltm_cut), levels=fact_level)

  # Create data frame
  plot_data <- data.frame(
    x = x,
    y = y,
    value = values,
    category = category,
    stringsAsFactors = FALSE
  )

  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$category), ]

  # Create the plot
  p <- ggplot(plot_data, aes(x = x, y = y, fill = category)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = pcolr,
                      name = if (level == 0.05) "p-value" else paste("At", round(level, 4), "level")) +
    scale_x_continuous(
      breaks = 1:length(rnpltm),
      labels = rnpltm,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = 1:length(rnpltm),
      labels = rev(rnpltm),
      expand = c(0, 0)
    ) +
    labs(
      title = mtitle,
      x = "",
      y = ""
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x = element_text(angle = if (max(nchar(rnpltm))/6 > 0.5) 90 else 0,
                                 hjust = if (max(nchar(rnpltm))/6 > 0.5) 1 else 0.5,
                                 vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(legendx, 0.99),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),
      legend.margin = margin(6, 6, 6, 6),
      aspect.ratio = 1
    ) +
    coord_fixed()

  # Display plot
  # if (newwd && interactive()) {
  # # Open new graphics device if requested and in interactive session
  # if (.Platform$OS.type == "windows") {
  # windows()
  # } else if (Sys.info()["sysname"] == "Darwin") {
  # quartz()
  # } else {
  # x11()
  # }
  # }

  if (newwd) {
    dev.new()
  }
 # print(p)
  return(invisible(p))
}
