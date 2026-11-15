#' Level Plot of a Matrix of p-values.
#'
#' @description
#' Plots the p-values of pairwise comparisons as a square grid of tiles, one
#' for each pair of levels, like a correlation plot: the lower triangle only,
#' the tiles coloured by the size of the p-value.
#' @usage
#' PMplot(pmatrix, level=0.05, mtitle=NULL, xylabel=NULL, margin=5, legendx=0.73,
#'   newwd=FALSE, values=FALSE, basesz=12)
#' @param pmatrix A matrix with p-values from pairwise comparisons (the lower
#' triangle is used), or a list of such matrices (e.g. \code{p_valueMatrix} of
#' \code{\link{predictmeans}} with \code{atvar}), which are drawn as blocks
#' down the diagonal of one plot.
#' @param level The level of p-value to be highlighted. With the default 0.05
#' the tiles show four bands (up to 0.01, 0.05, 0.1 and above); with any other
#' value two (up to \code{level} and above).
#' @param mtitle The main title in the graph.
#' @param xylabel The x and y labels in the graph: one for each row of the
#' matrix (or of the list's matrices together). The default uses the row
#' names of the matrix, else the row numbers.
#' @param margin Not used since 2.0 (ggplot2 sizes the margins from the
#' labels); kept so that 1.x calls still run.
#' @param legendx The x coordinate of the left side of the legend, as a
#' proportion of the plot's width. The default value is 0.73; the legend sits
#' in the empty upper triangle.
#' @param newwd A logical variable to indicate whether to print graph in a new window. The default is FALSE.
#' @param values Whether to print each p-value on its tile (the default is
#' \code{FALSE}). Useful for up to about 15 levels.
#' @param basesz The base font size. The default is 12.
#' @return A \code{ggplot} object (invisibly \code{NULL} when the matrix has
#' no pair of levels). Printing it draws the plot; it can be changed with the
#' usual ggplot2 functions. Before 2.0 the plot was drawn with base graphics
#' and nothing was returned: inside a function or a loop, wrap the call in
#' \code{print()}.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   set.seed(2013)
#'   pvalues <- runif(28)
#'   pmatrix <- matrix(0,8,8)
#'   pmatrix[lower.tri(pmatrix)] <- pvalues
#'   round(pmatrix, 4)
#'   PMplot(pmatrix)
#'   PMplot(pmatrix, values = TRUE)
#'
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH")
#'   PMplot(predictout$p_valueMatrix)
#' @name PMplot
NULL

PMplot <- function(pmatrix, level=0.05, mtitle=NULL, xylabel=NULL, margin=5, legendx=0.73, newwd=FALSE,
                   values=FALSE, basesz=12) {
  blocks <- NULL
  if (is.data.frame(pmatrix)) pmatrix <- as.matrix(pmatrix)
  if (is.list(pmatrix)) {
    mats <- lapply(pmatrix, as.matrix)
    blocks <- vapply(mats, nrow, integer(1))
    labs <- unlist(lapply(mats, function(m) if (is.null(rownames(m))) rep(NA_character_, nrow(m)) else rownames(m)))
    for (i in seq_along(mats)) mats[[i]][upper.tri(mats[[i]], diag = TRUE)] <- NA
    P <- do.call(adiag, c(unname(mats), pad = NA))
    if (anyNA(labs)) labs <- as.character(seq_len(nrow(P)))
  } else if (is.matrix(pmatrix)) {
    P <- pmatrix
    if (nrow(P) != ncol(P)) stop("'pmatrix' must be a square matrix (or a list of square matrices).", call. = FALSE)
    P[upper.tri(P, diag = TRUE)] <- NA
    labs <- if (is.null(rownames(P))) as.character(seq_len(nrow(P))) else rownames(P)
  } else stop("'pmatrix' must be a matrix of p-values or a list of such matrices.", call. = FALSE)
  if (!is.null(xylabel)) {
    if (length(xylabel) != nrow(P)) stop("'xylabel' must have one label for each of the ", nrow(P), " rows.", call. = FALSE)
    labs <- as.character(xylabel)
  }
  if (!any(is.finite(P))) {
    message("There are no p-values to plot.")
    return(invisible(NULL))
  }
  if (is.null(mtitle) || identical(mtitle, "NULL") || identical(mtitle, "")) mtitle <- "Level Plot of p-value Matrix"
  if (newwd) grDevices::dev.new()
  pm_ptile_plot(P, labs, level = level, title = mtitle, blocks = blocks, values = values, legendx = legendx,
                basesz = basesz)
}

## The bands of p-values shown by the level plots: four with level 0.05, else
## two. Returns list(breaks, labels (smallest p first), colours, title).
pm_p_bands <- function(level) {
  if (isTRUE(all.equal(level, 0.05))) {
    list(breaks = c(0.01, 0.05, 0.1), labels = c("p <= 0.01", "0.01 < p <= 0.05", "0.05 < p <= 0.1", "p > 0.1"),
         colours = c("#0D0DFF", "#5D5DFF", "#A1A1FF", "#E4E4FF"), title = "p-value")
  } else {
    lv <- format(round(level, 4), scientific = FALSE, drop0trailing = TRUE)
    list(breaks = level, labels = c(paste("p <=", lv), paste("p >", lv)), colours = c("#0D0DFF", "#C9C9FF"),
         title = "p-value")
  }
}

## Level plot of the p-values in the lower triangle of the square matrix P
## (row i, column j: the pair of levels i and j; other cells NA), the one
## drawing code of PMplot() and autoplot(type = "pvalue"). Square tiles, the
## first level at the top left, the legend in the empty upper triangle.
## labs: one label per row; blocks: the sizes of the groups down the diagonal
## (comparisons within the levels of 'atvar'), separated by lines.
pm_ptile_plot <- function(P, labs, level = 0.05, title = NULL, legend_title = NULL, blocks = NULL,
                          values = FALSE, legendx = 0.73, basesz = 12) {
  n <- nrow(P)
  k <- which(is.finite(P))
  d <- data.frame(row = (k - 1L) %% n + 1L, col = (k - 1L) %/% n + 1L, p = P[k])
  b <- pm_p_bands(level)
  d$band <- factor(b$labels[findInterval(d$p, b$breaks, left.open = TRUE) + 1L], levels = b$labels)
  # rows and columns without a pair at the edges (the first row and the last
  # column of a triangle) are left out, so the tiles fill a square
  rr <- range(d$row)
  cr <- range(d$col)
  long <- max(nchar(labs[cr[1]:cr[2]])) > 3
  lw <- if (n > 60) 0 else if (n > 30) 0.2 else 0.5         # no grid lines between very many tiles
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$band)) +
    ggplot2::geom_tile(colour = if (lw > 0) "white" else NA, linewidth = lw) +
    ggplot2::scale_fill_manual(values = stats::setNames(b$colours, b$labels), drop = FALSE,
                               name = if (is.null(legend_title)) b$title else legend_title) +
    ggplot2::scale_x_continuous(breaks = cr[1]:cr[2], labels = labs[cr[1]:cr[2]]) +
    ggplot2::scale_y_reverse(breaks = rr[1]:rr[2], labels = labs[rr[1]:rr[2]]) +
    ggplot2::coord_fixed(xlim = cr + c(-0.5, 0.5), ylim = rev(rr + c(-0.5, 0.5)), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = pm_wrap_title(title, basesz)) +
    ggplot2::theme_bw(base_size = basesz) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = if (long) ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5) else
                     ggplot2::element_text(),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.justification = c(0, 1),
                   legend.background = ggplot2::element_rect(fill = "white", colour = "grey40", linewidth = 0.3)) +
    pm_legend_inside(c(legendx, 0.99))
  if (length(blocks) > 1) {                # lines between the groups
    at <- cumsum(blocks)[-length(blocks)] + 0.5
    g <- g + ggplot2::geom_hline(yintercept = at, colour = "grey60", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = at, colour = "grey60", linewidth = 0.3)
  }
  if (isTRUE(values)) {
    d$txt <- ifelse(d$p < 0.001, "<.001", sub("^0\\.", ".", formatC(d$p, format = "f", digits = 3)))
    d$dark <- as.integer(d$band) <= length(b$labels)/2          # white text on the dark bands
    side <- max(diff(rr), diff(cr)) + 1
    g <- g + ggplot2::geom_text(data = d, ggplot2::aes(label = .data$txt, colour = .data$dark),
                                size = basesz/3.5*min(1, 9/side), show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none")
  }
  g
}

## A legend inside the panel at 'pos' (ggplot2 3.5.0 moved this setting).
pm_legend_inside <- function(pos) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::theme(legend.position = "inside", legend.position.inside = pos)
  } else ggplot2::theme(legend.position = pos)
}
