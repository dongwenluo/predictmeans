#' Level Plot of a Matrix of p-values
#'
#' @description
#' Plots the p-values of pairwise comparisons as a square grid of tiles, one
#' for each pair of levels, like a correlation plot: the lower triangle only,
#' the tiles coloured by the size of the p-value.
#' @usage
#' PMplot(pmatrix, conf_level=0.95, mtitle=NULL, xylabel=NULL, values=FALSE, basesz=12, ...)
#' @param pmatrix A matrix with p-values from pairwise comparisons (the lower
#' triangle is used), or a list of such matrices (e.g. \code{p_valueMatrix} of
#' \code{\link{predictmeans}} with \code{atvar}), which are drawn as one square
#' panel each, titled by the list's names (the levels of \code{atvar}).
#' @param conf_level The p-values below \code{1 - conf_level} are highlighted (\code{level = 0.05} before 2.0).
#' With the default 0.95 the tiles show four bands (up to 0.01, 0.05, 0.1 and above); with any other value two
#' (up to \code{1 - conf_level} and above).
#' @param mtitle The main title in the graph.
#' @param xylabel The x and y labels in the graph: one for each row of the
#' matrix (or of the list's matrices together). The default uses the row
#' names of the matrix, else the row numbers.
#' @param ... Not used: \code{level}, \code{margin}, \code{legendx}, \code{newwd} (before 2.0) stop with their 2.0
#'   spelling. (ggplot2 sizes the margins. The legend always lists every band, in the order of the p-values; it
#'   sits in the empty upper triangle of a single matrix of 6 or more rows, otherwise to the right of the plot.)
#' @param values Whether to print each p-value on its tile (the default is
#' \code{FALSE}). Useful for up to about 15 levels.
#' @param basesz The base font size. The default is 12.
#' @return A \code{ggplot} object (invisibly \code{NULL} when the matrix has
#' no pair of levels). Printing it draws the plot; it can be changed with the
#' usual ggplot2 functions. Before 2.0 the plot was drawn with base graphics
#' and nothing was returned: inside a function or a loop, wrap the call in
#' \code{print()}.
#' @seealso \code{\link{predictmeans}}, \code{\link{autoplot.predictmeans}}
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

#' @export
PMplot <- function(pmatrix, conf_level=0.95, mtitle=NULL, xylabel=NULL, values=FALSE, basesz=12, ...) {
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("PMplot", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  level <- pm_conf_alpha(conf_level)
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
  pm_ptile_plot(P, labs, level = level, title = mtitle, blocks = blocks, values = values, basesz = basesz,
                block_names = if (is.list(pmatrix)) names(pmatrix))
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
## first level at the top left. labs: one label per row; blocks: the sizes of
## the groups down the diagonal (comparisons within the levels of 'atvar'),
## drawn as one panel each, titled by 'block_names'. The legend sits in the
## empty upper triangle, anchored at its top-right corner, when there is one
## panel large enough; otherwise on the right, outside the panels.
pm_ptile_plot <- function(P, labs, level = 0.05, title = NULL, legend_title = NULL, blocks = NULL,
                          values = FALSE, basesz = 12, block_names = NULL) {
  n <- nrow(P)
  if (is.null(blocks)) blocks <- n
  grp <- rep(seq_along(blocks), blocks)
  k <- which(is.finite(P))
  d <- data.frame(row = (k - 1L) %% n + 1L, col = (k - 1L) %/% n + 1L, p = P[k])
  d$block <- grp[d$row]
  b <- pm_p_bands(level)
  d$band <- factor(b$labels[findInterval(d$p, b$breaks, left.open = TRUE) + 1L], levels = b$labels)
  # each panel shows its levels in their order; the first row and the last
  # column of a triangle have no pair and are left out, so the tiles fill a square
  ids <- paste(grp, labs, sep = "\r")                     # labels may repeat across blocks
  d$x <- factor(ids[d$col], levels = ids)
  d$y <- factor(ids[d$row], levels = rev(ids))
  nb <- length(unique(d$block))
  if (is.null(block_names)) block_names <- paste("Group", seq_along(blocks))
  d$panel <- factor(block_names[d$block], levels = unique(block_names))
  side <- max(tapply(d$row, d$block, function(r) length(unique(r))))
  long <- max(nchar(labs)) > 3
  lw <- if (side > 60) 0 else if (side > 30) 0.2 else 0.5   # no grid lines between very many tiles
  unid <- function(x) sub("^[^\r]*\r", "", x)
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$band)) +
    ggplot2::geom_tile(colour = if (lw > 0) "white" else NA, linewidth = lw, show.legend = TRUE) +
    ggplot2::scale_fill_manual(values = stats::setNames(b$colours, b$labels), breaks = b$labels, drop = FALSE,
                               name = if (is.null(legend_title)) b$title else legend_title) +
    ggplot2::scale_x_discrete(labels = unid, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(labels = unid, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL, title = pm_wrap_title(title, basesz)) +
    ggplot2::theme_bw(base_size = basesz) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   aspect.ratio = 1,
                   axis.text.x = if (long) ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5) else
                     ggplot2::element_text(),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.title.position = "plot",
                   legend.background = ggplot2::element_rect(fill = "white", colour = "grey40", linewidth = 0.3))
  if (nb > 1) {
    g <- g + ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free")
  }
  # inside the empty triangle only with one panel of 6 or more rows (the
  # legend takes about 4-5 rows); anchored at the corner, so it cannot run off
  g <- g + if (nb == 1 && side >= 6) {
    pm_legend_inside(c(0.99, 0.99)) + ggplot2::theme(legend.justification = c(1, 1))
  } else ggplot2::theme(legend.position = "right")
  if (isTRUE(values)) {
    d$txt <- ifelse(d$p < 0.001, "<.001", sub("^0\\.", ".", formatC(d$p, format = "f", digits = 3)))
    d$dark <- as.integer(d$band) <= length(b$labels)/2          # white text on the dark bands
    g <- g + ggplot2::geom_text(data = d, ggplot2::aes(label = .data$txt, colour = .data$dark),
                                size = basesz/3.5*min(1, 9/(side * min(nb, 3))), show.legend = FALSE) +
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
