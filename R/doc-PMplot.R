# Help page PMplot (roxygen2 source; man/PMplot.Rd is generated).
#' Level Plot of a Matrix of p-values.
#'
#' @description
#' Creates a plot of p-values of pairwise comparisons.
#' @usage
#' PMplot(pmatrix, level=0.05, mtitle=NULL, xylabel=NULL, margin=5, legendx=0.73, 
#'   newwd=FALSE)
#' @param pmatrix A matrix with p-values from pairwise comparisons. (This is a lower triangle matrix.)
#' @param level The level of p-value to be highlighted. Default is 0.05.
#' @param mtitle The main title in the graph.
#' @param xylabel The x and y labels in the graph.
#' @param margin A value for specifying x and y margins in the graph. The default value is 5.
#' @param legendx A value for specifying x coordinate of legend. The default value is 0.73.
#' @param newwd A logical variable to indicate whether to print graph in a new window. The default is FALSE.
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
#'
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", barplot=TRUE)
#'   PMplot(predictout$p_valueMatrix)
#' @name PMplot
NULL
