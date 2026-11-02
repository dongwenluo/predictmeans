# Help page CookD (roxygen2 source; man/CookD.Rd is generated).
#' Calculates and plots Cook's distances for a Linear (Mixed) Model
#'
#' @description
#' This function produces Cook's distance plots for a linear model obtained from functions \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, 
#'   or \code{lmer}.
#' @usage
#' CookD(model, group=NULL, plot=TRUE, idn=3, newwd=FALSE)
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @param group Name (in "quotes") for indicating how observations are deleted for Cook's distance calculation. If \code{group!=NULL} then deletions will be along levels of \code{group} variable, otherwise, will be along individual observations.
#' @param plot A logical variable; if it is true, a plot of Cook's distance will be presented. The default is TRUE.
#' @param idn An integer indicating the number of top Cook's distances to be labelled in the plot. The default value is 3.
#' @param newwd A logical variable to indicate whether to print graph in a new window. The default value is FALSE.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # library(lme4)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' CookD(fm)
#' @name CookD
NULL
