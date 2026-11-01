# Help page autoplot.predictmeans (roxygen2 source; man/autoplot.predictmeans.Rd is generated).
#' Plots and Summary of Predicted Means
#'
#' @description
#' \code{autoplot()} draws the predicted means, their pairwise p-values or the back-transformed means of a
#' \code{\link{predictmeans}} result as a \pkg{ggplot2} object. \code{summary()} prints the means with a footer that
#' states the degrees of freedom, the SED and LSD ranges, the p-value adjustment and the covariate values used.
#' Both are built from the \code{means} and \code{pairs} tables of the result: they do not compute statistics, open
#' graphics devices or print on their own, so they can be used in reports and Shiny apps
#' (e.g. \code{renderPlot(autoplot(pm))}) or made interactive with \code{plotly::ggplotly()}.
#' @usage
#' \method{autoplot}{predictmeans}(object, type = c("means", "pvalue", "response"),
#'          interval = c("ci", "lsd", "none"), ...)
#'
#' \method{summary}{predictmeans}(object, digits = 4, ...)
#'
#' \method{print}{summary.predictmeans}(x, ...)
#' @param object A result of \code{\link{predictmeans}}.
#' @param type \code{"means"} (the default) plots the predicted means; \code{"pvalue"} a heatmap of the adjusted
#'   pairwise p-values (within each level of \code{atvar} when it was given); \code{"response"} the back-transformed
#'   means.
#' @param interval Bars around the means: confidence limits (\code{"ci"}, the default), half the average LSD
#'   (\code{"lsd"}; means whose bars do not overlap differ significantly when the LSDs are equal), or none.
#' @param digits Number of significant digits printed (plus two).
#' @param x A result of \code{summary()}.
#' @param \dots Not used.
#' @return
#' \code{autoplot()} returns a \code{ggplot} object; \code{summary()} returns an object whose print method shows the
#' table and footer.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
#'   summary(pm)
#'   ggplot2::autoplot(pm)
#'   ggplot2::autoplot(pm, interval = "lsd")
#'   ggplot2::autoplot(pm, type = "pvalue")
#' @name autoplot.predictmeans
#' @aliases summary.predictmeans print.summary.predictmeans
NULL
