# Help page tidy.predictmeans (roxygen2 source; man/tidy.predictmeans.Rd is generated).
#' Predicted Means and Pairwise Comparisons as Data Frames
#'
#' @description
#' Return the \code{means} or \code{pairs} table of a \code{\link{predictmeans}} result. The tables have fixed column
#' names (see the Value section of \code{\link{predictmeans}}) and unrounded values, which makes them easy to use in
#' further analysis, reports, Shiny apps or AI tools. \code{tidy()} is the generic from package \pkg{generics} (also
#' used by \pkg{broom}); \code{as.data.frame()} needs no other package.
#' @usage
#' \method{tidy}{predictmeans}(x, what = c("means", "pairs"), ...)
#'
#' \method{as.data.frame}{predictmeans}(x, row.names = NULL, optional = FALSE,
#'               what = c("means", "pairs"), ...)
#' @param x A result of \code{\link{predictmeans}}.
#' @param what Which table: \code{"means"} (the default) or \code{"pairs"}.
#' @param row.names Optional row names for the result.
#' @param optional Not used.
#' @param \dots Not used.
#' @return
#' A data frame.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
#'   as.data.frame(pm)
#'   head(as.data.frame(pm, what = "pairs"))
#' @name tidy.predictmeans
#' @aliases as.data.frame.predictmeans
NULL
