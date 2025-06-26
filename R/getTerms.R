#' Get the expanded set of terms from a model
#'
#' Gets the expanded set of terms from a model specified using
#' Wilkinson-Rogers (1973) notation. For example, if your formula
#' is \code{y ~ A*B}, then this function will expand the formula to
#' \code{y ~ A + B + A:B} and return a vector containing \code{"A"},
#' \code{"B"}, and \code{"A:B"}.
#' @param model A fitted object from one of \code{aov}, \code{glm},
#' \code{glmm_TMB}, \code{gls}, \code{lm}, \code{lmer}, or \code{nlme}.
#'
#' @examples
#' ## A simple example without a model
#' f <- formula(y ~ A * B)
#' predictmeans:::getTerms(f)
#'
#' ## A more complex example with a fitted \code{lmerMod} object
#' Oats$nitro <- factor(Oats$nitro)
#' model <- lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = Oats)
#' predictmeans:::getTerms(model)
#'
#' @keywords internal
getTerms <- function(model) {
  attr(terms(model), "term.labels")
}
