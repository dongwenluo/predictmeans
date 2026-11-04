# Help page Body_length (roxygen2 source; man/Body_length.Rd is generated).
#' \code{Body_length} data
#'
#' @description
#' Simulated body lengths of beetles reared in containers within populations, under two
#' treatments in two habitats, of the kind used by Nakagawa and Schielzeth (2013) to illustrate
#' \eqn{R^2} for mixed models. Used in the \code{R2_glmm} vignette.
#' @usage
#' data(Body_length)
#' @format
#' A data frame with 960 rows and 6 columns:
#' \describe{
#' \item{Population}{factor, 12 populations.}
#' \item{Container}{factor, 120 containers (nested in populations).}
#' \item{Sex}{factor, \code{Female} or \code{Male}.}
#' \item{Habitat}{factor, \code{A} or \code{B}.}
#' \item{Treatment}{factor, \code{Cont} or \code{Exp}.}
#' \item{BodyL}{numeric, body length.}
#' }
#' @references
#' Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for obtaining
#' \eqn{R^2} from generalized linear mixed-effects models. \emph{Methods in Ecology and
#' Evolution}, 4, 133--142.
#' @examples
#' data(Body_length)
#' fm <- lmer(BodyL ~ Treatment + Habitat + Sex + (1|Population/Container), data = Body_length)
#' R2_glmm(fm)
#' @name Body_length
#' @docType data
NULL
