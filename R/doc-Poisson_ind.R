# Help page Poisson_ind (roxygen2 source; man/Poisson_ind.Rd is generated).
#' \code{Poisson_ind} data
#'
#' @description
#' Simulated egg counts of individual female beetles in containers within populations, under
#' two treatments in two habitats, of the kind used by Nakagawa and Schielzeth (2013). Used
#' in the \code{R2_glmm} vignette for a Poisson GLMM with an individual-level (overdispersion)
#' random effect.
#' @usage
#' data(Poisson_ind)
#' @format
#' A data frame with 480 rows and 7 columns:
#' \describe{
#' \item{Population}{factor, 12 populations.}
#' \item{Container}{factor, 60 containers (nested in populations).}
#' \item{Sex}{factor, \code{Female} only.}
#' \item{Habitat}{factor, \code{B} or \code{A}.}
#' \item{Treatment}{factor, \code{Exp} or \code{Cont}.}
#' \item{Egg}{integer, number of eggs.}
#' \item{individual}{factor, individual within container.}
#' }
#' @references
#' Nakagawa, S. and Schielzeth, H. (2013), A general and simple method for obtaining
#' \eqn{R^2} from generalized linear mixed-effects models. \emph{Methods in Ecology and
#' Evolution}, 4, 133--142.
#' @name Poisson_ind
#' @docType data
NULL
