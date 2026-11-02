# Help page se_ranef (roxygen2 source; man/se_ranef.Rd is generated).
#' Extract Standard Errors of Model Random Effects
#'
#' @description
#' These functions extract standard errors of model random effects
#'    from objects returned by modeling functions.
#' @details
#' \code{se_ranef} extracts standard errors of the random effects
#'     from objects returned by lmer, glmer and glmmTMB functions.
#' @usage
#' se_ranef(object, rand_term=NULL)
#' @param object object of \code{merMod} and \code{glmmTMB} fit
#' @param rand_term a name of random term in the model.
#' @return
#' \code{se_ranef} gives a list of standard errors for \code{ranef}.
#' @references
#' This function is modified from function 'se.ranef' at package 'arm'.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @name se_ranef
NULL
