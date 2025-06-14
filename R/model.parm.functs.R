#' model.frame method for objects of class 'gls'
#'
#' This is a specific \code{model.frame} method for objects of the class
#' \code{gls}. It simply calls \code{nlme::getData}. See \link[nlme]{getData}.
#'
#' @param formula An object of class \code{gls}.
#' @param ... Further arguments passed to \code{model.frame.formula}.
#' @importFrom stats coef model.frame model.matrix terms
#' @importFrom stats vcov
#' @method model.frame gls
#' @export
model.frame.gls <- function(formula, ...) {
  model.frame(formula(formula), nlme::getData(formula), ...)
}

#' model.frame method for objects of class 'lme'
#'
#' This is a specific \code{model.frame} method for objects of the class
#' \code{lme}.
#'
#' @param formula An object of class \code{lme}.
#' @param ... Further arguments passed to \code{model.frame.formula}.
#' @importFrom stats coef model.frame model.matrix terms
#' @importFrom stats vcov
#' @method model.frame lme
#' @export
model.frame.lme <- function(formula, ...) {
  formula$data
}

#' model.matrix method for objects of class 'aovlist'
#'
#' This is a specific \code{model.matrix} method for objects of the class
#' \code{aovlist}. It does not do anything except give an error because
#' \code{aovlist} is not supported.
#'
#' @param object An object of class \code{aovlist}.
#' @param ... Not used.
#' @method model.matrix aovlist
#' @export
model.matrix.aovlist <- function(object, ...) {
  stop(
    sQuote("predicted.means"),
    " does not support objects of class ",
    sQuote("aovlist")
  )
}

#' model.matrix method for objects of class 'gls'
#'
#' This is a specific \code{model.matrix} method for objects of the class
#' \code{gls}. See \link[nlme]{getData}.
#'
#' @param object An object of class \code{gls}.
#' @param ... Further arguments passed to \code{model.matrix.terms}.
#' @method model.matrix gls
#' @export
model.matrix.gls <- function(object, ...) {
  model.matrix(
    object = terms(object),
    data = nlme::getData(object),
    ...
  )
}

#' model.matrix method for objects of class 'lme'
#'
#' This is a specific \code{model.frame} method for objects of the class
#' \code{lme}. See \link[nlme]{lme}.
#'
#' @param object An object of class \code{lme}.
#' @param ... Further arguments passed to \code{model.matrix.terms}.
#' @method model.matrix lme
#' @export
model.matrix.lme <- function(object, ...) {
  ## DONGWEN: ... appears to already contain data, so I have changed the
  ## definition. The reason for the if statement is that varcomp seems to use
  ## one definition, and Kmatrix, another. TODO: We probably need a deeper dive
  ## here.
  dot_args <- list(...)

  # Only add data if not already passed in ...
  if (!"data" %in% names(dot_args)) {
    model.matrix(object = terms(object), data = model.frame(object), ...)
  } else {
    model.matrix(object = terms(object), ...)
  }
}

#' terms method for objects of class 'gls'
#'
#' This is a specific \code{terms} method for objects of the class \code{gls}.
#' See \link[nlme]{getData}.
#'
#' @param x An object of class \code{gls}.
#' @param ... Additional arguments passed to \code{terms}.
#' @method terms gls
#' @export
terms.gls <- function(x, ...) {
  terms(model.frame(x), ...)
}

#' terms method for objects of class 'lme'
#'
#' This is a specific \code{terms} method for objects of the class \code{lme}.
#' See \link[nlme]{lme}.
#'
#' @param x An object of class \code{lme}.
#' @param ... Not used.
#' @method terms lme
#' @export
terms.lme <- function(x, ...) {
  v <- x$terms

  if (is.null(v)) {
    stop("no terms component")
  }

  return(v)
}


# terms.merMod <- function (object, ...) {
#  v <- terms(object)
#  attr(v, "dataClasses") <- sapply(all.vars(formula(object, fixed.only=TRUE)),
# function(x) class(model.frame(object)[[x]])[1])
#  return(v)
# }
#
