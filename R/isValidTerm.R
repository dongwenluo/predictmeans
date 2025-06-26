#' Checks to see if the user has selected a valid model term
#'
#' A more user friendly check to see if the user has selected a valid model to
#' compute means for. As well as doing a simple check on presence in the model
#' the function also checks to see if the term is actually a legal model term in
#' terms of syntactical correctness, case, ordering of factors within
#' interactions, and both case and ordering. If any of these are true, then the
#' function will return \code{FALSE} but the user is provided with a *potentially*
#' valid option to call the function again. We say *potentially* because ultimately
#' this is a programmatic check and it cannot cover all possible cases.
#'
#' @param model A fitted object from one of \code{aov}, \code{glm},
#'   \code{glmm_TMB}, \code{gls}, \code{lm}, \code{lmer}, or \code{nlme}.
#' @param modelterm a character string containing one of the terms in the fitted
#'   model. Note that this means the possibly expanded model in Wilkinson Rogers
#'   (1973) notation. For example, \code{Y ~ A * B} expands to \code{Y ~ A + B +
#'   A:B}, so legitimate choices for \code{modelterm} in this model would be
#'   \code{"A"}, \code{"B"}, and \code{"A:B"}.
#'
#' @examples
#' f <- formula(y ~ A * B)
#' predictmeans:::isValidTerm(f, "A")
#' \dontrun{
#' ## This causes a warning so not is not run during R CHECK
#' predictmeans:::isValidTerm(f, "a")
#' }
#' predictmeans:::isValidTerm(f, "A:B")
#' \dontrun{
#' ## These cause warnings so not is not run during R CHECK
#' predictmeans:::isValidTerm(f, "B:A")
#' predictmeans:::isValidTerm(f, "a:b")
#' predictmeans:::isValidTerm(f, "b:a")
#' }
#'
#' ## A more complex example with a fitted \code{lmerMod} object
#' Oats$nitro <- factor(Oats$nitro)
#' model <- lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = Oats)
#' predictmeans:::isValidTerm(model, "nitro:Variety")
#' \dontrun{
#' ## This causes a warning so not is not run during R CHECK
#' ## Random effects are not valid terms for predictmeans
#' predictmeans:::isValidTerm(model, "Block")
#' }
#'
#' @keywords internal#'
#'
#' @importFrom combinat permn
isValidTerm <- function(model, modelterm) {
  ## First check that modelterm is scalar.
  if (length(modelterm) > 1) {
    warning("modelterm must be a scalar. Only using first element.")
    modelterm <- modelterm[1]
  }

  ## Now test to see if the user is asking for a valid
  ## model term
  isValidModelTermName <- function(modelterm) {
    tryCatch(
      {
        terms(as.formula(paste("y ~", modelterm)))
        TRUE
      },
      error = function(e) {
        FALSE
      }
    )
  }

  if (!isValidModelTermName(modelterm)) {
    warning(
      paste0(
        modelterm,
        " is not a syntactically correct model term.",
        " Perhaps you made a typo?"
      )
    )
    return(FALSE)
  }

  ModelTerms <- getTerms(model)

  ## First simple check
  if (modelterm %in% ModelTerms) {
    return(TRUE)
  } else {
    cat(paste0(
      "I can't find ",
      modelterm,
      ". I am going to try some other checks.\n"
    ))
  }

  ## Check case sensitivity
  if (tolower(modelterm) %in% tolower(ModelTerms)) {
    correctCase <- ModelTerms[match(tolower(modelterm), tolower(ModelTerms))]
    msg <- paste0(
      "I found your model term, but the case is not the same.",
      " Perhaps you meant ",
      correctCase,
      "?\n",
      "If you did call the function again with \"",
      correctCase,
      "\""
    )
    warning(msg)
    return(FALSE)
  }

  ## Is it an interaction?
  if (grepl(":", modelterm)) {
    ## check against the permuations
    Factors <- unlist(strsplit(modelterm, ":"))
    Perms <- permn(Factors)
    Perms <- lapply(Perms, paste0, collapse = ":")
    Perms <- as.vector(do.call(rbind, Perms))
    Interactions <- ModelTerms[grepl(":", ModelTerms)]

    if (any(Perms %in% Interactions)) {
      correctOrder <- Perms[match(Interactions, Perms)]
      msg <- paste0(
        "I found your model term, but the order is not the same as",
        " the model. Perhaps you meant ",
        correctOrder,
        "?\n",
        "If you did, then call the function again with \"",
        correctOrder,
        "\"."
      )
      warning(msg)
      return(FALSE)
      ## check case as well
    } else if (any(tolower(Perms) %in% tolower(Interactions))) {
      correctOrder <- Perms[match(tolower(Interactions), tolower(Perms))]
      correctCase <- Interactions[match(tolower(correctOrder), tolower(Interactions))]

      msg <- paste0(
        "I found your model term, but the case is not the same and",
        " the order is incorrect. Perhaps you meant ",
        correctCase,
        "?\n",
        "If you did, then call the function again with \"",
        correctCase,
        "\"."
      )
      warning(msg)
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}
