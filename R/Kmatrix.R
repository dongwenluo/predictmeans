#' Matrix of Coefficients in a Linear Model
#'
#' This function obtains a matrix of coefficients for parametric models such as
#' \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#'
#'
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which model term's
#' predicted mean to be calculated.  The \code{modelterm} must be given exactly
#' as it appears in the printed model, e.g. "A" or "A:B".
#' @param covariate A numerical vector to specify values of covariates for
#' calculating predicted means, default values are the means of the associated
#' covariates. It also can be the name of one covariate in the model.
#' @param covariateV A numeric vector or list of numeric vector, then
#' covariatemeans will produce the result for \code{covariate} at value of
#' \code{covariateV}.
#' @param data In some cases, you need to provide the data set used in model
#' fitting, especially when you have applied some variable trnasformation in
#' the model.
#' @param prtnum An option for printing covariate info on the screen or not.
#' The default is FALSE.
#' @return \item{K}{Coefficients matrix} \item{fctnames}{A model frame contains
#' factor(s) info in the model.} \item{response}{The name of response variable
#' in the model.}
#' @author This function heavily depends on the codes from package "lsmeans".
#' @references Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R.
#' (2004), \emph{Prediction in linear mixed models}, Australian and New Zealand
#' Journal of Statistics, 46(3), 325-347.
#' @examples
#'
#'   library(predictmeans)
#'   data(Oats, package="nlme")
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'   Kmatrix(fm, "Variety", prtnum=TRUE)$K
#'   Kmatrix(fm, "Variety", 0.5, prtnum=TRUE)$K
#'  # Kmatrix(fm, "Variety", "nitro")$K
#'   Kmatrix(fm, "Variety", "nitro", covariateV=seq(0, 0.6, 0.1))$K
#'
#' @importFrom methods slot
#' @importFrom stats as.formula delete.response formula get_all_vars
#' @importFrom stats na.omit na.pass
#' @export

Kmatrix <- function(model,
                    modelterm,
                    covariate = NULL,
                    covariateV = NULL,
                    data = NULL,
                    prtnum = FALSE) {
  if (inherits(model, "mer") || inherits(model, "merMod")) {
    if (!lme4::isLMM(model) && !lme4::isGLMM(model)) {
      stop("Can't handle a nonlinear mixed model")
    }
    thecall <- slot(model, "call")
    contrasts <- attr(model.matrix(model), "contrasts")
  } else if (inherits(model, "lme")) {
    thecall <- model$call
    contrasts <- model$contrasts
  } else if (inherits(model, "gls")) {
    thecall <- model$call
    contrasts <- model$contrasts
  } else if (inherits(model, "lm")) {
    thecall <- model$call
    contrasts <- attr(model.matrix(model), "contrasts")
  } else if (inherits(model, "glmmTMB")) {
    thecall <- model$call
    contrasts <- attr(model.matrix(model), "contrasts")
  } else {
    stop(paste("Can't handle an model of class", class(model)[1]))
  }

  cov.reduce <- function(x, name) {
    mean(x, na.rm = TRUE)
  }

  fac.reduce <- function(coefs, lev) {
    apply(coefs, 2, mean, na.rm = TRUE)
  }

  # Get model formula and index of response
  Terms <- terms(model)
  if (!is.null(attr(Terms, "offset"))) {
    offset_term <- deparse(attr(Terms, "variables")[attr(Terms, "offset") +
                                                      1])
    offset_n <- unlist(strsplit(trimws(gsub(
      "[\\(\\)]", " ", offset_term
    )), " "))
    offset_n <- offset_n[length(offset_n)]
  }
  yname <- as.character(attr(Terms, "variables"))[[2]]
  Terms <- delete.response(Terms)

  # get the pure formula w/o extra stuff
  formrhs <- formula(Terms)
  nm <- all.vars(formrhs)
  nm <- nm[nm != "pi"]

  # Figure out if any are coerced to factor or ordered
  anm <- all.names(formrhs)
  coerced <- anm[1 + grep("factor|as.factor|ordered|as.ordered", anm)]

  # Obtain a simplified formula -- needed to recover the data in the model
  form <- as.formula(paste("~", paste(nm, collapse = "+")))
  envir <- attr(Terms, ".Environment")
  #eval(thecall$data, envir=envir)

  if (inherits_any(model, c("mer", "merMod", "lme", "gls"))) {
    ## NOTE DONGWEN: have made getData specific to lme4
    model_data <- get_all_vars(formula(terms(model)), lme4::getData(model))
    model_data <- droplevels(na.omit(model_data))
    model_data <- get_all_vars(formrhs, model_data)
  } else {
    if (!is.null(data)) {
      if (grepl("[,]", yname)) {
        yname <- unlist(strsplit(yname, "[,] "))[2]
        yname <- gsub("\\)", "", unlist(strsplit(yname, " - "))) # ynames
      }
      model_data <- as.data.frame(data[, c(yname, nm)])
    } else {
      model_data <- model.frame(model)
    }
  }

  if (is.null(covariate) || !setequal(covariate, modelterm)) {
    vars <- unlist(strsplit(modelterm, "\\:"))
    if (any(!base::is.element(vars, names(model_data)[sapply(model_data, is.factor)]))) {
      stop(paste(vars, "must be factor(s) or not in the fixed effects!"))
    }
  }

  X <- model.frame(
    form,
    model_data,
    subset = eval(thecall$subset, enclos = envir),
    na.action = na.omit,
    drop.unused.levels = TRUE
  )
  preddf <- X
  baselevs <- xlev <- matdat <- list()
  all.var.names <- names(X)

  for (xname in all.var.names) {
    obj <- X[[xname]]
    if (is.factor(obj)) {
      xlev[[xname]] <- levels(obj)
      baselevs[[xname]] <- levels(obj)
    } else if (is.matrix(obj)) {
      # Matrices -- reduce columns thereof, but don't add to baselevs
      matdat[[xname]] <- apply(obj, 2, cov.reduce, xname)
    } else {
      # single numeric pred but coerced to a factor - use unique values
      if (length(grep(xname, coerced)) > 0) {
        baselevs[[xname]] <- sort(unique(obj))
        # Ordinary covariates - summarize if not in 'at' arg
      } else {
        baselevs[[xname]] <- cov.reduce(obj, xname)
      }
    }
  }

  factor_names <- c(names(xlev), coerced)
  covlevname <- setdiff(names(baselevs), factor_names)

  if (!is.null(attr(Terms, "offset"))) {
    covlevname <- covlevname[covlevname != offset_n]
  }
  if (length(factor_names) != 0) {
    n.table <- table(preddf[, factor_names, drop = FALSE])
    ndf <- data.frame(n.table)
  } else {
    ndf <- NULL
  }

  if ((!is.null(covariate) && !unique(covariate %in% c("NULL", ""))) &&
        all(is.numeric(covariate))) {
    baselevs[covlevname] <- as.list(covariate)
  }
  if ((!is.null(covariate) && !unique(covariate %in% c("NULL", ""))) &&
        is.character(covariate) && all(covariate %in% covlevname))  {
    if (length(covariate) == 1) {
      # if (as.is) {
      # baselevs[[covariate]] <- sort(unique(X[[covariate]]))
      # }else{
      if ((!is.null(covariateV) && !unique(covariateV %in% c("NULL", ""))) && is.vector(covariateV)) {
        baselevs[[covariate]] <- covariateV
      } else {
        baselevs[[covariate]] <- sort(unique(c(seq(
          min(X[[covariate]]), max(X[[covariate]]), length = 50
        ), X[[covariate]])))
      }
      # }
    } else {
      if ((!is.null(covariateV) && !unique(covariateV %in% c("NULL", ""))) && is.list(covariateV)) {
        for (i in seq_along(covariate)) {
          baselevs[[covariate[i]]] <- sort(unique(covariateV[[i]]))
        }
      } else {
        for (i in covariate) {
          baselevs[[i]] <- sort(unique(c(seq(min(X[[i]]), max(X[[i]]), length = 30))))
        }
      }
    }
  }

  if (all(length(covlevname) != 0, prtnum)) {
    cat("\n", "The predicted means are estimated at \n\n")
    print(round(unlist(baselevs[covlevname]), 4))
    cat("\n")
  }

  # OK. Now make a grid of the factor levels of interest, along w/ covariate "at" values
  grid <- do.call("expand.grid", baselevs)

  # add any matrices
  for (nm in names(matdat)){
    grid[[nm]] <- matrix(rep(matdat[[nm]], each = nrow(grid)), nrow = nrow(grid))
  }

  # Now make a new dataset with just the factor combs and covariate values we want for prediction
  # WARNING -- This will overwrite X, so get anything you need from X BEFORE we get here
  if (!is.null(ndf)) {
    ndf <- merge(grid, ndf, sort = FALSE)
  }
  m <- model.frame(Terms, grid, na.action = na.pass, xlev = xlev)
  X <- model.matrix(Terms, m, contrasts.arg = contrasts)
  max_ord <- max(attr(Terms, "order"))
  if (!setequal(covariate, modelterm) && max_ord > 1 && length(factor_names) > 1  && any(ndf$Freq == 0)) {
    warning("Missing treatments' combination appeared, predicted means maybe misleading!")
    X[ndf$Freq == 0, ] <- NA
    # X <- X[ndf$Freq!=0, ]
  }

  # All factors (excluding covariates)
  allFacs <- all.var.names

  ### Array of indexes for rows of X, organized by dimensions
  row.indexes <- array(seq_len(nrow(X)), sapply(baselevs, length))

  # convert a string to a formula
  form <- as.formula(paste("~", paste(modelterm, collapse = "+")))

  # These are the variables involved; and the label to use in the results
  facs <- all.vars(form)
  if ((!is.null(covariate) && !unique(covariate %in% c("NULL", ""))) &&
        all(is.character(covariate), all(!covariate %in% facs))) {
    facs <- c(facs, covariate)
  }
  if (any(sapply(facs,
                 function(nm) {
                   length(grep(nm, allFacs)) == 0
                 }))) {
    stop(paste(
      "Unknown factor(s) in specification:",
      paste(form, collapse = " ")
    ))
  }

  # create the grid of factor combinations
  levs <- list()
  ## TODO: Can we remove this loop? It looks like it
  for (f in facs) {
    levs[[f]] <- baselevs[[f]]
  }
  combs <- do.call("expand.grid", levs)

  fctnames <- do.call("expand.grid", levs[rev(names(levs))])
  fctnames <- fctnames[, rev(names(fctnames)), drop = FALSE]
  rnK <- do.call("paste", c(fctnames, sep = ":"))

  K <- plyr::alply(row.indexes, match(facs, names(baselevs)), function(idx) {
    fac.reduce(X[idx, , drop = FALSE], "")
  })
  K <- as.matrix(as.data.frame(K))
  dimnames(K)[[2]] <- do.call("paste", c(combs, sep = ":"))
  K <- t(K)
  K <- K[rnK, , drop = FALSE]
  if (dim(model.matrix(model))[2] != dim(K)[2]) {
    if (length(setdiff(colnames(model.matrix(model)), colnames(K))) != 0) {
      stop(paste("You may need drop empty levels in factor", modelterm))
    }
    K <- K[, colnames(model.matrix(model)), drop = FALSE]
  }

  # fctnames <- cbind(fctnames, K[, setdiff(covlevname, names(fctnames)), drop=FALSE])
  extra_covn <- setdiff(covlevname, names(fctnames))
  if (length(extra_covn) > 0) {
    ## TODO: Can we remove this loop? It looks like it
    for (i in extra_covn){
      fctnames[[i]] <- baselevs[[i]]
    }
  }

  return(list(
    K = K,
    fctnames = fctnames,
    response = yname,
    preddf = preddf
  ))
}
