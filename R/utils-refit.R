## Refit a model with changed arguments from inside the package.
##
## stats::update() evaluates the new call in the caller's frame, which here is
## package code: data fitted in a function (or under a name that also exists
## elsewhere on the search path, such as 'Oats') is then not found, or worse,
## a different data set with the same name is used. pm_refit() puts the data
## frame itself into the call and evaluates the call where the model was
## fitted. With check = TRUE it also verifies that the refit has the same
## coefficients as the original model.
pm_refit <- function(model, data = NULL, ..., check = FALSE) {
  # Replace arguments in the original call (formulas are never changed here).
  cl <- stats::getCall(model)
  args <- list(...)
  # a formula given as a variable (lmer(f, ...)): put the formula itself in
  # the call, as the variable may not exist where the call is evaluated
  fa <- intersect(c("formula", "fixed", "model"), names(cl))[1]
  if (!is.na(fa) && is.name(cl[[fa]]) && !fa %in% names(args)) {
    fo <- tryCatch(stats::formula(model), error = function(e) NULL)
    if (inherits(fo, "formula")) cl[[fa]] <- fo
  }
  # the same for a family given as a variable (glm(y ~ x, family = fam))
  if (is.name(cl$family) && !"family" %in% names(args)) {
    fam <- tryCatch(stats::family(model), error = function(e) NULL)
    if (inherits(fam, "family")) cl$family <- fam
  }
  for (nm in names(args)) cl[[nm]] <- args[[nm]]
  if (!is.null(data)) cl$data <- data
  env <- pm_fit_env(model)
  fit <- eval(cl, env)
  if (check) {
    old <- names(pm_fixef(model))
    new <- names(pm_fixef(fit))
    if (!identical(old, new)) {
      stop("Refitting the model gave different coefficients from the fitted model; ",
           "please check that the data used in fitting are available.", call. = FALSE)
    }
  }
  fit
}

## Environment the model was fitted in. Most model classes record it in their
## terms (lm, glm, lme, merMod, glmmTMB, gam). gls models don't, so look up
## the call stack, innermost first and skipping this package's own functions,
## for the frame in which the 'data' argument of the call exists.
pm_fit_env <- function(model) {
  env <- tryCatch(attr(model$terms, ".Environment"), error = function(e) NULL)
  if (!is.environment(env) && !inherits(model, "gls")) {
    env <- tryCatch(attr(terms(model), ".Environment"), error = function(e) NULL)
  }
  if (is.environment(env)) return(env)
  data_expr <- tryCatch(stats::getCall(model)$data, error = function(e) NULL)
  if (!is.null(data_expr)) {
    own <- asNamespace("predictmeans")
    for (i in rev(seq_len(sys.nframe()))) {
      fn_env <- environment(sys.function(i))
      if (is.environment(fn_env) && identical(topenv(fn_env), own)) next
      f <- sys.frame(i)
      d <- tryCatch(eval(data_expr, f), error = function(e) NULL)
      if (is.data.frame(d)) return(f)
    }
  }
  globalenv()
}

## The data used to fit a gls or lme model. nlme::getData() looks for the data
## in the wrong place for models fitted inside a function; fall back to the
## fitting environment found by pm_fit_env().
pm_get_data <- function(model) {
  cl <- stats::getCall(model)
  # first where the model was fitted: nlme::getData() can find another data
  # set of the same name (e.g. nlme's own Oats) when the model was fitted in
  # a function
  d <- if (is.null(cl$data)) NULL else tryCatch(eval(cl$data, pm_fit_env(model)), error = function(e) NULL)
  if (!is.data.frame(d)) {
    d <- tryCatch(nlme::getData(model), error = function(e) NULL)
    if (is.data.frame(d)) return(d)
    if (is.null(cl$data)) {
      stop("Can't find the data used to fit the model; please refit it with a 'data' argument.", call. = FALSE)
    }
    stop("Can't find the data used to fit the model; please supply them via 'data'.", call. = FALSE)
  }
  if (!is.null(cl$subset)) {
    keep <- eval(cl$subset, d, pm_fit_env(model))
    d <- d[keep, , drop = FALSE]
  }
  na_act <- model$na.action
  if (!is.null(na_act) && inherits(na_act, "omit")) d <- d[-na_act, , drop = FALSE]
  d
}

## Fixed-effect estimates for any supported model class.
pm_fixef <- function(model) {
  if (inherits(model, "glmmTMB")) return(glmmTMB::fixef(model)$cond)
  if (inherits(model, "merMod")) return(lme4::fixef(model))
  if (inherits(model, "lme")) return(nlme::fixef(model))
  coef(model)
}

## The data of an lmer/glmer fit: the data frame of the call if it can be
## found, else the model frame stored in the fit (enough when the formula's
## variables are plain columns; transformed variables then need the data).
pm_mer_data <- function(model) {
  d <- tryCatch(nlme::getData(model), error = function(e) NULL)
  if (is.data.frame(d)) return(d)
  fr <- stats::model.frame(model)
  fr <- fr[, !grepl("^\\(", names(fr)), drop = FALSE]
  as.data.frame(lapply(fr, function(v) if (is.matrix(v)) as.vector(v) else v))
}
