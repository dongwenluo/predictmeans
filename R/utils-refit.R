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
## terms (lm, glm, lme, merMod, glmmTMB). gls models don't, and mgcv's gam
## records the global environment whatever it was fitted in; so when the
## recorded environment does not hold the 'data' of the call, look up the
## call stack, innermost first and skipping this package's own functions, for
## the frame in which it exists.
pm_fit_env <- function(model) {
  env <- tryCatch(attr(model$terms, ".Environment"), error = function(e) NULL)
  if (!is.environment(env) && !inherits(model, "gls")) {
    env <- tryCatch(attr(terms(model), ".Environment"), error = function(e) NULL)
  }
  data_expr <- tryCatch(stats::getCall(model)$data, error = function(e) NULL)
  has_data <- function(e) is.data.frame(tryCatch(eval(data_expr, e), error = function(err) NULL))
  # a frame that holds the data itself (a name, not found by inheritance
  # from the global environment)
  has_own <- function(e) {
    if (is.name(data_expr)) exists(as.character(data_expr), envir = e, inherits = FALSE) && has_data(e) else has_data(e)
  }
  # gam records the global environment wherever it was fitted: a data set of
  # the same name there (another session object) is not evidence; look in the
  # calling frames first
  unreliable <- inherits(model, "gam") && identical(env, globalenv())
  if (!unreliable && is.environment(env) && (is.null(data_expr) || has_data(env))) return(env)
  if (!is.null(data_expr)) {
    own <- asNamespace("predictmeans")
    for (i in rev(seq_len(sys.nframe()))) {
      fn_env <- environment(sys.function(i))
      if (is.environment(fn_env) && identical(topenv(fn_env), own)) next
      f <- sys.frame(i)
      if (has_own(f)) return(f)
    }
    if (unreliable && has_data(globalenv())) return(globalenv())
  }
  if (is.environment(env)) env else globalenv()
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
  na_act <- if (isS4(model)) attr(stats::model.frame(model), "na.action") else model$na.action
  if (!is.null(na_act) && inherits(na_act, c("omit", "exclude"))) d <- d[-na_act, , drop = FALSE]   # (na.exclude too)
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

## lme4::refit() with a new response for the rows the fit used. lme4 1.1-35
## (and other versions before the fix) drops the fit's NA rows again from a
## response that has none ("replacement has 58 rows, data has 59"), unless
## the response carries an "na.action" attribute: so mark it as subset. All
## refits of lme4 fits with a new response go through here.
pm_mer_refit <- function(model, y, ...) {
  na <- attr(model@frame, "na.action")
  if (!is.null(na) && is.null(dim(y)) && is.null(attr(y, "na.action"))) attr(y, "na.action") <- na
  lme4::refit(model, newresp = y, ...)
}
