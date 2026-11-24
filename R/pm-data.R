## The observations of a fitted model, with its response on the scales the
## means are shown on. One function for every model class predictmeans
## handles (lm, aov, aovlist, glm, gls, lme, lmer, glmer, glmmTMB, gam,
## semireg); plots read it for the data points, and pm_response_info() /
## pm_scale_label() say what the response is for axis labels and titles.

#' Data of a Fitted Model, with the Response on Both Scales
#'
#' @description
#' \code{pm_data()} returns the rows used to fit a model (after \code{subset} and missing values), with the response
#' on the scale of the model and on the original scale: \code{.y} is the response as modelled (e.g.
#' \code{log(y + 1)}; for a GLM the observed response, as a proportion for \code{cbind(successes, failures)}), and
#' \code{.y_orig} the response on the original scale (\code{y}), taken from the data or, when only the transformed
#' response is available, computed with the inverse transformation. A transformation given to
#' \code{\link{predictmeans}} as \code{trans} (for a response transformed before the fit) is used the same way.
#' It is also the generic through which a new model class can say where its data are (see \code{\link{pm_coef}}).
#' @param model A fitted model, or a \code{\link{semireg}} fit.
#' @param data Optional: the data the model was fitted to, when they cannot be found from the model.
#' @param trans Optional: a function that back-transforms the response (as \code{trans} of
#'   \code{\link{predictmeans}}), used when the model formula has an untransformed response.
#' @param \dots Not used.
#' @return A data frame with the variables of the data, \code{.y}, \code{.y_orig} and \code{.weights} (prior weights,
#'   when there are any), and attribute \code{"response"}: a list with \code{expr} (the response as written, e.g.
#'   \code{"log(y + 1)"}), \code{name} (the response variable, e.g. \code{"y"}), \code{scale} of the model's means
#'   (\code{"identity"}, \code{"transformed"} or \code{"link"}), \code{link}, \code{family} and \code{back} (what
#'   back-transforms the means: \code{"link"}, \code{"formula"}, \code{"trans"} or \code{NA}).
#' @examples
#'   oats <- as.data.frame(nlme::Oats)
#'   fm <- lm(log(yield + 1) ~ Variety, data = oats)
#'   d <- pm_data(fm)
#'   head(d[c("yield", ".y", ".y_orig")])
#'   attr(d, "response")
#' @export
pm_data <- function(model, ...) UseMethod("pm_data")

#' @rdname pm_data
#' @exportS3Method pm_data default
pm_data.default <- function(model, data = NULL, trans = NULL, ...) {
  d <- pm_fit_rows(model, data)
  info <- pm_response_info(model, trans = trans)
  y <- pm_eval_response(model, d)
  if (is.matrix(y) && ncol(y) == 2) y <- y[, 1]/rowSums(y)       # cbind(successes, failures)
  d$.y <- if (is.null(y)) NA_real_ else as.vector(y)
  d$.y_orig <- d$.y
  if (identical(info$back, "formula") || identical(info$back, "trans")) {
    # the original variable when the formula transforms it; with a user's
    # 'trans' the response variable is the transformed one
    orig <- if (identical(info$back, "formula") && info$name %in% names(d)) d[[info$name]]
    d$.y_orig <- if (is.numeric(orig) && length(orig) == nrow(d)) orig else as.numeric(info$inverse(d$.y))
  }
  w <- tryCatch(stats::weights(model), error = function(e) NULL)
  if (is.numeric(w) && length(w) == nrow(d) && any(w != 1)) d$.weights <- w
  info$inverse <- NULL
  attr(d, "response") <- info
  d
}

## The rows used in the fit: the 'data' given, the call's data in the fitting
## environment (subset and na.action applied; pm_get_data()), the data a
## semireg fit keeps, or the model frame (matrix columns flattened).
pm_fit_rows <- function(model, data = NULL) {
  if (inherits(model, "semireg")) return(as.data.frame(model$data))
  base <- if (inherits(model, "aovlist")) model else model
  if (is.null(data)) {
    data <- if (inherits(model, c("glmerMod", "glmmTMB"))) {
      tryCatch(pm_model_rows(model), error = function(e) NULL)
    } else tryCatch(suppressWarnings(pm_get_data(base)), error = function(e) NULL)
  }
  fr <- tryCatch(stats::model.frame(base), error = function(e) NULL)
  if (!is.data.frame(data)) {
    if (is.null(fr)) stop("Can't find the data used to fit the model; please give them via 'data'.", call. = FALSE)
    fr <- fr[, !grepl("^\\(", names(fr)), drop = FALSE]
    return(as.data.frame(lapply(fr, function(v) if (is.matrix(v) && ncol(v) == 1) as.vector(v) else v),
                         check.names = FALSE))
  }
  data <- as.data.frame(data)
  if (!is.null(fr) && nrow(data) != nrow(fr) && all(rownames(fr) %in% rownames(data))) {
    data <- data[rownames(fr), , drop = FALSE]
  }
  data
}

## The response as modelled, evaluated in the rows 'd' (NULL if it cannot be).
pm_eval_response <- function(model, d) {
  tt <- tryCatch(stats::terms(pm_fit(model)), error = function(e) NULL)
  e <- if (!is.null(tt) && isTRUE(attr(tt, "response") > 0)) attr(tt, "variables")[[attr(tt, "response") + 1]]
  if (is.null(e)) return(NULL)
  env <- attr(tt, ".Environment")
  y <- tryCatch(eval(e, d, if (is.environment(env)) env else globalenv()), error = function(err) NULL)
  if (is.null(y) || NROW(y) != nrow(d)) {
    fr <- tryCatch(stats::model.frame(model), error = function(err) NULL)
    y <- if (!is.null(fr) && nrow(fr) == nrow(d)) stats::model.response(fr)
  }
  y
}

## What the response is: list(expr, name, scale, link, family, back,
## inverse). scale (of the model's means): "link" (a non-identity link),
## "transformed" (a transformed response in the formula), else "identity".
## back (what back-transforms the means): "link", "formula", "trans" (the
## user's), or NA (none).
pm_response_info <- function(model, trans = NULL, link_scale = NULL) {
  m <- pm_fit(model)
  tt <- tryCatch(stats::terms(m), error = function(e) NULL)
  e <- if (!is.null(tt) && isTRUE(attr(tt, "response") > 0)) attr(tt, "variables")[[attr(tt, "response") + 1]]
  expr <- if (is.null(e)) "response" else paste(deparse(e), collapse = "")
  vars <- if (is.null(e)) character(0) else all.vars(e)
  if (is.call(e) && identical(e[[1]], as.name("cbind"))) vars <- vars[1]
  name <- if (length(vars) == 1) vars else expr
  fam <- tryCatch(stats::family(m), error = function(err) NULL)
  link <- if (is.list(fam)) fam$link else NULL
  family <- if (is.list(fam)) fam$family else NULL
  linked <- !is.null(pm_link_trans(m))
  rt <- pm_response_trans(m)
  user <- is.function(trans) && !identical(trans, I)
  if (is.null(link_scale)) link_scale <- linked && !user
  scale <- if (linked) "link" else if (!is.null(rt)) "transformed" else "identity"
  back <- if (link_scale) "link" else if (user) "trans" else if (!is.null(rt)) "formula" else NA_character_
  inverse <- if (user) trans else if (!is.null(rt)) rt$fun else NULL
  list(expr = expr, name = name, scale = scale, link = if (linked) link else NULL, family = family, back = back,
       inverse = inverse)
}

## Axis labels and titles for means of a predictmeans result (its settings):
## 'scale' "model" (the scale of estimate) or "response" (back-transformed).
pm_scale_label <- function(settings, scale = c("model", "response"), what = c("axis", "title")) {
  scale <- match.arg(scale)
  what <- match.arg(what)
  r <- settings$response
  term <- settings$modelterm
  if (is.null(r)) {
    return(if (what == "axis") (if (scale == "model") "Predicted mean" else "Back-transformed mean") else
      paste0("Predicted means for '", term, "'"))
  }
  prob <- isTRUE(r$family %in% c("binomial", "quasibinomial", "betabinomial", "beta", "beta_family")) ||
    isTRUE(r$link %in% c("logit", "probit", "cloglog", "cauchit"))
  if (scale == "model" || is.na(r$back)) {
    axis <- if (r$scale == "link") paste0(r$link, "(", if (prob) "probability" else "mean", " of ", r$name, ")") else
      r$expr
    if (what == "axis") return(axis)
    return(paste0("Means of ", axis, " for '", term, "'"))
  }
  bt <- settings$back_transform
  est <- if (is.null(bt)) "" else bt$estimate
  if (r$back == "trans") {
    axis <- paste0("trans(", if (r$scale == "link") paste0(r$link, " scale") else r$expr, ")")
    if (what == "axis") return(axis)
    return(paste0("Means of '", term, "' back-transformed by 'trans'"))
  }
  noun <- if (identical(settings$component, "overall")) "overall mean" else if (grepl("bias adjusted", est)) "mean" else
    if (grepl("median", est)) "median" else if (prob) "probability" else "mean"
  axis <- paste0(noun, " of ", r$name)
  if (grepl("typical cluster", est)) axis <- paste0(axis, ", typical cluster")
  if (what == "axis") return(pm_ucfirst(axis))
  how <- if (r$back == "link") paste0("inverse ", r$link, " link") else paste0("back-transformed from ", r$expr)
  if (grepl("bias adjusted", est)) how <- paste0(how, ", bias-adjusted")
  plural <- sub("^(overall mean|mean|median|probability)", "\\1s", axis)
  plural <- sub("^probabilitys", "probabilities", plural)
  paste0(pm_ucfirst(plural), " for '", term, "' (", how, ")")
}

pm_ucfirst <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

## settings$response of a result: pm_response_info() without the function,
## 'back' NA when the means are not back-transformed.
pm_response_label_info <- function(model, trans, link_scale, back_transformed, from_formula = FALSE) {
  m <- if (inherits(model, "pm_aovlist")) attr(model, "aovlist") else model
  user_trans <- if (link_scale || from_formula) NULL else trans
  r <- tryCatch(pm_response_info(m, trans = user_trans, link_scale = link_scale),
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  if (!back_transformed) r$back <- NA_character_
  r$inverse <- NULL
  r
}
