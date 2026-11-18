#' Functions removed in predictmeans 2.0
#'
#' @description
#' These functions were removed in 2.0. Each stops with what to use instead; \code{vignette("migrating")} shows
#' the 1.1.1 code beside the 2.0 code.
#' \tabular{ll}{
#'   \strong{Removed} \tab \strong{Use} \cr
#'   \code{smZ()}, \code{smoothZ =} \tab smooths in the formula of \code{\link{semireg}} or \code{\link{semireg_tmb}},
#'     e.g. \code{smoothZ = list(sm = smZ(x))} is \code{s(x, bs = "os")} \cr
#'   \code{permglmer()}, \code{permglmmTMB()} \tab \code{\link{permlmer}()} (they were the same function) \cr
#'   \code{ci_mcp()} \tab \code{\link{predictmeans}(..., letterCI = TRUE)}: letters from comparison intervals,
#'     whose overlap agrees with the pairwise tests (overlap of confidence intervals does not) \cr
#'   \code{R2_glmm()} \tab \code{\link{R2_model}(model)}: the same values, as numbers \cr
#'   \code{permanova.lmer()} \tab \code{\link{permmodels}(model, type = )}: Freedman-Lane permutations, so the
#'     p-values differ from the refits of \code{permanova.lmer()} \cr
#' }
#' @param ... Not used.
#' @name predictmeans-defunct
#' @aliases smZ permglmer permglmmTMB permanova.lmer ci_mcp R2_glmm
#' @keywords internal
NULL

## Arguments and functions removed in 2.0 (dev/review-options.md, the
## maintainer's decisions D1-D10 of 2026-11-18). One table; every front end
## calls pm_defunct() on its call, and each removed function calls
## pm_defunct_fun(). A removed spelling always stops with the 2.0 spelling:
## never ignored through '...', never a warning (rule no. 1: no silent wrong
## answer). vignette("migrating") lists the same rows.
##
## Columns: fun (the function), arg (the 1.1.1 argument), use (the 2.0
## spelling; "%s" is filled with the value given, through 'value': "same" or
## "one_minus" for level -> conf.level).
pm_defunct_rows <- function() {
  r <- function(fun, arg, use, value = "") data.frame(fun = fun, arg = arg, use = use, value = value)
  plot_set <- function(fun, args) r(fun, args, paste0("plot_args = list(", args, " = ...) or autoplot()"))
  means <- c("predictmeans", "contrastmeans", "covariatemeans")
  print_dig <- "print(x, digits = ...): nothing is printed during the call"
  rbind(
    # one vocabulary for the means functions (D1, D2, D4)
    r(means, "Df", "df_method = %s", "same"),
    r(c(means, "pm_slopes", "semipred", "permstats", "PMplot"), "level", "conf.level = %s", "one_minus"),
    r("varcomp", "level", "conf.level = %s (it was a confidence level already)", "same"),
    r(c("covariatemeans", "semipred", "Kmatrix"), "covariateV", "at = %s", "same"),
    r("df_term", "type", "df_method = %s", "same"),
    r("df_term", "vcov_adj", "nothing: the adjusted covariance is computed from df_method"),
    # back-transformation: detected from the formula or the family
    r(c("predictmeans", "covariatemeans"), "transOff", "trans = function(x) exp(x) - offset (the offset inside trans)"),
    r("predictmeans", "count", "nothing: Poisson and negative binomial means are detected; otherwise trans ="),
    r(c("predictmeans", "covariatemeans"), "responsen", "nothing: the response is read from the model"),
    # printing (D10)
    r(c("predictmeans", "Kmatrix"), "prtnum", print_dig),
    r("predictmeans", "ndecimal", print_dig),
    r(c("semipred", "semireg", "semireg_tmb", "permmodels", "bootmodels"), "prt", "print(x)"),
    # plots (D10): analysis functions keep only 'plot'
    plot_set("predictmeans", c("lsd_bar", "plotord", "lineplot", "plottitle", "plotxlab", "plotylab", "jitterv", "basesz")),
    r("predictmeans", c("mplot", "barplot", "pplot", "bkplot", "which.plots"), "plot = c(\"mean\", \"bar\", \"pval\", \"back\") (the plots wanted)"),
    r("predictmeans", "prtplt", "plot = TRUE: plots are drawn and returned (autoplot())"),
    plot_set("covariatemeans", c("trellis", "plotord", "mtitle", "ci", "point", "jitterv")),
    plot_set("semipred", c("trellis", "scales", "plotord", "ci", "point", "jitterv", "threeD")),
    r(c("predictmeans", "covariatemeans", "residplot", "CookD", "PMplot"), "newwd", "nothing: ggplot2 draws in the current device"),
    r("residplot", "ask", "nothing: one ggplot2 figure"),
    r("PMplot", c("margin", "legendx"), "nothing: ggplot2 sizes the margins and the legend"),
    # renamed for one vocabulary
    r("semipred", "semireg", "model = %s (the first argument, as in every function)", "same"),
    r("semipred", "simultaneous", "band = \"simultaneous\""),
    r("residplot", "level", "re_level = %s", "same"),
    r("CookD", "idn", "id = %s", "same"),
    r("permlmer", "lmer0", "model0 = %s", "same"),
    r("permlmer", "lmer1", "model1 = %s", "same"),
    # semireg (D5, D9)
    r(c("semireg", "semireg_tmb"), "smoothZ", "smooths in the formula, e.g. y ~ s(x, bs = \"os\") for smoothZ = list(sm = smZ(x)) (see vignette(\"migrating\"))"),
    r("semireg", c("ncenter", "nscale", "resp_scale"), "nothing: the variables are used as they are (centring moved knots given in the data's units)"),
    r("semireg", "ngbinomial", "semireg_tmb(family = glmmTMB::nbinom2) (the maximum-likelihood fit)"),
    r(c("semireg", "semireg_tmb"), "predict_info", "nothing: the prediction information is always kept"),
    r("permmodels", "fo", "nothing: the formula is read from the model")
  )
}
pm_defunct_table <- pm_defunct_rows()

## Removed functions and their 2.0 replacement (D5-D9).
pm_defunct_funs <- c(
  smZ = "smooths in the formula of semireg(), e.g. y ~ s(x, bs = \"os\") for smoothZ = list(sm = smZ(x)) (see vignette(\"migrating\"))",
  permglmer = "permlmer()",
  permglmmTMB = "permlmer()",
  permanova.lmer = "permmodels(model, type = ) (Freedman-Lane permutations; p-values differ from permanova.lmer()'s refits)",
  ci_mcp = "predictmeans(..., letterCI = TRUE): letters from comparison intervals that agree with the tests",
  R2_glmm = "R2_model(model) (the same values, as numbers)"
)

## The message for a removed argument: its 2.0 spelling with the value given
## where the table says how.
pm_defunct_msg <- function(row, value) {
  use <- row$use
  if (nzchar(row$value) && grepl("%s", use, fixed = TRUE)) {
    v <- if (row$value == "one_minus" && is.numeric(value) && length(value) == 1 && value > 0 && value < 0.5) {
      format(1 - value)
    } else if (row$value == "one_minus") {
      "0.95"
    } else {
      paste(deparse(value, width.cutoff = 60L)[1], collapse = "")
    }
    use <- sprintf(use, v)
  }
  paste0("'", row$arg, "' was removed in predictmeans 2.0: use ", use, ".")
}

## Stop when a call names a removed argument of 'fun'. 'call' is
## match.call() of the front end (its '...' holds the removed names, which
## are no longer formals); 'env' evaluates a value for the message (never
## needed for the stop itself). With 'strict', any other name in '...' stops
## too ("unused argument"), as R does for functions without '...'.
pm_defunct <- function(fun, call, env = parent.frame(2), strict = FALSE, formals_of = NULL) {
  given <- names(as.list(call))[-1]
  given <- given[!is.na(given) & nzchar(given)]
  if (!length(given)) return(invisible(NULL))
  tab <- pm_defunct_table[pm_defunct_table$fun == fun, , drop = FALSE]
  # match.call() gives the function's own arguments in full; the other names
  # (in '...') may abbreviate a removed one, as R allowed (bar = for barplot =)
  own <- if (is.null(formals_of)) character(0) else names(formals_of)
  k <- ifelse(given %in% own, NA_integer_, pmatch(given, tab$arg, duplicates.ok = TRUE))
  hit <- which(!is.na(k))
  if (length(hit)) {
    msgs <- vapply(hit, function(h) {
      row <- tab[k[h], , drop = FALSE]
      expr <- call[[given[h]]]
      value <- if (nzchar(row$value)) tryCatch(eval(expr, env), error = function(e) expr) else NULL
      # the value when it is short (level = a -> conf.level = 0.95), else what was written (model = fit)
      if (!(is.atomic(value) && length(value) <= 10)) value <- expr
      pm_defunct_msg(row, value)
    }, character(1))
    stop(paste0(fun, "(): ", paste(msgs, collapse = "\n  "),
                "\nSee vignette(\"migrating\") for every change from 1.1.1."), call. = FALSE)
  }
  if (strict && !is.null(formals_of)) {
    extra <- setdiff(given, names(formals_of))
    if (length(extra)) stop(fun, "(): unused argument", if (length(extra) > 1) "s", ": ",
                            paste(extra, collapse = ", "), call. = FALSE)
  }
  invisible(NULL)
}

## The body of a removed function.
pm_defunct_fun <- function(fun) {
  .Defunct(msg = paste0(fun, "() was removed in predictmeans 2.0: use ", pm_defunct_funs[[fun]],
                        ".\nSee vignette(\"migrating\") for every change from 1.1.1."))
}
