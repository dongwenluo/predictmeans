#' Fitting Semi Parametric Models Using lme4 Ecosystem
#'
#' @description
#' Fit a semi parametric model based on lme4 ecosystem including \code{lmer}, \code{glmer} and \code{glmer.nb}.
#' @details
#' A semi parametric model can be parameterized as a linear (or generalized linear) mixed model in which its random effects are smooth functions of some covariates (named <e2><80><98>smooth term<e2><80><99>). \code{semireg} follows the approach suggested by Wand and Ormerod (2008) and represents the 'smooth term' using O'Sullivan-type of Z.
#' @section Smooth terms in the formula:
#' Smooths can be written in the model formula, as for \code{mgcv::gam}, next to fixed and lme4 random terms:
#' \code{y ~ A + s(x, by = A) + (1 | Block)}. Each smooth is set up by mgcv (bases, knots, identifiability
#' constraints) and split, as in \pkg{gamm4}, into an unpenalised part, which is added to the fixed effects as
#' ordinary terms, and penalised parts, each a block of random effects with one variance. The fixed terms are
#' added automatically (terms already written are not doubled), so a smooth by a factor cannot lose its group
#' slopes:
#' \tabular{lll}{
#'   \strong{Written} \tab \strong{Meaning} \tab \strong{Fixed terms added} \cr
#'   \code{s(x)}, \code{s(x, bs = "cr", k = 12)} \tab one curve \tab \code{x} \cr
#'   \code{s(x, by = f)} \tab a curve per level of \code{f}, a variance each (\code{id = 1}: one shared) \tab \code{f + f:x} \cr
#'   \code{s(x, by = z)}, \code{z} numeric \tab varying coefficient \eqn{f(x) z} \tab \code{z + z:x} \cr
#'   \code{s(x, f, bs = "fs")} \tab random curves per level, one variance \tab none \cr
#'   \code{s(x1, x2)} \tab isotropic surface (thin plate) \tab \code{x1 + x2} \cr
#'   \code{s(x, bs = "os", k = 10)} \tab O'Sullivan spline with \code{k} interior knots, as \code{smZ(x, k = 10)} \tab \code{x} \cr
#'   \code{s(x1, x2, bs = "ztps", k = 6)} \tab thin plate of \code{Ztps()}, as \code{smZ(cbind(x1, x2), type = "Ztps")} \tab \code{x1 + x2} \cr
#'   \code{s(x, bs = "cs")}, \code{"ts"} \tab shrinkage smooths (no unpenalised part) \tab none \cr
#' }
#' The fits equal \code{gam(..., method = "REML")} for gaussian models (with \code{s(g, bs = "re")} for a random
#' intercept), and \code{semipred()} and \code{predictmeans()} give the curves, means and their standard errors of
#' \code{predict.gam()}. The variables of formula smooths are not centred or scaled (\code{ncenter} and
#' \code{nscale} apply to the other numeric covariates). Tensor products (\code{te()}, \code{t2()}, \code{ti()}) are
#' not available yet. Smooths are given either in the formula or in \code{smoothZ}, not both.
#' @usage
#' semireg(formula, data, family = NULL, ngbinomial=FALSE, REML = TRUE, 
#'         smoothZ = list(), ncenter=TRUE, nscale=FALSE, resp_scale=FALSE, 
#'         control = lmerControl(optimizer="bobyqa"), start = NULL, 
#'         verbose = FALSE, drop.unused.levels=TRUE, subset, weights, 
#'         offset, contrasts = NULL,  prt=TRUE, predict_info=TRUE, ...)
#' @param formula A two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors.
#' @param data A data frame or list containing the model response variable and covariates required by the formula. By default the variables are taken from environment(formula and smoothZ), typically the environment from which semireg is called.
#' @param family A GLM family, see glm and family.
#' @param ngbinomial Logical scalar - Should a negative binomial GLMMs be used?                  .
#' @param REML Logical scalar - Should the estimates be chosen to optimize the REML criterion (as opposed to the log-likelihood)?
#' @param smoothZ The older form of the smooths (see \sQuote{Smooth terms in the formula} for the newer one): a list includes a set of smooth Z matrixs (called 'smooth term') used in the mixed effects model, the name of 'smooth term' should be different any variables in the model, each 'smooth term' is the result of function \code{smZ}. e.g.
#' smoothZ=list(sm1=smZ(x1), sm2=smZ(x2, by=f1), sm3=smZ(x3, by=f2, group=TRUE), ...) where 'sm1' to 'sm3' should be new variable names in the \code{data}, and x1 to x3 are covariates, and f1, f2 are factors.
#' @param ncenter Logical scalar - Should the numeric predictors to be centered or not?
#' @param nscale Logical scalar - Should the numeric predictors to be scaled or not?
#' @param resp_scale Logical scalar - Should the response be involved in the scaling action or not?
#' @param control A list (of correct class, resulting from lmerControl() or glmerControl() respectively) containing control parameters, including the nonlinear optimizer to be used and parameters to be passed through to the nonlinear optimizer, see the *lmerControl documentation for details.
#' @param start Starting value list as used by lmer or glmer.
#' @param verbose Passed on to fitting lme4 fitting routines.
#' @param drop.unused.levels By default unused levels are dropped from factors before fitting. For some smooths involving factor variables you might want to turn this off. Only do so if you know what you are doing.
#' @param subset An optional expression indicating the subset of the rows of data that should be used in the fit. This can be a logical vector, or a numeric vector indicating which observation numbers are to be included, or a character vector of the row names to be included. All observations are included by default.
#' @param weights An optional vector of <e2><80><98>prior weights<e2><80><99> to be used in the fitting process. Should be NULL or a numeric vector.
#' @param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See model.offset.
#' @param contrasts An optional list. See the contrasts.arg of model.matrix.default.
#' @param prt Logical scalar - Should the info to be print on screen in the middle of the process or not?
#' @param predict_info Logical scalar - Should provide the info for function semipred or not?
#' @param ... Further arguments for passing on to model setup routines.
#' @return
#' \item{semer}{A mer model used in the fitting.}
#' \item{data}{A data.frame with generated variables in the fitting.}
#' \item{fomul_vars}{Name of variables in the formula of semireg model.}
#' \item{sm_vars}{Name of variables in the smoothZ list.}
#' \item{smoothZ_call}{A call used to produce smooth terms in the fitting.}
#' \item{knots_lst}{Knots used in each smooth term in the fitting.}
#' \item{range_lst}{Range of covariate used in each smooth term in the fitting.}
#' \item{cov_lst}{Covariance matrix list for each smooth term.}
#' \item{u_lst}{Random effects list for each smooth term.} 
#' \item{type_lst}{Smooth type list of smooth terms.} 
#' \item{CovMat}{Covariance matrix for all smooth terms.}
#' \item{Cov_ind}{Covariance matrix index for each smooth term.}
#' \item{Cov_indN}{Covariance matrix index for each smooth term when \code{group=TRUE} in \code{smoothZ} argument.}
#' \item{df}{Degree of freedom of all random terms.}
#' \item{lmerc}{Call used in the mer model in the fitting.}
#' \item{smooths}{For smooths written in the formula: the mgcv smooth objects and the maps that give the random-effect
#'   design of each block for new data (used by \code{semipred()} and \code{predictmeans()}).}
#' @references
#' Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with O'Sullivan penalized splines. \emph{Australian and New Zealand Journal of Statistics.} \bold{50}, 179-198.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(WarsawApts, package = "HRW")
#'   WarsawApts$district <- factor(WarsawApts$district)
#'   # smooths in the formula: a curve per district; district and its slopes are added as fixed terms
#'   fit0 <- semireg(areaPerMzloty ~ s(construction.date, by = district), data = WarsawApts, prt = FALSE)
#'   semipred(fit0, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"), prt = FALSE)$plt
#'   predictmeans(fit0, "district", covariate = c(construction.date = 1980))$means
#' }
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(WarsawApts, package = "HRW")
#'   fit1 <- semireg(areaPerMzloty ~ construction.date,
#'                   smoothZ=list(
#'                     grp=smZ(construction.date, k=25)
#'                   ),
#'                   data = WarsawApts)
#'   sp_out1 <- semipred(fit1, "construction.date", "construction.date")
#'
#'   WarsawApts$district <- factor(WarsawApts$district)
#'   fit2 <- semireg(areaPerMzloty ~ construction.date*district, resp_scale = TRUE,
#'                   smoothZ=list(group=smZ(construction.date, k=15,
#'                                          by = district, group=TRUE)),
#'                   data=WarsawApts)
#'   sp_out2_1 <- semipred(fit2, "district", "construction.date")
#'   sp_out2_2 <- semipred(fit2, "district", "construction.date", contr=c(2,1))
#' }
#' \donttest{
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(indonRespir, package = "HRW")
#'   fit3 <- semireg(respirInfec ~ age+vitAdefic + female + height
#'                   + stunted + visit2 + visit3 + visit4  + visit5 + visit6+(1|idnum),
#'                   smoothZ=list(
#'                     grp=smZ(age)
#'                   ),
#'                   family = binomial,
#'                   data = indonRespir)
#'   sp_out3 <- semipred(fit3, "age", "age")
#'   library(ggplot2)
#'   print(sp_out3$plt +
#'     geom_rug(data = subset(indonRespir, respirInfec==0), sides = "b", col="deeppink") +
#'     geom_rug(data = subset(indonRespir, respirInfec==1), sides = "t", col="deeppink") +
#'     coord_cartesian(ylim = c(0, 0.2)))
#' }
#' }
#' @name semireg
NULL

semireg <- function(formula, data, family = NULL, ngbinomial=FALSE, REML = TRUE, 
                     smoothZ = list(), ncenter=TRUE, nscale=FALSE, resp_scale=FALSE, 
                     control = lmerControl(optimizer="bobyqa"), start = NULL, 
                     verbose = FALSE, drop.unused.levels=TRUE, subset, weights, 
                     offset, contrasts = NULL, prt=TRUE, predict_info=TRUE, ...)
{
  mc <- match.call()
  env <- parent.frame()
  environment(formula) <- env
  mc$formula <- formula
  # data, smooth designs and the formula with placeholders (R/semireg-engine.R)
  st <- pm_sr_setup(mc, data, substitute(smoothZ), env, ncenter = ncenter, nscale = nscale, resp_scale = resp_scale)

  gaus <- FALSE
  if (is.null(family)) {
    gaus <- TRUE
  }else{
    ## copied from glm()
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = env)
    if (is.function(family)) 
      family <- family()
    if (!inherits(family, "family")) stop("unknown family type")
    gaus <- family$family == "gaussian" && family$link == "identity"
  }  
  if (ngbinomial) gaus <- FALSE       # a negative binomial GLMM (theta by ML, R/semireg-engine.R)
  
  lmerc <- mc                         # the call to lmer / glmer
  lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
  if (ngbinomial) {
    lmerc[[1]] <- as.name("glmer")
    lmerc$family <- poisson()
  }
  lmerc <- lmerc[!names(lmerc) %in% c("smoothZ", "ncenter", "nscale", "resp_scale", "ngbinomial", "prt", "predict_info")]
  if (!gaus) lmerc$REML <- NULL
  if (!length(st$Zt))                 # no smooth: lmer / glmer itself
    return(eval.parent(lmerc))
  lmerc$formula <- st$formula

  # the structure without fitting, the smooth designs in Z, one fit
  fit <- pm_sr_fit_lme4(lmerc, st, gaus, family, ngbinomial, control, start, verbose, drop.unused.levels, mc, env, ...)
  semer <- fit$semer
  if (gaus) {
    # lmerTest's class (Satterthwaite df), with the deviance function of the
    # fitted model (with the smooth designs, not the placeholder terms)
    semer <- suppressWarnings(as_lmerModLT(semer, fit$devfun))
    semer@call <- lmerc
    # what a refit needs, so that update() (used by lmerTest::ranova) refits
    # through semireg() with the same smooth bases
    semer <- methods::new("semerLmerTest", semer,
                          semireg_info = list(data = st$data, smoothZ = st$smoothZN, formula = st$formula,
                                              REML = isTRUE(REML), control = control, devfun = fit$devfun))
  }
  if (prt) print(semer)

  info <- NULL
  if (predict_info) {
    Z <- getME(semer, "Z")
    w <- if (gaus) semer@resp$weights else semer@resp$sqrtWrkWt()^2
    u_ind <- vapply(getME(semer, "Ztlist"), function(x) x@Dim[1], numeric(1))
    names(u_ind) <- sub("\\.\\(Intercept\\)$", "", names(u_ind))
    info <- pm_sr_info(getME(semer, "X"), Z, getME(semer, "Lambda"), w, unlist(ranef(semer)), u_ind,
                       names(st$Zt), st$smoothZm)
  }
  pm_sr_result(semer, st, mc, info, list(lmerc = lmerc))
}

## A semireg() lmer fit that remembers its smooth bases (see update below).
methods::setClass("semerLmerTest", contains = "lmerModLmerTest",
                  slots = c(semireg_info = "list"))

## Refit of a semireg() fit with a new formula (as lmerTest::ranova() does):
## the smooth terms (1 | name) still in the formula keep their spline bases,
## instead of being refitted as ordinary random intercepts.
update.semerLmerTest <- function(object, formula = NULL, ...) {
  info <- object@semireg_info
  if (is.null(formula) || length(list(...))) {
    upd <- utils::getS3method("update", "lmerModLmerTest", envir = asNamespace("lmerTest"))
    return(upd(object, formula = formula, ...))
  }
  # ranova() only drops random terms: take the response and the fixed part
  # from the original formula (newer R/lmerTest can pass a reduced formula
  # whose response is damaged, e.g. 'log' for 'log(grade)'), and only the
  # random terms from the new one
  if (is.character(formula))
    formula <- stats::as.formula(if (length(formula) == 3) paste(formula[2], formula[1], formula[3])
                                 else paste(formula, collapse = " "))
  f <- info$formula
  bars <- pm_findbars(formula)
  bar_grp <- vapply(bars, function(b) paste(all.vars(b[[3]]), collapse = ":"), character(1))
  sm_names <- names(info$smoothZ)
  keep <- intersect(sm_names, bar_grp)
  other <- bars[!bar_grp %in% sm_names]
  fixed <- if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(f) else lme4::nobars(f)
  rhs <- paste(c(paste(deparse(fixed[[3]], width.cutoff = 500L), collapse = " "),
                 vapply(other, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"),
                        character(1))), collapse = " + ")
  f_new <- stats::as.formula(paste(paste(deparse(f[[2]]), collapse = ""), "~", rhs),
                             env = environment(stats::formula(object)))
  e <- new.env(parent = asNamespace("predictmeans"))
  e$sr_data <- info$data
  e$sr_Z <- info$smoothZ[keep]
  cl <- call("semireg", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z),
             REML = info$REML, ncenter = FALSE, nscale = FALSE, prt = FALSE, predict_info = FALSE)
  fit <- suppressMessages(suppressWarnings(eval(cl, e)))
  if (is.list(fit) && !is.null(fit$semer)) fit$semer else fit
}
