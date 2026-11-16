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
#'   \code{(s(x, k = 10) | g)} \tab random curve per level of \code{g}: \code{(x | g)} (correlated intercepts and slopes) plus a spline deviation per level, one variance \tab none \cr
#'   \code{s(x1, x2)} \tab isotropic surface (thin plate) \tab \code{x1 + x2} \cr
#'   \code{s(x, bs = "os", k = 10)} \tab O'Sullivan spline with \code{k} interior knots, as \code{smZ(x, k = 10)}; \code{xt = list(knots = , range = )} sets the interior knots and the range \tab \code{x} \cr
#'   \code{s(x1, x2, bs = "ztps", k = 6)} \tab thin plate of \code{Ztps()}, as \code{smZ(cbind(x1, x2), type = "Ztps")}; \code{xt = list(knots = )} sets the knots (a two-column matrix) \tab \code{x1 + x2} \cr
#'   \code{s(x, bs = "cs")}, \code{"ts"} \tab shrinkage smooths (no unpenalised part) \tab none \cr
#'   \code{s(x, f)} \tab common curve + a deviation curve per level of \code{f}, one variance for the deviations \tab \code{x + f + f:x} \cr
#'   \code{s(x, f1:f2)} \tab the same with the cells of \code{f1:f2} \tab \code{x + f1:f2 + f1:f2:x} \cr
#'   \code{s(x, f1*f2)} \tab common curve + deviation curves for \code{f1}, \code{f2} and \code{f1:f2}, a variance each \tab \code{x + f1*f2 + f1:x + f2:x + f1:f2:x} \cr
#'   \code{s(x1, x2, f)} \tab common surface + a deviation surface per level \tab \code{x1 + x2 + f + f:x1 + f:x2} \cr
#' }
#' A factor term among the arguments of \code{s()} gives this ANOVA-like split (not an mgcv form): the deviation
#' curves of a factor term are random effects with one variance, so \code{varcomp()} shows how much the curves differ
#' by \code{f1}, by \code{f2} and by their interaction, and a variance of zero means the curves of that term have the
#' same shape (their levels still differ by the fixed terms). The pieces are set up without constraints between them:
#' the mixed model separates them by their variances, as nested random effects. Each level's curve is the fixed part
#' plus the common curve plus its deviations.
#'
#' A random curve \code{(s(x) | g)} (O'Sullivan basis unless \code{bs} is given) is the group-specific curve model of
#' Durban et al. (2005) and Harezlak, Ruppert and Wand (2018, section 9.2): each level of \code{g} (a subject) has its
#' own straight line, \code{(x | g)}, and its own spline deviation, all levels sharing one variance; e.g. heights of
#' boys over age, \code{height ~ s(age, bs = "os", k = 20) + (s(age, k = 10) | id)}. Unlike \code{s(x, g, bs = "fs")},
#' the intercepts and slopes are correlated, as \code{pdSymm(~ x)} in \code{nlme::lme}. The random curves are random
#' effects: \code{semipred()} and \code{predictmeans()} give the population curves without them, and
#' \code{anova()} tests whether they are needed (their variance is zero: each subject's curve is a straight line).
#' The fits equal \code{gam(..., method = "REML")} for gaussian models (with \code{s(g, bs = "re")} for a random
#' intercept), and \code{semipred()} and \code{predictmeans()} give the curves, means and their standard errors of
#' \code{predict.gam()}. The variables of formula smooths are not centred or scaled (\code{ncenter} and
#' \code{nscale} apply to the other numeric covariates). Tensor products (\code{te()}, \code{t2()}, \code{ti()}) are
#' not available yet. Smooths are given either in the formula or in \code{smoothZ}, not both. Objects named in
#' \code{s()} (e.g. \code{xt = list(knots = kn)}) are taken from the environment of the formula.
#'
#' The penalised-spline examples fitted by \code{nlme::lme} in Wand's course scripts and in \pkg{HRW} (additive
#' models, additive mixed models, thin-plate surfaces with given knots, group-specific curves with a contrast curve)
#' give the same REML log-likelihood and fitted values with \code{semireg()} (data in
#' \code{system.file("extdata", "semireg", package = "predictmeans")}; \code{tests/testthat/test-semireg-golden.R}).
#' The optimizer is \code{nloptwrap} (lme4's default for \code{lmer}) since 2.0: on these examples it reached the
#' same optimum as \code{bobyqa} (1.x's default; \code{control = lmerControl(optimizer = "bobyqa")}) faster, also
#' for binomial and Poisson fits. Each smooth block of a formula smooth is scaled as a random intercept (sum of
#' squares equal to the number of observations), so that the optimizer starts near the optimum; the variance of a
#' smooth block reported by \code{VarCorr()} or \code{varcomp()} is on that scale.
#'
#' Moving from the \code{smoothZ} form: \code{y ~ x, smoothZ = list(sm = smZ(x, k = 10))} is
#' \code{y ~ s(x, bs = "os", k = 10)}; \code{smZ(x, by = f, group = TRUE)} is \code{s(x, by = f, bs = "os")};
#' \code{smZ(x, by = f)} (one variance) is \code{s(x, by = f, bs = "os", id = 1)}; \code{smZ(cbind(x1, x2),
#' type = "Ztps")} is \code{s(x1, x2, bs = "ztps")}; \code{smZ(x, intKnots = kn, range.x = r)} is
#' \code{s(x, bs = "os", xt = list(knots = kn, range = r))}. The fits are the same (\code{vignette("smooth")}).
#' @usage
#' semireg(formula, data, family = NULL, ngbinomial=FALSE, REML = TRUE, 
#'         smoothZ = list(), ncenter=TRUE, nscale=FALSE, resp_scale=FALSE, 
#'         control = lmerControl(optimizer="nloptwrap"), start = NULL, 
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
#' @param control A list (of correct class, resulting from lmerControl() or glmerControl() respectively) containing control parameters, including the nonlinear optimizer to be used and parameters to be passed through to the nonlinear optimizer, see the *lmerControl documentation for details. The default is \code{lmerControl(optimizer = "nloptwrap")} (\code{"bobyqa"} in 1.x), also used for generalized fits.
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
#' Durban, M., Harezlak, J., Wand, M.P. and Carroll, R.J. (2005). Simple fitting of subject-specific curves for
#' longitudinal data. \emph{Statistics in Medicine}, \bold{24}, 1153-1167.
#'
#' Harezlak, J., Ruppert, D. and Wand, M.P. (2018). \emph{Semiparametric Regression with R}. Springer, New York.
#'
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
                     control = lmerControl(optimizer="nloptwrap"), start = NULL, 
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
    # an lmerMod that keeps the deviance function of the fitted model (with
    # the smooth designs, not the placeholder terms): Satterthwaite df and
    # varcomp() need it, a refit from the call would not have the smooths
    semer@call <- lmerc
    semer <- methods::new("semerMod", semer,
                          semireg_info = list(devfun = fit$devfun, data = st$data, smoothZ = st$smoothZN,
                                              formula = st$formula, REML = isTRUE(REML), control = control))
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
  # what a refit needs (tests of the smooth terms, anova.semireg())
  refit <- list(engine = "lme4", formula = st$formula, data = st$data, smoothZ = st$smoothZN,
                family = if (gaus || ngbinomial) NULL else family, ngbinomial = ngbinomial, REML = REML,
                control = control, weights = semer@resp$weights,
                offset = if (any(semer@resp$offset != 0)) semer@resp$offset)
  pm_sr_result(semer, st, mc, info, list(lmerc = lmerc, devfun0 = fit$devfun0, refit = refit))
}

## A gaussian semireg() fit: an lmerMod that keeps the deviance function of
## the fitted model (lmerTest is no longer needed since 2.0).
methods::setClass("semerMod", contains = "lmerMod", slots = c(semireg_info = "list"))

## Refit of a gaussian semireg() fit with a new random part (as
## lmerTest::ranova() does): the smooth terms (1 | name) still in the formula
## keep their spline bases, instead of being refitted as plain random
## intercepts. anova(fit) tests the smooth terms without it.
update.semerMod <- function(object, formula = NULL, ...) {
  info <- object@semireg_info
  if (is.null(formula) || length(list(...))) {
    stop("update() of a semireg() fit only takes a new random part (as lmerTest::ranova()); refit with semireg().",
         call. = FALSE)
  }
  if (is.character(formula))
    formula <- stats::as.formula(if (length(formula) == 3) paste(formula[2], formula[1], formula[3])
                                 else paste(formula, collapse = " "))
  # the response and the fixed part from the original formula, only the
  # random terms from the new one
  f <- info$formula
  bars <- pm_findbars(formula)
  bar_grp <- vapply(bars, function(b) paste(all.vars(b[[3]]), collapse = ":"), character(1))
  sm_names <- names(info$smoothZ)
  keep <- intersect(sm_names, bar_grp)
  other <- bars[!bar_grp %in% sm_names]
  fixed <- pm_nobars(f)
  e <- new.env(parent = asNamespace("predictmeans"))
  rhs <- paste(c(paste(deparse(fixed[[3]], width.cutoff = 500L), collapse = " "),
                 vapply(other, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"),
                        character(1))), collapse = " + ")
  f_new <- stats::as.formula(paste(paste(deparse(f[[2]]), collapse = ""), "~", rhs), env = e)
  e$sr_data <- info$data
  e$sr_Z <- info$smoothZ[keep]
  cl <- call("semireg", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z), REML = info$REML,
             ncenter = FALSE, nscale = FALSE, prt = FALSE, predict_info = FALSE)
  fit <- suppressMessages(suppressWarnings(eval(cl, e)))
  if (is.list(fit) && !is.null(fit$semer)) fit$semer else fit
}
