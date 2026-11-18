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
#'   \code{s(x1, x2, bs = "ztps", k = 6)} \tab thin plate of \code{Ztps()}, as \code{smZ(cbind(x1, x2), type = "Ztps", k = 6)} with \code{k = 6}; \code{xt = list(knots = )} sets the knots (a two-column matrix) \tab \code{x1 + x2} \cr
#'   \code{s(x, bs = "user", xt = list(basis = f))} \tab a basis of your own: \code{f(x)} (\code{f(x1, x2)} for two covariates) returns a matrix with a row per value, whose columns get iid random coefficients with one variance (as a matrix in \code{smoothZ}); checked before the fit (see below); \code{xt = list(basis = f, null = "constant")} or \code{"none"} for a smaller fixed part \tab \code{x} (\code{x1 + x2}) \cr
#'   \code{t2(x1, x2)}, \code{t2(x1, x2, bs = c("cr", "ps"), k = c(6, 8))} \tab tensor-product surface (covariates on different scales), a variance per penalty; \code{te()} is fitted as \code{t2()} \tab \code{x1 + x2 + x1:x2} \cr
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
#' \code{predict.gam()}. The variables are used as they are (not centred or scaled). Tensor products: \code{t2(x1, x2)} (Wood, Scheipl and Faraway 2013), a surface for covariates on different scales,
#' whose penalties (one per product of the marginal penalised and unpenalised parts; \code{full = TRUE} for more)
#' are independent random effects, each with its own variance; the fixed terms \code{x1 + x2 + x1:x2} are added,
#' and \code{by = f} gives a surface per level. It is set up without mgcv's sum-to-zero constraint on the whole
#' product (which leaves a little of the penalised space in the unpenalised one), so its fixed part is exactly those
#' terms; the fit equals \code{gam()} with that basis and its penalties, and \code{gam()}'s own \code{t2()} when there
#' are no other random effects. \code{te()} is fitted as \code{t2()} (with a message): its overlapping penalties
#' cannot be written as random effects; \code{ti()} is not available. Objects named in
#' \code{s()} (e.g. \code{xt = list(knots = kn)}) are taken from the environment of the formula.
#'
#' A basis of your own (\code{bs = "user"}) must give the same row for the same covariate value whatever data it is
#' given, as prediction evaluates it on new values: knots and range belong in the function, e.g.
#' \code{f <- function(x) ZOSull(x, intKnots = c(2, 4, 6, 8), range.x = c(0, 10))}. It is checked before the fit, and
#' the fit stops with what is wrong and how to put it right when it fails on the data, does not return a numeric
#' matrix with a row per value (also for a single value), gives missing or infinite values, gives different values on
#' two calls, gives other rows for the same values on part of the data (knots or range taken from the data), or has
#' linearly dependent columns or columns in the span of its fixed part; a warning says when it fails just beyond the
#' ends of the data (the slopes of \code{semipred(deriv = )} there need it). Its values are checked again at
#' prediction.
#'
#' The penalised-spline examples fitted by \code{nlme::lme} in Wand's course scripts and in \pkg{HRW} (additive
#' models, additive mixed models, thin-plate surfaces with given knots, group-specific curves with a contrast curve)
#' give the same REML log-likelihood and fitted values with \code{semireg()} (\code{tests/testthat/test-semireg-golden.R});
#' their data are in the package: \code{\link{growthIndiana}}, \code{\link{SBMDfem}}, \code{\link{ozoneCalif}},
#' \code{\link{scallop}}, \code{\link{wageTUnion}}, \code{\link{copper}}, whose help pages show the analyses.
#' The optimizer is \code{nloptwrap} (lme4's default for \code{lmer}) since 2.0: on these examples it reached the
#' same optimum as \code{bobyqa} (1.x's default; \code{control = lmerControl(optimizer = "bobyqa")}) faster, also
#' for binomial and Poisson fits. Each smooth block of a formula smooth is scaled as a random intercept (sum of
#' squares equal to the number of observations), so that the optimizer starts near the optimum; the variance of a
#' smooth block reported by \code{VarCorr()} or \code{varcomp()} is on that scale.
#'
#' Moving from 1.x: \code{smZ()} and \code{smoothZ} were removed in 2.0; \code{vignette("migrating")} shows each
#' \code{smoothZ = list(sm = smZ(...))} form beside its formula form, which gives the same fit (checked in
#' \code{dev/check-smz-migration.R}). The variables are no longer centred (1.x's \code{ncenter = TRUE}), which put knots
#' given in the data's units (\code{smZ(x, intKnots = )}) in the wrong place.
#' @usage
#' semireg(formula, data, family = NULL, REML = TRUE, control = lmerControl(optimizer = "nloptwrap"),
#'         start = NULL, verbose = FALSE, drop.unused.levels = TRUE, subset, weights,
#'         offset, contrasts = NULL, ...)
#' @param formula A two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors.
#' @param data A data frame or list containing the model response variable and covariates required by the formula. By default the variables are taken from environment(formula), typically the environment from which semireg is called.
#' @param family A GLM family, see glm and family.
#' @param REML Logical scalar - Should the estimates be chosen to optimize the REML criterion (as opposed to the log-likelihood)?
#' @param control A list (of correct class, resulting from lmerControl() or glmerControl() respectively) containing control parameters, including the nonlinear optimizer to be used and parameters to be passed through to the nonlinear optimizer, see the *lmerControl documentation for details. The default is \code{lmerControl(optimizer = "nloptwrap")} (\code{"bobyqa"} in 1.x), also used for generalized fits.
#' @param start Starting value list as used by lmer or glmer.
#' @param verbose Passed on to fitting lme4 fitting routines.
#' @param drop.unused.levels By default unused levels are dropped from factors before fitting. For some smooths involving factor variables you might want to turn this off. Only do so if you know what you are doing.
#' @param subset An optional expression indicating the subset of the rows of data that should be used in the fit. This can be a logical vector, or a numeric vector indicating which observation numbers are to be included, or a character vector of the row names to be included. All observations are included by default.
#' @param weights An optional vector of <e2><80><98>prior weights<e2><80><99> to be used in the fitting process. Should be NULL or a numeric vector.
#' @param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See model.offset.
#' @param contrasts An optional list. See the contrasts.arg of model.matrix.default.
#' @param ... Further arguments for passing on to model setup routines. Arguments removed in 2.0 (\code{smoothZ},
#'   \code{ncenter}, \code{nscale}, \code{resp_scale}, \code{prt}, \code{predict_info}, \code{ngbinomial}) stop with their 2.0 spelling; nothing is printed during the fit
#'   (\code{print(fit$semer)}). For a negative binomial model use \code{\link{semireg_tmb}(family =
#'   glmmTMB::nbinom2)}: the 1.x \code{ngbinomial = TRUE} estimated theta from a Poisson fit, which is not the
#'   maximum-likelihood fit.
#' @return
#' A list of class \code{semireg}; \code{print()} shows the fitted model and the random-effect term of each
#' smooth, \code{summary()} the summary of the fitted model, \code{anova()} tests of its terms. Its elements:
#' \item{semer}{A mer model used in the fitting.}
#' \item{data}{A data.frame with generated variables in the fitting.}
#' \item{fomul_vars}{Name of variables in the formula of semireg model.}
#' \item{sm_vars}{Names of the variables of the smooths.}
#' \item{cov_lst}{Covariance matrix list for each smooth term.}
#' \item{u_lst}{Random effects list for each smooth term.} 
#' \item{CovMat}{Covariance matrix for all smooth terms.}
#' \item{Cov_ind}{Covariance matrix index for each smooth term.}
#' \item{Cov_indN}{Covariance matrix index for each smooth term of several blocks (e.g. \code{s(x, by = f)}).}
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
#'
#' Wood, S.N., Scheipl, F. and Faraway, J.J. (2013). Straightforward intermediate rank tensor product smoothing in
#' mixed models. \emph{Statistics and Computing}, \bold{23}, 341-360.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' ## 1. Curves in a mixed model (nlme::Dialyzer): a curve of the ultrafiltration rate over
#' ##    pressure for each blood flow rate, a random intercept and slope per dialyzer; the
#' ##    fixed terms QB + QB:pressure are added automatically
#' data(Dialyzer, package = "nlme")
#' fm <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
#'               data = Dialyzer)
#' anova(fm)       # fixed effects (Kenward-Roger df) and a likelihood ratio test of each smooth
#' varcomp(fm)     # variance components, the smooths' among them
#' predictmeans(fm, "QB", covariate = c(pressure = 1))$means
#' semipred(fm, "QB", "pressure", contr = c("300", "200"), band = "simultaneous")$plt
#'
#' ## 2. A designed experiment: treatments on the whole plots of six blocks, curves over time;
#' ##    treatment means at a time get Satterthwaite df from the whole-plot stratum
#' set.seed(11)
#' sp <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
#' sp$y <- 2 + sin(sp$t/2) + 0.08*(sp$A == "a2")*sp$t + rnorm(6, 0, 0.5)[sp$Block] +
#'   rnorm(18, 0, 0.3)[interaction(sp$Block, sp$A)] + rnorm(nrow(sp), 0, 0.3)
#' fsp <- semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A),
#'                data = sp)
#' predictmeans(fsp, "A", covariate = c(t = 5), pairwise = TRUE)$pairs
#'
#' ## 3. Two crossed factors: a common curve and deviation curves for f1, f2 and f1:f2, each
#' ##    with its own variance; which deviations are needed?
#' set.seed(1)
#' dx <- expand.grid(x = seq(0, 1, length = 15), f1 = factor(c("p", "q")), f2 = factor(c("u", "v")),
#'                   rep = 1:3)
#' dx$y <- sin(2*pi*dx$x) + (dx$f1 == "q")*sin(pi*dx$x)^2 + rnorm(nrow(dx), 0, 0.2)
#' fx <- semireg(y ~ s(x, f1*f2, bs = "os", k = 6), data = dx)
#' varcomp(fx, test = TRUE)
#'
#' ## 4. An additive model with a varying coefficient: a smooth in x1, and the effect of z
#' ##    changing smoothly with x2 (numeric 'by')
#' set.seed(2)
#' da <- data.frame(x1 = runif(200), x2 = runif(200), z = rnorm(200))
#' da$y <- sin(3*da$x1) + da$z*sin(2*pi*da$x2) + rnorm(200, 0, 0.3)
#' fa <- semireg(y ~ s(x1, bs = "os") + s(x2, by = z, bs = "os"), data = da)
#' anova(fa)$smooth
#'
#' ## 5. A basis of your own (checked before the fit): a seasonal cycle on a straight trend
#' set.seed(3)
#' ds <- data.frame(day = 1:730)
#' ds$y <- 10 + 0.002*ds$day + 2*sin(2*pi*ds$day/365) + rnorm(730)
#' season <- function(d) cbind(sin(2*pi*d/365), cos(2*pi*d/365), sin(4*pi*d/365), cos(4*pi*d/365))
#' fs <- semireg(y ~ day + s(day, bs = "user", xt = list(basis = season, null = "none")),
#'               data = ds)
#' semipred(fs, covariate = "day")$plt
#' ## 6. A tensor-product surface for two covariates on different scales (data(ozoneCalif)),
#' ##    and a curve in a third
#' data(ozoneCalif)
#' fo <- semireg(ozone ~ t2(inversion.base.height, inversion.base.temp) +
#'                 s(daggett.pressure.gradient, bs = "os", k = 15), data = ozoneCalif)
#' anova(fo)$smooth          # a variance per penalty of the surface
#' semipred(fo, covariate = c("inversion.base.height", "inversion.base.temp"))$plt
#' \donttest{
#' ## 7. Binary responses with a random intercept per child (HRW::indonRespir): the probability
#' ##    of a respiratory infection over age
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(indonRespir, package = "HRW")
#'   fb <- semireg(respirInfec ~ s(age, bs = "os", k = 10) + vitAdefic + female + height + stunted +
#'                   visit2 + visit3 + visit4 + visit5 + visit6 + (1 | idnum),
#'                 family = binomial, data = indonRespir)
#'   semipred(fb, "age", "age")$plt
#' }
#'
#' ## 8. Subject-specific curves: each subject's own straight line, (Days | Subject), and its
#' ##    own spline deviation, one variance for all; are the subjects' curves straight lines?
#' fr <- semireg(Reaction ~ s(Days, bs = "os", k = 4) + (s(Days, k = 4) | Subject),
#'               data = lme4::sleepstudy)
#' anova(fr)$smooth
#'
#' ## 9. Curves by group, the price of apartments in Warsaw by district (HRW::WarsawApts)
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(WarsawApts, package = "HRW")
#'   WarsawApts$district <- factor(WarsawApts$district)
#'   fw <- semireg(areaPerMzloty ~ s(construction.date, by = district, bs = "os", k = 15),
#'                 data = WarsawApts)
#'   semipred(fw, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"))$plt
#'   predictmeans(fw, "district", covariate = c(construction.date = 1980))$means
#' }
#' }
#' @name semireg
NULL

semireg <- function(formula, data, family = NULL, REML = TRUE, control = lmerControl(optimizer = "nloptwrap"),
                    start = NULL, verbose = FALSE, drop.unused.levels = TRUE, subset, weights, offset, contrasts = NULL,
                    ...)
{
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("semireg", match.call(), formals_of = formals(sys.function()))
  mc <- match.call()
  mc[[1]] <- pm_semireg_lme4
  eval(mc, parent.frame())
}

## The fit of semireg(). The package's own refits (update(), anova(),
## variance-component tests) call it with the smooth blocks of the fit
## already built ('smoothZ') and predict_info = FALSE: the covariance of the
## smooth coefficients is not needed there.
pm_semireg_lme4 <- function(formula, data, family = NULL, REML = TRUE,
                            smoothZ = list(),
                            control = lmerControl(optimizer="nloptwrap"), start = NULL,
                            verbose = FALSE, drop.unused.levels=TRUE, subset, weights,
                            offset, contrasts = NULL, predict_info=TRUE, ...)
{
  mc <- match.call()
  mc[[1]] <- as.name("semireg")
  env <- parent.frame()
  environment(formula) <- env
  mc$formula <- formula
  # lme4 has no covariance structures: lmer would fit ar1(t + 0 | id) as an
  # unstructured term
  if (length(pm_special_terms(formula))) {
    stop("Covariance structures such as ar1(), ou() or cs() are glmmTMB terms; fit the model with semireg_tmb().",
         call. = FALSE)
  }
  # data, smooth designs and the formula with placeholders (R/semireg-engine.R)
  st <- pm_sr_setup(mc, data, substitute(smoothZ), env)

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

  lmerc <- mc                         # the call to lmer / glmer
  lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
  lmerc <- lmerc[!names(lmerc) %in% c("smoothZ", "predict_info")]
  if (!gaus) lmerc$REML <- NULL
  if (!length(st$Zt))                 # no smooth: lmer / glmer itself
    return(eval.parent(lmerc))
  lmerc$formula <- st$formula

  # the structure without fitting, the smooth designs in Z, one fit
  fit <- pm_sr_fit_lme4(lmerc, st, gaus, family, control, start, verbose, drop.unused.levels, mc, env, ...)
  semer <- fit$semer
  if (gaus) {
    # an lmerMod that keeps the deviance function of the fitted model (with
    # the smooth designs, not the placeholder terms): Satterthwaite df and
    # varcomp() need it, a refit from the call would not have the smooths
    semer@call <- lmerc
    # the model's theta is the vector the deviance function writes into (lme4
    # maps it in C++): a copy of its own, so that evaluating the deviance
    # function (Satterthwaite df, varcomp()) cannot move the fitted model
    semer@theta <- semer@theta + 0
    semer <- methods::new("semerMod", semer,
                          semireg_info = list(devfun = fit$devfun, cache = new.env(parent = emptyenv()),
                                              data = st$data, smoothZ = st$smoothZN,
                                              formula = st$formula, REML = isTRUE(REML), control = control))
  }
  info <- NULL
  if (predict_info) {
    Z <- getME(semer, "Z")
    w <- if (gaus) semer@resp$weights else semer@resp$sqrtWrkWt()^2
    u_ind <- vapply(getME(semer, "Ztlist"), function(x) x@Dim[1], numeric(1))
    names(u_ind) <- sub("\\.\\(Intercept\\)$", "", names(u_ind))
    info <- pm_sr_info(getME(semer, "X"), Z, getME(semer, "Lambda"), w, unlist(ranef(semer)), u_ind,
                       names(st$Zt))
  }
  # what a refit needs (tests of the smooth terms, anova.semireg())
  refit <- list(engine = "lme4", formula = st$formula, data = st$data, smoothZ = st$smoothZN,
                family = if (gaus) NULL else family, REML = REML,
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
  cl <- call("pm_semireg_lme4", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z), REML = info$REML,
             predict_info = FALSE)
  fit <- suppressMessages(suppressWarnings(eval(cl, e)))
  if (is.list(fit) && !is.null(fit$semer)) fit$semer else fit
}
