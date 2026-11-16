#' Fitting Semi Parametric Models Using glmmTMB
#'
#' @description
#' Fit a semi parametric model based on glmmTMB.
#' @details
#' A semi parametric model can be parameterized as a linear (or generalized linear) mixed model in which its random effects are smooth functions of some covariates (named <e2><80><98>smooth term<e2><80><99>). \code{semireg_tmb} follows the approach suggested by Wand and Ormerod (2008) and represents the 'smooth term' using O'Sullivan-type of Z.
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
#' semireg_tmb(formula, data, family = gaussian(), smoothZ = list(), ziformula = ~0, 
#'         dispformula = ~1, weights = NULL, offset = NULL, contrasts = NULL, na.action, 
#' 		se = TRUE, verbose = FALSE, doFit = TRUE, control = glmmTMBControl(), 
#' 		REML = TRUE, start = NULL, map = NULL, sparseX = NULL, prt=TRUE, 
#' 		predict_info=TRUE)
#' @param formula A two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors.
#' @param data data frame (tibbles are OK) containing model variables. Not required, but strongly recommended; if \code{data} is not specified, downstream methods such as prediction with new data (\code{predict(fitted_model, newdata = ...)}) will fail. If it is necessary to call \code{glmmTMB} with model variables taken from the environment rather than from a data frame, specifying \code{data=NULL} will suppress the warning message.
#' @param family a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. See \code{\link{family}} for a generic discussion of families or \code{\link{family_glmmTMB}} for details of \code{glmmTMB}-specific families.
#' @param smoothZ The older form of the smooths (see \sQuote{Smooth terms in the formula} for the newer one): a list includes a set of smooth Z matrixs (called 'smooth term') used in the mixed effects model, the name of 'smooth term' should be different any variables in the model, each 'smooth term' is the result of function \code{smZ}. e.g.
#' smoothZ=list(sm1=smZ(x1), sm2=smZ(x2, by=f1), sm3=smZ(x3, by=f2, group=TRUE), ...) where 'sm1' to 'sm3' should be new variable names in the \code{data}, and x1 to x3 are covariates, and f1, f2 are factors.
#' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.
#' @param dispformula a \emph{one-sided} formula for dispersion containing only fixed effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value), forcing variance into the random effects. The precise value can be controlled via \code{control=glmmTMBControl(zero_dispval=...)}; the default value is \code{sqrt(.Machine$double.eps)}.
#' @param weights weights, as in \code{glm}. Not automatically scaled to have sum 1.
#' @param offset offset for conditional model (only).
#' @param contrasts an optional list, e.g., \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
#' @param na.action a function that specifies how to handle observations
#' containing \code{NA}s.  The default action (\code{na.omit},
#' inherited from the 'factory fresh' value of
#' \code{getOption("na.action")}) strips any observations with any
#' missing values in any variables. Using \code{na.action = na.exclude}
#' will similarly drop observations with missing values while fitting the model,
#' but will fill in \code{NA} values for the predicted and residual
#' values for cases that were excluded during the fitting process
#' because of missingness.
#' @param se whether to return standard errors.
#' @param verbose whether progress indication should be printed to the console.
#' @param doFit whether to fit the full model, or (if FALSE) return the preprocessed data and parameter objects, without fitting the model.
#' @param control control parameters, see \code{\link{glmmTMBControl}}.
#' @param REML whether to use REML estimation rather than maximum likelihood.
#' @param start starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi} (conditional modes for conditional and zero-inflation models); \code{theta}, \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional and z-i models); \code{psi} (extra family parameters, e.g., shape for Tweedie models).
#' @param map a list specifying which parameter values should be fixed to a constant value rather than estimated. \code{map} should be a named list containing factors corresponding to a subset of the internal parameter names (see \code{start} parameter). Distinct factor values are fitted as separate parameter values, \code{NA} values are held fixed: e.g., \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three fixed-effect parameters of the conditional model and fix the fourth parameter to its starting value. In general, users will probably want to use \code{start} to specify non-default starting values for fixed parameters. See \code{\link[TMB]{MakeADFun}} for more details.
#' @param sparseX a named logical vector containing (possibly) elements named "cond", "zi", "disp" to indicate whether fixed-effect model matrices for particular model components should be generated as sparse matrices, e.g. \code{c(cond=TRUE)}. Default is all \code{FALSE}
#' @param prt Logical scalar - Should the info to be print on screen in the middle of the process or not?
#' @param predict_info Logical scalar - Should provide the info for function semipred or not? In case of there is a correlation theta parameter appearing, you may set predict=FALSE.
#' @return
#' \item{semer}{A glmmTMB model used in the fitting.}
#' \item{data}{A data.frame with generated variables in the fitting.}
#' \item{fomul_vars}{Name of variables in the formula of semireg_tmb model.}
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
#' \item{tmbf}{The glmmTMB model result using doFit=FALSE.}
#' @references
#' Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with O'Sullivan penalized splines. \emph{Australian and New Zealand Journal of Statistics.} \bold{50}, 179-198.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(WarsawApts, package = "HRW")
#'   fit1 <- semireg_tmb(areaPerMzloty ~ construction.date,
#'                       smoothZ=list(
#'                         grp=smZ(construction.date, k=25)
#'                       ),
#'                       data = WarsawApts)
#'   sp_out1 <- semipred(fit1, "construction.date", "construction.date")
#'
#'   WarsawApts$district <- factor(WarsawApts$district)
#'   fit2 <- semireg_tmb(areaPerMzloty ~ construction.date*district,
#'                       smoothZ=list(group=smZ(construction.date, k=15,
#'                                              by = district, group=TRUE)),
#'                       data=WarsawApts)
#'   sp_out2_1 <- semipred(fit2, "district", "construction.date")
#'   sp_out2_2 <- semipred(fit2, "district", "construction.date", contr=c(2,1))
#' }
#' if (requireNamespace("HRW", quietly = TRUE)) {
#'   data(indonRespir, package = "HRW")
#'   fit3 <- semireg_tmb(respirInfec ~ age+vitAdefic + female + height
#'                       + stunted + visit2 + visit3 + visit4  + visit5 + visit6+(1|idnum),
#'                       smoothZ=list(
#'                         grp=smZ(age)
#'                       ),
#'                       family = binomial,
#'                       data = indonRespir)
#'   sp_out3 <- semipred(fit3, "age", "age")
#'   library(ggplot2)
#'   print(sp_out3$plt +
#'     geom_rug(data = subset(indonRespir, respirInfec==0), sides = "b", col="deeppink") +
#'     geom_rug(data = subset(indonRespir, respirInfec==1), sides = "t", col="deeppink") +
#'     coord_cartesian(ylim = c(0, 0.2)))
#' }
#' @name semireg_tmb
NULL

semireg_tmb <- function(formula, data, family = gaussian(), smoothZ = list(), ziformula = ~0, dispformula = ~1, 
                        weights = NULL, offset = NULL, contrasts = NULL, na.action, se = TRUE, 
                        verbose = FALSE, doFit = TRUE, control = glmmTMBControl(), REML = TRUE, 
                        start = NULL, map = NULL, sparseX = NULL, prt=TRUE, predict_info=TRUE)
{
  mc <- match.call() 
  env <- parent.frame()
  environment(formula) <- env
  mc$formula <- formula
  # data, smooth designs and the formula with placeholders (R/semireg-engine.R);
  # no centring for glmmTMB fits
  st <- pm_sr_setup(mc, data, substitute(smoothZ), env)

  glmmTMBcall <- mc 
  glmmTMBcall[[1]] <- quote(glmmTMB::glmmTMB)
  glmmTMBcall <- glmmTMBcall[!names(glmmTMBcall) %in% c("smoothZ", "prt", "predict_info")]
  if (!length(st$Zt)) return(eval.parent(glmmTMBcall))

  # the structure without fitting, the smooth designs in Z, one fit
  fit <- pm_sr_fit_tmb(glmmTMBcall, st, env)
  semer <- fit$semer
  tmbf <- fit$tmbf
  if (prt) print(semer)  
  
  info <- NULL
  if (predict_info) {
    re <- pm_re_terms(semer)
    Z <- getME(semer, "Z")
    # random effects in the column order of Z (getME(, "b") also holds the
    # zero-inflation ones, after these)
    uHat <- getME(semer, "b")[seq_len(ncol(Z))]
    if (all(vapply(re, function(tm) length(tm$comps) == 1, logical(1)))) {
      names(uHat) <- names(unlist(ranef(semer)$cond))   # same order then
    }
    X <- getME(semer, "X")
    family_n <- semer$modelInfo$family$family
    w <- if (family_n == "gaussian") tmbf$data.tmb$weights else {
      # the information weights 1/Var(working response | b) (pm_working_tmb(),
      # checked against simulate()), times sigma^2: pm_sr_info() works
      # relative to sigma^2 and semipred() multiplies by it, and for families
      # whose sigma() is not a residual SD (negative binomial: theta; beta:
      # phi) the weights alone gave covariances off by that factor
      wk <- tryCatch(pm_working_tmb(semer), error = function(e) NULL)
      w_info <- if (!is.null(wk)) 1/wk$resvar else {
        fam <- semer$modelInfo$family
        eta <- as.vector(cbind(X, Z) %*% c(fixef(semer)$cond, uHat))
        w0 <- as.vector(fam$mu.eta(eta)^2/fam$variance(fam$linkinv(eta))) * tmbf$data.tmb$weights
        if (length(tmbf$data.tmb$size) > 0) w0 <- w0 * tmbf$data.tmb$size
        w0
      }
      stats::sigma(semer)^2 * w_info
    }
    u_ind <- vapply(re, function(tm) length(tm$cols), numeric(1))
    names(u_ind) <- vapply(re, `[[`, "", "group")
    info <- pm_sr_info(X, Z, pm_tmb_lambda(semer, re), w, uHat, u_ind, names(st$Zt), st$smoothZm)
  }
  pm_sr_result(semer, st, mc, info, list(tmbf = tmbf))
}
