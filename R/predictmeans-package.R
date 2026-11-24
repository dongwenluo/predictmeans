#' Predicted Means for Linear and Semiparametric Models
#'
#' @description
#' \if{html}{\figure{logo.png}{options: style='float: right' alt='predictmeans logo' width='120'}}
#'
#' This package provides functions to diagnose and make inferences from linear, generalized linear, mixed and semiparametric models, such as those obtained from 'aov', 'lm', 'glm', 'gls', 'lme', 'lmer', 'glmer', 'glmmTMB', 'gam' and 'semireg'. Inferences include predicted means and standard errors, contrasts, multiple comparisons, permutation tests, adjusted R-square and graphs (\code{autoplot()}). Vignettes (\code{browseVignettes("predictmeans")}) show worked examples. For AI assistants and scripted use, \code{system.file("llms-context.md", package = "predictmeans")} gives a compact description of the result tables, defaults and common recipes.
#' @section Workflow:
#' \if{html}{\figure{workflow.png}{options: width='720' alt='From a fitted model to the report: diagnostics, tests, means and contrasts, variance components, then tables and plots'}}
#' @section What predictmeans aims for:
#' \itemize{
#'   \item Defaults chosen on evidence (simulation, literature): Kenward-Roger df and adjusted SEs for REML
#'     \code{lmer} fits; permutation p-values \eqn{(b + 1)/(B + 1)}, never 0; Genstat-style stratum SEs for
#'     \code{aov} models with \code{Error()} strata.
#'   \item Checked against outside references: SAS PROC MIXED and Genstat runs (\code{dev/sas},
#'     \code{dev/genstat} in the package sources), \pkg{pbkrtest}, \pkg{nlme} and \pkg{emmeans}, in tests with
#'     stated tolerances.
#'   \item Works across \code{lm}, \code{glm}, \code{aov}, \code{lmer}, \code{glmer}, \code{lme}, \code{gls},
#'     \code{glmmTMB}, \code{gam} and \code{semireg} fits; when a result cannot be computed correctly it stops with
#'     a clear message.
#'   \item Advice triggered by this fit and data, read with \code{\link{advice}}.
#'   \item Results are the same for any \code{ncore}.
#' }
#' @seealso \code{\link{predictmeans}}, \code{\link{autoplot.predictmeans}}, \code{\link{tidy.predictmeans}},
#'   \code{\link{advice}}, \code{\link{pm_check}}
#' @details
#' \tabular{ll}{
#' Package: \tab predictmeans\cr
#' Type: \tab Package\cr
#' License: \tab GPL (>= 2)\cr
#' }
#' @references
#' Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
#' 	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347.
#' @author
#' Dongwen Luo, Siva Ganesh, John Koolaard and James Curran
#'
#' Maintainer: Dongwen Luo <dongwen.luo@agresearch.co.nz>
#' @name predictmeans-package
#' @keywords package
#' @import parallel lme4 glmmTMB ggplot2
#' @importMethodsFrom Matrix * solve t which diag crossprod tcrossprod chol
#' @importFrom numDeriv hessian jacobian
#' @importFrom MASS theta.ml negative.binomial
#' @importFrom splines ns bs
#' @importFrom nlme getData intervals getResponse
#' @importFrom grDevices dev.new
#' @importFrom Matrix .bdiag Diagonal forceSymmetric Matrix nearPD
#' @importFrom graphics abline axis box hist identify image legend lines par plot
#'   points text title
#' @importFrom methods is slot as
#' @importFrom stats anova as.formula aov coef complete.cases confint
#'   delete.response density dist family fitted formula ftable gaussian lm
#'   na.action lm.influence loess.smooth logLik make.link median model.response
#'   model.frame model.matrix na.omit na.pass p.adjust pf poisson ppoints predict
#'   pt ptukey qnorm qqline qqnorm qt quantile resid residuals sd setNames sigma
#'   terms update vcov xtabs aggregate get_all_vars
#' @importFrom utils as.roman combn edit
"_PACKAGE"
