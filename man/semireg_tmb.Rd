\name{semireg_tmb}
\alias{semireg_tmb}
\title{
  Fitting Semi Parametric Models Using glmmTMB
}
\description{
   Fit a semi parametric model based on glmmTMB.
}
\usage{
semireg_tmb(formula, data, family = gaussian(), smoothZ = list(), ziformula = ~0, 
        dispformula = ~1, weights = NULL, offset = NULL, contrasts = NULL, na.action, 
		se = TRUE, verbose = FALSE, doFit = TRUE, control = glmmTMBControl(), 
		REML = TRUE, start = NULL, map = NULL, sparseX = NULL, prt=TRUE, 
		predict_info=TRUE)
}

\arguments{	
\item{formula}{
A two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors.
}	

\item{data}{data frame (tibbles are OK) containing model variables. Not required, but strongly recommended; if \code{data} is not specified, downstream methods such as prediction with new data (\code{predict(fitted_model, newdata = ...)}) will fail. If it is necessary to call \code{glmmTMB} with model variables taken from the environment rather than from a data frame, specifying \code{data=NULL} will suppress the warning message.}

\item{family}{a family function, a character string naming a family function, or the result of a call to a family function (variance/link function) information. See \code{\link{family}} for a generic discussion of families or \code{\link{family_glmmTMB}} for details of \code{glmmTMB}-specific families.}

\item{smoothZ}{
A list includes a set of smooth Z matrixs (called 'smooth term') used in the mixed effects model, the name of 'smooth term' should be different any variables in the model, each 'smooth term' is the result of function \code{smZ}. e.g.
smoothZ=list(sm1=smZ(x1), sm2=smZ(x2, by=f1), sm3=smZ(x3, by=f2, group=TRUE), ...) where 'sm1' to 'sm3' should be new variable names in the \code{data}, and x1 to x3 are covariates, and f1, f2 are factors.
}

\item{ziformula}{a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.}

\item{dispformula}{a \emph{one-sided} formula for dispersion containing only fixed effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value), forcing variance into the random effects. The precise value can be controlled via \code{control=glmmTMBControl(zero_dispval=...)}; the default value is \code{sqrt(.Machine$double.eps)}.}

\item{weights}{weights, as in \code{glm}. Not automatically scaled to have sum 1.}

\item{offset}{offset for conditional model (only).}

\item{contrasts}{an optional list, e.g., \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.}

\item{na.action}{a function that specifies how to handle observations
containing \code{NA}s.  The default action (\code{na.omit},
inherited from the 'factory fresh' value of
\code{getOption("na.action")}) strips any observations with any
missing values in any variables. Using \code{na.action = na.exclude}
will similarly drop observations with missing values while fitting the model,
but will fill in \code{NA} values for the predicted and residual
values for cases that were excluded during the fitting process
because of missingness.}

\item{se}{whether to return standard errors.}

\item{verbose}{whether progress indication should be printed to the console.}

\item{doFit}{whether to fit the full model, or (if FALSE) return the preprocessed data and parameter objects, without fitting the model.}

\item{control}{control parameters, see \code{\link{glmmTMBControl}}.}

\item{REML}{whether to use REML estimation rather than maximum likelihood.}

\item{start}{starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi} (conditional modes for conditional and zero-inflation models); \code{theta}, \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional and z-i models); \code{psi} (extra family parameters, e.g., shape for Tweedie models).}

\item{map}{a list specifying which parameter values should be fixed to a constant value rather than estimated. \code{map} should be a named list containing factors corresponding to a subset of the internal parameter names (see \code{start} parameter). Distinct factor values are fitted as separate parameter values, \code{NA} values are held fixed: e.g., \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three fixed-effect parameters of the conditional model and fix the fourth parameter to its starting value. In general, users will probably want to use \code{start} to specify non-default starting values for fixed parameters. See \code{\link[TMB]{MakeADFun}} for more details.}

\item{sparseX}{a named logical vector containing (possibly) elements named "cond", "zi", "disp" to indicate whether fixed-effect model matrices for particular model components should be generated as sparse matrices, e.g. \code{c(cond=TRUE)}. Default is all \code{FALSE}}

\item{prt}{
Logical scalar - Should the info to be print on screen in the middle of the process or not?
}

\item{predict_info}{
Logical scalar - Should provide the info for function semipred or not? In case of there is a correlation theta parameter appearing, you may set predict=FALSE.  
}
}

\details{
A semi parametric model can be parameterized as a linear (or generalized linear) mixed model in which its random effects are smooth functions of some covariates (named ‘smooth term’). \code{semireg_tmb} follows the approach suggested by Wand and Ormerod (2008) and represents the 'smooth term' using O'Sullivan-type of Z. 
}

\value{
\item{semer}{A glmmTMB model used in the fitting.}
\item{data}{A data.frame with generated variables in the fitting.}
\item{fomul_vars}{Name of variables in the formula of semireg_tmb model.}
\item{sm_vars}{Name of variables in the smoothZ list.}
\item{smoothZ_call}{A call used to produce smooth terms in the fitting.}
\item{knots_lst}{Knots used in each smooth term in the fitting.}
\item{range_lst}{Range of covariate used in each smooth term in the fitting.}
\item{cov_lst}{Covariance matrix list for each smooth term.}
\item{u_lst}{Random effects list for each smooth term.} 
\item{type_lst}{Smooth type list of smooth terms.} 
\item{CovMat}{Covariance matrix for all smooth terms.}
\item{Cov_ind}{Covariance matrix index for each smooth term.}
\item{Cov_indN}{Covariance matrix index for each smooth term when \code{group=TRUE} in \code{smoothZ} argument.}
\item{df}{Degree of freedom of all random terms.}
\item{tmbf}{The glmmTMB model result using doFit=FALSE.}
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\references{
Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with O'Sullivan penalized splines. \emph{Australian and New Zealand Journal of Statistics.} \bold{50}, 179-198. 
}

\examples{
## Not run
# library(predictmeans)
# library(HRW) 
# data(WarsawApts)
# help(WarsawApts)
# str(WarsawApts)
# fit1 <- semireg_tmb(areaPerMzloty ~ construction.date,
#                smoothZ=list(
#                  grp=smZ(construction.date, k=25)
#                ),
#                data = WarsawApts)
# sp_out1 <- semipred(fit1, "construction.date", "construction.date")
# 
# WarsawApts$district <- factor(WarsawApts$district)
# fit2 <- semireg_tmb(areaPerMzloty ~ construction.date*district, resp_scale = TRUE,
#                 smoothZ=list(group=smZ(construction.date, k=15,
#                                        by = district, group=TRUE)), 
#                 data=WarsawApts)
# sp_out2_1 <- semipred(fit2, "district", "construction.date")
# sp_out2_2 <- semipred(fit2, "district", "construction.date", contr=c(2,1))
# 
# data(indonRespir)
# help(indonRespir)
# str(indonRespir)
# fit3 <- semireg_tmb(respirInfec ~ age+vitAdefic + female + height
#                + stunted + visit2 + visit3 + visit4  + visit5 + visit6+(1|idnum),
#                smoothZ=list(
#                  grp=smZ(age)
#                ),
#                family = binomial,
#                data = indonRespir)
# sp_out3 <- semipred(fit3, "age", "age")
# library(ggplot2)
# sp_out3$plt+
#   geom_rug(data = subset(indonRespir, respirInfec==0), sides = "b", col="deeppink") +
#   geom_rug(data = subset(indonRespir, respirInfec==1), sides = "t", col="deeppink")+
#   ylim(0, 0.2)                     
}