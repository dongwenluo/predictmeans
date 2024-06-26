\name{Kmatrix}
\alias{Kmatrix}
\title{
Matrix of Coefficients in a Linear Model
}
\description{
  This function obtains a matrix of coefficients for parametric models such as \code{aov}, \code{lm}, 
  \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
}
\usage{
Kmatrix(model, modelterm, covariate=NULL, covariateV=NULL, data=NULL, prtnum=FALSE)
}
\arguments{
  \item{model}{
Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
}
  \item{modelterm}{
  Name (in "quotes") for indicating which model term's predicted mean to be calculated. 
  The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B". 
}
  \item{covariate}{
A numerical vector to specify values of covariates for calculating predicted means, 
default values are the means of the associated covariates. It also can be the name of one covariate in the model.
}
  \item{covariateV}{
  A numeric vector or list of numeric vector, then covariatemeans will produce the result for \code{covariate} at value of \code{covariateV}.
}
  \item{data}{
In some cases, you need to provide the data set used in model fitting, especially when you have applied some variable trnasformation in the model.}
  \item{prtnum}{
An option for printing covariate info on the screen or not. The default is FALSE.
}
}

\value{
 \item{K}{Coefficients matrix}
 \item{fctnames}{A model frame contains factor(s) info in the model.}
 \item{response}{The name of response variable in the model.}
}
\references{
  Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347. 
}
\author{
This function heavily depends on the codes from package "lsmeans".
}

\examples{
  library(predictmeans)
  data(Oats, package="nlme")
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  Kmatrix(fm, "Variety", prtnum=TRUE)$K
  Kmatrix(fm, "Variety", 0.5, prtnum=TRUE)$K
 # Kmatrix(fm, "Variety", "nitro")$K
  Kmatrix(fm, "Variety", "nitro", covariateV=seq(0, 0.6, 0.1))$K
}

