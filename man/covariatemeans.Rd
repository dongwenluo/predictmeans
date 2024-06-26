\name{covariatemeans}
\alias{covariatemeans}
\title{
  Predicted Means of a Linear Model with Covariate Variable(s)
}
\description{
  This function obtains predicted means with graph for a new set of covariate values.
}
\usage{
  covariatemeans(model, modelterm=NULL, covariate, as.is=FALSE, covariateV=NULL, 
  data=NULL, level=0.05, Df=NULL, trans=NULL, transOff=0, responsen=NULL, trellis=TRUE, 
  plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=FALSE)
}
\arguments{
  \item{model}{
  Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}. 
}
  \item{modelterm}{
  Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
  The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B". 
}
  \item{covariate}{
  Name (in "quotes") of one the covariate variables in the \code{model}.
}
  \item{as.is}{
A logic value to specify wheather or not using original covariate values' rage for graph, the default is FALSE.
}
  \item{covariateV}{
  A numeric vector when as.is is FALSE, then covariatemeans will produce the result for \code{covariate} at value of \code{covariateV}.
}
  \item{data}{
In some cases, you need to provide the data set used in model fitting, especially when you have applied some variable trnasformation in the model.}
  \item{level}{
  A significant level for calculating confident interval. The default value is 0.05.
}

\item{Df}{
A degree of freedom for calculating CI of predicted means (you can manually specified ddf here). For the above models, ddf is obtained from the function automatically.
}
  \item{trans}{
A function object for calculating the back transformed means, e.g. \code{trans=exp}.
}
  \item{transOff}{
When you use \code{trans=exp(x+1)}, then \code{transOff=1}, the default is 0.
}
  \item{responsen}{
  Name (in "quotes") of the back transformed response variable in the \code{model}.
}
  \item{trellis}{
  A logical variable. If set to TRUE (default), a trellis plots of predicted means with CI will be drawn.
}
 \item{plotord}{
  A numeric (or character) vector specifying the order of plotting for two or three way interaction (e.g. 
  \code{plotord = c(2, 1, 3)} or \code{plotord = c("B", "A", "C")} will put the second variable 'B' in 
  \code{modelterm} on the \code{X} axis, the first variable 'A' as the grouping variable, and the third one 
  'C' as the panel variable). The defaults are \code{c(1, 2)} and \code{c(1, 2, 3)} for two and three way 
  interactions. 
}
  \item{mtitle}{
  The main title in the graph.
}
\item{ci}{
A logical variable to indicate whether to print confidence interval. The default value is TRUE.
} 
\item{point}{
A logical variable to indicate whether to print raw data points. The default value is TRUE.
} 
  \item{jitterv}{
  A degree of jitter in x and y direction in the graph. The default is zero.
}
\item{newwd}{
A logical variable to indicate whether to print graph in a new window. The default value is FALSE.
} 
}

\value{
 \item{Predicted Means}{A table of predicted means.}
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
  library(predictmeans)
  data(Oats, package="nlme")
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  covariatemeans(fm, "Variety", covariate="nitro")
  covariatemeans(fm, "Variety", covariate="nitro", covariateV=seq(0, 0.6, 0.1))$pltdf
}
