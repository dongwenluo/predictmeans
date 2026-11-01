\name{varcomp}
\alias{varcomp}
\title{
Calculate SE and CI of variance components for \code{lmer}, \code{glmer}, \code{lme}, \code{glmmTMB} model
}
\description{
This function calculates SE and CI of variance components for \code{lmer}, \code{glmer}, \code{lme}, \code{glmmTMB} model. For \code{lme} models it needs package \pkg{lmeInfo}.
}
\usage{
varcomp(model, ci=TRUE, level=0.95) 
}
\arguments{
  \item{model}{
Model object returned by \code{lmer}, \code{glmer}, \code{lme}, \code{glmmTMB}. 
}
  \item{ci}{
a logical value to indicates wheather or not to simulate a confidence interval for \code{lmer} model, the default value is TRUE.
}
  \item{level}{
level of confidence of CI, the default value is 0.95.
}
}
\value{
Variance components table. 
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
\dontrun{varcomp(fm)}
fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
varcomp(fm1)

data(Orthodont, package="nlme")
mod <- lmer(distance ~ age + (age|Subject), data=Orthodont)
\dontrun{varcomp(mod)}
mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
varcomp(mod1)
}
