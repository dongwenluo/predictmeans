\name{permlmer}
\alias{permlmer}
\title{
Permutation Test of random or fixed effects for \code{lmer} model.  
}
\description{
This function provides permutation tests for the terms in a linear mixed model of \code{lmer}.
}
\usage{
permlmer(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed)
}
\arguments{
  \item{lmer0}{
\code{lmer} model under H0, note that \code{lmer0} model must nest within \code{lmer1} model. 
}
  \item{lmer1}{
\code{lmer} model under H1, note that \code{lmer0} model must nest within \code{lmer1} model. 
}
  \item{nperm}{
Number of permutation, the default value is 999.
}
  \item{ncore}{
Number of core for parallel computing, the default value is 3.
}
  \item{plot}{
Plot permutation distribution or not, the default value is FALSE.
}
  \item{seed}{
Specify a random number generator seed, for reproducible results.
}
}

\value{
Permutation p-value. 
}
\references{
    Oliver E. Lee and Thomas M. Braun (2012),
    \emph{Permutation Tests for Random Effects in Linear Mixed Models. Biometrics}, Journal 68(2).
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
# library(predictmeans)
## Test random effects
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
# fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
# fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
# anova(fm1, fm2, fm3)
# permlmer(fm3, fm2)
# permlmer(fm2, fm1)

## Test fixed effects
# Oats$nitro <- factor(Oats$nitro)
# fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# permlmer(fm0, fm)
}
