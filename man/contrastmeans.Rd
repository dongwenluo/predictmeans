\name{contrastmeans}
\alias{contrastmeans}
\title{Linear Contrast Tests for a Linear Model}
\description{
  Performs t-tests (or permuted t-tests) of specified contrasts for linear models 
  obtained from functions \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, 
  or \code{lmer}. 
}
\usage{
  contrastmeans(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", Df, permlist)
}
\arguments{
  \item{model}{
  Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}. 
}
  \item{modelterm}{
  Name (in "quotes") for indicating which factor term's contrast to be calculated. 
  The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B". 
}
  \item{ctrmatrix}{
  A specified contrast matrix. If \code{ctrmatrix} is missing, the programe will ask user to enter it. 
}
  \item{ctrnames}{ Names of the specified contrasts, e.g. c("A vs D", "C vs B", ...) 
}
  \item{adj}{   Name (in "quote") for indicating a method for adjusting p-values of pairwise comparisons. 
  The choices are "none", "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr". 
  The default method is "none".
}  
  \item{Df}{ A denominator degree of freedom for \code{modelterm}. (For \code{glmer} models 
    the \code{Df} needs to be specified, while for the other models, \code{Df} is obtained from the fitted model automatically). 
}

  \item{permlist}{ A model parameter list containing \code{nsim} parameters produced by the function \code{permmodels}. When \code{permlist != NULL}, the option \code{Df} will be non-functional. This is a key 
	option for the permutation test. 
}
}
\value{
There are two components in the output which are
\item{Table}{A table showing t-test results for the specified linear contrasts.}
\item{K}{A contrast matrix.}
}
\references{
    Torsten Hothorn, Frank Bretz and Peter Westfall (2008),
    \emph{Simultaneous Inference in General Parametric Models. Biometrical}, Journal 50(3), 346--363.
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\examples{
library(predictmeans)
# ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## Not run: 
## The contrast has a contrast matrix as follows:
#     0:Golden Rain 0:Marvellous 0:Victory 
#[1,]            -1            0         1 
#[2,]             0            0         1 
#     0.2:Golden Rain 0.2:Marvellous 0.2:Victory 
#[1,]               0              0           0 
#[2,]               0              0           0 
#     0.4:Golden Rain  0.4:Marvellous 0.4:Victory
#[1,]               0               0           0
#[2,]               0              -1           0
#      0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#[1,]                0              0           0
#[2,]                0              0           0

# 1. Enter above contrast matrix into a pop up window, then close the window
# contrastmeans(fm, "nitro:Variety")
 
# 2. Construct the contrast matrix directly
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)
}
