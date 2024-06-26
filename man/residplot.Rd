\name{residplot}
\alias{residplot}
\title{
Diagnostic Plots for a Linear (Mixed) Model
}
\description{
This function produces diagnostic plots for linear models including 'aov', 'lm', 'glm', 'gls', 'lme' and 'lmer'. 
}
\usage{
  residplot(model, group = "none", level = 1, slope = FALSE, id = FALSE, newwd=FALSE,
    ask=FALSE)
}
\arguments{
  \item{model}{
Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}. 
}
  \item{group}{
Name (in "quotes") for indicating the variable used to show grouping in the residual vs predicted plot. If variable is a term in the model, then group will be a name of the variable such as \code{group="A"}, otherwise group will be the actual variable such as \code{group=data$A}.
}
  \item{level}{
An integer 1, 2, etc. used to specify a level of the random effect for plotting. The default value is 1.
}
  \item{slope}{
  A logical variable. If set to TRUE, a Q-Q plot of random slope will be drawn.
}
  \item{id}{
  A logical variable. If set to TRUE, outliers in the residual vs fitted plot can be identified interactively.
}
\item{newwd}{
A logical variable to indicate whether to print graph in a new window. The default is FALSE.
}
  \item{ask}{
  logical. If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn.
}
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}
\examples{
## Note that the order of levels of nested random effects is oposite 
## between lme and lmer objects.

library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
residplot(fm, level=2)    #lme: level=2 for random effect "Block:Variety"

#  Not Run
#  library(lme4)
#  fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#  residplot(fm) # lmer: By default level=1 for random effect "Block:Variety"
}
