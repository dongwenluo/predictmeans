\name{ci_mcp}
\alias{ci_mcp}
\title{
Multiple Comparisons Based on the Confidence Intervals
}
\description{
This function produces letter representations for a multiple comparison test by analyzing the confidence intervals associated with the mean values of different treatments. In particular, if the confidence intervals of two treatments overlap, it indicates that there is no significant difference between them. Conversely, if the confidence intervals do not overlap, it indicates that the treatments are significantly different from each other.
}
\usage{
ci_mcp(LL, UL, trt_n=NULL)
}
\arguments{
  \item{LL}{
  Lower limits of treatments' confidence interval.
}
  \item{UL}{
  Upper limits of treatments' confidence interval. 
}
  \item{trt_n}{
  Treatments' names. 
}
}

\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\references{
	Vanessa, C. (05 October 2022),
	\emph{Confidence tricks: the 83.4\% confidence interval for comparing means}, https://vsni.co.uk/blogs/confidence_trick. 
}

\examples{
  library(predictmeans)
  ci_mcp(LL=c(68.2566,  87.7566, 103.0899, 112.2566), UL=c(90.5212, 110.0212, 125.3545, 134.5212))
  
  data("Oats", package="nlme")
  Oats$nitro <- factor(Oats$nitro)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  predictmeans(fm, "nitro", adj="BH", plot=FALSE)$mean_table
  predictmeans(fm, "nitro", pair=TRUE, level=0.166, letterCI = TRUE, plot=FALSE)$mean_table
}