\name{permindex}
\alias{permindex}
\title{
Permutation Index
}
\description{
This function obtains permutation index for a dataset.
}
\usage{
permindex(data, block=NULL, group=NULL, nsim=4999, seed)
}
\arguments{
  \item{data}{
Data object used in the \code{model} fitting.
}
  \item{block}{
Name (in "quotes") for the blocking factor in the \code{data}.
}
  \item{group}{
Name (in "quotes") for the group factor in the \code{data}.
}
  \item{nsim}{
The number of permutations. The default is 4999.  
}
  \item{seed}{
Specify a random number generator seed, for reproducible results.
}
}
\value{
  A matrix has 'nsim' columns of permuted index.
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}
\examples{
  library(predictmeans)
  block <- rep(1:3, each=12)
  group <- rep(rep(1:3, each=4), 3)
  data <- data.frame(block, group)
  cbind(data, permindex(data, block="block", group="group", nsim=5))  
                        # Permute group as a whole within each block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, block="block",  nsim=5)) 
                        # Permute obs within each block only.
  cbind(data, permindex(data, group="group", nsim=5)) 
                        # Permute groups as a whole block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, nsim=5))  # Free permutation.
}
