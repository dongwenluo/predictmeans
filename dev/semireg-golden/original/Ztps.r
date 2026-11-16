########## R function: Ztps ##########

# For creation of thin plate spline-type Z matrices.

# Last changed: 26 FEB 2008

Ztps <- function(x,knots)
{
   # Set up thin plate spline generalised
   # covariance function:
  
   tps.cov <- function(r,m=2,d=1)
   {     
       r <- as.matrix(r)
       num.row <- nrow(r)
       num.col <- ncol(r)
       r <- as.vector(r)
       nzi <- (1:length(r))[r!=0]
       ans <- rep(0,length(r))
       if ((d+1)%%2!=0)    
          ans[nzi] <- (abs(r[nzi]))^(2*m-d)*log(abs(r[nzi])) # d is even
       else
          ans[nzi] <- (abs(r[nzi]))^(2*m-d)

       if (num.col>1) ans <- matrix(ans,num.row,num.col)     # d is odd
       return(ans)
   }

   # Set up function for matrix square-roots:

   matrix.sqrt <- function(A)
   {
      sva <- svd(A)
      if (min(sva$d)>=0)
         Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
      else
         stop("Matrix square root is not defined")
      return(Asqrt)
   }
   
   # Obtain  matrix of inter-knot distances:

   numKnots <- nrow(knots)
   
   dist.mat <- matrix(0,numKnots,numKnots)
   dist.mat[lower.tri(dist.mat)] <- dist(as.matrix(knots))
   dist.mat <- dist.mat + t(dist.mat)

   Omega <- tps.cov(dist.mat,d=2)

   # Obtain preliminary Z matrix of knot to data covariances:

   x.knot.diffs.1 <- outer(x[,1],knots[,1],"-")
   x.knot.diffs.2 <- outer(x[,2],knots[,2],"-")
   x.knot.dists <- sqrt(x.knot.diffs.1^2+x.knot.diffs.2^2)

   prelim.Z <- tps.cov(x.knot.dists,m=2,d=2)

   # Transform to canonical form: 
  
   sqrt.Omega <- matrix.sqrt(Omega)
   Z <- t(solve(sqrt.Omega,t(prelim.Z)))
   
   return(Z)
}

########## End of Ztps ##########

