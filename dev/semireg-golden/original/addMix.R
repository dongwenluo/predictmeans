########## R script: addMix ##########

# For illustration of additive mixed models
# in a simple binary offset situation,
# using simulated data.
#
# Fitting is performed via the lme() function
# within the linear mixed model framework.
#
# Variance calculations are performed using the
# streamlined methodology  described in the paper:
# Smith, A.D.A.C. & Wand, M.P. (2008).
# Streamlined variance calculations for semiparametric
# mixed models. Statistics in Medicine, 27, 435-448.

# Last changed: 20 MAY 2009

# Load required packages:

library(nlme) ; library(lattice)

# Load function for O'Sullivan spline basis computation:

source("ZOSull.r")

# Generate synthetic data according to (5) of Smith & Wand (2008).

set.seed(39402) ; m <- 250 ; nVals <- sample(1:4,m,replace=TRUE)
betaVal <- 0.3 ; sigU <- 0.5 ; sigEps <- 0.2
f <- function(x) return(-sin(2*pi*x))
U <- rnorm(m,0,sigU)
sVals <- NULL ; xVals <- NULL
Uvals <- NULL ; idNum <- NULL
for (i in 1:m)
{
   idNum <- c(idNum,rep(i,nVals[i]))
   stt <- runif(1,0,1-0.05*(nVals[i]-1))
   sVals <- c(sVals,seq(stt,by=0.05,length=nVals[i]))
   xCurr <- sample(c(0,1),1)
   xVals <- c(xVals,rep(xCurr,nVals[i]))
   Uvals <- c(Uvals,rep(U[i],nVals[i]))
}
epsVals <- rnorm(sum(nVals),0,sigEps)
yVals <- f(sVals) + betaVal*xVals + Uvals + epsVals

# Fit an additive mixed model:

numIntKnots <- 15 ; p <- 1
intKnots <- quantile(unique(sVals),seq(0,1,length=
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.s <- c(0,1)
X <- cbind(rep(1,length(yVals)),sVals,xVals)
ZSpline <- ZOSull(sVals,range.s,intKnots)
numKnots <- ncol(ZSpline)
#========================================================================== 
groupVec <- factor(rep(1,length(yVals)))
ZBlock <- list(groupVec=pdIdent(~ZSpline-1),idNum=pdIdent(~1))
dataFr_lme <- groupedData(yVals~xVals|groupVec,
                      data=data.frame(yVals,X,ZSpline,idNum))
fit <- lme(yVals~-1+X,data=dataFr_lme,random=ZBlock)
betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)
uSplineHat <- uHat[1:ncol(ZSpline)]

# Perform streamlined variance estimation for variability
# bars and standard error of binary predictor:

sig.eps.hat <- fit$sigma
sig.u.hat <- intervals(fit)$reStruct$groupVec$est
sig.U.hat <- intervals(fit)$reStruct$idNum$est
M11inv <- diag(c(rep(0,p+2),rep(sig.eps.hat^2/sig.u.hat^2,numKnots)))
M11inv <- M11inv + crossprod(cbind(X,ZSpline))
M.21 <- rowsum(cbind(X,ZSpline),idNum)
for (i in 1:m)
{
   h.i <- M.21[i,] ; n.i <- M.21[i,1]
   M11inv <- (M11inv - sig.U.hat^2*tcrossprod(h.i)
                     /(sig.U.hat^2*n.i+sig.eps.hat^2))
}
covMat <- (sig.eps.hat^2)*solve(M11inv)

# Report results:

cat("\n Estimated effect of x=1 and standard error is:\n\n")
outp <- c(betaHat[3],sqrt(covMat[3,3]))
names(outp) <- c("est.","s.e.")
print(signif(outp,3))  ; cat("\n\n")
ng <- 101 ; sValsg <- seq(0,1,length=ng)
ZgSpline <- ZOSull(sValsg,range.s,intKnots)
plotMatrix0 <- cbind(rep(1,ng),sValsg,rep(0,ng),ZgSpline)  
fhatgREML <- plotMatrix0 %*% c(betaHat, uSplineHat)
Var0 <- diag(plotMatrix0%*%covMat%*%t(plotMatrix0))
Var1 <- 2*plotMatrix0%*%covMat[,3] + covMat[3,3]
xLabs <- paste("x =",as.character(xVals))
pobj <- xyplot(yVals~sVals|xLabs,groups=idNum,xlab="value of s",
        ylab="value of y",subscripts=TRUE,
        panel=function(x,y,subscripts,groups)
        {
           panel.grid()  ; panel.superpose(x,y,subscripts,groups,
                                   type="b",pch=16)
           panelInd <- any(xVals[subscripts]==1)
           panel.xyplot(sValsg,fhatgREML+panelInd*betaHat[3],
                        lwd=3,type="l",col="navy")
           for (multFac in c(-2,2))
              panel.xyplot(sValsg,fhatgREML+panelInd*betaHat[3]+multFac
                *sqrt(Var0+panelInd*Var1),lty=2,lwd=2,type="l",col="navy")
        })
print(pobj)

########## End of addMix ##########
