########## R script: wageTUnionLMM.Rs ##########

# For illustrating penalised spline logistic 
# regression using the generalised linear mixed 
# model approach.
#
# The R function lme() for linear mixed model 
# fitting is used.

# Last changed: 12 AUG 2008 by M.P. Wand.

# Load required packages:

library(nlme) ; library(MASS)

# Load function for O'Sullivan spline basis computation:

source("ZOSull.r")

# Read in and extract data:

wageTUnion <- read.table("wageTUnion.txt",header=T)
x <- wageTUnion[,1]
y <- wageTUnion[,2]

# Contruct X and Z matrix:

X <- cbind(rep(1,length(x)),x)
numIntKnots <- 15
intKnots <- quantile(unique(x),seq(0,1,length=
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x <- c(min(x),max(x))
Z <- ZOSull(x,range.x,intKnots)

# Obtain fit using the linear mixed model function lme() :

group <- rep(1,length(x))
data.fr <- groupedData(y~x|group,data=data.frame(x,y))
fit <- glmmPQL(y~-1+X,random=list(group=pdIdent(~-1+Z)),
               data=data.fr,family=binomial)

# Extract the main parameters:

betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)
sig.u.hat <- fit$sigma*exp(unlist(fit$modelStruct))

# Plot the fit over a fine grid:

ng <- 101
xg <- seq(min(x),max(x),length=ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,range.x,intKnots)
etaHatg <- Xg%*%betaHat + Zg%*%uHat
muHatg <- 1/(1+exp(-etaHatg))
plot(x,y,bty="l",xlab="wage",ylab="est. P(trade union membership)",ylim=c(0,0.5),
     type="n",main="fit using logistic linear mixed model approach and PQL")
lines(xg,muHatg,col="blue",lwd=2)
rug(x[y==0],col="DeepPink")
rug(x[y==1],col="DeepPink",side=3)

# Compute the degrees of freedom and add to the plot:

Cmat <- cbind(X,Z)
etaHat <- as.vector(X%*%betaHat + Z%*%uHat)
wVec <- exp(etaHat)/((1+exp(etaHat))^2)
CTWC <- crossprod(Cmat*wVec,Cmat)
lamVal <- 1/(sig.u.hat^2)
Dmat <- diag(c(0,0,rep(1,ncol(Z))))
dfVal <- sum(diag(solve(CTWC+lamVal*Dmat,CTWC)))
text(10,0.45,adj=0,paste("Fit has",signif(dfVal,3),"degrees of freedom."),col="blue")

# Compute and add variablity bars:

Cg <- cbind(Xg,Zg)
expEtag <- exp(etaHatg)
sdg <- sqrt(diag(Cg%*%solve(CTWC+lamVal*Dmat,t(Cg))))*expEtag/((1+expEtag)^2)
lower <- muHatg - 2*sdg ; upper <- muHatg + 2*sdg
lower[lower<0.001] <- 0.001
polygon(c(xg,rev(xg)),c(lower,rev(upper)),col="cyan",border=FALSE)
lines(xg,muHatg,col="blue",lwd=2)

########## End of wageTUnionLMM ##########
