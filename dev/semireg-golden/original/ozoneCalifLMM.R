########## R script: ozoneCalifLMM.Rs ##########

# For illustrating penalised spline fitting
# of a generalised additive model using the 
# linear mixed model approach.
#
# The R function lme() for linear mixed model 
# fitting is used.

# Last changed: 05 AUG 2010.

# Load required package and function:

library(nlme) ; source("ZOsull.r")

# Read in and extract data:

ozoneCalif <- read.table("ozoneCalif.txt",header=T)
pairs(ozoneCalif)
cat("Hit Enter to continue.\n")
ans <- readline()

y <- ozoneCalif$ozone
x1 <- ozoneCalif$inversion.base.height
x2 <- ozoneCalif$daggett.pressure.gradient
x3 <- ozoneCalif$inversion.base.temp

# Set up matrices for additive model:

X <- cbind(rep(1,length(y)),x1,x2,x3)

numIntKnots <- 15
intKnots1 <- quantile(unique(x1),seq(0,1,length=
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x1 <- c(min(x1),max(x1))
Z1 <- ZOSull(x1,range.x1,intKnots1)
K1 <- ncol(Z1)

intKnots2 <- quantile(unique(x2),seq(0,1,length=
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x2 <- c(min(x2),max(x2))
Z2 <- ZOSull(x2,range.x2,intKnots2)
K2 <- ncol(Z2)

intKnots3 <- quantile(unique(x3),seq(0,1,length=
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x3 <- c(min(x3),max(x3))
Z3 <- ZOSull(x3,range.x3,intKnots3)
K3 <- ncol(Z3)
Z <- cbind(Z1,Z2,Z3)

re.block.inds <- list(1:K1,(K1+1):(K1+K2),(K1+K2+1):(K1+K2+K3))         

# Fit model using lme() and extract coefficient estimates:

Z.block  <-  list()
for (i in 1:length(re.block.inds))
   Z.block[[i]] <- as.formula(paste("~Z[,c(",paste(
             re.block.inds[[i]],collapse=","),")]-1"))

group <- rep(1,length(y))
data.fr <- groupedData(y~-1+X|group,data=data.frame(y,X,Z))
                  
fit <- lme(y~-1+X,random=list(group=pdBlocked(Z.block,
           pdClass=rep("pdIdent",length(Z.block)))),data=data.fr)

betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)

sig.epsHat <- fit$sigma
sig.uHat <- sig.epsHat*exp(unlist(fit$modelStruct))

# Plot estimated functions:

par(mfrow=c(1,3))
ng <- 101
ylim.val <- c(0,30)

# Plot first function fit (with other variables
# set at their averages):

x1g <- seq(min(x1),max(x1),length=ng)
x2g <- rep(mean(x2),length=ng)
x3g <- rep(mean(x3),length=ng)
Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)
Z3g <- ZOSull(x3g,range.x3,intKnots3)
Xg <- cbind(rep(1,ng),x1g,x2g,x3g)
Zg <- cbind(Z1g,Z2g,Z3g)
f1g <- Xg%*%betaHat + Zg%*%uHat

plot(x1g,f1g,col="blue",lwd=4,type="l",bty="l",
     xlab="inversion base height",ylab="ozone",ylim=ylim.val)
rug(x1)

# Plot second function fit (with other variables
# set at their averages):

x1g <- rep(mean(x1),length=ng)
x2g <- seq(min(x2),max(x2),length=ng)
x3g <- rep(mean(x3),length=ng)
Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)
Z3g <- ZOSull(x3g,range.x3,intKnots3)
Xg <- cbind(rep(1,ng),x1g,x2g,x3g)
Zg <- cbind(Z1g,Z2g,Z3g)
f2g <- Xg%*%betaHat + Zg%*%uHat

plot(x2g,f2g,col="red",lwd=4,type="l",bty="l",
     xlab="Daggett pressure gradient",ylab="ozone",ylim=ylim.val)
rug(x2)

# Plot third function fit (with other variables
# set at their averages):

x1g <- rep(mean(x1),length=ng)
x2g <- rep(mean(x2),length=ng)
x3g <- seq(min(x3),max(x3),length=ng)
Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)
Z3g <- ZOSull(x3g,range.x3,intKnots3)
Xg <- cbind(rep(1,ng),x1g,x2g,x3g)
Zg <- cbind(Z1g,Z2g,Z3g)
f3g <- Xg%*%betaHat + Zg%*%uHat

plot(x3g,f3g,col="green",lwd=4,type="l",bty="l",
     xlab="inversion base temperature",ylab="ozone",ylim=ylim.val)
rug(x3)

########## End of ozoneCalifLMM ##########


