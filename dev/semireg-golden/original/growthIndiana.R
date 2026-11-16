########## R script: growthIndiana ##########

# For fitting a subject specific curves model
# to growth curve data on Indiana children.

# Last changed: 25 MAR 2013

# Load in required packages and R functions:

library(lattice) ; library(nlme)
source("ZOSull.r") ; source("wait.r")

# Read in and plot the data:

growthIndiana <- read.table("growthIndiana.txt",header=TRUE)

# Extract data for African-American subjects:

growthIndBlkMal <- growthIndiana[(growthIndiana$male==1)&(growthIndiana$black==1),]
age <- growthIndBlkMal$age
height <- growthIndBlkMal$height
idnumOrig <-  growthIndBlkMal$idnum

# Create new (ordered) ID numbers:

idnum <- rep(NA,length(idnumOrig))
uqID <- unique(idnumOrig)
for (i in 1:length(uqID))
      idnum[idnumOrig==uqID[i]] <- i
growthIndBlkMal$idnum <- idnum

growthIndBlkMalGD <-  groupedData(height~age|idnum,data=growthIndBlkMal)
tmp <- trellis.par.get("layout.heights")
tmp$strip <- 0
trellis.par.set("layout.heights",tmp)
figRaw <- xyplot(height~age|idnum,groups=idnum,data=growthIndBlkMalGD,
                 layout=c(4,7),main="Indiana growth data (black male subset)",
                 xlab="age (years)",ylab="height (centimetres)",
                 panel=function(x,y,subscripts,groups)
                 {
                    panel.grid() 
                    panel.superpose(x,y,subscripts,groups,
                             type="b",col="DeepPink",pch=1,cex=0.5)
                 })

print(figRaw)
wait()

# Set up X matrix, and Z matrix for mean curve:

cat("Starting fitting...\n")
numObs <- length(height)
numSubj <- length(unique(idnum))
uqID <- unique(idnum)

X <- model.matrix(height~age)
numIntKnots <- 20
intKnots <- quantile(unique(age),
                     seq(0,1,length=numIntKnots+2))[-c(1,numIntKnots+2)]

range.age <- c(5.5,20)
Z <- ZOSull(age,range.age,intKnots)

# Set up Z matrix for subject specific curves:

numIntKnotsSubj <- 10	
intKnotsSubj <- quantile(unique(age),
        seq(0,1,length=numIntKnotsSubj+2))[-c(1,numIntKnotsSubj+2)]

ZSubj <- ZOSull(age,range.age,intKnotsSubj)

# Fit using linear mixed model software:

dummyId <- factor(rep(1,numObs))
Z.block<-list(dummyId=pdIdent(~-1+Z),
              idnum=pdSymm(~age),
              idnum=pdIdent(~-1+ZSubj))
data.fr <- groupedData(height~X[,-1]|rep(1,length=numObs),
                    data=data.frame(height,X,Z,ZSubj,idnum))
fit <- lme(height~-1+X,data=data.fr,random=Z.block)
cat("Finished fitting.\n")

# Plot the fits:

ng <- 101
ageg <- seq(range.age[1],range.age[2],length=ng)
Xg <- cbind(rep(1,ng),ageg)
Zg <- ZOSull(ageg,range.age,intKnots)
Zsubjg <- ZOSull(ageg,range.age,intKnotsSubj)
betaHat <- as.vector(fit$coef$fixed)
uHat <- as.vector(fit$coef$random[[1]])
fhatg <- Xg%*%betaHat + Zg%*%uHat

curvEsts <- matrix(NA,ng,numSubj)
for (i in 1:numSubj)
{
   indsi <- (1:numObs)[idnum==i]
   uLinHati <- as.vector(fit$coef$random[[2]][i,])
   uSplHati <- as.vector(fit$coef$random[[3]][i,])
   ghati <- Xg%*%uLinHati + Zsubjg%*%uSplHati
   curvEsts[,i] <- fhatg + ghati
}

figFit <- xyplot(height~age|idnum,groups=idnum,data=growthIndBlkMalGD,
                 layout=c(4,8),main="fitted subject specific curves",
                 xlab="age (years)",ylab="height (centimetres)",
                 panel=function(x,y,subscripts,groups)
                 {
                    panel.grid() 
                    colNum <- idnum[subscripts][1]
                    panel.superpose(x,y,subscripts,groups,
                             type="b",col="DeepPink",pch=1,cex=0.5)
                    panel.xyplot(ageg,curvEsts[,colNum],col="blue",type="l")

                 })
print(figFit)
wait()

# Plot the residuals:

print(plot(fit))

############ End of growthIndiana ###########








