########## R script: growthIndianaCont ##########

# For fitting a subject specific curves model to growth 
# curve for Indiana children, that includes a contrast 
# for comparison of black and white children growth rates.

# Last changed: 11 JUL 2014

# Load in required packages and R functions:

library(lattice) ; library(nlme)
source("ZOSull.r") ; source("wait.r")

# Read in and plot the data:

growthIndiana <- read.table("growthIndiana.txt", header=TRUE)
growthIndMalOrig <- growthIndiana[growthIndiana$male==1,]
idnum <-  growthIndMalOrig$idnum
ageOrig <- growthIndMalOrig$age
heightOrig <- growthIndMalOrig$height

# Obtain unit interval-transformed versions of the main 
# variables to enhance numerical stability:

minAge <- min(ageOrig) ; maxAge <- max(ageOrig)
minHeight <- min(heightOrig) ; maxHeight <- max(heightOrig)
rangeAge <- maxAge - minAge
age <- (ageOrig - minAge)/rangeAge
rangeHeight <- maxHeight - minHeight
height <- (heightOrig - minHeight)/rangeHeight
growthIndMal <- growthIndMalOrig
growthIndMal$age <- age
growthIndMal$height <- height
black <- growthIndMal$black
growthIndMalGDOrig <-  groupedData(heightOrig~ageOrig|idnum,data=growthIndMalOrig)
growthIndMalGD <-  groupedData(height~age|idnum,data=growthIndMal)

# Create new idnum with natural ordering:

idnumNew <- rep(NA,length(idnum))
for (i in 1:length(unique(idnum)))
{
   inds <- (1:length(idnum))[idnum==unique(idnum)[i]]
   idnumNew[inds] <- i
}
idnum <- idnumNew

# Display data:

colourVec <- c("dodgerblue","deeppink","blue","darkmagenta")

tmp <- trellis.par.get("layout.heights")
tmp$strip <- 0
trellis.par.set("layout.heights",tmp)
tmp <- trellis.par.get("par.main.text")
tmp$col <- colourVec[2]
trellis.par.set("par.main.text",tmp)
tmp <- trellis.par.get("par.main.text")
tmp$cex <- 2
trellis.par.set("par.main.text",tmp)
figRaw <- xyplot(heightOrig~ageOrig|idnum,groups=idnum,data=growthIndMalGDOrig,
                    layout=c(9,13),main="black males this colour",
                    xlab="age (years)",ylab="height (centimetres)",
                    panel=function(x,y,subscripts,groups)
                    {
                       colInd <- black[subscripts[1]] + 1
                       panel.grid() 
                       panel.superpose(x,y,subscripts,groups,
                                type="b",col=colourVec[colInd],pch=1,cex=0.25)
                    })
print(figRaw)
wait()

# Set up X matrix, and Z matrix for mean baseline and contrast curves:

numObs <- length(height)
numSubj <- length(unique(idnum))
   
X <- cbind(rep(1,numObs),age,black,age*black)

numIntKnots <- 20
intKnots <- quantile(unique(age),seq(0,1,length=numIntKnots+2))[-c(1,numIntKnots+2)]
range.age <- c(1.01*min(age) - 0.01*max(age),1.01*max(age) - 0.01*min(age))

Z <- ZOSull(age,range.age,intKnots)
Zbase <- Z
Zcont <- black*Z

# Set up Z matrix for subject specific curves:

numIntKnotsSubj <- 10	
intKnotsSubj <- quantile(unique(age),
        seq(0,1,length=numIntKnotsSubj+2))[-c(1,numIntKnotsSubj+2)]

ZSubj <- ZOSull(age,range.age,intKnotsSubj)

# Set up structure for fitting via linear mixed model software:

ncZ <- ncol(Z)
re.block.inds <- list(1:ncZ,(ncZ+1):(2*ncZ),(2*ncZ+1):(3*ncZ))
dummyId <- factor(rep(1,numObs))
Z.block <- list(dummyId=pdIdent(~-1+Zbase),
                dummyId=pdIdent(~-1+Zcont),
                idnum=pdSymm(~age),
                idnum=pdIdent(~-1+ZSubj))

data.fr <- groupedData(height~X[,-1]|rep(1,length=numObs),
                       data=data.frame(height,X,Z,idnum))

# Fit using lme():

cat("Starting fitting (takes several minutes)...\n")
print(system("date"))
fit <- lme(height~-1+X,data=data.fr,random=Z.block)
cat("Finished fitting.\n")
print(system("date"))
#save(fit, file="growthIndianaCont_fit.RData")
summary(fit)
(sig.epsHat <- fit$sigma)
(sig.uHat <- sig.epsHat*exp(unlist(fit$modelStruct)))

(lamVal <- (sig.epsHat/sig.uHat)^2)

# Set up plotting grids:

ng <- 101
ageg <- seq(range.age[1],range.age[2],length=ng)
Xg <- cbind(rep(1,ng),ageg)
Zg <- ZOSull(ageg,range.age,intKnots)
betaHat <- as.vector(fit$coef$fixed)
uHatBase <- as.vector(fit$coef$random[[1]])
uHatCont <- as.vector(fit$coef$random[[2]])

fhatBaseg <- Xg%*%betaHat[1:2] + Zg%*%uHatBase
Contg <-  Xg%*%betaHat[3:4] + Zg%*%uHatCont
fhatBlackg <- fhatBaseg + Contg

# Plot fitted curves:

Xsubjg <- cbind(rep(1,ng),ageg)
ZSubjg <- ZOSull(ageg,range.age,intKnotsSubj)
curvEsts <- matrix(NA,ng,numSubj)
for (i in 1:numSubj)
{
   indsi <- (1:numObs)[idnum==i]
   uLinHati <- as.vector(fit$coef$random[[3]][i,])
   uSplHati <- as.vector(fit$coef$random[[4]][i,])
   ghati <- Xsubjg%*%uLinHati + ZSubjg%*%uSplHati
   curvEsts[,i] <- fhatBaseg + black[indsi[1]]*Contg + ghati
}

tmp <- trellis.par.get("layout.heights")
tmp$strip <- 0
trellis.par.set("layout.heights",tmp)
tmp <- trellis.par.get("par.main.text")
tmp$col <- "black"
trellis.par.set("par.main.text",tmp)
figFit <- xyplot(heightOrig~ageOrig|idnum,groups=idnum,data=growthIndMalGDOrig,
              layout=c(9,13),main="fitted subject specific curves",
              xlab="age (years)",ylab="height (centimetres)",
              panel=function(x,y,subscripts,groups)
              {
                 panel.grid() 
                 columnNum <- idnum[subscripts][1]
                 colourInd <- black[subscripts[1]] + 1
                 panel.superpose(x,y,subscripts,groups,
                          type="b",col=colourVec[colourInd],pch=1,cex=0.5)
                 panel.xyplot((rangeAge*ageg+minAge),(rangeHeight*curvEsts[,columnNum]+minHeight),
                              type="l",col=colourVec[colourInd+2])

              })
print(figFit)
wait()

# Obtain penalty matrix:

sigsq.epsHat <- fit$sigma^2
sig.uHat <- as.numeric(sqrt(sigsq.epsHat)*exp(unlist(fit$modelStruct)))

sigsq.Gbl <- sig.uHat[6]^2
sigsq.Gblcont <- sig.uHat[5]^2
cholSigma <- rbind(c(sig.uHat[2],sig.uHat[3]),c(0,sig.uHat[4]))
SigmaHat <- crossprod(cholSigma)
sigsq.Sbj <- sig.uHat[1]^2

DmatGbl <- (sigsq.epsHat/sigsq.Gbl)*diag(ncol(Z))
DmatGblcont <- (sigsq.epsHat/sigsq.Gblcont)*diag(ncol(Zg))
DmatLinSbj <- sigsq.epsHat*kronecker(diag(numSubj),solve(SigmaHat))
DmatSplSbj <- (sigsq.epsHat/sigsq.Sbj)*diag(ncol(ZSubj)*numSubj)

dimVec <- c(4,nrow(DmatGbl),nrow(DmatGblcont),nrow(DmatLinSbj),
             nrow(DmatSplSbj))

lamMat <- matrix(0,sum(dimVec),sum(dimVec))
csdV <- cumsum(dimVec)
lamMat[(csdV[1]+1):csdV[2],(csdV[1]+1):csdV[2]] <- DmatGbl
lamMat[(csdV[2]+1):csdV[3],
       (csdV[2]+1):csdV[3]] <- DmatGblcont
lamMat[(csdV[3]+1):csdV[4],
       (csdV[3]+1):csdV[4]] <- DmatLinSbj
lamMat[(csdV[4]+1):csdV[5],
       (csdV[4]+1):csdV[5]] <- DmatSplSbj

# Obtain C matrix:

uqID <- unique(idnum)
Cmat <- cbind(X,Z,Zcont)
for (iSbj in 1:numSubj)
{
   newCols <- matrix(0,numObs,2)
   indsCurr <- (1:numObs)[idnum==uqID[iSbj]]
   newCols[indsCurr,] <- X[indsCurr,1:2]
   Cmat <- cbind(Cmat,newCols)
}
for (iSbj in 1:numSubj)
{
   newCols <- matrix(0,numObs,ncol(ZSubj))
   indsCurr <- (1:numObs)[idnum==uqID[iSbj]]
   newCols[indsCurr,] <- ZSubj[indsCurr,]
   Cmat <- cbind(Cmat,newCols)
}

# Obtain full covariance matrix:

CTC <- crossprod(Cmat)
fullCovMat <- solve(CTC + lamMat)

# Find subset of covariance corresponding to the contrast curve:

contInds <- c(3,4,(ncZ+5):(2*ncZ+4))
contCovMat <- fullCovMat[contInds,contInds]

# Obtain approximate pointwise 95% confidence limits:

Cg <- cbind(Xg,Zg)
sdg <- sqrt(sigsq.epsHat)*sqrt(diag(Cg%*%contCovMat%*%t(Cg)))
lowerg <- Contg - qnorm(0.975)*sdg ; upperg <- Contg + qnorm(0.975)*sdg

# Plot contrast with variability band:

plot(rangeAge*ageg+minAge,rangeHeight*Contg,bty="l",
     type="n",xlab="age (years)",ylab="height (cm)",
     ylim=rangeHeight*range(c(lowerg,upperg)), 
     main="contrast curve (mean height difference)",
     cex.main=2,col.main="darkgreen")

polygon(c(rangeAge*ageg+minAge,rev(rangeAge*ageg+minAge)),
        (rangeHeight*c(lowerg,rev(upperg))),col="PaleGreen",
        border=FALSE)

lines(rangeAge*ageg+minAge,rangeHeight*Contg,col="DarkGreen",lwd=2)
abline(0,0,col="red")

wait()

# Check residuals:

print(plot(fit))

############ End of growthIndianaCont ###########
# 
# > summary(fit)
# Linear mixed-effects model fit by REML
# Data: data.fr 
# AIC       BIC   logLik
# -12750.2 -12687.28 6386.102
# 
# Random effects:
#   Formula: ~-1 + Zbase | dummyId
# Structure: Multiple of an Identity
# Zbase1    Zbase2    Zbase3    Zbase4    Zbase5    Zbase6    Zbase7    Zbase8
# StdDev: 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676
# Zbase9   Zbase10   Zbase11   Zbase12   Zbase13   Zbase14   Zbase15   Zbase16
# StdDev: 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676
# Zbase17   Zbase18   Zbase19   Zbase20   Zbase21   Zbase22
# StdDev: 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676 0.7100676
# 
# Formula: ~-1 + Zcont | dummyId %in% dummyId
# Structure: Multiple of an Identity
# Zcont1    Zcont2    Zcont3    Zcont4    Zcont5    Zcont6    Zcont7    Zcont8
# StdDev: 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805
# Zcont9   Zcont10   Zcont11   Zcont12   Zcont13   Zcont14   Zcont15   Zcont16
# StdDev: 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805
# Zcont17   Zcont18   Zcont19   Zcont20   Zcont21   Zcont22
# StdDev: 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805 0.6660805
# 
# Formula: ~age | idnum %in% dummyId %in% dummyId
# Structure: General positive-definite
#             StdDev     Corr  
# (Intercept) 0.05582902 (Intr)
# age         0.04559385 0.133 
# 
# Formula: ~-1 + ZSubj | idnum %in% idnum %in% dummyId %in% dummyId
# Structure: Multiple of an Identity
# ZSubj1   ZSubj2   ZSubj3   ZSubj4   ZSubj5   ZSubj6   ZSubj7   ZSubj8   ZSubj9
# StdDev: 1.000212 1.000212 1.000212 1.000212 1.000212 1.000212 1.000212 1.000212 1.000212
# ZSubj10  ZSubj11  ZSubj12    Residual
# StdDev: 1.000212 1.000212 1.000212 0.006830225
# 
# Fixed effects:  height ~ -1 + X 
# Error in `.rowNamesDF<-`(x, value = value) : 
#   duplicate 'row.names' are not allowed
# In addition: Warning message:
#   non-unique value when setting 'row.names': ‘X’ 
# > (sig.epsHat <- fit$sigma)
# [1] 0.006830225
# > (sig.uHat <- sig.epsHat*exp(unlist(fit$modelStruct)))
# reStruct.idnum  reStruct.idnum1  reStruct.idnum2  reStruct.idnum3 reStruct.dummyId 
# 1.00021159       0.05561525       0.00729784       0.04536514       0.66608051 
# reStruct.dummyId 
# 0.71006755 
# > (lamVal <- (sig.epsHat/sig.uHat)^2)
# reStruct.idnum  reStruct.idnum1  reStruct.idnum2  reStruct.idnum3 reStruct.dummyId 
# 4.663223e-05     1.508281e-02     8.759541e-01     2.266864e-02     1.051518e-04 
# reStruct.dummyId 
# 9.252746e-05 
# 




