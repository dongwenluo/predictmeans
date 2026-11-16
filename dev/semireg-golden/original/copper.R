########## R script: copper ##########

# For illustration of geoadditive models
# using the "copper" data.
#
# Fitting is performed via the lme() function
# within the linear mixed model framework.

# Last changed: 03 FEB 2008

# Load in required packages:

library(nlme)

# Load in required functions:

source("ZOSull.r") ; source("Ztps.r")
source("wait.r")

# Read in copper data

copper <- read.table("copper.txt",header=TRUE)

grade <- copper$grade
xcoord <- copper$xcoord
ycoord <- copper$ycoord
zcoord <- copper$zcoord

n <- length(grade)

# Look at the data:

pairs(copper) ; wait()

# Do a graphical check of the knots:

copperKnots <- read.table("copperKnots.txt")

plot(xcoord,ycoord,col="red")
points(copperKnots[,1],copperKnots[,2],cex=2,col="green3",lwd=2)
wait()

# Set up Z matrix for geographical component:

Zgeog <- Ztps(cbind(xcoord,ycoord),copperKnots)

# Set up Z matrix for `zcoord' component:

intKnots <- seq(80,120,by=5)
range.zcoord <- c(min(zcoord),max(zcoord))
Zzcoord <- ZOSull(zcoord,range.zcoord,intKnots)

X <- cbind(rep(1,n),xcoord,ycoord,zcoord)
Z <- cbind(Zgeog,Zzcoord)

# Set up various auxiliary structures required to trick
# lme() into fitting a geoadditive model:

re.block.inds <- list(1:ncol(Zgeog),(ncol(Zgeog)+1):(ncol(Zgeog)+ncol(Zzcoord)))
Z.block  <-  list()
for (i in 1:length(re.block.inds))
   Z.block[[i]] <- as.formula(paste("~Z[,c(",paste(
             re.block.inds[[i]],collapse=","),")]-1"))
group <- rep(1,n)
dataFr <- groupedData(grade~-1+X|group,data=data.frame(grade,X,Z))

# Fit model using lme() and extract coefficient estimates:
           
fit <- lme(grade~-1+X,random=list(group=pdBlocked(Z.block,
           pdClass=rep("pdIdent",length(Z.block)))),data=dataFr)

betaHat <- as.vector(fit$coef$fixed)
uHat <- unlist(fit$coef$random)

# Display the estimated additive components, starting
# with that for (xcoord,ycoord) (with zcoord set at its average):

meshSize <- 51 ; numPix <- meshSize^2

xGrid <- seq(min(xcoord),max(xcoord),length=meshSize)
yGrid <- seq(min(ycoord),max(ycoord),length=meshSize)
xyMesh <- expand.grid(xGrid,yGrid)
ZgeogMesh <- Ztps(xyMesh,copperKnots)

zAve <- rep(mean(zcoord),numPix)
ZzAve <- ZOSull(zAve,range.zcoord,intKnots)
ZMesh <- cbind(ZgeogMesh,ZzAve)
XMesh <- as.matrix(cbind(rep(1,numPix),xyMesh,zAve))

fHatxy <- matrix(XMesh%*%betaHat + ZMesh%*%uHat,meshSize,meshSize)
filled.contour(xGrid,yGrid,fHatxy,col=terrain.colors(20))
wait()

# Now do plot for zcoord with xcoord and ycoord set at their averages:

gridSize <- 101
xAve <- seq(mean(xcoord),length=gridSize)
yAve <- seq(mean(ycoord),length=gridSize)
ZgeogAve <- Ztps(cbind(xAve,yAve),copperKnots)

zGrid <- seq(min(zcoord),max(zcoord),length=gridSize)

ZzGrid <- ZOSull(zGrid,range.zcoord,intKnots)
ZGrid <- cbind(ZgeogAve,ZzGrid)
XGrid <- as.matrix(cbind(rep(1,gridSize),xAve,yAve,zGrid))

fHatz <- XGrid%*%betaHat + ZGrid%*%uHat

plot(zGrid,fHatz,bty="l",xlab="zcoord",ylab="mean grade",type="n")
lines(zGrid,fHatz,col="blue",lwd=2)

############ End of copper ############
