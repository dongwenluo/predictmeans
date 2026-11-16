
########## R script: scallop.Rs (Updated for modern rgl) ##########

# For bivariate smoothing (low-rank kriging) of 
# some scallop catch data. Visualisation using 
# 3D spin graphics is also provided.

# Original: 17 NOV 2008 by M.P. Wand.
# Updated: Modern rgl (*3d API) compatibility.

# Load in required packages:
library(nlme)
library(rgl)

# Load in required functions:
source("pointsInPoly.r")


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
#----------------------------------------------------------------------------

clean <- function(x)
  return(x[is.na(x)==0])
#----------------------------------------------------------------------------
# Read in and extract data:
scallop <- read.table("scallop.txt", header = TRUE)

n <- nrow(scallop)
x <- scallop$longitude
y <- scallop$latitude
z <- asinh(scallop$tot.catch)

# Read in knots (obtained elsewhere using the CLARA algorithm):
scallopKnots <- read.table("scallopKnots.txt")

# Read in boundary polygon, corresponding to the
# location of the data (obtained earlier by eye):
scallopBdry <- read.table("scallopBdry.txt")

# Plot the knots and boundary together with the data:
par(mfrow = c(1, 1))
plot(x, y, bty = "l", xlab = "longitude", ylab = "latitude", col = "red",
     xlim = c(-73.8, -71.3), ylim = c(38.5, 41.2),
     main = "knot locations for bivariate spline smoothing and\nboundary used for plotting purposes.")
points(scallopKnots[, 1], scallopKnots[, 2], pch = 1, cex = 2, col = "green3", lwd = 2)
lines(scallopBdry[, 1], scallopBdry[, 2], col = "blue", lwd = 2)

# Perform low-rank kriging using linear mixed model representation and software:
X <- cbind(rep(1, n), x, y)
Z <- Ztps(cbind(x, y), scallopKnots)

group <- rep(1, n)
data.fr <- groupedData(y ~ x | group, data = data.frame(x, y))
fit <- lme(z ~ -1 + X, random = pdIdent(~ -1 + Z), data = data.fr)

# Display the estimated surface via an image plot:
betaHat <- as.vector(fit$coef$fixed)
uHat <- unlist(fit$coef$random)

ng <- 201
xgrid <- seq(-73.8, -71.3, length = ng)
ygrid <- seq(38.5, 41.2, length = ng)
xymesh <- as.matrix(expand.grid(xgrid, ygrid))

Xmesh <- as.matrix(cbind(rep(1, nrow(xymesh)), xymesh))
Zmesh <- Ztps(xymesh, scallopKnots)

# Switch off pixels outside the boundary of the data:
numPix <- ng^2
outInds <- (1:numPix)[pointsInPoly(xymesh, scallopBdry) == FALSE]
zmesh <- Xmesh %*% betaHat + Zmesh %*% uHat
zmesh[outInds] <- NA

# Generate image plot:
zmesh <- matrix(zmesh, ng, ng)
image(xgrid, ygrid, zmesh, bty = "l", xlab = "longitude", ylab = "latitude",
      col = terrain.colors(1000),
      main = "image plot of bivariate smooth of arcsinh(total catch)")
lines(scallopBdry, lwd = 2)
points(x, y, col = "red", cex = 0.7)

# Add legend:
legx <- c(-72, -71.5)
legy <- c(39.0, 39.2)
legPolyx <- c(legx, rev(legx))
legPolyy <- c(rep(legy[1], 2), rep(legy[2], 2))
zRange <- range(clean(zmesh))
legxgrid <- seq(min(legx), max(legx), length = ng)
legygrid <- seq(min(legy), max(legy), length = ng)
legz <- matrix(rep(seq(0, 1, length = ng), ng), ng, ng)
image(legxgrid, legygrid, legz, col = terrain.colors(1000), add = TRUE)
polygon(legPolyx, legPolyy, lwd = 2)
text(legx[1], 1.25 * legy[2] - 0.25 * legy[1], signif(zRange[1], 3))
text(legx[2], 1.25 * legy[2] - 0.25 * legy[1], signif(zRange[2], 3))
text((legx[1] + legx[2]) / 2, 1.25 * legy[1] - 0.25 * legy[2], "mean arcsinh(total catch)")

# -------------------------------------------------------------------------
# 3D graphical display of fit (Updated to modern *3d functions)
# -------------------------------------------------------------------------

# Transformation of the data to a cube anchored at the origin:
x <- 2 * (x - min(x)) / (max(x) - min(x))
y <- 2 * (y - min(y)) / (max(y) - min(y))
z <- 2 * (z - min(z)) / (max(z) - min(z))

# Set up 'Right-Hand Coordinate System' versions using modern *3d functions
# Note: Keeping the (x, z, -y) mapping preserves M.P. Wand's exact camera 
# orientation where Y (z) is vertical in OpenGL space.
rhcs.spheres <- function(x, y, z, radius, ...) {
  spheres3d(x, z, -y, radius = radius, ...)
}

rhcs.texts <- function(x, y, z, text, ...) {
  texts3d(x, z, -y, texts = text, ...)
}

rhcs.triangles <- function(x, y, z, ...) {
  triangles3d(x, z, -y, ...)
}

rhcs.quads <- function(x, y, z, ...) {
  quads3d(x, z, -y, ...)
}

rhcs.surface <- function(x, y, z, ...) {
  # In modern rgl, surface3d(x, y, z) treats z as the height matrix by default.
  # To match legacy rgl.surface(x, z, -y) where the 2nd coordinate is height,
  # we pass coords = c(1, 3, 2) so z remains the vertical axis in (x, z, -y) space.
  surface3d(x, -y, z, coords = c(1, 3, 2), ...)
}

rhcs.clear <- function(type = "shapes") {
  view3d(theta = 135, phi = 15)
  clear3d(type = type)
}

# Set parameters:
bg.col <- "lightskyblue1"
bs.col <- "blue1"
ph.col <- "green3"
ax.col <- "DeepPink"
tx.col <- "maroon"      
sp.rad <- 0.03
xtx.pos <- 2.2 ; ytx.pos <- 2.2 ; ztx.pos <- 2.2
xax.low <- -0.2; yax.low <- -0.2; zax.low <- -0.2
xax.upp <- 2.1 ; yax.upp <- 2.1 ; zax.upp <- 2.1
al.val <- 0.7
sz.val <- 15

rhcs.lines3d <- function(x, y, z, add = FALSE, size = 3, ...) {
  if (length(x) != length(y) || length(x) != length(z)) 
    stop("Lengths of x, y, z do not match")
  if (!add) 
    clear3d()
  # Modern lines3d connects sequential vertices automatically and uses lwd instead of size
  lines3d(x, z, -y, lwd = size, ...)
}

createAxesandBase <- function() {
  rhcs.clear()
  bg3d(color = bg.col)
  rhcs.lines3d(c(xax.low, xax.upp), rep(0, 2), rep(0, 2),
               size = 3, col = ax.col, add = TRUE)
  rhcs.lines3d(rep(0, 2), c(yax.low, yax.upp), rep(0, 2),
               size = 3, col = ax.col, add = TRUE)
  rhcs.lines3d(rep(0, 2), rep(0, 2), c(zax.low, zax.upp),
               size = 3, col = ax.col, add = TRUE)
  rhcs.texts(xtx.pos, 0, 0, "longitude", col = tx.col)
  rhcs.texts(0, ytx.pos, 0, "latitude", col = tx.col)
  rhcs.texts(0, 0, ztx.pos, "arcsinh(catch)", col = tx.col)
  rhcs.quads(c(xax.low, xax.low, xax.upp, xax.upp),
             c(yax.low, yax.upp, yax.upp, yax.low),
             rep(0, 4), col = bs.col, alpha = al.val)
}

# Open a fresh 3D device window
open3d(windowRect = c(50, 50, 800, 800))

createAxesandBase()
rhcs.spheres(x, y, z, col = "red", radius = sp.rad)

createAxesandBase()
rhcs.spheres(x, y, z, col = "red", radius = sp.rad)

ngrid <- 40
xgrid <- seq(0, 2, length.out = ngrid)
ygrid <- seq(0, 2, length.out = ngrid)
xymesh <- as.matrix(expand.grid(xgrid, ygrid))

scallopKnots <- rbind(
  c(1.969469, 1.683451), c(1.816794, 1.654677), c(1.755731, 1.482015),
  c(1.404579, 1.266188), c(1.267174, 1.294962), c(1.038166, 1.122301),
  c(0.839699, 1.079135), c(0.9618335, 0.8489168), c(0.9465633, 0.7482032),
  c(0.7786271, 0.762586), c(0.8091585, 0.5755416), c(0.6564926, 0.3165492),
  c(0.4732862, 0.2302184), c(0.2900798, 0), c(0.0916032, 0.1438876),
  c(0.1679361, 0.1151049), c(0.2442782, 0.2446011), c(0.229008, 0.4604281),
  c(0.3511425, 0.6043157), c(0.5038176, 0.4604281), c(0.7022942, 0.618707),
  c(0.5343489, 0.7913686), c(0.5954208, 0.9640303), c(0.7480958, 0.8633081),
  c(0.6564926, 1.079135), c(0.381683, 1.151075), c(0.3053409, 1.237405),
  c(0.1984766, 1.251797), c(0.1374048, 1.338128), c(0.5343489, 1.122301),
  c(0.8091585, 1.19424), c(0.8091585, 1.294962), c(1.053437, 1.338128),
  c(0.9313022, 1.410067), c(0.9313022, 1.553955), c(1.129770, 1.611511),
  c(1.160310, 1.482015), c(1.236643, 1.410067), c(1.389318, 1.410067),
  c(1.694659, 1.956835), c(1.511453, 1.841721), c(1.282445, 1.913669),
  c(1.206112, 1.712225), c(1.389318, 1.697842), c(1.297715, 1.582729),
  c(1.465651, 1.525181), c(1.572525, 1.697842), c(1.664128, 1.841721)
)

scallopBdry <- rbind(
  c(0.1581562, 1.569938), c(0.4091438, 1.490103), c(0.6083402, 1.423573),
  c(0.731842, 1.578809), c(0.8672956, 1.756222), c(1.046572, 1.973552),
  c(1.381222, 2.084435), c(1.684001, 2.142094), c(1.946940, 2.044517),
  c(1.931005, 1.822752), c(2.102313, 1.609856), c(1.847342, 1.370349),
  c(1.472853, 1.002218), c(1.221865, 0.6606982), c(0.9788456, 0.3324846),
  c(0.5963884, -0.0312115), c(0.1940116, -0.1820123), c(-0.08486342, -0.01790554),
  c(-0.1924295, 0.4123203), c(0.02668658, 0.6606982), c(0.2497866, 0.7937577),
  c(0.1223009, 1.002218), c(-0.1167349, 1.277207), c(-0.09283128, 1.530021),
  c(0.1581562, 1.569938)
)

innerInds <- pointsInPoly(xymesh, scallopBdry)

n <- length(z)

X <- cbind(rep(1, n), x, y)
Z <- Ztps(cbind(x, y), scallopKnots)

group <- rep(1, n)
data.fr <- groupedData(y ~ x | group, data = data.frame(x, y))
fit <- lme(z ~ -1 + X, random = pdIdent(~ -1 + Z), data = data.fr)

betaHat <- as.vector(fit$coef$fixed)
uHat <- unlist(fit$coef$random)

Xmesh <- as.matrix(cbind(rep(1, nrow(xymesh)), xymesh))

Zmesh <- Ztps(xymesh, scallopKnots)
zmesh <- as.vector(Xmesh %*% betaHat + Zmesh %*% uHat)

lowTh <- -0.1
zmesh[zmesh < lowTh] <- lowTh
al.val <- 0.75
alvec <- rep(al.val, length(zmesh))
alvec[innerInds == FALSE] <- 0

numCol <- 501
zcolgrid <- seq(min(zmesh), max(zmesh), length.out = numCol)
tz <- round((numCol - 1) * (zmesh - min(zmesh)) / (max(zmesh) - min(zmesh))) + 1
colorlut <- terrain.colors(numCol)
colvec <- colorlut[tz] # assign colors to heights for each point
colvec[innerInds == FALSE] <- bg.col

# Reshape zmesh into an ngrid x ngrid matrix for modern surface3d
zmesh_mat <- matrix(zmesh, nrow = ngrid, ncol = ngrid)

# Render the smoothed kriging surface using modern surface3d via rhcs.surface
rhcs.surface(xgrid, ygrid, zmesh_mat, color = colvec, alpha = alvec)

########## End of scallop.Rs ##########