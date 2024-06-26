\name{smZ}
\alias{smZ}
\title{
  Generate Sparse Matrix Z for penalized spline smoothing
}
\description{Constructs a sparse matrix (Z) of a spline function with for a covariate with(out) group. 
}
\usage{
smZ(x, k=6, intKnots=NULL, range.x=NULL, degree=3, 
    type=c("ZOSull", "Ztps", "ns", "bs", "bernstein", "bSpline",
         "nSpline", "cSpline", "iSpline", "mSpline", "smspline"), 
    by=NULL, group=FALSE, intercept=FALSE, pred=FALSE, ...)
}

\arguments{	
\item{x}{x covariate for the smooth function. Missing values are allowed and will be returned as they are.}
\item{k}{Degree of freedom that equals to the column number of the returned matrix. One can specify df rather than knots, then the function chooses df - degree - as.integer(intercept) internal knots at suitable quantiles of x ignoring missing values and those x outside of the boundary. If internal knots are specified via knots, the specified df will be ignored.}
\item{intKnots}{Ordered array of length smaller than that of x and consisting of unique numbers between min(x) and max(x) that specifies the positions of internal knots, that define the spline basis (see the Wand and Ormerod (2008) reference below for full mathematical details).}
\item{range.x}{Array of length 2 such that range.x[1] >= min(x) and range.x[2] <= max(x).}
\item{degree}{Integer: degree of (truncated) polynomial.}
\item{type}{Type of splines including "ZOSull", "Ztps", "ns", "bs", "bernstein", "bSpline", "nSpline", "cSpline", "iSpline", "mSpline" and "smspline", the default is "ZOSull".}
\item{by}{Factor for group wise splines.}
\item{group}{When \code{by != NULL}, producing group wise splines with radom effects separately.}
\item{intercept}{If TRUE, all of the spline basis functions are returned. Notice that when using I-Spline for monotonic regression, intercept = TRUE should be set even when an intercept term is considered additional to the spline basis functions.}
\item{pred}{If TRUE, the function \code{smZ} will be applied for prediction purpose, this option mainly used by function \code{semipred} internally.}
\item{...}{Further arguments for passing on to model setup routines, such as
drv: either 0,1 or 2 with a default value of 0. If drv = 1 then the first derivatives of the O'Sullivan spline basis functions are computed instead. Similarly, if drv = 2 then the second derivatives are computed.
}
}
\value{
\item{Z}{A (or a list of) spline design matrix used in the list \code{smoothZ}.}
}
\author{
Dongwen Luo, Siva Ganesh and John Koolaard
}

\references{
O'Sullivan, F. (1986). A statistical perspective on ill-posed inverse problems (with discussion). \emph{Statistical Science}, \bold{1}, 505-527.
}

\examples{
x <- seq.int(0, 1, by = 0.01)
knots <- c(0.3, 0.5, 0.6)

zosuMat <- smZ(x, intKnots = knots)
bsMat <- smZ(x, intKnots = knots, degree = 2, type="bs")
isMat <- smZ(x, intKnots = knots, degree = 2, type="iSpline")

splst <- list(zosuMat, bsMat, isMat)
for (i in splst) {
op <- par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
matplot(x, i, type = "l", ylab = "I-spline basis")
abline(v = knots, lty = 2, col = "gray")
## reset to previous plotting settings
par(op)
}

f <- gl(4, 25, length=length(x))
zosuMat_by <- smZ(x, intKnots = knots, by=f) # one sparse matrix
str(zosuMat_by)

zosuMat_by <- smZ(x, intKnots = knots, by=f, group=TRUE) # a list of sparse matrix
str(zosuMat_by)
}