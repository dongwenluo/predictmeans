#' Generate Sparse Matrix Z for penalized spline smoothing
#'
#' Constructs a sparse matrix (Z) of a spline function with for a covariate
#' with(out) group.
#'
#'
#' @param x x covariate for the smooth function. Missing values are allowed and
#' will be returned as they are.
#' @param k Degree of freedom that equals to the column number of the returned
#' matrix. One can specify df rather than knots, then the function chooses df -
#' degree - as.integer(intercept) internal knots at suitable quantiles of x
#' ignoring missing values and those x outside of the boundary. If internal
#' knots are specified via knots, the specified df will be ignored.
#' @param intKnots Ordered array of length smaller than that of x and
#' consisting of unique numbers between min(x) and max(x) that specifies the
#' positions of internal knots, that define the spline basis (see the Wand and
#' Ormerod (2008) reference below for full mathematical details).
#' @param range.x Array of length 2 such that range.x[1] >= min(x) and
#' range.x[2] <= max(x).
#' @param degree Integer: degree of (truncated) polynomial.
#' @param type Type of splines including "ZOSull", "Ztps", "ns", "bs",
#' "bernstein", "bSpline", "nSpline", "cSpline", "iSpline", "mSpline" and
#' "smspline", the default is "ZOSull".
#' @param by Factor for group wise splines.
#' @param group When \code{by != NULL}, producing group wise splines with radom
#' effects separately.
#' @param intercept If TRUE, all of the spline basis functions are returned.
#' Notice that when using I-Spline for monotonic regression, intercept = TRUE
#' should be set even when an intercept term is considered additional to the
#' spline basis functions.
#' @param pred If TRUE, the function \code{smZ} will be applied for prediction
#' purpose, this option mainly used by function \code{semipred} internally.
#' @param ... Further arguments for passing on to model setup routines, such as
#' drv: either 0,1 or 2 with a default value of 0. If drv = 1 then the first
#' derivatives of the O'Sullivan spline basis functions are computed instead.
#' Similarly, if drv = 2 then the second derivatives are computed.
#' @return \item{Z}{A (or a list of) spline design matrix used in the list
#' \code{smoothZ}.}
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @references O'Sullivan, F. (1986). A statistical perspective on ill-posed
#' inverse problems (with discussion). \emph{Statistical Science}, \bold{1},
#' 505-527.
#' @examples
#'
#' x <- seq.int(0, 1, by = 0.01)
#' knots <- c(0.3, 0.5, 0.6)
#'
#' zosuMat <- smZ(x, intKnots = knots)
#' bsMat <- smZ(x, intKnots = knots, degree = 2, type="bs")
#' isMat <- smZ(x, intKnots = knots, degree = 2, type="iSpline")
#'
#' splst <- list(zosuMat, bsMat, isMat)
#' for (i in splst) {
#' op <- par(mar = c(2.5, 2.5, 0.2, 0.1), mgp = c(1.5, 0.5, 0))
#' matplot(x, i, type = "l", ylab = "I-spline basis")
#' abline(v = knots, lty = 2, col = "gray")
#' ## reset to previous plotting settings
#' par(op)
#' }
#'
#' f <- gl(4, 25, length=length(x))
#' zosuMat_by <- smZ(x, intKnots = knots, by=f) # one sparse matrix
#' str(zosuMat_by)
#'
#' zosuMat_by <- smZ(x, intKnots = knots, by=f, group=TRUE) # a list of sparse matrix
#' str(zosuMat_by)
#' @importFrom methods as
#' @importFrom stats quantile
#' @export
smZ <- function(x, k=6, intKnots=NULL, range.x=NULL, degree=3,
                type=c("ZOSull", "Ztps", "ns", "bs", "bernstein", "bSpline",
                       "nSpline", "cSpline", "iSpline", "mSpline", "smspline"),
                by=NULL, group=FALSE, intercept=FALSE, pred=FALSE, ...){

  type <- as.character(type)
  type <- match.arg(type)

  if (type != "Ztps") {
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) {
      x <- x[!nax]
    }
    if (is.null(range.x)) {
      range.x <- c(1.01*min(x) - 0.01*max(x), 1.01*max(x) - 0.01*min(x))
    }
  } else {
    x <- as.matrix(x)
    nax <- rowSums(is.na(x)) > 0
    if (nas <- any(nax)) {
      x <- x[!nax,]
    }
  }

  if (type=="ZOSull") {
    if (is.factor(by)) {
      if (is.null(intKnots)) {
        k <- min(k, floor(length(x)/nlevels(by))-2)
        intKnots <-  quantile(unique(x), seq(0,1,length = k+2)[-c(1,k+2)])
      }
      ZZ <- ZOSull(x, range.x=range.x, intKnots=intKnots, ...)
      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(ZZ))
        nmat[!nax, ] <- ZZ
        ZZ <- nmat
      }
      by_lst <- lapply(split(data.frame(ZZ), by), as.matrix)[levels(by)]
      Z <- lapply(by_lst, function(x) {
        Z_i <- matrix(0, nrow(ZZ), ncol(x))
        Z_i[as.numeric(rownames(x)),] <- x
        as(Z_i, "sparseMatrix")
      })
      if (!group) Z <- do.call("cbind", Z)
      attr(Z,"range.x") <- range.x
      attr(Z,"knots") <- intKnots
    } else {
      if (is.null(intKnots)) {
        k <- min(k, length(x)-2)
        intKnots <-  quantile(unique(x), seq(0,1,length = k+2)[-c(1,k+2)])
      }
      Z <- ZOSull(x, intKnots=intKnots, range.x=range.x, ...)
      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(Z))
        nmat[!nax, ] <- Z
        Z <- nmat
      }
      Z <- as(Z, "sparseMatrix")
      attr(Z,"range.x") <- range.x
      attr(Z,"knots") <- intKnots
    }
  } else if (type=="Ztps") {
    if (is.factor(by)) {
      ZZ <- Ztps(x, k=k, knots=intKnots, range.x=range.x)
      if (is.null(intKnots)) {
        Z_knots <- attr(ZZ,"knots")
      }

      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(ZZ))
        nmat[!nax, ] <- ZZ
        ZZ <- nmat
      }
      by_lst <- lapply(split(data.frame(ZZ), by), as.matrix)[levels(by)]
      Z <- lapply(by_lst, function(x) {
        Z_i <- matrix(0, nrow(ZZ), ncol(x))
        Z_i[as.numeric(rownames(x)),] <- x
        as(Z_i, "sparseMatrix")
      })
      if (!group) {
        Z <- do.call("cbind", Z)
      }
    } else {
      Z <- Ztps(x, k=k, knots=intKnots, range.x=range.x)
      if (is.null(intKnots)) {
        Z_knots <- attr(Z,"knots")
      }
      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(Z))
        nmat[!nax, ] <- Z
        Z <- nmat
      }
      Z <- as(Z, "sparseMatrix")
    }
    if (is.null(intKnots)) {
      attr(Z,"knots") <- Z_knots
    } else {
      attr(Z,"knots") <- intKnots
    }
    if (!is.null(range.x)) {
      attr(Z,"range.x") <- range.x
    }
  } else if (type %in% c("bs", "ns", "bernstein", "bSpline", "nSpline", "cSpline", "iSpline", "mSpline")) {
    ZZ <- switch(type,
                 "bs"=as.matrix(splines::bs(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept)),
                 "ns"=as.matrix(splines::ns(x, df=k, knots = intKnots, Boundary.knots = range.x, intercept=intercept)),
                 "bernstein"=as.matrix(splines2::bernsteinPoly(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept)),
                 "bSpline"=as.matrix(splines2::bSpline(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept)),
                 "nSpline"=as.matrix(splines2::naturalSpline(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept)),
                 "cSpline"=as.matrix(splines2::cSpline(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept)),
                 "iSpline"=as.matrix(splines2::iSpline(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept)),
                 "mSpline"=as.matrix(splines2::mSpline(x, df=k, knots = intKnots, degree = degree, Boundary.knots = range.x, intercept=intercept))
    )
    Z_attr <- attributes(ZZ)[3:6]
    attr(ZZ,"class") <- "matrix"

    if (nas) {
      nmat <- matrix(NA, length(nax), ncol(ZZ))
      nmat[!nax, ] <- ZZ
      ZZ <- nmat
    }
    if (is.factor(by)) {
      by_lst <- lapply(split(data.frame(ZZ), by), as.matrix)[levels(by)]
      Z <- lapply(by_lst, function(x) {
        Z_i <- matrix(0, nrow(ZZ), ncol(x))
        Z_i[as.numeric(rownames(x)),] <- x
        as(Z_i, "sparseMatrix")
      })
      if (!group) {
        Z <- do.call("cbind", Z)
      }
    } else {
      Z <- as(ZZ, "sparseMatrix")
    }
    attr(Z,"degree") <- Z_attr$degree
    attr(Z,"knots") <- Z_attr$knots
    attr(Z,"range.x") <- Z_attr$Boundary.knots
    attr(Z,"intercept") <- Z_attr$intercept
  } else if (type=="smspline") {
    if (is.factor(by)) {
      ZZ <- lmeSplines::smspline(x)
      if (pred) {
        ZZ <- lmeSplines::approx.Z(Z=range.x, oldtimes=intKnots, newtimes=x)
      }
      if (ncol(ZZ)*nlevels(by) > nrow(ZZ) - 4) {
        stop("'smspline' is not suitable for this data!")
      }
      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(ZZ))
        nmat[!nax, ] <- ZZ
        ZZ <- nmat
      }
      by_lst <- lapply(split(data.frame(ZZ), by), as.matrix)[levels(by)]
      Z <- lapply(by_lst, function(x) {
        Z_i <- matrix(0, nrow(ZZ), ncol(x))
        Z_i[as.numeric(rownames(x)),] <- x
        as(Z_i, "sparseMatrix")
      })
      if (!group) {
        Z <- do.call("cbind", Z)
      }
    } else {
      Z <- lmeSplines::smspline(x)
      if (pred) {
        Z <- lmeSplines::approx.Z(Z=range.x, oldtimes=intKnots, newtimes=x)
      }
      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(Z))
        nmat[!nax, ] <- Z
        Z <- nmat
      }
      Z <- as(Z, "sparseMatrix")
    }
    attr(Z,"knots") <- lmeSplines::smspline.v(x)$Xs[,2]
    attr(Z,"range.x") <- lmeSplines::smspline.v(x)$Zs
  }

  attr(Z,"type") <- type
  return(Z=Z)
}




