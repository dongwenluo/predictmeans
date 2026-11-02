smZ <- function(x, k=6, intKnots=NULL, range.x=NULL, degree=3, 
                type=c("ZOSull", "Ztps", "ns", "bs", "bernstein", "bSpline",
                      "nSpline", "cSpline", "iSpline", "mSpline", "smspline"), 
                by=NULL, group=FALSE, intercept=FALSE, pred=FALSE, ...){
  
  type <- as.character(type)
  type <- match.arg(type)
  
  if (type!="Ztps") {
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) x <- x[!nax]
	if (is.null(range.x)) range.x <- c(1.01*min(x) - 0.01*max(x), 
                                       1.01*max(x) - 0.01*min(x))
  }else{
    x <- as.matrix(x)
    nax <- rowSums(is.na(x)) > 0
    if (nas <- any(nax)) x <- x[!nax,]  
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
    }else{
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
  }else if(type=="Ztps"){
    if (is.factor(by)) {
      ZZ <- Ztps(x, k=k, knots=intKnots, range.x=range.x)
      if (is.null(intKnots)) Z_knots <- attr(ZZ,"knots")
	  
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
    }else{
      Z <- Ztps(x, k=k, knots=intKnots, range.x=range.x)
	  if (is.null(intKnots)) Z_knots <- attr(Z,"knots")
      if (nas) {
        nmat <- matrix(NA, length(nax), ncol(Z))
        nmat[!nax, ] <- Z
        Z <- nmat
      }  
      Z <- as(Z, "sparseMatrix")
    }
    if (is.null(intKnots)) attr(Z,"knots") <- Z_knots else attr(Z,"knots") <- intKnots
	if (!is.null(range.x)) attr(Z,"range.x") <- range.x
  }else if(type %in% c("bs", "ns", "bernstein", "bSpline", "nSpline", "cSpline", "iSpline", "mSpline")){ 
    if (!type %in% c("bs", "ns")) pm_need("splines2", paste0("type = \"", type, "\""))
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
      if (!group) Z <- do.call("cbind", Z)
    }else{      
      Z <- as(ZZ, "sparseMatrix")
    }
    attr(Z,"degree") <- Z_attr$degree
    attr(Z,"knots") <- Z_attr$knots
    attr(Z,"range.x") <- Z_attr$Boundary.knots
    attr(Z,"intercept") <- Z_attr$intercept
  }else if(type=="smspline"){
    pm_need("lmeSplines", "type = \"smspline\"")
    if (is.factor(by)) {
      ZZ <- lmeSplines::smspline(x)	
	  if (pred) ZZ <- lmeSplines::approx.Z(Z=range.x, oldtimes=intKnots, newtimes=x)
      if (ncol(ZZ)*nlevels(by) > nrow(ZZ) - 4) stop("'smspline' is not suitable for this data!")	  
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
    }else{
      Z <- lmeSplines::smspline(x)
	  if (pred) Z <- lmeSplines::approx.Z(Z=range.x, oldtimes=intKnots, newtimes=x)
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




