semireg <- function(formula, data, family = NULL, ngbinomial=FALSE, REML = TRUE, 
                     smoothZ = list(), ncenter=TRUE, nscale=FALSE, resp_scale=FALSE, 
                     control = lmerControl(optimizer="bobyqa"), start = NULL, 
                     verbose = FALSE, drop.unused.levels=TRUE, subset, weights, 
                     offset, contrasts = NULL, prt=TRUE, predict_info=TRUE, ...)
{
  mc <- match.call()
  environment(formula) <- parent.frame()
  mc$formula <- formula
  fomul_vars <- all.vars(terms(mc$formula))
 # fomul_vars <- all.vars(terms(as.formula(mc$formula)))
  sm_vars <- all.vars(as.list(mc)$smoothZ)
  sm_vars <- intersect(sm_vars, names(data))
  if (any(unlist(lapply(data[, setdiff(sm_vars, fomul_vars)], is.factor)))) stop("Any factor in 'smoothZ' list must be in the formula!")
  response_n <- fomul_vars[1] 
  fomul_vars <- fomul_vars[-1]
  
  data <- droplevels(as.data.frame(na.omit(data[, intersect(all.vars(mc), names(data))])))
  
  if (resp_scale) numeric_n <- names(data)[sapply(data, is.numeric)] else numeric_n <- setdiff(names(data)[sapply(data, is.numeric)], response_n)
  if (any(ncenter, nscale)) {
    for (i in numeric_n) {
      data_i <- data[,i]
      if (!is.matrix(data_i)){  
	    scaled_i <- scale(data[i], center=ncenter, scale=nscale)
        data[i] <- as.numeric(scaled_i)
        if (ncenter) attr(data, paste(i, "mean", sep="_")) <- attr(scaled_i,"scaled:center")
        if (nscale) attr(data, paste(i, "sd", sep="_")) <- attr(scaled_i,"scaled:scale")
      }
    }
	attr(data, "numeric_var") <- numeric_n
  }
  smoothZm <-  eval(substitute(smoothZ), data, parent.frame())
  Z_namem <- names(smoothZm)
  
  # if (any(!sapply(smoothZm, is, class2 = "sparseMatrix"))) {
  # smoothZm <- lapply(smoothZm, function(x) {
  # attr(x, "class") <- "matrix"
  # x_N_attr <- attributes(x)
  # x_N <- as(x, "sparseMatrix")
  
  # for (i in setdiff(names(x_N_attr), c("dim", "class"))) attr(x_N,i) <- x_N_attr[[i]]
  # if (!is.null(attr(x_N,"Boundary.knots"))) {
  # attr(x_N,"range.x") <- attr(x_N,"Boundary.knots")
  # attr(x_N,"Boundary.knots") <- NULL
  # }
  # x_N
  # })
  # }
  
  smoothZN <- list()
  for (i in Z_namem){
    if (is.list(smoothZm[[i]])) smoothZN <- c(smoothZN, unlist(smoothZm[i], recursive=FALSE)) else smoothZN <- c(smoothZN, smoothZm[i])
  }
  smoothZt <- lapply(smoothZN, t)
  
  gaus <- FALSE
  if (is.null(family)) {
    gaus <- TRUE
  }else{
    ## copied from glm()
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    if (!inherits(family, "family")) stop("unknown family type")
    gaus <- family$family == "gaussian" && family$link == "identity"
  }  
  
  lmerc <- mc                         # create a call to lmer
  
  lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
  if (ngbinomial) {
    lmerc[[1]] <- as.name("glmer")
    lmerc$family <- poisson()
  }
  lmerc$smoothZ <- NULL
  lmerc$ncenter <- NULL
  lmerc$nscale <- NULL
  lmerc$resp_scale <- NULL
  lmerc$ngbinomial <- NULL
  lmerc$prt <- NULL
  lmerc$predict_info <- NULL
  
  if (!gaus) lmerc$REML <- NULL
  
  if (!length(smoothZt))              # call [g]lmer instead
    return(eval.parent(lmerc))
  Zt_name <- names(smoothZt)	
  
  stopifnot(is.list(smoothZt),        # check the smoothZt argument
            length(Zt_name) == length(smoothZt),
            all(sapply(smoothZt, is, class2 = "sparseMatrix")))
  
  # adding the constructed variables to the model frame avoiding name duplication		  
  for (i in Zt_name) { 
    data[[i]] <- factor(rep(1:nrow(smoothZt[[i]]), length=ncol(smoothZt[[i]])))
    lmerc$formula <- paste(paste(format(lmerc$formula), collapse = ""), "+ (1|",i,")")
  }
  
  lmerc$formula <- as.formula(lmerc$formula)
  assign(as.character(lmerc$data), data)
  lmf <- eval(lmerc)
  
  pnms <- names(smoothZt)
  pp <- lmf@pp
  resp <- lmf@resp
  fl <- lmf@flist
  stopifnot(all(pnms %in% names(fl)))
  asgn <- attr(fl, "assign")
  Zt <- pp$Zt
  for (i in seq_along(smoothZt)) {
    tn <- which(match(pnms[i], names(fl)) == asgn)
    if (length(tn) > 1)
      stop("a smoothZt factor must be associated with only one r.e. term")
    ind <- (lmf@Gp)[tn:(tn+1L)]
    rowsi <- (ind[1]+1L):ind[2]
    stopifnot(all(dim(Zt[rowsi,])==dim(smoothZt[[i]])))
    Zt[rowsi,] <- smoothZt[[i]]
  }
  reTrms <- list(Zt=Zt,theta=lmf@theta,Lambdat=pp$Lambdat,Lind=pp$Lind,
                 lower=lmf@lower,flist=lmf@flist,cnms=lmf@cnms, Gp=lmf@Gp)
  dfl <- list(fr=lmf@frame, X=pp$X, reTrms=reTrms, start=lmf@theta, drop.unused.levels = drop.unused.levels)
  if (ngbinomial) {
    dfl$family <- poisson()
    devfun <- do.call(mkGlmerDevfun,dfl)
    opt <- optimizeGlmer(devfun,start=start,verbose=verbose,optimizer=control$optimizer, control=control$optCtrl, ...)
    devfun <- updateGlmerDevfun(devfun, dfl$reTrms)
    opt <- optimizeGlmer(devfun, stage=2,start=start,verbose=verbose,optimizer=control$optimizer, control=control$optCtrl, ...)
    semer1 <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)
    
    Y <- model.response(model.frame(semer1))
    mu <- na.omit(fitted(semer1))
    theta_v <- MASS::theta.ml(Y, mu, weights = semer1@resp$weights, limit = 20)
    
    dfl$family <- bquote(MASS::negative.binomial(theta = .(theta_v)))
    devfun <- do.call(mkGlmerDevfun,dfl)
    opt <- optimizeGlmer(devfun,start=start,verbose=verbose,optimizer=control$optimizer, control=control$optCtrl, ...)
    devfun <- updateGlmerDevfun(devfun, dfl$reTrms)
    opt <- optimizeGlmer(devfun, stage=2,start=start,verbose=verbose,optimizer=control$optimizer, control=control$optCtrl, ...)
    semer <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)
    
  }else{
    if (gaus) {
      dfl$REML = resp$REML > 0L
      devfun <- do.call(mkLmerDevfun,dfl)
      opt <- optimizeLmer(devfun, start=start, verbose=verbose, optimizer=control$optimizer, control=control$optCtrl, ...) 
      semer0 <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)
      
      orig_call <- lmerc
      args <- as.list(lmerc)
      args$devFunOnly <- TRUE
      if(!"control" %in% names(as.list(lmerc)))
        args$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
      Call <- as.call(c(list(quote(lme4::lmer)), args[-1]))
      
      options(warn = -1)
      assign(as.character(Call$data), data)
      devfun <- suppressMessages(eval(Call))
      semer <- as_lmerModLT(semer0, devfun)
      # Restore the right 'call' in model:
      semer@call <- orig_call
      options(warn = 0)
    } else {
      dfl$family <- family
      devfun <- do.call(mkGlmerDevfun,dfl)
      opt <- optimizeGlmer(devfun,start=start,verbose=verbose,optimizer=control$optimizer, control=control$optCtrl, ...)
      devfun <- updateGlmerDevfun(devfun, dfl$reTrms)
      opt <- optimizeGlmer(devfun, stage=2,start=start,verbose=verbose,optimizer=control$optimizer, control=control$optCtrl, ...)
      semer <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc) 
    }
  }
  # Gaussian fits keep what a refit needs, so that update() (used by
  # lmerTest::ranova) refits through semireg() with the same smooth bases.
  if (methods::is(semer, "lmerModLmerTest")) {
    semer <- methods::new("semerLmerTest", semer,
                          semireg_info = list(data = data, smoothZ = smoothZN,
                                              REML = isTRUE(REML), control = control))
  }
  if (prt) print(semer)
  
  if (predict_info){
  n_beta <- getME(semer, "p")
  C <- cbind(getME(semer, "X"), getME(semer, "Z"))

  if (gaus) {
    CTC <- Matrix::tcrossprod(t(C))	
  }else{
 	family_n <- semer@resp$family$family
	link_n <- semer@resp$family$link 
    yobs <- unique(semer@resp$y)
     if (family_n =="binomial" && all(yobs%in%c(0,1))){
      wVec <- switch(link_n,
                     "logit" = 3/(pi^2),
                     "probit" = 1,
                     "cloglog" = 6/(pi^2))	
    }else wVec <- semer@resp$sqrtWrkWt()^2 
  
    CTC <- crossprod(C*wVec,C)
  }
  
  beta_cov <- Matrix(0, n_beta, n_beta)
  #  b_cov <- Matrix::solve(Matrix::tcrossprod(getME(semer, "Lambda")))
  
  b_cov_inv <- Matrix::tcrossprod(getME(semer, "Lambda")) 
  if (any(eigen(b_cov_inv)$values <= 0)) {
    Amatrixpd <- nearPD(b_cov_inv)$mat
    b_cov_inv <- Amatrixpd@x
    dim(b_cov_inv) <- Amatrixpd@Dim
  }
  b_cov <- Matrix::solve(b_cov_inv)
  
  # b_cov <- Matrix::tcrossprod(getME(semer, "Lambda"))*(sigma(semer)^2)
  
  D <- .bdiag(list(beta_cov, b_cov))
  fullCovMat <- Matrix::solve(CTC + D)
  # Compute the degrees of freedom  formula
  df <- sum(diag(fullCovMat%*%CTC)) 		
  
  u_ind <- sapply(getME(semer, "Ztlist"), function(x) x@Dim[1])
  u_indn <- sub("\\.\\(Intercept\\)$", "",names(u_ind))
  vcov_ind <- cumsum(c(n_beta, u_ind))
  names(vcov_ind) <- c("beta", u_indn)
  Zt_name <- intersect(u_indn, Zt_name)
  uHat <- unlist(ranef(semer))
  uHat_ind <- vcov_ind-n_beta
  cov_lst <- u_lst <- vector("list", length(Zt_name))
  names(cov_lst) <- names(u_lst) <- Zt_name
  for (i in Zt_name) {
    mch_id <- match(i, names(vcov_ind))
    contInds <- c(1:n_beta, (vcov_ind[mch_id-1]+1):vcov_ind[mch_id])
    cov_lst[[i]] <- fullCovMat[contInds, contInds]
    u_lst[[i]] <- uHat[(uHat_ind[mch_id-1]+1):uHat_ind[mch_id]]  
  }
  vcov_indN <- vcov_ind
  names(vcov_indN) <- gsub("\\..*","",names(vcov_ind))
  vcov_indN <- rev(vcov_indN)[unique(names(vcov_indN))]
  uHat_indN <- vcov_indN-n_beta
  
  for (i in setdiff(names(vcov_indN), names(vcov_ind))) {
    mch_idN <- match(i, names(vcov_indN))
    contIndsN <- c(1:n_beta, (vcov_indN[mch_idN-1]+1):vcov_indN[mch_idN])
    cov_lst[[i]] <- fullCovMat[contIndsN, contIndsN]
    #  u_lst[[i]] <- uHat[(uHat_indN[mch_idN-1]+1):uHat_indN[mch_idN]]  
  }  
  
  sm_contInds <- 1:n_beta
  for (i in Zt_name) {
    mch_id <- match(i, names(vcov_ind))
    sm_contInds <- c(sm_contInds, (vcov_ind[mch_id-1]+1):vcov_ind[mch_id])
  }  
  cov_lst$sm_Cov <- fullCovMat[sm_contInds, sm_contInds]
  
  Z_namem <- intersect(names(vcov_indN), Z_namem)
  knots_lst <- range_lst <- type_lst <- vector("list", length(Z_namem))
  names(knots_lst) <- names(range_lst) <- names(type_lst) <- Z_namem
  for (i in Z_namem) {
    knots_lst[[i]] <- attr(smoothZm[[i]], "knots")
    range_lst[[i]] <- attr(smoothZm[[i]], "range.x")
    type_lst[[i]] <- attr(smoothZm[[i]], "type")
  }
  }else knots_lst <- range_lst <- type_lst <- cov_lst <- u_lst <- fullCovMat <- vcov_ind <- vcov_indN <- df <- NULL
  
  ans <- list(semer=semer, data=data, fomul_vars=fomul_vars, sm_vars=sm_vars, smoothZ_call=mc$smoothZ, knots_lst=knots_lst, range_lst=range_lst, type_lst=type_lst, cov_lst=cov_lst, u_lst=u_lst, CovMat=fullCovMat, Cov_ind=vcov_ind, Cov_indN=vcov_indN, df=df, lmerc=lmerc)
  class(ans) <- c("semireg", "list")
  return(ans)
}




## A semireg() lmer fit that remembers its smooth bases (see update below).
methods::setClass("semerLmerTest", contains = "lmerModLmerTest",
                  slots = c(semireg_info = "list"))

## Refit of a semireg() fit with a new formula (as lmerTest::ranova() does):
## the smooth terms (1 | name) still in the formula keep their spline bases,
## instead of being refitted as ordinary random intercepts.
update.semerLmerTest <- function(object, formula = NULL, ...) {
  info <- object@semireg_info
  if (is.null(formula) || length(list(...))) {
    upd <- utils::getS3method("update", "lmerModLmerTest", envir = asNamespace("lmerTest"))
    return(upd(object, formula = formula, ...))
  }
  f <- stats::as.formula(paste(deparse(formula, width.cutoff = 500L), collapse = " "))
  f <- stats::update(stats::formula(object), f)
  bars <- pm_findbars(f)
  bar_grp <- vapply(bars, function(b) paste(all.vars(b[[3]]), collapse = ":"), character(1))
  sm_names <- names(info$smoothZ)
  keep <- intersect(sm_names, bar_grp)
  other <- bars[!bar_grp %in% sm_names]
  fixed <- if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(f) else lme4::nobars(f)
  rhs <- paste(c(paste(deparse(fixed[[3]], width.cutoff = 500L), collapse = " "),
                 vapply(other, function(b) paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"),
                        character(1))), collapse = " + ")
  f_new <- stats::as.formula(paste(paste(deparse(f[[2]]), collapse = ""), "~", rhs),
                             env = environment(stats::formula(object)))
  e <- new.env(parent = asNamespace("predictmeans"))
  e$sr_data <- info$data
  e$sr_Z <- info$smoothZ[keep]
  cl <- call("semireg", formula = f_new, data = quote(sr_data), smoothZ = quote(sr_Z),
             REML = info$REML, ncenter = FALSE, nscale = FALSE, prt = FALSE, predict_info = FALSE)
  fit <- suppressMessages(suppressWarnings(eval(cl, e)))
  if (is.list(fit) && !is.null(fit$semer)) fit$semer else fit
}
