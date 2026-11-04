semireg_tmb <- function(formula, data, family = gaussian(), smoothZ = list(), ziformula = ~0, dispformula = ~1, 
                        weights = NULL, offset = NULL, contrasts = NULL, na.action, se = TRUE, 
                        verbose = FALSE, doFit = TRUE, control = glmmTMBControl(), REML = TRUE, 
                        start = NULL, map = NULL, sparseX = NULL, prt=TRUE, predict_info=TRUE)
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
  
  smoothZm <-  eval(substitute(smoothZ), data) 
  Z_namem <- names(smoothZm)
  smoothZN <- list()
  for (i in Z_namem){
    if (is.list(smoothZm[[i]])) smoothZN <- c(smoothZN, unlist(smoothZm[i], recursive=FALSE)) else smoothZN <- c(smoothZN, smoothZm[i])
  }
  
  smoothZt <- lapply(smoothZN, t)
  
  glmmTMBcall <- mc 
  
  glmmTMBcall[[1]] <- quote(glmmTMB::glmmTMB)
  glmmTMBcall$smoothZ <- NULL
  glmmTMBcall$prt <- NULL
  glmmTMBcall$predict_info <- NULL
  
  if (!length(smoothZt))             
    return(eval.parent(glmmTMBcall))
  Zt_name <- names(smoothZt)
  
  stopifnot(is.list(smoothZt),        # check the smoothZt argument
            length(Zt_name) == length(smoothZt),
            all(sapply(smoothZt, is, class2 = "sparseMatrix")))
  
  # adding the constructed variables to the semer frame avoiding name duplication		  
  for (i in Zt_name) { 
    data[[i]] <- factor(rep(1:nrow(smoothZt[[i]]), length=ncol(smoothZt[[i]])))
    glmmTMBcall$formula <- paste(paste(format(glmmTMBcall$formula), collapse = ""), "+ (1|",i,")")
  }
  
  glmmTMBcall$formula <- as.formula(glmmTMBcall$formula)
  assign(as.character(glmmTMBcall$data), data)
  glmmTMBcall$doFit <- FALSE
  tmbf <- eval(glmmTMBcall)
  
  cond_reTrms <- tmbf$condList$reTrms
  fl <- cond_reTrms$flist
  stopifnot(all(Zt_name %in% names(fl)))
  asgn <- attr(fl, "assign")
  
  # Zt_list <- cond_reTrms$Ztlist
  Zt <- cond_reTrms$Zt
  for (i in seq_along(smoothZt)) {
    tn <- which(match(Zt_name[i], names(fl)) == asgn)
    if (length(tn) > 1)
      stop("a smoothZt factor must be associated with only one r.e. term")
    ind <- (cond_reTrms$Gp)[tn:(tn+1L)]
    rowsi <- (ind[1]+1L):ind[2]
    stopifnot(all(dim(Zt[rowsi,])==dim(smoothZt[[i]])))
    Zt[rowsi,] <- smoothZt[[i]]
  }
  
  tmbf$data.tmb$Z <- t(Zt)
  semer <- fitTMB(tmbf)  
  if (prt) print(semer)  
  
  if (predict_info){
    beta <- fixef(semer)$cond
    n_beta <- length(beta)
    s <- sigma(semer) 
    
    cond_reTrms <- tmbf$condList$reTrms
    if (is.null(cond_reTrms$Lambdat) || is.null(cond_reTrms$Lind) || is.null(cond_reTrms$lower)) {
      # glmmTMB >= 1.1.12 no longer stores Lambdat, Lind and lower; rebuild
      # them (they depend only on the terms and factor levels) in glmmTMB's
      # own term order
      mk <- if (requireNamespace("reformulas", quietly = TRUE)) reformulas::mkReTrms else lme4::mkReTrms
      rt <- mk(pm_findbars(glmmTMBcall$formula), tmbf$fr, reorder.terms = FALSE)
      cond_reTrms[c("Lambdat", "Lind", "lower")] <- rt[c("Lambdat", "Lind", "lower")]
    }
    n_diag <- sum(cond_reTrms$lower==0)
    theta_ini <- getME(semer, "theta")
    if (length(theta_ini)!=n_diag) {
      theta <- cond_reTrms$lower
      theta[theta==0] <- exp(theta_ini[1:n_diag])/s
      theta[theta!=0] <- theta_ini[-1:-n_diag]/s
    }else theta <- exp(theta_ini)/s
    
    Lambdat <- cond_reTrms$Lambdat
    Lambdat@x <- theta[cond_reTrms$Lind]
    b_cov_inv <- Matrix::tcrossprod(t(Lambdat)) 
    
    if (any(eigen(b_cov_inv)$values <= 0)) {
      Amatrixpd <- nearPD(b_cov_inv)$mat
      b_cov_inv <- Amatrixpd@x
      dim(b_cov_inv) <- Amatrixpd@Dim
    }
    b_cov <- Matrix::solve(b_cov_inv)
    beta_cov <- Matrix(0, n_beta, n_beta)
    D <- .bdiag(list(beta_cov, b_cov))
    
    weights <- tmbf$data.tmb$weights
    family_n <- semer$modelInfo$family$family
    if (family_n=="binomial") yobs <- unique(tmbf$data.tmb$yobs)
    link_n <- semer$modelInfo$family$link
    uHat <- unlist(ranef(semer)$cond) # cbind(uHat, getME(semer, "b"))
    C <- cbind(getME(semer, "X"), getME(semer, "Z"))
    
    if (family_n != "gaussian"){
      if (family_n =="binomial" && all(yobs%in%c(0,1))){
        wVec <- switch(link_n,
                       "logit" = 3/(pi^2),
                       "probit" = 1,
                       "cloglog" = 6/(pi^2))	
      }else{
        if (inherits(semer, "glmmTMB")) {
          weight_fun <- function(eta, mod){
            family_fun <- mod$modelInfo$family
            family_fun$mu.eta(eta)^2/family_fun$variance(family_fun$linkinv(eta))
          }
          y_hat <- as.vector(C%*%c(beta, uHat))
          wVec <- as.vector(weight_fun(y_hat, semer))*tmbf$data.tmb$weights 
          if (length(tmbf$data.tmb$size) > 0) wVec <- wVec*tmbf$data.tmb$size
        }else wVec <- semer@resp$sqrtWrkWt()^2 
      }
    }else wVec <- weights 
    
    if (isLMM(semer)) { 
      CTC <- Matrix::tcrossprod(t(C))	
    }else{
      CTC <- crossprod(C*wVec,C)
    }
    
    fullCovMat <- Matrix::solve(CTC + D)
    # Compute the degrees of freedom  formula
    df <- sum(diag(fullCovMat%*%CTC)) 
    
    u_ind <- cond_reTrms$nl
    u_indn <- names(u_ind) 
    vcov_ind <- cumsum(c(n_beta, u_ind))
    names(vcov_ind) <- c("beta", u_indn)  
    
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
  
  ans <- list(semer=semer, data=data, fomul_vars=fomul_vars, sm_vars=sm_vars, smoothZ_call=mc$smoothZ, knots_lst=knots_lst, range_lst=range_lst, type_lst=type_lst, cov_lst=cov_lst, u_lst=u_lst, CovMat=fullCovMat, Cov_ind=vcov_ind, Cov_indN=vcov_indN, df=df, tmbf=tmbf)
  class(ans) <- c("semireg", "list")
  return(ans) 
}

