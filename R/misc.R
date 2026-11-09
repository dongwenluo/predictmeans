adiag <- function (..., pad = as.integer(0), do.dimnames = TRUE) # function from package 'magic'
{
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1], list(pad = pad)))
    return(do.call("Recall", c(list(args[[1]]), list(jj),
                               list(pad = pad))))
  }
  a <- args[[1]]
  b <- args[[2]]
  if (is.null(b)) {
    return(a)
  }
  if (is.null(dim(a)) & is.null(dim(b))) {
    dim(a) <- rep(1, 2)
    dim(b) <- rep(1, 2)
  }
  if (is.null(dim(a)) & length(a) == 1) {
    dim(a) <- rep(1, length(dim(b)))
  }
  if (is.null(dim(b)) & length(b) == 1) {
    dim(b) <- rep(1, length(dim(a)))
  }
  if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
    stop("a and b must have identical number of dimensions")
  }
  s <- array(pad, dim.a + dim.b)
  s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
  ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) +
                  dim.a[[i]])
  out <- do.call("[<-", c(list(s), ind, list(b)))
  n.a <- dimnames(a)
  n.b <- dimnames(b)
  if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
    dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
    names(dimnames(out)) <- names(n.a)
  }
  return(out)
}

vec2mat2 <- function (x, sep = "-")
{
  splits <- strsplit(x, sep)
  n.spl <- sapply(splits, length)
  if (any(n.spl != 2)) {
    stop("Names must contain exactly one '", sep, "' each;  instead got ",
         paste(x, collapse = ", "))
  }
  x2 <- t(as.matrix(as.data.frame(splits)))
  dimnames(x2) <- list(x, NULL)
  x2
}

multcompLetters <- function (x, compare = "<", threshold = 0.05,   # function from package 'multcompView'
                             Letters = c(letters, LETTERS, "."), reversed = FALSE)
{
  x.is <- deparse(substitute(x))
  if (inherits(x, "dist")) {
    x <- as.matrix(x)
  }
  if (!is.logical(x)) {
    x <- do.call(compare, list(x, threshold))
  }
  dimx <- dim(x)
  {
    if ((length(dimx) == 2) && (dimx[1] == dimx[2])) {
      Lvls <- dimnames(x)[[1]]
      if (length(Lvls) != dimx[1]) {
        stop("Names requred for ", x.is)
      } else {
        x2. <- t(outer(Lvls, Lvls, paste, sep = ""))
        x2.n <- outer(Lvls, Lvls, function(x1, x2) nchar(x2))
        x2.2 <- x2.[lower.tri(x2.)]
        x2.2n <- x2.n[lower.tri(x2.n)]
        x2a <- substring(x2.2, 1, x2.2n)
        x2b <- substring(x2.2, x2.2n + 1)
        x2 <- cbind(x2a, x2b)
        x <- x[lower.tri(x)]
      }
    } else {
      namx <- names(x)
      if (length(namx) != length(x)) {
        stop("Names required for ", x.is)
      }
      x2 <- vec2mat2(namx)
      Lvls <- unique(as.vector(x2))
    }
  }
  n <- length(Lvls)
  LetMat <- array(TRUE, dim = c(n, 1), dimnames = list(Lvls, NULL))
  k2 <- sum(x)
  if (k2 == 0) {
    Ltrs <- rep(Letters[1], n)
    names(Ltrs) <- Lvls
    dimnames(LetMat)[[2]] <- Letters[1]
    return(Ltrs)
  }
  distinct.pairs <- x2[x, , drop = FALSE]
  absorb <- function(A.) {
    k. <- dim(A.)[2]
    if (k. > 1) {
      for (i. in 1:(k. - 1)) for (j. in (i. + 1):k.) {
        if (all(A.[A.[, j.], i.])) {
          A. <- A.[, -j., drop = FALSE]
          return(absorb(A.))
        }
        else {
          if (all(A.[A.[, i.], j.])) {
            A. <- A.[, -i., drop = FALSE]
            return(absorb(A.))
          }
        }
      }
    }
    A.
  }
  for (i in 1:k2) {
    dpi <- distinct.pairs[i, ]
    ijCols <- (LetMat[dpi[1], ] & LetMat[dpi[2], ])
    if (any(ijCols)) {
      A1 <- LetMat[, ijCols, drop = FALSE]
      A1[dpi[1], ] <- FALSE
      LetMat[dpi[2], ijCols] <- FALSE
      LetMat <- cbind(LetMat, A1)
      LetMat <- absorb(LetMat)
    }
  }
  sortCols <- function(B) {
    firstRow <- apply(B, 2, function(x) which(x)[1])
    B <- B[, order(firstRow)]
    firstRow <- apply(B, 2, function(x) which(x)[1])
    reps <- (diff(firstRow) == 0)
    if (any(reps)) {
      nrep <- table(which(reps))
      irep <- as.numeric(names(nrep))
      k <- dim(B)[1]
      for (i in irep) {
        i. <- i:(i + nrep[as.character(i)])
        j. <- (firstRow[i] + 1):k
        B[j., i.] <- sortCols(B[j., i., drop = FALSE])
      }
    }
    B
  }
  LetMat. <- sortCols(LetMat)
  if (reversed) {
    LetMat. <- LetMat.[, rev(1:ncol(LetMat.))]
  }
  k.ltrs <- dim(LetMat.)[2]
  makeLtrs <- function(kl, ltrs = Letters) {
    kL <- length(ltrs)
    if (kl < kL) {
      return(ltrs[1:kl])
    }
    ltrecurse <- c(paste(ltrs[kL], ltrs[-kL], sep = ""),
                   ltrs[kL])
    c(ltrs[-kL], makeLtrs(kl - kL + 1, ltrecurse))
  }
  Ltrs <- makeLtrs(k.ltrs, Letters)
  dimnames(LetMat.)[[2]] <- Ltrs
  LetVec <- rep(NA, n)
  names(LetVec) <- Lvls
  for (i in 1:n) LetVec[i] <- paste(Ltrs[LetMat.[i, ]], collapse = "")
  nch.L <- nchar(Ltrs)
  blk.L <- rep(NA, k.ltrs)
  for (i in 1:k.ltrs) blk.L[i] <- paste(rep(" ", nch.L[i]),
                                        collapse = "")
  monoVec <- rep(NA, n)
  names(monoVec) <- Lvls
  for (j in 1:n) {
    ch2 <- blk.L
    if (any(LetMat.[j, ])) {
      ch2[LetMat.[j, ]] <- Ltrs[LetMat.[j, ]]
    }
    monoVec[j] <- paste(ch2, collapse = "")
  }
  return(monoVec)
}

######################## Functions from lmerTest begin #################

get_contrasts_type1 <- function(model) {
  terms <- terms(model)
  X <- model.matrix(model)
  p <- ncol(X)
  if(p == 0L) {
    return(list(matrix(numeric(0L), nrow=0L))) # no fixef
  }
  if(p == 1L && attr(terms, "intercept")) { # intercept-only model
    return(list(matrix(numeric(0L), ncol=1L)))
  }
  # Compute 'normalized' doolittle factorization of XtX:
  L <- if(p == 1L) {
    matrix(1L)
  } else {
    t(doolittle(crossprod(X))$L)
  }
  dimnames(L) <- list(colnames(X), colnames(X))
  # Determine which rows of L belong to which term:
  ind.list <- term2colX(terms, X)[attr(terms, "term.labels")]
  lapply(ind.list, function(rows) L[rows, , drop=FALSE])
}

term2colX <- function(terms, X) {
  # Compute map from terms to columns in X using the assign attribute of X.
  # Returns a list with one element for each term containing indices of columns
  #   in X belonging to that term.
  if(is.null(asgn <- attr(X, "assign"))) {
    stop("Invalid design matrix:",
         "design matrix 'X' should have a non-null 'assign' attribute",
         call. = FALSE)
  }
  term_names <- attr(terms, "term.labels")
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if(has_intercept) {
    c("(Intercept)", term_names)[asgn + 1]
  } else {
    term_names[asgn[asgn > 0]]
  }
  if(!length(col_terms) == ncol(X)) {
    # should never happen.
    stop("An error happended when mapping terms to columns of X")
  }
  # get names of terms (including aliased terms)
  nm <- union(unique(col_terms), term_names)
  res <- lapply(setNames(as.list(nm), nm), function(x) numeric(0L))
  map <- split(seq_along(col_terms), col_terms)
  res[names(map)] <- map
  res[nm] # order appropriately
}

doolittle <- function(x, eps = 1e-6) {
  if(!is.matrix(x) || ncol(x) != nrow(x) || !is.numeric(x)){
    stop("argument 'x' should be a numeric square matrix")
  }
  stopifnot(ncol(x) > 1L)
  n <- nrow(x)
  L <- U <- matrix(0, nrow=n, ncol=n)
  diag(L) <- rep(1, n)
  for(i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for(j in 1:n) {
      U[i,j] <- x[i,j]
      if (im1 > 0) {
        for(k in 1:im1) {
          U[i,j] <- U[i,j] - L[i,k] * U[k,j]
        }
      }
    }
    if ( ip1 <= n ) {
      for ( j in ip1:n ) {
        L[j,i] <- x[j,i]
        if ( im1 > 0 ) {
          for ( k in 1:im1 ) {
            L[j,i] <- L[j,i] - L[j,k] * U[k,i]
          }
        }
        L[j, i] <- if(abs(U[i, i]) < eps){
          0
          } else {
            L[j,i] / U[i,i]
          }
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list( L=L, U=U )
}

# df_term <- function (model, term) {
# if(!getME(model, "is_REML"))
# stop("Kenward-Roger's method is only available for REML model fits")
# if(!requireNamespace("pbkrtest", quietly = TRUE))
# stop("pbkrtest package required for Kenward-Roger's method")
# Lc <- get_contrasts_type1(model)[[term]]
# Lc <- Kmatrix(model, term)$K
# vcov_beta_adj <- try(pbkrtest::vcovAdj(model), silent=TRUE) # Adjusted vcov(beta)
# ddf <- try(pbkrtest::Lb_ddf(L=Lc, V0=vcov(model),
# Vadj=vcov_beta_adj), silent=TRUE) # vcov_beta_adj need to be dgeMatrix!
# if(any(inherits(vcov_beta_adj, "try-error"), inherits(ddf, "try-error"))) {
# warning("Unable to compute Kenward-Roger Df: using Satterthwaite instead")
# if(!inherits(model, "lmerModLmerTest")) model <- as_lmerModLmerTest(model)
# beta <- model@beta
# # Compute Var(L beta) and eigen-decompose:
# VLbeta <- Lc %*% model@vcov_beta %*% t(Lc) # Var(contrast) = Var(Lbeta)
# eig_VLbeta <- eigen(VLbeta)
# P <- eig_VLbeta$vectors
# d <- eig_VLbeta$values
# tol <- max(sqrt(.Machine$double.eps) * d[1], 0)
# pos <- d > tol
# q <- sum(pos) # rank(VLbeta)

# PtL <- crossprod(P, Lc)[1:q, ]
# if(q == 1) { # 1D case:
# ddf <- contest1D(model, PtL, rhs=0, confint=FALSE)$df
# return(ddf)
# } # multi-D case proceeds:

# # Compute q-list of gradients of (PtL)' cov(beta) (PtL) wrt. varpar vector:
# grad_PLcov <- lapply(1:q, function(m) {
# vapply(model@Jac_list, function(J) qform(PtL[m, ], J), numeric(1L))
# })
# # Compute degrees of freedom for the q t-statistics:
# nu_m <- vapply(1:q, function(m) {
# 2*(d[m])^2 / qform(grad_PLcov[[m]], model@vcov_varpar) }, numeric(1L)) # 2D_m^2 / g'Ag
# # Compute ddf for the F-value:
# ddf <- get_Fstat_ddf(nu_m, tol=1e-8)
# }
# return(ddf)
# }

# ##########
# get_Fstat_ddf <- function(nu, tol=1e-8) {
# fun <- function(nu) {
# if(any(nu <= 2)) 2 else {
# E <- sum(nu / (nu - 2))
# 2 * E / (E - (length(nu))) # q = length(nu) : number of t-statistics
# }
# }
# stopifnot(length(nu) >= 1,
# # all(nu > 0), # returns 2 if any(nu < 2)
# all(sapply(nu, is.numeric)))
# if(length(nu) == 1L) return(nu)
# if(all(abs(diff(nu)) < tol)) return(mean(nu))
# if(!is.list(nu)) fun(nu) else vapply(nu, fun, numeric(1L))
# }

########################

df_term <- function(model, modelterm, covariate=NULL, ctrmatrix=NULL, ctrnames=NULL, type=c("Kenward-Roger", "Satterthwaite", "containment"), vcov_adj=NULL) {

  stopifnot(inherits(model, "lmerMod"))

  type <- as.character(type)
  type <- match.arg(type)
  if (type == "Kenward-Roger" && !getME(model, "is_REML")) {
    # Kenward-Roger is defined for REML fits; Satterthwaite also works for ML.
    warning("Kenward-Roger df need a REML fit: using Satterthwaite instead")
    type <- "Satterthwaite"
  }

  if (type == "containment") {
    # one df per term (SAS DDFM=CONTAIN): every row gets the df of 'modelterm'
    if (missing(modelterm)) stop("containment df need 'modelterm'")
    cdf <- pm_containment_df(model)
    if (!modelterm %in% names(cdf)) stop("'", modelterm, "' is not a term of the model")
    n_row <- if (is.null(ctrmatrix)) nrow(Kmatrix(model, modelterm, covariate)$K) else NROW(rbind(ctrmatrix))
    return(rep(unname(cdf[modelterm]), n_row))
  }

  if (!is.null(ctrmatrix)) {
    if (is.vector(ctrmatrix)) {
      Lc <- matrix(ctrmatrix, nrow=1)
    } else {
      Lc <- ctrmatrix
    }
    stopifnot(is.numeric(Lc), ncol(Lc)==length(fixef(model)))

    if (!is.null(ctrnames)) {
      rownames(Lc) <- ctrnames
    }
  } else {
    Lc <- Kmatrix(model, modelterm, covariate)$K
  }

  if (type=="Kenward-Roger") {
    # Adjusted vcov(beta); reuse it when the caller has already computed it.
    vcov_beta_adj <- if (is.null(vcov_adj) || inherits(vcov_adj, "pm_satt")) try(pm_vcov_kr(model), silent=TRUE) else vcov_adj
    # All rows at once; same values as pbkrtest::Lb_ddf() row by row.
    ddf <- try(pm_kr_ddf(Lc, V0=vcov(model), Vadj=vcov_beta_adj), silent=TRUE)

    if (any(inherits(vcov_beta_adj, "try-error"),
	        inherits(ddf, "try-error"),
			ddf >= nrow(model.frame(model)),
			ddf <= 0)) {
      warning("Unable to compute Kenward-Roger Df: using Satterthwaite instead")
      type <- "Satterthwaite"
    }
  }

  if (type == "Satterthwaite") {
    # all rows at once; same values as lmerTest. 'vcov_adj' may carry the
    # Satterthwaite ingredients already computed for this model.
    parts <- if (inherits(vcov_adj, "pm_satt")) vcov_adj else pm_satt_parts(model)
    ddf <- pm_satt_ddf(Lc, parts)
  }
  return(ddf)
}

##########
as_lmerModLT <- function(model, devfun, tol=1e-8) {
  is_reml <- getME(model, "is_REML")
  # Coerce 'lme4-model' to 'lmerModLmerTest':
  res <- as(model, "lmerModLmerTest")
  # Set relevant slots of the new model object:
  res@sigma <- sigma(model)
  res@vcov_beta <- as.matrix(vcov(model))
  varpar_opt <- unname(c(res@theta, res@sigma))
  # Compute Hessian:
  h <- numDeriv::hessian(func=devfun_vp, x=varpar_opt, devfun=devfun,
                         reml=is_reml)
  # Eigen decompose the Hessian:
  eig_h <- eigen(h, symmetric=TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  zero <- evals > -tol & evals < tol
  if(sum(neg) > 0) { # negative eigenvalues
    eval_chr <- if(sum(neg) > 1) {
      "eigenvalues"
    } else {
      "eigenvalue"
    }
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf("Model failed to converge with %d negative %s: %s",
                    sum(neg), eval_chr, evals_num), call.=FALSE)
  }
  # Note: we warn about negative AND zero eigenvalues:
  if(sum(zero) > 0) { # some eigenvalues are zero
    eval_chr <- if(sum(zero) > 1) {
      "eigenvalues"
    } else {
      "eigenvalue"
    }
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf("Model may not have converged with %d %s close to zero: %s",
                    sum(zero), eval_chr, evals_num))
  }
  # Compute vcov(varpar):
  pos <- eig_h$values > tol
  q <- sum(pos)
  # Using the Moore-Penrose generalized inverse for h:
  h_inv <- with(eig_h, {
    vectors[, pos, drop=FALSE] %*% diag(1/values[pos], nrow=q) %*%
      t(vectors[, pos, drop=FALSE]) })
  res@vcov_varpar <- 2 * h_inv # vcov(varpar)
  # Compute Jacobian of cov(beta) for each varpar and save in list:
  Jac <- numDeriv::jacobian(func=get_covbeta, x=varpar_opt, devfun=devfun)
  res@Jac_list <- lapply(1:ncol(Jac), function(i)
    array(Jac[, i], dim=rep(length(res@beta), 2))) # k-list of jacobian matrices
  res
}

##########
devfun_vp <- function(varpar, devfun, reml) {
  nvarpar <- length(varpar)
  sigma2 <- varpar[nvarpar]^2
  theta <- varpar[-nvarpar]
  df_envir <- environment(devfun)
  devfun(theta) # Evaluate deviance function at varpar
  n <- nrow(df_envir$pp$V)
  # Compute deviance for ML:
  dev <- df_envir$pp$ldL2() + (df_envir$resp$wrss() + df_envir$pp$sqrL(1))/sigma2 +
    n * log(2 * pi * sigma2)
  if(!reml) {
    return(dev)
  }
  # Adjust if REML is used:
  RX <- df_envir$pp$RX() # X'V^{-1}X ~ crossprod(RX^{-1}) = cov(beta)^{-1} / sigma^2
  dev + 2*c(determinant(RX)$modulus) - ncol(RX) * log(2 * pi * sigma2)
}

##########
get_covbeta <- function(varpar, devfun) {
  nvarpar <- length(varpar)
  sigma <- varpar[nvarpar] # residual std.dev.
  theta <- varpar[-nvarpar] # ranef var-par
  devfun(theta) # evaluate REML or ML deviance 'criterion'
  df_envir <- environment(devfun) # extract model environment
  sigma^2 * tcrossprod(df_envir$pp$RXi()) # vcov(beta)
}

######################## Functions from lmerTest end #################


#######################
# https://rpubs.com/bbolker/waldvar

# waldVar2 <- function(object) {
# ## test for/warn if ML fit?
# dd <- lme4::devfun2(object,useSc=TRUE,signames=FALSE)
# nvp <- length(attr(dd,"thopt"))+1 ## variance parameters (+1 for sigma)
# pars <- attr(dd,"optimum")[seq(nvp)] ## var params come first
# hh <- numDeriv::hessian(dd,pars)
# ## factor of 2: deviance -> negative log-likelihood
# vv <- 2*solve(hh)
# nn <- tn(object)
# dimnames(vv) <- list(nn,nn)
# return(vv)
# }

# tn <- function(object) {
# c(names(getME(object,"theta")),"sigma")
# }

# confintlmer <- function (object, parm, level = 0.95, ...)
# {
# cf <- coef(object)
# pnames <- names(cf)
# if (missing(parm))
# parm <- pnames
# else if (is.numeric(parm))
# parm <- pnames[parm]
# a <- (1 - level)/2
# a <- c(a, 1 - a)
# pct <- format.perc(a, 3)
# fac <- qnorm(a)
# ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
# ses <- sqrt(diag(object$vcov))[parm]
# ci[] <- cf[parm] + ses %o% fac
# ci
# }

# format.perc <- function (probs, digits) {
# paste(format(100 * probs, trim = TRUE, scientific = FALSE,
# digits = digits), "%")
# }

###################### for print
print.pdmlist = function(x, ...){
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  pos = grep('predictmeansPlot|predictmeansBKPlot|predictmeansBarPlot|p_valueMatrix', names(x))
  # the stable tables (means, pairs, settings) are not printed, so the printed
  # output is the same as before; use tidy() or as.data.frame() to see them
  pos = union(pos, which(names(x) %in% c("means", "pairs", "settings")))
  if (length(pos)) x = x[names(x)[-pos]]
  NextMethod()
}
######################
# vcov.VarCorr.merMod <- function(object,fit,...) {
# if (isREML(fit)) {
# warning("refitting model with ML")
# fit <- refitML(fit)
# }
# if (!require("numDeriv")) stop("numDeriv package required")
# useSc <- attr(object,"useSc")
# dd <- lme4:::devfun2(fit,useSc=useSc,signames=FALSE)
# vdd <- as.data.frame(object,order="lower.tri")
# pars <- vdd[,"sdcor"]
# npar0 <- length(pars)
# if (isGLMM(fit)) {
# pars <- c(pars,fixef(fit))
# }
# hh1 <- hessian(dd,pars)
# vv2 <- 2*solve(hh1)
# if (isGLMM(fit)) {
# vv2 <- vv2[1:npar0,1:npar0,drop=FALSE]
# }
# nms <- apply(vdd[,1:3],1,
# function(x) paste(na.omit(x),collapse="."))
# dimnames(vv2) <- list(nms,nms)
# return(vv2)
# }

# http://rstudio-pubs-static.s3.amazonaws.com/28864_dd1f084207d54f5ea67c9d1a9c845d01.html

#######################################################
# from package merDeriv vcov.lmerMod.R

vcov_lmerMod <- function (object, ...) {
  if (!is(object, "lmerMod")) {
    stop("vcov.lmerMod() only works for lmer() models.")
  }
  dotdotdot <- list(...)
  if ("full" %in% names(dotdotdot)) {
    full <- dotdotdot$full
  }
  else {
    full <- FALSE
  }
  if ("information" %in% names(dotdotdot)) {
    information <- dotdotdot$information
  }
  else {
    information <- "expected"
  }
  if (!(full %in% c("TRUE", "FALSE"))) {
    stop("invalid 'full' argument supplied")
  }
  if (!(information %in% c("expected", "observed"))) {
    stop("invalid 'information' argument supplied")
  }
  if ("ranpar" %in% names(dotdotdot)) {
    ranpar <- dotdotdot$ranpar
  }
  else {
    ranpar <- "var"
  }
  parts <- getME(object, "ALL")
  yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
  uluti <- length(parts$theta)
  Zlam <- tcrossprod(parts$Z, parts$Lambdat)
  V <- (tcrossprod(Zlam, Zlam) + Matrix::Diagonal(parts$n, 1)) * (parts$sigma)^2
  M <- solve(chol(V))
  invV <- tcrossprod(M, M)
  LambdaInd <- parts$Lambda
  LambdaInd@x[] <- parts$Lind
  invVX <- crossprod(parts$X, invV)
  Pmid <- solve(crossprod(parts$X, t(invVX)))
  P <- invV - tcrossprod(crossprod(invVX, Pmid), t(invVX))
  fixvar <- solve(tcrossprod(crossprod(parts$X, invV), t(parts$X)))
  if (full == FALSE) {
    fixvar
  }
  else {
    fixhes <- tcrossprod(crossprod(parts$X, invV), t(parts$X))
    uluti <- length(parts$theta)
    devV <- vector("list", (uluti + 1))
    devLambda <- vector("list", uluti)
    score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)
    for (i in 1:uluti) {
      devLambda[[i]] <- Matrix::forceSymmetric(LambdaInd == i,
                                               uplo = "L")
      devV[[i]] <- tcrossprod(tcrossprod(parts$Z, t(devLambda[[i]])),
                              parts$Z)
    }
    devV[[(uluti + 1)]] <- Matrix::Diagonal(nrow(parts$X), 1)
    ranhes <- matrix(NA, nrow = (uluti + 1), ncol = (uluti +
                                                       1))
    entries <- rbind(matrix(rep(1:(uluti + 1), each = 2),
                            (uluti + 1), 2, byrow = TRUE), t(combn((uluti + 1),
                                                                   2)))
    entries <- entries[order(entries[, 1], entries[, 2]),
    ]
    if (parts$devcomp$dims[["REML"]] == 0) {
      if (information == "expected") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries,
                                                        1, function(x) as.numeric((1/2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(invV,
                                                                                                                                           devV[[x[1]]]), invV), t(devV[[x[2]]])))))
      }
      if (information == "observed") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- unlist(apply(entries,
                                                               1, function(x) as.vector(-as.numeric((1/2) *
                                                                                                      lav_matrix_trace(tcrossprod(tcrossprod(crossprod(invV,
                                                                                                                                                       devV[[x[1]]]), invV), t(devV[[x[2]]])))) +
                                                                                          tcrossprod((tcrossprod((crossprod(yXbe, tcrossprod(tcrossprod(crossprod(invV,
                                                                                                                                                                  devV[[x[1]]]), invV), t(devV[[x[2]]])))),
                                                                                                                 invV)), t(yXbe)))))
      }
    }
    if (parts$devcomp$dims[["REML"]] > 0) {
      if (information == "expected") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries,
                                                        1, function(x) as.numeric((1/2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(P,
                                                                                                                                           devV[[x[1]]]), P), t(devV[[x[2]]])))))
      }
      if (information == "observed") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(entries,
                                                        1, function(x) -as.numeric((1/2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(P,
                                                                                                                                            devV[[x[1]]]), P), t(devV[[x[2]]])))) + tcrossprod((tcrossprod((crossprod(yXbe,
                                                                                                                                                                                                                      tcrossprod(tcrossprod(crossprod(invV, devV[[x[1]]]),
                                                                                                                                                                                                                                            P), t(devV[[x[2]]])))), invV)), t(yXbe)))
      }
    }
    ranhes <- Matrix::forceSymmetric(ranhes, uplo = "L")
    if (information == "expected") {
      varcov_beta <- matrix(0, length(devV), length(parts$beta))
    }
    if (information == "observed") {
      varcov_beta <- matrix(NA, length(devV), length(parts$beta))
      for (j in 1:(length(devV))) {
        varcov_beta[j, ] <- as.vector(tcrossprod(crossprod(parts$X,
                                                           (tcrossprod(crossprod(invV, devV[[j]]), invV))),
                                                 t(yXbe)))
      }
    }
    if (ranpar == "var") {
      ranhes <- ranhes
      varcov_beta <- varcov_beta
    }
    else if (ranpar == "sd") {
      sdcormat <- as.data.frame(VarCorr(object, comp = "Std.Dev"),
                                order = "lower.tri")
      sdcormat$sdcor2[which(is.na(sdcormat$var2))] <- sdcormat$sdcor[which(is.na(sdcormat$var2))] *
        2
      sdcormat$sdcor2[which(!is.na(sdcormat$var2))] <- sdcormat$vcov[which(!is.na(sdcormat$var2))]/sdcormat$sdcor[which(!is.na(sdcormat$var2))]
      varcov_beta <- sweep(varcov_beta, MARGIN = 1, sdcormat$sdcor2,
                           `*`)
      weight <- apply(entries, 1, function(x) sdcormat$sdcor2[x[1]] *
                        sdcormat$sdcor2[x[2]])
      ranhes[lower.tri(ranhes, diag = TRUE)] <- weight *
        ranhes[lower.tri(ranhes, diag = TRUE)]
      ranhes <- Matrix::forceSymmetric(ranhes, uplo = "L")
    }
    else {
      stop("ranpar needs to be var or sd for lmerMod object.")
    }
    full_varcov <- solve(rbind(cbind(fixhes, t(varcov_beta)),
                               cbind(varcov_beta, ranhes)))
    colnames(full_varcov) <- c(names(parts$fixef), paste("cov",
                                                         names(parts$theta), sep = "_"), "residual")
    callingFun <- try(deparse(sys.call(-2)), silent = TRUE)
    if (length(callingFun) > 1) {
      callingFun <- paste(callingFun, collapse = "")
    }
    if (!inherits(callingFun, "try-error") & grepl("summary.merMod",
                                                   callingFun)) {
      return(fixvar)
    }
    else {
      return(full_varcov)
    }
  }
}

vcov_glmerMod <- function(object, ...) {
  if (!(family(object)$family %in% c("binomial", "poisson"))) {
    stop("family has to be binomial or poisson")
  }

  dotdotdot <- list(...)
  if("full" %in% names(dotdotdot)){
    full <- dotdotdot$full
  } else {
    full <- FALSE
  }

  if("ranpar" %in% names(dotdotdot)){
    ranpar <- dotdotdot$ranpar
  } else {
    ranpar <- "var"
  }

  if (full == FALSE) {
    full_vcov <- vcov.merMod(object)
  } else {
    if (length(getME(object, "l_i")) > 1L) {
      stop("Multiple cluster variables detected. This type of model is currently not supported for full vcov.")
    }

    ## Hessian was based on deviance function, which is the
    ## -2*LogLik. That's why divided by -2
    if (object@devcomp$dims[['nAGQ']] == 0L) {
      stop("For full vcov, nAGQ of at least 1 is required.")
    }
    full_vcov_noorder <- -solve(object@optinfo$derivs$Hessian/(-2))

    ## Block order in Hessian was theta, beta. Reorganize to
    ## put fixed parameter block first to match with score
    ## matrix order.

    ## count parameter numbers
    p <- nrow(full_vcov_noorder)
    pfix <- length(object@beta)
    pran <- length(object@theta)
    ## reorder four blocks
    full_vcov <- matrix(NA, nrow(full_vcov_noorder), ncol(full_vcov_noorder))
    full_vcov[1:pfix, 1:pfix] <- full_vcov_noorder[(pran + 1):p, (pran + 1):p]
    full_vcov[(pfix + 1):p, (pfix + 1): p] <- full_vcov_noorder[1:pran, 1:pran]
    full_vcov[(pfix + 1):p, 1:pfix] <- full_vcov_noorder[1:pran, (pran + 1): p]
    full_vcov[1:pfix, (pfix + 1): p] <- full_vcov_noorder[(pran + 1): p, 1:pran]


    ## reparameterize for sd and var for random variance/covariance parameters.
    if (ranpar == "theta") {
       full_vcov <- full_vcov
    }

    if (ranpar == "sd") {
        dd <- devfun2(object,useSc=FALSE,signames=TRUE)
        nvp <- length(attr(dd,"thopt"))
        pars <- attr(dd,"optimum")
        pars <- pars[!is.na(names(pars))]
        hh <- hessian(dd, pars)/(-2)

        full_vcov_noorder <- -solve(hh)
        full_vcov <- matrix(NA, nrow(full_vcov_noorder),
          ncol(full_vcov_noorder))
        full_vcov[1:pfix, 1:pfix] <-
          full_vcov_noorder[(pran + 1):p, (pran + 1):p]
        full_vcov[(pfix + 1):p, (pfix + 1): p] <-
          full_vcov_noorder[1:pran, 1:pran]
        full_vcov[(pfix + 1):p, 1:pfix] <-
          full_vcov_noorder[1:pran, (pran + 1): p]
        full_vcov[1:pfix, (pfix + 1): p] <-
            full_vcov_noorder[(pran + 1): p, 1:pran]
     }

    if (ranpar == "var"){
        dd <- devfun2(object,useSc=FALSE,signames=TRUE)
        nvp <- length(attr(dd,"thopt"))
        pars <- attr(dd,"optimum")
        pars <- pars[!is.na(names(pars))]
        hh <- hessian(dd, pars)/(-2)

        # Covariance on the scale of devfun2() (SDs and correlations, then
        # the fixed effects), then the delta method to variances and
        # covariances: var = sd^2, cov = cor*sd_a*sd_b. (A diagonal
        # rescaling of the Hessian with d(sd)/d(var) = 0.5*sd^(-1/2) was
        # used before; the derivative is 0.5/sd, and cov also depends on
        # the two SDs.)
        vc_sd <- -solve(hh)
        sdcormat <- as.data.frame(VarCorr(object, comp = "Std.Dev"),
          order = "lower.tri")
        Kr <- diag(2*sdcormat$sdcor, pran)
        is_sd <- is.na(sdcormat$var2)
        for (k in which(!is_sd)) {
          same <- is_sd & sdcormat$grp == sdcormat$grp[k]
          ia <- max(which(same & sdcormat$var1 == sdcormat$var1[k] & seq_len(pran) < k))
          ib <- min(which(same & sdcormat$var1 == sdcormat$var2[k] & seq_len(pran) > k))
          r <- sdcormat$sdcor[k]
          Kr[k, ] <- 0
          Kr[k, k] <- sdcormat$sdcor[ia]*sdcormat$sdcor[ib]
          Kr[k, ia] <- r*sdcormat$sdcor[ib]
          Kr[k, ib] <- r*sdcormat$sdcor[ia]
        }
        K <- as.matrix(Matrix::bdiag(Kr, diag(pfix)))
        full_vcov_noorder <- K %*% vc_sd %*% t(K)
        full_vcov <- matrix(NA, nrow(full_vcov_noorder),
          ncol(full_vcov_noorder))
        full_vcov[1:pfix, 1:pfix] <-
          full_vcov_noorder[(pran + 1):p, (pran + 1):p]
        full_vcov[(pfix + 1):p, (pfix + 1): p] <-
          full_vcov_noorder[1:pran, 1:pran]
        full_vcov[(pfix + 1):p, 1:pfix] <-
          full_vcov_noorder[1:pran, (pran + 1): p]
        full_vcov[1:pfix, (pfix + 1): p] <-
          full_vcov_noorder[(pran + 1): p, 1:pran]
    }
    if (!(ranpar %in% c("sd", "theta", "var"))){
       stop("ranpar needs to be sd, theta or var for glmerMod object.")
    }

    ## name the matrix
    parts <- getME(object, c("fixef", "theta"))
    colnames(full_vcov) <- c(names(parts$fixef), paste("cov",
      names(parts$theta), sep="_"))
  }
  return(full_vcov)
}

lav_matrix_trace <- function (..., check = TRUE)
{
  if (nargs() == 0L) {
    return(as.numeric(NA))
  }
  dots <- list(...)
  if (is.list(dots[[1]])) {
    mlist <- dots[[1]]
  }
  else {
    mlist <- dots
  }
  nMat <- length(mlist)
  if (nMat == 1L) {
    S <- mlist[[1]]
    if (check) {
      stopifnot(NROW(S) == NCOL(S))
    }
    out <- sum(S[lav_matrix_diag_idx(n = NROW(S))])
  }
  else if (nMat == 2L) {
    out <- sum(mlist[[1]] * t(mlist[[2]]))
  }
  else if (nMat == 3L) {
    A <- mlist[[1]]
    B <- mlist[[2]]
    C <- mlist[[3]]
    B2 <- B %*% C
    out <- sum(A * t(B2))
  }
  else {
    M1 <- mlist[[1]]
    M2 <- mlist[[2]]
    for (m in 3L:nMat) {
      M2 <- M2 %*% mlist[[m]]
    }
    out <- sum(M1 * t(M2))
  }
  out
}

lav_matrix_diag_idx <- function (n = 1L)
{
  1L + (seq_len(n) - 1L) * (n + 1L)
}

########## R-function: ZOSull ##########
# For creation of O'Sullivan-type Z matrices.
# Last changed: 04 OCT 2021 by M.P.Wand.

ZOSull <- function(x,intKnots, range.x, drv=0) {
  if (drv>2) {
    stop("splines not smooth enough for more than 2 derivatives")
  }

  # Set defaults for `range.x' and `intKnots'
  if (missing(range.x)) {
    range.x <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))
  }

  if (missing(intKnots))
  {
    numIntKnots <- min(length(unique(x)),35)
    intKnots <- quantile(unique(x),seq(0,1,length=
                                         (numIntKnots+2))[-c(1,(numIntKnots+2))])
  }
  numIntKnots <- length(intKnots)

  # Obtain the penalty matrix.
  allKnots <- c(rep(range.x[1],4),intKnots,rep(range.x[2],4))
  K <- length(intKnots) ; L <- 3*(K+8)
  xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+
               rep(allKnots,each=3)[-c(1,2,L)])/2
  wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
  Bdd <- splines::spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                             outer.ok=TRUE)$design
  Omega     <- t(Bdd*wts)%*%Bdd

  # Use the spectral decomposition of Omega to obtain Z.
  eigOmega <- eigen(Omega)
  indsZ <- 1:(numIntKnots+2)
  UZ <- eigOmega$vectors[,indsZ]
  LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))

  # Perform stability check.
  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX <- eigOmega$vectors[,indsX]
  Lmat <- cbind(UX,LZ)
  stabCheck <- t(crossprod(Lmat,t(crossprod(Lmat,Omega))))
  if (sum(stabCheck^2) > 1.0001*(numIntKnots+2)) {
    print("WARNING: NUMERICAL INSTABILITY ARISING\\
              FROM SPECTRAL DECOMPOSITION")
  }

  # Obtain B and post-multiply by LZ matrix to get Z.
  B <- splines::spline.des(allKnots,x,derivs=rep(drv,length(x)),
                           outer.ok=TRUE)$design

  Z <- crossprod(t(B),LZ)

  # Order the columns of Z with respect to the eigenvalues of "Omega":
  Z <- Z[,order(eigOmega$values[indsZ])]

  # Add the `range.x' and 'intKnots' as attributes
  # of the return object.
  attr(Z,"range.x") <- range.x
  attr(Z,"knots") <- intKnots

  # Return Z matrix with 2 attributes.
  return(Z)
}

############ R-function: Ztps ############
Ztps <- function(x, k, knots=NULL, range.x=NULL) {
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
    if ((d+1)%%2!=0) {
      ans[nzi] <- (abs(r[nzi]))^(2*m-d)*log(abs(r[nzi])) # d is even
    } else {
      ans[nzi] <- (abs(r[nzi]))^(2*m-d)
    }

    if (num.col>1) {
      ans <- matrix(ans,num.row,num.col)     # d is odd
    }
    return(ans)
  }

  # Set up function for matrix square-roots:
  matrix.sqrt <- function(A)
  {
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
    } else {
      stop("Matrix square root is not defined")
    }
    return(Asqrt)
  }

  if(is.null(knots)) {
    x1_grid <- seq(range(x[,1])[1], range(x[,1])[2], length = k)
    x2_grid <- seq(range(x[,2])[1], range(x[,2])[2], length = k)
    knots <- expand.grid(x1_grid, x2_grid)
    names(knots) <- colnames(x)
  }

  if (!is.null(range.x)) {
    inBdry <- pm_points_in_poly(knots, range.x)
    knots <- knots[inBdry, ]
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
  attr(Z,"knots") <- knots
  if (!is.null(range.x)) {
    attr(Z,"range.x") <- range.x
  }
  return(Z)
}
########################################################################
#==========================================================================
#  https://www.r-bloggers.com/2021/08/r-dataframe-merge-while-keeping-orders-of-row-and-column/
#—————————————————————–
# Function : f_loj_krc
#—————————————————————–
# Left outer join while keeping orders of input rows and columns
# Meaning of input arguments are the same as those of merge()
#—————————————————————–
f_loj_krc <- function(x, y, by.x, by.y) {

  # save row id
  x.temp <- x; x.temp$temp.id <- 1:nrow(x.temp);

  # each column names
  x.cn <- colnames(x); y.cn <- colnames(y)

  # replace column names of y with same names of x
  # to avoid duplicate fields
  for(i in 1:length(by.y)) {
    colnames(y)[which(y.cn == by.y[i])] <- by.x[i]
  }
  by.y <- by.x # since two fields are the same now

  # new column names of y
  y.cn <- colnames(y)

  # remove only joining key fields which are redundant
  # and keep only new informative fields
  y.cn.not.key <- setdiff(y.cn, by.y)

  # left outer join
  df <- merge(x = x.temp, y = y, by.x=by.x, by.y=by.y, all.x = TRUE)

  # recover the original rows and columns orders
  df <- df[order(df$temp.id),c(x.cn, y.cn.not.key)]; rownames(df) <- NULL

  return(df)
}

########################################################
# To perform a multiple comparison test based on the confidence intervals for each treatment's mean value. Specifically, if the confidence intervals for two treatments overlap, then they are not significantly different from each other, and if the confidence intervals do not overlap, then the treatments are significantly different from each other.

## LL -- Lower Limit of CI
## UL -- Upper Limit of CI
## trt_n -- names of treatment

ci_mcp <- function(LL, UL, trt_n=NULL) {

  stopifnot("Check your LL and UL input!"={
    is.numeric(LL)
	is.numeric(UL)
	length(LL)==length(UL)
	all(LL <= UL)   # a mean with SE 0 has LL = UL
  })
  trt_len <- length(LL)
  if(is.null(trt_n) || length(unique(trt_n))!=trt_len) {
    trt_n <- as.character(1:trt_len)
  }

  ci_mcp_letters_0 <- rep("A", trt_len)
  names(ci_mcp_letters_0) <- trt_n

  results <- matrix(NA_real_, nrow = trt_len, ncol = trt_len)

  for (i in 1:trt_len) {
    for (j in (i+1):trt_len) {
      if (j > trt_len) break
      ci1 <- c(LL[i],  UL[i])
      ci2 <-  c(LL[j],  UL[j])
      if (max(ci1) < min(ci2) || max(ci2) < min(ci1)) {
        results[i,j] <- 0.01
      } else {
        results[i,j] <- 0.08
      }
    }
  }

  if (all(unique(na.omit(as.vector(results)))==0.08)) {
    ci_mcp_letters <- ci_mcp_letters_0
  } else {
    rownames(results) <- colnames(results) <- trt_n
    results[lower.tri(results)] <- t(results)[lower.tri(results)]
    ci_mcp_letters <- multcompLetters(results, Letters=LETTERS)
  }
  return(ci_mcp_letters)
}

########################################################
# To handle aovlist object refit by lmer function
aovlist_lmer <- function(object) {
  stopifnot(inherits(object, "aovlist"))
  mod_df <- model.frame(object)
  lmer_call = match.call(aov, attr(object, "call"))

  trms = terms(object)
  response <- as.character(attr(trms, "variables"))[[2]]

  # Find the Error terms
  trms_label = attr(trms, "term.labels")
  err.idx = grep("^Error\\(", trms_label)

  # Original Error term
  error_term <- trms_label[err.idx]
  error_term <- gsub("Error\\(|\\)", "", error_term)
  error_terms <- strsplit(error_term, "\\+\\s+")[[1]]

  # Reformat into lmer format
  lmer_random_parts <- paste0("(1|", error_terms, ")", collapse = " + ")
  lmer_fix_parts <- paste0(trms_label[-err.idx], collapse = " + ")

  lmer_call$formula <- as.formula(paste(response, " ~ ", lmer_fix_parts, "+", lmer_random_parts, sep=""))
  lmer_call[[1]] <- as.name("lmer")
  assign(as.character(lmer_call$data), mod_df)
  lmer_object <- eval(lmer_call)
  return(lmer_object)
}



## Resolve the df method for an lmer model: "auto" is Kenward-Roger for REML
## fits and Satterthwaite for ML fits. Returns list(df_type, vcov_adj).
pm_df_method <- function(model, df_method) {
  if (!inherits(model, "lmerMod")) {
    if (df_method != "auto") warning("'df_method' applies to lmer models only; ignored.")
    return(list(df_type = NULL, vcov_adj = NULL))
  }
  reml <- lme4::getME(model, "is_REML")
  df_type <- switch(df_method, auto = if (reml) "Kenward-Roger" else "Satterthwaite",
                    df_method)
  if (df_type == "Kenward-Roger" && !reml) {
    warning("Kenward-Roger df need a REML fit: using Satterthwaite instead")
    df_type <- "Satterthwaite"
  }
  if (df_method == "auto" && !reml) message("The model was fitted by ML: using Satterthwaite df.")
  vcov_adj <- switch(df_type,
    "Kenward-Roger" = tryCatch(pm_vcov_kr(model), error = function(e) NULL),
    "Satterthwaite" = pm_satt_parts(model),
    NULL)
  list(df_type = df_type, vcov_adj = vcov_adj)
}

## findbars() moved from lme4 to reformulas; use reformulas when installed.
pm_findbars <- function(f) {
  if (requireNamespace("reformulas", quietly = TRUE)) reformulas::findbars(f) else lme4::findbars(f)
}
