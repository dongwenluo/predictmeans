#' Fitting Semi Parametric Models Using lme4 Ecosystem
#'
#' Fit a semi parametric model based on lme4 ecosystem including \code{lmer},
#' \code{glmer} and \code{glmer.nb}.
#'
#' A semi parametric model can be parameterized as a linear (or generalized
#' linear) mixed model in which its random effects are smooth functions of some
#' covariates (named ‘smooth term’). \code{semireg} follows the approach
#' suggested by Wand and Ormerod (2008) and represents the 'smooth term' using
#' O'Sullivan-type of Z.
#'
#' @param formula A two-sided linear formula object describing both the
#' fixed-effects and random-effects part of the model, with the response on the
#' left of a ~ operator and the terms, separated by + operators, on the right.
#' Random-effects terms are distinguished by vertical bars ("|") separating
#' expressions for design matrices from grouping factors.
#' @param data A data frame or list containing the model response variable and
#' covariates required by the formula. By default the variables are taken from
#' environment(formula and smoothZ), typically the environment from which
#' semireg is called.
#' @param family A GLM family, see glm and family.
#' @param ngbinomial Logical scalar - Should a negative binomial GLMMs be used?
#' .
#' @param REML Logical scalar - Should the estimates be chosen to optimize the
#' REML criterion (as opposed to the log-likelihood)?
#' @param smoothZ A list includes a set of smooth Z matrixs (called 'smooth
#' term') used in the mixed effects model, the name of 'smooth term' should be
#' different any variables in the model, each 'smooth term' is the result of
#' function \code{smZ}. e.g. smoothZ=list(sm1=smZ(x1), sm2=smZ(x2, by=f1),
#' sm3=smZ(x3, by=f2, group=TRUE), ...) where 'sm1' to 'sm3' should be new
#' variable names in the \code{data}, and x1 to x3 are covariates, and f1, f2
#' are factors.
#' @param ncenter Logical scalar - Should the numeric predictors to be centered
#' or not?
#' @param nscale Logical scalar - Should the numeric predictors to be scaled or
#' not?
#' @param resp_scale Logical scalar - Should the response be involved in the
#' scaling action or not?
#' @param control A list (of correct class, resulting from lmerControl() or
#' glmerControl() respectively) containing control parameters, including the
#' nonlinear optimizer to be used and parameters to be passed through to the
#' nonlinear optimizer, see the *lmerControl documentation for details.
#' @param start Starting value list as used by lmer or glmer.
#' @param verbose Passed on to fitting lme4 fitting routines.
#' @param drop.unused.levels By default unused levels are dropped from factors
#' before fitting. For some smooths involving factor variables you might want
#' to turn this off. Only do so if you know what you are doing.
#' @param subset An optional expression indicating the subset of the rows of
#' data that should be used in the fit. This can be a logical vector, or a
#' numeric vector indicating which observation numbers are to be included, or a
#' character vector of the row names to be included. All observations are
#' included by default.
#' @param weights An optional vector of ‘prior weights’ to be used in the
#' fitting process. Should be NULL or a numeric vector.
#' @param offset This can be used to specify an a priori known component to be
#' included in the linear predictor during fitting. This should be NULL or a
#' numeric vector of length equal to the number of cases. One or more offset
#' terms can be included in the formula instead or as well, and if more than
#' one is specified their sum is used. See model.offset.
#' @param contrasts An optional list. See the contrasts.arg of
#' model.matrix.default.
#' @param prt Logical scalar - Should the info to be print on screen in the
#' middle of the process or not?
#' @param predict_info Logical scalar - Should provide the info for function
#' semipred or not?
#' @param ...  Further arguments for passing on to model setup routines.
#' @return \item{semer}{A mer model used in the fitting.} \item{data}{A
#' data.frame with generated variables in the fitting.} \item{fomul_vars}{Name
#' of variables in the formula of semireg model.} \item{sm_vars}{Name of
#' variables in the smoothZ list.} \item{smoothZ_call}{A call used to produce
#' smooth terms in the fitting.} \item{knots_lst}{Knots used in each smooth
#' term in the fitting.} \item{range_lst}{Range of covariate used in each
#' smooth term in the fitting.} \item{cov_lst}{Covariance matrix list for each
#' smooth term.} \item{u_lst}{Random effects list for each smooth term.}
#' \item{type_lst}{Smooth type list of smooth terms.} \item{CovMat}{Covariance
#' matrix for all smooth terms.} \item{Cov_ind}{Covariance matrix index for
#' each smooth term.} \item{Cov_indN}{Covariance matrix index for each smooth
#' term when \code{group=TRUE} in \code{smoothZ} argument.} \item{df}{Degree of
#' freedom of all random terms.} \item{lmerc}{Call used in the mer model in the
#' fitting.}
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @references Wand, M.P. and Ormerod, J.T. (2008). On semiparametric
#' regression with O'Sullivan penalized splines. \emph{Australian and New
#' Zealand Journal of Statistics.} \bold{50}, 179-198.
#' @examples
#'
#' \dontrun{
#'  library(predictmeans)
#'  library(HRW)
#'  data(WarsawApts)
#'  help(WarsawApts)
#'  str(WarsawApts)
#'  fit1 <- semireg(areaPerMzloty ~ construction.date,
#'                  smoothZ=list(
#'                    grp=smZ(construction.date, k=25)
#'                  ),
#'                  data = WarsawApts)
#'  sp_out1 <- semipred(fit1, "construction.date", "construction.date")
#'
#'  WarsawApts$district <- factor(WarsawApts$district)
#'  fit2 <- semireg(areaPerMzloty ~ construction.date*district, resp_scale = TRUE,
#'                  smoothZ=list(group=smZ(construction.date, k=15,
#'                                         by = district, group=TRUE)),
#'                  data=WarsawApts)
#'  sp_out2_1 <- semipred(fit2, "district", "construction.date")
#'  sp_out2_2 <- semipred(fit2, "district", "construction.date", contr=c(2,1))
#'
#'  data(indonRespir)
#'  help(indonRespir)
#'  str(indonRespir)
#'  fit3 <- semireg(respirInfec ~ age+vitAdefic + female + height
#'                  + stunted + visit2 + visit3 + visit4  + visit5 + visit6+(1|idnum),
#'                  smoothZ=list(
#'                    grp=smZ(age)
#'                  ),
#'                  family = binomial,
#'                  data = indonRespir)
#'  sp_out3 <- semipred(fit3, "age", "age")
#'  library(ggplot2)
#'  sp_out3$plt+
#'    geom_rug(data = subset(indonRespir, respirInfec==0), sides = "b", col="deeppink") +
#'    geom_rug(data = subset(indonRespir, respirInfec==1), sides = "t", col="deeppink")+
#'    ylim(0, 0.2)
#' }
#'
#' @importFrom lme4 mkGlmerDevfun mkLmerDevfun mkMerMod optimizeGlmer
#' @importFrom lme4 optimizeLmer updateGlmerDevfun
#' @importFrom Matrix .bdiag Matrix nearPD
#' @importFrom stats fitted model.response
#' @export
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
  if (any(unlist(lapply(data[, setdiff(sm_vars, fomul_vars)], is.factor)))) {
    stop("Any factor in 'smoothZ' list must be in the formula!")
  }
  response_n <- fomul_vars[1]
  fomul_vars <- fomul_vars[-1]

  data <- droplevels(as.data.frame(na.omit(data[, intersect(all.vars(mc), names(data))])))

  if (resp_scale) {
    numeric_n <- names(data)[sapply(data, is.numeric)]
  } else {
    numeric_n <- setdiff(names(data)[sapply(data, is.numeric)], response_n)
  }
  if (any(ncenter, nscale)) {
    for (i in numeric_n) {
      data_i <- data[,i]
      if (!is.matrix(data_i)) {
        scaled_i <- scale(data[i], center=ncenter, scale=nscale)
        data[i] <- as.numeric(scaled_i)
        if (ncenter) {
          attr(data, paste(i, "mean", sep="_")) <- attr(scaled_i,"scaled:center")
        }
        if (nscale) {
          attr(data, paste(i, "sd", sep="_")) <- attr(scaled_i,"scaled:scale")
        }
      }
    }
    attr(data, "numeric_var") <- numeric_n
  }
  smoothZm <-  eval(substitute(smoothZ), data)
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
    if (is.list(smoothZm[[i]])) {
      smoothZN <- c(smoothZN, unlist(smoothZm[i], recursive=FALSE))
    } else {
      smoothZN <- c(smoothZN, smoothZm[i])
    }
  }
  smoothZt <- lapply(smoothZN, t)

  gaus <- FALSE
  if (is.null(family)) {
    gaus <- TRUE
  } else {
    ## copied from glm()
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) {
      family <- family()
    }
    if (!inherits(family, "family")) {
      stop("unknown family type")
    }
    gaus <- family$family == "gaussian" && family$link == "identity"
  }

  lmerc <- mc                         # create a call to lmer

  lmerc[[1]] <- ifelse(gaus, as.name("lmer"), as.name("glmer"))
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

  if (!gaus) {
    lmerc$REML <- NULL
  }

  if (!length(smoothZt)) {              # call [g]lmer instead
    return(eval.parent(lmerc))
  }
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
    if (length(tn) > 1) {
      stop("a smoothZt factor must be associated with only one r.e. term")
    }
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

  } else {
    if (gaus) {
      dfl$REML = resp$REML > 0L
      devfun <- do.call(mkLmerDevfun,dfl)
      opt <- optimizeLmer(devfun, start=start, verbose=verbose, optimizer=control$optimizer, control=control$optCtrl, ...)
      semer0 <- mkMerMod(environment(devfun), opt, reTrms, lmf@frame, mc)

      orig_call <- lmerc
      args <- as.list(lmerc)
      args$devFunOnly <- TRUE
      if (!"control" %in% names(as.list(lmerc))) {
        args$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
      }
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
  if (prt) {
    print(semer)
  }

  if (predict_info) {
    n_beta <- getME(semer, "p")
    C <- cbind(getME(semer, "X"), getME(semer, "Z"))

    if (gaus) {
      CTC <- Matrix::tcrossprod(t(C))
    } else {
      family_n <- semer@resp$family$family
      link_n <- semer@resp$family$link
      yobs <- unique(semer@resp$y)
      if (family_n =="binomial" && all(yobs%in%c(0,1))) {
        wVec <- switch(link_n,
                       "logit" = 3/(pi^2),
                       "probit" = 1,
                       "cloglog" = 6/(pi^2))
      } else {
        wVec <- semer@resp$sqrtWrkWt()^2
      }

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


