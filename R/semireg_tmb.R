#' Fitting Semi Parametric Models Using glmmTMB
#'
#' Fit a semi parametric model based on glmmTMB.
#'
#' A semi parametric model can be parameterized as a linear (or generalized
#' linear) mixed model in which its random effects are smooth functions of some
#' covariates (named ‘smooth term’). \code{semireg_tmb} follows the approach
#' suggested by Wand and Ormerod (2008) and represents the 'smooth term' using
#' O'Sullivan-type of Z.
#'
#' @param formula A two-sided linear formula object describing both the
#' fixed-effects and random-effects part of the model, with the response on the
#' left of a ~ operator and the terms, separated by + operators, on the right.
#' Random-effects terms are distinguished by vertical bars ("|") separating
#' expressions for design matrices from grouping factors.
#' @param data data frame (tibbles are OK) containing model variables. Not
#' required, but strongly recommended; if \code{data} is not specified,
#' downstream methods such as prediction with new data
#' (\code{predict(fitted_model, newdata = ...)}) will fail. If it is necessary
#' to call \code{glmmTMB} with model variables taken from the environment
#' rather than from a data frame, specifying \code{data=NULL} will suppress the
#' warning message.
#' @param family a family function, a character string naming a family
#' function, or the result of a call to a family function (variance/link
#' function) information. See \code{family} for a generic discussion of
#' families or \code{family_glmmTMB} for details of \code{glmmTMB}-specific
#' families.
#' @param smoothZ A list includes a set of smooth Z matrixs (called 'smooth
#' term') used in the mixed effects model, the name of 'smooth term' should be
#' different any variables in the model, each 'smooth term' is the result of
#' function \code{smZ}. e.g. smoothZ=list(sm1=smZ(x1), sm2=smZ(x2, by=f1),
#' sm3=smZ(x3, by=f2, group=TRUE), ...) where 'sm1' to 'sm3' should be new
#' variable names in the \code{data}, and x1 to x3 are covariates, and f1, f2
#' are factors.
#' @param ziformula a \emph{one-sided} (i.e., no response variable) formula for
#' zero-inflation combining fixed and random effects: the default \code{~0}
#' specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation
#' formula identical to the right-hand side of \code{formula} (i.e., the
#' conditional effects formula); terms can also be added or subtracted.
#' \strong{When using \code{~.} as the zero-inflation formula in models where
#' the conditional effects formula contains an offset term, the offset term
#' will automatically be dropped}. The zero-inflation model uses a logit link.
#' @param dispformula a \emph{one-sided} formula for dispersion containing only
#' fixed effects: the default \code{~1} specifies the standard dispersion given
#' any family. The argument is ignored for families that do not have a
#' dispersion parameter. For an explanation of the dispersion parameter for
#' each family, see \code{sigma}. The dispersion model uses a log link. In
#' Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to
#' be 0 (actually a small non-zero value), forcing variance into the random
#' effects. The precise value can be controlled via
#' \code{control=glmmTMBControl(zero_dispval=...)}; the default value is
#' \code{sqrt(.Machine$double.eps)}.
#' @param weights weights, as in \code{glm}. Not automatically scaled to have
#' sum 1.
#' @param offset offset for conditional model (only).
#' @param contrasts an optional list, e.g., \code{list(fac1="contr.sum")}. See
#' the \code{contrasts.arg} of \code{model.matrix.default}.
#' @param na.action a function that specifies how to handle observations
#' containing \code{NA}s.  The default action (\code{na.omit}, inherited from
#' the 'factory fresh' value of \code{getOption("na.action")}) strips any
#' observations with any missing values in any variables. Using \code{na.action
#' = na.exclude} will similarly drop observations with missing values while
#' fitting the model, but will fill in \code{NA} values for the predicted and
#' residual values for cases that were excluded during the fitting process
#' because of missingness.
#' @param se whether to return standard errors.
#' @param verbose whether progress indication should be printed to the console.
#' @param doFit whether to fit the full model, or (if FALSE) return the
#' preprocessed data and parameter objects, without fitting the model.
#' @param control control parameters, see \code{glmmTMBControl}.
#' @param REML whether to use REML estimation rather than maximum likelihood.
#' @param start starting values, expressed as a list with possible components
#' \code{beta}, \code{betazi}, \code{betad} (fixed-effect parameters for
#' conditional, zero-inflation, dispersion models); \code{b}, \code{bzi}
#' (conditional modes for conditional and zero-inflation models); \code{theta},
#' \code{thetazi} (random-effect parameters, on the standard deviation/Cholesky
#' scale, for conditional and z-i models); \code{psi} (extra family parameters,
#' e.g., shape for Tweedie models).
#' @param map a list specifying which parameter values should be fixed to a
#' constant value rather than estimated. \code{map} should be a named list
#' containing factors corresponding to a subset of the internal parameter names
#' (see \code{start} parameter). Distinct factor values are fitted as separate
#' parameter values, \code{NA} values are held fixed: e.g.,
#' \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three
#' fixed-effect parameters of the conditional model and fix the fourth
#' parameter to its starting value. In general, users will probably want to use
#' \code{start} to specify non-default starting values for fixed parameters.
#' See \code{MakeADFun} for more details.
#' @param sparseX a named logical vector containing (possibly) elements named
#' "cond", "zi", "disp" to indicate whether fixed-effect model matrices for
#' particular model components should be generated as sparse matrices, e.g.
#' \code{c(cond=TRUE)}. Default is all \code{FALSE}
#' @param prt Logical scalar - Should the info to be print on screen in the
#' middle of the process or not?
#' @param predict_info Logical scalar - Should provide the info for function
#' semipred or not? In case of there is a correlation theta parameter
#' appearing, you may set predict=FALSE.
#' @return \item{semer}{A glmmTMB model used in the fitting.} \item{data}{A
#' data.frame with generated variables in the fitting.} \item{fomul_vars}{Name
#' of variables in the formula of semireg_tmb model.} \item{sm_vars}{Name of
#' variables in the smoothZ list.} \item{smoothZ_call}{A call used to produce
#' smooth terms in the fitting.} \item{knots_lst}{Knots used in each smooth
#' term in the fitting.} \item{range_lst}{Range of covariate used in each
#' smooth term in the fitting.} \item{cov_lst}{Covariance matrix list for each
#' smooth term.} \item{u_lst}{Random effects list for each smooth term.}
#' \item{type_lst}{Smooth type list of smooth terms.} \item{CovMat}{Covariance
#' matrix for all smooth terms.} \item{Cov_ind}{Covariance matrix index for
#' each smooth term.} \item{Cov_indN}{Covariance matrix index for each smooth
#' term when \code{group=TRUE} in \code{smoothZ} argument.} \item{df}{Degree of
#' freedom of all random terms.} \item{tmbf}{The glmmTMB model result using
#' doFit=FALSE.}
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
#'  fit1 <- semireg_tmb(areaPerMzloty ~ construction.date,
#'                      smoothZ=list(
#'                        grp=smZ(construction.date, k=25)
#'                      ),
#'                      data = WarsawApts)
#'  sp_out1 <- semipred(fit1, "construction.date", "construction.date")
#'
#'  WarsawApts$district <- factor(WarsawApts$district)
#'  fit2 <- semireg_tmb(areaPerMzloty ~ construction.date*district, resp_scale = TRUE,
#'                      smoothZ=list(group=smZ(construction.date, k=15,
#'                                             by = district, group=TRUE)),
#'                      data=WarsawApts)
#'  sp_out2_1 <- semipred(fit2, "district", "construction.date")
#'  sp_out2_2 <- semipred(fit2, "district", "construction.date", contr=c(2,1))
#'
#'  data(indonRespir)
#'  help(indonRespir)
#'  str(indonRespir)
#'  fit3 <- semireg_tmb(respirInfec ~ age+vitAdefic + female + height
#'                      + stunted + visit2 + visit3 + visit4  + visit5 + visit6+(1|idnum),
#'                      smoothZ=list(
#'                        grp=smZ(age)
#'                      ),
#'                      family = binomial,
#'                      data = indonRespir)
#'  sp_out3 <- semipred(fit3, "age", "age")
#'  library(ggplot2)
#'  sp_out3$plt+
#'    geom_rug(data = subset(indonRespir, respirInfec==0), sides = "b", col="deeppink") +
#'    geom_rug(data = subset(indonRespir, respirInfec==1), sides = "t", col="deeppink")+
#'    ylim(0, 0.2)
#' }
#' @importFrom glmmTMB fitTMB glmmTMBControl
#' @importFrom stats gaussian na.action setNames
#' @export
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
  if (any(unlist(lapply(data[, setdiff(sm_vars, fomul_vars)], is.factor)))) {
    stop("Any factor in 'smoothZ' list must be in the formula!")
  }
  response_n <- fomul_vars[1]
  fomul_vars <- fomul_vars[-1]

  data <- droplevels(as.data.frame(na.omit(data[, intersect(all.vars(mc), names(data))])))

  smoothZm <-  eval(substitute(smoothZ), data)
  Z_namem <- names(smoothZm)
  smoothZN <- list()
  for (i in Z_namem){
    if (is.list(smoothZm[[i]])) {
      smoothZN <- c(smoothZN, unlist(smoothZm[i], recursive=FALSE))
    } else {
      smoothZN <- c(smoothZN, smoothZm[i])
    }
  }

  smoothZt <- lapply(smoothZN, t)

  glmmTMBcall <- mc

  glmmTMBcall[[1]] <- as.name("glmmTMB")
  glmmTMBcall$smoothZ <- NULL
  glmmTMBcall$prt <- NULL
  glmmTMBcall$predict_info <- NULL

  if (!length(smoothZt)) {
    return(eval.parent(glmmTMBcall))
  }
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
    if (length(tn) > 1) {
      stop("a smoothZt factor must be associated with only one r.e. term")
    }
    ind <- (cond_reTrms$Gp)[tn:(tn+1L)]
    rowsi <- (ind[1]+1L):ind[2]
    stopifnot(all(dim(Zt[rowsi,])==dim(smoothZt[[i]])))
    Zt[rowsi,] <- smoothZt[[i]]
  }

  tmbf$data.tmb$Z <- t(Zt)
  semer <- fitTMB(tmbf)
  if (prt) {
    print(semer)
  }

  if (predict_info) {
    beta <- fixef(semer)$cond
    n_beta <- length(beta)
    s <- sigma(semer)
    Lambdat <- reTrms_tmb(semer)$Lambdat
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
    if (family_n=="binomial") {
      yobs <- unique(tmbf$data.tmb$yobs)
    }
    link_n <- semer$modelInfo$family$link
    uHat <- unlist(ranef(semer)$cond) # cbind(uHat, getME(semer, "b"))
    C <- cbind(getME(semer, "X"), getME(semer, "Z"))

    if (family_n != "gaussian") {
      if (family_n =="binomial" && all(yobs%in%c(0,1))) {
        wVec <- switch(link_n,
                       "logit" = 3/(pi^2),
                       "probit" = 1,
                       "cloglog" = 6/(pi^2))
      } else {
        if (inherits(semer, "glmmTMB")) {
          weight_fun <- function(eta, mod){
            family_fun <- mod$modelInfo$family
            family_fun$mu.eta(eta)^2/family_fun$variance(family_fun$linkinv(eta))
          }
          y_hat <- as.vector(C%*%c(beta, uHat))
          wVec <- as.vector(weight_fun(y_hat, semer))*tmbf$data.tmb$weights
          if (length(tmbf$data.tmb$size) > 0) {
            wVec <- wVec*tmbf$data.tmb$size
          }
        } else {
          wVec <- semer@resp$sqrtWrkWt()^2
        }
      }
    } else {
      wVec <- weights
    }

    if (isLMM(semer)) {
      CTC <- Matrix::tcrossprod(t(C))
    } else {
      CTC <- crossprod(C*wVec,C)
    }

    # fullCovMat <- Matrix::solve(CTC + D)

    AA <- CTC + D
    if (any(eigen(AA)$values <= 0)) {
      AAmatrixpd <- nearPD(AA, doSym =TRUE)$mat
      AA <- AAmatrixpd@x
      dim(AA) <- AAmatrixpd@Dim
    }
    fullCovMat <- Matrix::solve(AA)

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
  } else {
    knots_lst <- range_lst <- type_lst <- cov_lst <- u_lst <- fullCovMat <- vcov_ind <- vcov_indN <- df <- NULL
  }

  ans <- list(semer=semer, data=data, fomul_vars=fomul_vars, sm_vars=sm_vars, smoothZ_call=mc$smoothZ, knots_lst=knots_lst, range_lst=range_lst, type_lst=type_lst, cov_lst=cov_lst, u_lst=u_lst, CovMat=fullCovMat, Cov_ind=vcov_ind, Cov_indN=vcov_indN, df=df, tmbf=tmbf)
  class(ans) <- c("semireg", "list")
  return(ans)
}

