#' Permutation Test of random or fixed effects for \code{lmer} model.
#'
#' This function provides permutation tests for the terms in a linear mixed
#' model of \code{lmer}.
#'
#'
#' @param lmer0 \code{lmer} model under H0, note that \code{lmer0} model must
#' nest within \code{lmer1} model.
#' @param lmer1 \code{lmer} model under H1, note that \code{lmer0} model must
#' nest within \code{lmer1} model.
#' @param nperm Number of permutation, the default value is 999.
#' @param ncore Number of core for parallel computing, the default value is 3.
#' @param plot Plot permutation distribution or not, the default value is
#' FALSE.
#' @param seed Specify a random number generator seed, for reproducible
#' results.
#' @return Permutation p-value.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @references Oliver E. Lee and Thomas M. Braun (2012), \emph{Permutation
#' Tests for Random Effects in Linear Mixed Models. Biometrics}, Journal 68(2).
#' @examples
#' \dontrun{
#' library(predictmeans)
#' # Test random effects
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
#' fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
#' anova(fm1, fm2, fm3)
#' permlmer(fm3, fm2)
#' permlmer(fm2, fm1)
#'
#' # Test fixed effects
#' Oats$nitro <- factor(Oats$nitro)
#' fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' permlmer(fm0, fm)
#' }
#'
#' @importFrom parallel clusterEvalQ
#' @importFrom stats density logLik
#' @export
permlmer <- function(lmer0, lmer1, nperm = 999, ncore=3L, plot=FALSE, seed){

  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) {
    stop("The model must be a lmer object!")
  }
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) {
    stop("Please check the response in your model!")
  }

  c_deparse <- function (...)
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")

  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1)
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1)

  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }

  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))

  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) {
    stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  }

  ref <- ifelse(setequal(fixef0name, fixef1name), TRUE, FALSE)

  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  Lambda1n <- getME(lmer1, "Lambda")
  Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X") %*% fixef0)
  errors <- getME(lmer0, "y") - xbeta

  #Weighting the residuals.
  wterrors <- wt %*% errors

  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors))
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) {
      set.seed(seed)
    }
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut %*% sample(wterrors))))
  } else {
    if (!missing(seed)) {
      set.seed(seed)
    }
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      as.data.frame(xbeta[rowindex]+as.vector(Ut %*% wterrors[rowindex]))
    })
  }

  # Calculating the likelihood ratio test statistic for each permutation.
  lrtest1 <- 2*(logLik(lmer1, REML=ref)-logLik(lmer0, REML=ref))
  lrtest1 <- ifelse(lrtest1 < 0, 0, lrtest1)

  if (.Platform$OS.type=="windows") {

    ## DONGWEN: I have added this function because the code was failing under
    ## devtools::check(). There's a long explanation, but I am trying something
    ## suggested by ChatGPT

    get_safe_cores <- function(default = 2) {
      limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_")
      if (limit == "TRUE") {
        return(min(default, 2))
      }
      return(min(default, parallel::detectCores()))
    }

    ## this tries to respect the request from the function call to permlmer
    ncore = get_safe_cores(ncore)


    cl <- makeCluster(ncore)
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment())
    lrtest2 <- parLapplyLB(cl, permy, function(x) {
      LRT <- try(2*(logLik(refit(lmer1, x), REML=ref) - logLik(refit(lmer0, x), REML=ref)), TRUE)
      LRT <- ifelse(is.numeric(LRT), LRT, NA)
    })
    stopCluster(cl)
  } else {
    lrtest2 <- mclapply(permy, function(x) {
      LRT <- try(2*(logLik(suppressMessages(refit(lmer1, x)), REML=ref) - logLik(suppressMessages(refit(lmer0, x)), REML=ref)), TRUE)
      LRT <- ifelse(is.numeric(LRT), LRT, NA)
    }, mc.cores=ncore)
  }

  #Calculating the p-values.
  lrtest <- na.omit(unlist(lrtest2))
  lrtest <- ifelse(lrtest < 0, 0, lrtest)
  perm_p <- (sum(lrtest >= lrtest1) +1)/(length(lrtest) + 1)
  aod <- anova(lmer0, lmer1, refit=!ref)
  aod$'Perm-p' <- c(NA, perm_p)

  if (plot) {
    # dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}
