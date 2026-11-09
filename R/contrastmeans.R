#' Linear Contrast Tests for a Linear Model
#'
#' @description
#' Performs t-tests (or permuted t-tests) of specified contrasts for linear models 
#'   obtained from functions \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, 
#'   or \code{lmer}.
#' @usage
#' contrastmeans(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", Df, permlist,
#'                 vcov_adjust=TRUE, weights=c("equal", "proportional", "outer", "cells", "flat"),
#'                 df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment"))
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's contrast to be calculated. 
#'   The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B".
#' @param ctrmatrix A specified contrast matrix. If \code{ctrmatrix} is missing, the programe will ask user to enter it. The columns can be named by level (e.g. from \code{\link{specify_contrast}}); levels that are not given get zero weight.
#' @param ctrnames Names of the specified contrasts, e.g. c("A vs D", "C vs B", ...)
#' @param adj Name (in "quote") for indicating a method for adjusting p-values of pairwise comparisons. 
#'   The choices are "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "sidak", "scheffe" and
#'   "mvt", in any letter case. "sidak", "scheffe" (rank of \code{ctrmatrix}) and "mvt" (single-step multivariate t, the exact counterpart of SAS
#'   \code{ADJUST=SIMULATE}; needs package \pkg{mvtnorm}) are as in \pkg{emmeans}.
#'   The default method is "none".
#' @param Df A denominator degree of freedom for \code{modelterm}. (For \code{glmer} models 
#'     the \code{Df} needs to be specified, while for the other models, \code{Df} is obtained from the fitted model automatically).
#' @param permlist A model parameter list containing \code{nsim} parameters produced by the function \code{permmodels}. When \code{permlist != NULL}, the option \code{Df} will be non-functional. This is a key 
#' 	option for the permutation test.
#' @param vcov_adjust For \code{lmer} models fitted by REML: if \code{TRUE} (the default), standard errors use the Kenward-Roger
#'   adjusted covariance matrix of the fixed effects, as SAS \code{DDFM=KR} and \pkg{emmeans} do, to go with the
#'   Kenward-Roger degrees of freedom. \code{FALSE} uses the unadjusted covariance and reproduces results of
#'   predictmeans 1.x. It has no effect for other models.
#' @param weights Passed to \code{\link{Kmatrix}}: how the cells of factors not in \code{modelterm} are weighted when averaging
#'   over them (\code{"equal"}, the default, \code{"proportional"}, \code{"outer"}, \code{"cells"} or \code{"flat"}).
#' @param df_method Denominator df for \code{lmer} models: \code{"auto"} (Kenward-Roger for REML fits, Satterthwaite for ML
#'   fits), \code{"Kenward-Roger"}, \code{"Satterthwaite"} or \code{"containment"} (SAS \code{DDFM=CONTAIN}; see
#'   \code{\link{predictmeans}}).
#' @return
#' There are two components in the output which are
#' \item{Table}{A table showing t-test results for the specified linear contrasts.}
#' \item{K}{A contrast matrix.}
#' @references
#' Torsten Hothorn, Frank Bretz and Peter Westfall (2008),
#'     \emph{Simultaneous Inference in General Parametric Models. Biometrical}, Journal 50(3), 346--363.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' # ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # library(lme4)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'
#' ## The contrast has a contrast matrix as follows:
#' #     0:Golden Rain 0:Marvellous 0:Victory 
#' #[1,]            -1            0         1 
#' #[2,]             0            0         1 
#' #     0.2:Golden Rain 0.2:Marvellous 0.2:Victory 
#' #[1,]               0              0           0 
#' #[2,]               0              0           0 
#' #     0.4:Golden Rain  0.4:Marvellous 0.4:Victory
#' #[1,]               0               0           0
#' #[2,]               0              -1           0
#' #      0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#' #[1,]                0              0           0
#' #[2,]                0              0           0
#'
#' # 1. Enter above contrast matrix into a pop up window, then close the window
#' if (interactive()) contrastmeans(fm, "nitro:Variety")
#'  
#' # 2. Construct the contrast matrix directly
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)
#' @name contrastmeans
NULL

contrastmeans <- function(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", Df, permlist, vcov_adjust=TRUE,
                          weights=c("equal", "proportional", "outer", "cells", "flat"),
                          df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment")) {
  vcov_adjust_default <- missing(vcov_adjust)
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)

  if (inherits(model, "aovlist")) {
    model <- aovlist_lmer(model)
  }

  if (is.null(ctrnames) || all(ctrnames%in%c("NULL", ""))) {
    ctrnames <- NULL
  }
  weights <- match.arg(weights)
  adj <- pm_match_adj(adj, tukey = FALSE)
  KK <- Kmatrix(model, modelterm, weights=weights)$K
  termsLabel <- rownames(KK)
  if (missing(ctrmatrix)) {
    if (!interactive()) {
      stop("Please supply 'ctrmatrix'; the contrast editor needs an interactive session.")
    }
    ques <- paste("\nHow many contrasts you want set up for ", sQuote(modelterm), "? ", sep="")
    nctr <- as.integer(readline(ques))
    nrK <- length(termsLabel)
    xnew <- matrix(rep(0, nctr*nrK), ncol=nrK)
    colnames(xnew) <- termsLabel
    xnew <- edit(xnew)
    ctrmatrix <- xnew[1:nctr, 1:nrK, drop=FALSE]
    rownames(ctrmatrix) <- ctrnames
    if (!(all(rowSums(ctrmatrix)==0)))  {
      cat("\n", "The contrast matrix is:\n\n")
      print(ctrmatrix)
      cat("\n\n")
      stop("\n", "Please check the row",  sQuote(which(rowSums(ctrmatrix)!=0)), "of the contrast matrix!\n\n")
    }
  } else {
    if (is.null(dim(ctrmatrix))) {
      ctrmatrix <- matrix(ctrmatrix, nrow=1, dimnames=list(NULL, names(ctrmatrix)))
    }
    if (!is.null(colnames(ctrmatrix)) && all(colnames(ctrmatrix) %in% termsLabel)) {
      # columns named by level (e.g. from specify_contrast()): put them in the
      # order of the means, with zero weight for levels not mentioned
      full <- matrix(0, nrow(ctrmatrix), length(termsLabel), dimnames=list(rownames(ctrmatrix), termsLabel))
      full[, colnames(ctrmatrix)] <- ctrmatrix
      ctrmatrix <- full
    } else {
      if (ncol(ctrmatrix) != length(termsLabel)) {
        stop("'ctrmatrix' must have one column for each of the ", length(termsLabel), " levels of ",
             sQuote(modelterm), ": ", paste(termsLabel, collapse=", "))
      }
      colnames(ctrmatrix) <- termsLabel
    }
    if (!is.null(ctrnames)) rownames(ctrmatrix) <- ctrnames
    nctr <- nrow(ctrmatrix)
  }
  if (all((is.null(dim(ctrmatrix)) | dim(ctrmatrix)[1]==1), missing(permlist)))  {
    adj <- "none"
  }

  # Means for missing cells (NaN rows of K) or other non-estimable means can
  # only be left out of a contrast, not used in it.
  est_flag <- attr(KK, "estimable")
  if (is.null(est_flag)) est_flag <- rep(TRUE, nrow(KK))
  bad_row <- !est_flag | apply(KK, 1, function(r) any(is.na(r)))
  used_bad <- colSums(abs(ctrmatrix[, bad_row, drop=FALSE])) > 0
  if (any(used_bad)) {
    stop("The contrast(s) use mean(s) that cannot be estimated: ",
         paste(termsLabel[bad_row][used_bad], collapse=", "))
  }
  KK[bad_row, ] <- 0
  df_method <- match.arg(df_method)
  dfm <- pm_df_method(model, df_method)
  df_type <- dfm$df_type
  vcov_adj <- dfm$vcov_adj
  if (identical(df_type, "containment") && vcov_adjust_default) vcov_adjust <- FALSE
  # containment df with the adjusted covariance only when asked for explicitly
  if (identical(df_type, "containment") && vcov_adjust && getME(model, "is_REML"))
    vcov_adj <- tryCatch(pm_vcov_kr(model), error = function(e) NULL)
  mp <- mymodelparm(model, df_type = df_type, vcov_adj = vcov_adj)
  if (vcov_adjust) {
    if (is.null(vcov_adj) || inherits(vcov_adj, "pm_satt")) {
      # the default applies to REML lmer fits only; say so only when asked for
      if (!vcov_adjust_default) warning("'vcov_adjust = TRUE' needs an lmer model fitted by REML; ignored.")
    } else {
      mp$vcov <- as.matrix(vcov_adj)
    }
  }
  K <- KK[, mp$estimable, drop=FALSE]
  rK <- ctrmatrix%*%K


  if (all(missing(permlist), missing(Df))){
    df_custom <- pm_df(model, rK)
	  if (!is.null(df_custom)) {
	    Df <- rep_len(as.numeric(df_custom), nrow(rK))
	  } else if (length(mp$df)==1 && mp$df!=0) {
		  Df <- mp$df
	  } else {
		  if (inherits(model, "lme")) {
			vars <- c(unlist(strsplit(modelterm, "\\:")), modelterm)
			Df <- min(terms(model$fixDF)[vars], na.rm=TRUE)
			mDf <- max(terms(model$fixDF)[vars], na.rm=TRUE)
			if (length(vars) > 2) cat("\n", "Denominator degree of freedom for",
			  sQuote(modelterm), "and its marginal terms vary between", sQuote(Df), "and",
			  sQuote(mDf), ".\n","Probabilities will be calculated using", sQuote(Df), "Df.",  "\n")
		  } else if (inherits(model, "lmerMod")) {
			  Df <- df_term(model, modelterm, ctrmatrix = rK, type = df_type, vcov_adj = vcov_adj)
      } else {
        stop("You need provide Df for the model!")
      }
	  }
  }

  cm <- rK%*%mp$coef
  vcovm <- mp$vcov
  vcov.contr <- rK %*% tcrossprod(vcovm, rK)
  ses <- sqrt(diag(vcov.contr))
  t.v <- cm/ses
  nr <- nrow(rK)
  dv <- t(1/ses)
  cor.contr <- as.matrix(vcov.contr * (t(dv) %*% dv))

  if (missing(permlist)) {
   # t.p.value <- 2*pt(-abs(t.v), Df)
	t.p.value <- apply(cbind(t.v, Df), 1, function(x) 2*pt(-abs(x[1]), x[2]))
    t.p.value <- if (adj %in% pm_adjust_extra_methods) {
      pm_adjust_extra(adj, as.vector(t.v), t.p.value, Df, rank = qr(ctrmatrix)$rank, R = cor.contr)
    } else p.adjust(t.p.value, adj)
    out.put <- cbind(cm, ses, t.v, Df, t.p.value)
    colnames(out.put) <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
    rownames(out.put) <- rownames(ctrmatrix)
    attr(out.put,"Note") <- paste("The p-value is adjusted by", sQuote(adj), "method, if p-value = 0 means p-value < 0.0001.")
    }else{
  	nsim <- length(permlist[[1]])
    tValue <- function(x, rK){
    cm <- rK%*%x$coef
    vcovm <- x$vcov
    vcov.contr <- rK %*% tcrossprod(vcovm, rK)
    ses <- sqrt(diag(vcov.contr))
    t.v <- cm/ses
    return(t.v)
  }

  if (nctr==1) {
    per.p <- (sum(sapply(permlist[[1]], function(x) round(abs(tValue(x, rK)), 6) >= round(abs(t.v), 6)))+1)/(nsim + 1)
  } else {
    per.p <- (rowSums(sapply(permlist[[1]], function(x) round(abs(tValue(x, rK)), 6) >= round(abs(t.v), 6)))+1)/(nsim + 1)
  }
  out.put <- cbind(cm, ses, t.v, per.p)
  colnames(out.put) <- c("Estimate", "Std. Error", "t value", "Permuted Pr(>|t|)")
  rownames(out.put) <- rownames(ctrmatrix)
  attr(out.put,"Note") <- paste("The permuted p-value is obtained using", sQuote(nsim), "permutations.")
  }
  return(list("The t tests of the specified contrasts"=round(out.put, 4), "K"=ctrmatrix))
}
