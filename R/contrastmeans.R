#' Linear Contrast Tests for a Linear Model
#'
#' Performs t-tests (or permuted t-tests) of specified contrasts for linear
#' models obtained from functions \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, or \code{lmer}.
#'
#'
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's
#' contrast to be calculated.  The \code{modelterm} must be given exactly as it
#' appears in the printed model, e.g. "A" or "A:B".
#' @param ctrmatrix A specified contrast matrix. If \code{ctrmatrix} is
#' missing, the programe will ask user to enter it.
#' @param ctrnames Names of the specified contrasts, e.g. c("A vs D", "C vs B",
#' ...)
#' @param adj Name (in "quote") for indicating a method for adjusting p-values
#' of pairwise comparisons.  The choices are "none", "tukey", "holm",
#' "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr".  The default
#' method is "none".
#' @param Df A denominator degree of freedom for \code{modelterm}. (For
#' \code{glmer} models the \code{Df} needs to be specified, while for the other
#' models, \code{Df} is obtained from the fitted model automatically).
#' @param permlist A model parameter list containing \code{nsim} parameters
#' produced by the function \code{permmodels}. When \code{permlist != NULL},
#' the option \code{Df} will be non-functional. This is a key option for the
#' permutation test.
#' @return There are two components in the output which are \item{Table}{A
#' table showing t-test results for the specified linear contrasts.} \item{K}{A
#' contrast matrix.}
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @references Torsten Hothorn, Frank Bretz and Peter Westfall (2008),
#' \emph{Simultaneous Inference in General Parametric Models. Biometrical},
#' Journal 50(3), 346--363.
#' @examples
#'
#' library(predictmeans)
#' # ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # library(lme4)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'
#' ## Not run:
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
#' # contrastmeans(fm, "nitro:Variety")
#'
#' # 2. Construct the contrast matrix directly
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)
#' @export
contrastmeans <- function(model, modelterm, ctrmatrix, ctrnames=NULL, adj="none", Df, permlist) {
  options(scipen=6)

  if (inherits(model, "aovlist")) {
    model <- aovlist_lmer(model)
  }

  if (is.null(ctrnames) || all(ctrnames%in%c("NULL", ""))) {
    ctrnames <- NULL
  }
  K <- Kmatrix(model, modelterm)$K
  termsLabel <- rownames(K)
  if (missing(ctrmatrix)) {
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
    colnames(ctrmatrix) <- termsLabel
    rownames(ctrmatrix) <- ctrnames
    nctr <- nrow(ctrmatrix)
  }
  if (all((is.null(dim(ctrmatrix)) | dim(ctrmatrix)[1]==1), missing(permlist)))  {
    adj <- "none"
  }
  rK <- ctrmatrix%*%K

  mp <- mymodelparm(model)
  if (all(missing(permlist), missing(Df))) {
    if (length(mp$df)==1 && mp$df!=0) {
      Df <- mp$df
    } else {
      if (inherits(model, "lme")) {
        vars <- c(unlist(strsplit(modelterm, "\\:")), modelterm)
        Df <- min(terms(model$fixDF)[vars], na.rm=TRUE)
        mDf <- max(terms(model$fixDF)[vars], na.rm=TRUE)
        if (length(vars) > 2) {
          cat("\n", "Denominator degree of freedom for",
              sQuote(modelterm), "and its marginal terms vary between", sQuote(Df), "and",
              sQuote(mDf), ".\n","Probabilities will be calculated using", sQuote(Df), "Df.",  "\n")
        }
      } else if (inherits(model, "lmerMod")) {
        Df <- df_term(model, ctrmatrix = rK)
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
    t.p.value <- p.adjust(t.p.value, adj)
    out.put <- cbind(cm, ses, t.v, Df, t.p.value)
    colnames(out.put) <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
    rownames(out.put) <- ctrnames
    attr(out.put,"Note") <- paste("The p-value is adjusted by", sQuote(adj), "method, if p-value = 0 means p-value < 0.0001.")
  } else {
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
      per.p <- (sum(sapply(permlist[[1]], function(x) abs(tValue(x, rK)) > abs(t.v)))+1)/(nsim + 1)
    } else {
      per.p <- (rowSums(sapply(permlist[[1]], function(x) abs(tValue(x, rK)) > abs(t.v)))+1)/(nsim + 1)
    }
    out.put <- cbind(cm, ses, t.v, per.p)
    colnames(out.put) <- c("Estimate", "Std. Error", "t value", "Permuted Pr(>|t|)")
    rownames(out.put) <- ctrnames
    attr(out.put,"Note") <- paste("The permuted p-value is obtained using", sQuote(nsim), "permutations.")
  }
  return(list("The t tests of the specified contrasts"=round(out.put, 4), "K"=ctrmatrix))
}
