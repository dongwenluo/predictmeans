#' Predicted Means of a Linear Model
#'
#' This function obtains predicted means, SE of means, SED of means, LSDs and
#' plots of means with SE bar or LSD bar for parametric models such as
#' \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' The function also perfomrs pairwise comparisons and permutation tests.
#'
#'
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's
#' predicted mean to be calculated.  The \code{modelterm} must be factors and
#' given exactly as it appears in the printed model, e.g. "A" or "A:B".
#' @param data In some cases, you need to provide the data set used in model
#' fitting, especially when you have applied some variable trnasformation in
#' the model.
#' @param pairwise An option for showing pair-wise LSDs and p-values, or not.
#' The default is FALSE.
#' @param atvar When \code{pairwise = TRUE}, a quoted name indicating within
#' levels of which variable in \code{modelterm} the multiple comparison will be
#' performed.
#' @param adj Name (in "quote") for indicating a method for adjusting p-values
#' of pairwise comparisons.  The choices are "none", "tukey", "holm",
#' "hochberg", "hommel", "bonferroni", "BH", "BY" and "fdr".  The default
#' method is "none". Note that LSD can't be adjusted except for "bonferroni"
#' method.
#' @param Df A degree of freedom for calculating LSD. For the above models, Df
#' is obtained from the function automatically.
#' @param level A significant level for calculating LSD, CI etc. The default
#' value is 0.05.
#' @param covariate A numerical vector to specify values of covariates for
#' calculating predicted means. The default values are the means of the
#' associated covariates.
#' @param meandecr A logical variable to indicate whether to print letters for
#' multiple comparisons by decreasing order of means in the mean_table.  The
#' default is NULL which indicates the mean order follows the associated factor
#' levels.
#' @param letterCI A logical variable to indicate printed letters for multiple
#' comparisons by whether or not CI overlap in the mean_table.  The default is
#' FALSE. Note that the method of examining overlap is more conservative (i.e.,
#' rejects the null hypothesis less often) than the standard method when the
#' null hypothesis is true.
#' @param trans A function object for calculating the back transformed means,
#' e.g. \code{trans=exp}.
#' @param transOff When you use \code{trans=exp(x+1)}, then \code{transOff=1},
#' the default is 0.
#' @param responsen Name (in "quotes") of the back transformed response
#' variable in the \code{model}.
#' @param count An option for indicating the back transformed mean values are
#' counts or not. The default is FALSE.
#' @param prtnum An option for printing covariate information on the screen, or
#' not. The default is TRUE.
#' @param permlist A model parameter list produced by the function
#' \code{permmodels}. When \code{permlist != NULL}, the option \code{Df} will
#' be non-functional. This is a key option for pairwise comparisons via
#' permutation tests.
#' @param ncore Number of core for parallel computing when \code{permlist !=
#' NULL}, the default value is 3.
#' @param ndecimal An option for specifying number of decimal point to be print
#' at predicted means table. The default is 4.
#' @return \item{Predicted Means}{A table of predicted means.} \item{Standard
#' Error of Means}{A table of standard errors of predicted means.}
#' \item{Standard Error of Differences}{Standard errors of differences between
#' predicted means.} \item{LSD}{Least significant differences between predicted
#' means.} \item{Pairwise p-value}{A matrix with t-values above the diagonal
#' and p-values below the diagonal, or matrix of pairwise comparison p-values
#' for each level of \code{atvar}.} \item{mean_table}{A summary of predicted
#' means result including 'Predicted means', 'Standard error', 'Df' and 'CIs'.
#' When \code{trans!=NULL} or \code{trans!=I}, a table of back transformed
#' means with CIs are also shown.} \item{p_valueMatrix}{p_value matrix for
#' pairwise comparison.}
#' @note The \code{predictmeans} function becomes confused if a factor or
#' covariate is changed to the other in a model formula. Consequently, formulae
#' that include calls \code{as.factor}, \code{factor}, or \code{numeric} (e.g.
#' \code{as.factor(income)}) will cause errors. Instead, create the modified
#' variables outside of the model formula (e.g., \code{fincome <-
#' as.factor(income)}) and then use them in the model formula.
#'
#' Factors cannot have colons in level names (e.g., \code{"level:A"}); the
#' \code{predictmeans} function will confuse the colons with interactions;
#' rename levels to avoid colons.
#'
#' For \code{predictmeans} function, it is assumed that methods \code{coef},
#' \code{vcov}, \code{model.matrix}, \code{model.frame} and \code{terms} are
#' available for \code{model}.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @references Maghsoodloo Saeed, Ching-Ying Huang (2010), \emph{Comparing the
#' overlapping of two independent confidence intervals with a single confidence
#' interval for two normal population parameters}, Journal of Statistical
#' Planning and Inference, 140(11), 3295-3305.
#' https://www.sciencedirect.com/science/article/pii/S0378375810002405.
#'
#' Torsten Hothorn, Frank Bretz and Peter Westfall (2008), \emph{Simultaneous
#' Inference in General Parametric Models. Biometrical}, Journal 50(3),
#' 346-363.
#'
#' Welham S., Cullis B., Gogel B., Gilmour A., & Thompson R. (2004),
#' \emph{Prediction in linear mixed models}, Australian and New Zealand Journal
#' of Statistics, 46(3), 325-347.
#' @examples
#'
#'   library(predictmeans)
#'   ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   predictmeans(fm, "nitro", adj="BH")
#'   predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", line=FALSE)
#'   predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH",
#'     barplot=TRUE, line=FALSE)
#'   names(predictout)
#'   print(predictout$predictmeansPlot)
#'   print(predictout$predictmeansBarPlot)
#' @importFrom lme4 devfun2 vcov.merMod
#' @importFrom Matrix forceSymmetric
#' @importFrom parallel clusterExport makeCluster mclapply parLapplyLB stopCluster
#' @importFrom numDeriv hessian
#' @importFrom stats family ftable make.link p.adjust pnorm pt ptukey qt
#' @importFrom stats xtabs

#' @export
predictmeansN <- function (model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL,
                           level=0.05, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = I, transOff=0,
                           responsen=NULL, count=FALSE, prtnum=TRUE, permlist=NULL, ncore=3L, ndecimal=4L) {
  options(scipen=6)
  if (any(missing(model), missing(modelterm))) {
    stop("The arguments 'model', and 'modelterm' must be provided!")
  }
  if (!(modelterm  %in%  attr(terms(model), "term.labels"))) {
    stop(paste("The", modelterm, "must be exactly a term in the model (especially check the order of interaction)."))
  }

  # if (inherits(model, "aovlist")) stop("Plese use model 'lme' instead of ???'aov'!")
  if (inherits(model, "aovlist")) {
    model <- aovlist_lmer(model)
  }
  if (identical(trans, I) && inherits(model, "glm")) {
    trans <- model$family$linkinv  # identical(trans, make.link("log")$linkinv)
    if (model$family$family  %in%  c("poisson", "quasipoisson")) {
      count=TRUE
    }
  }
  if (identical(trans, I) && inherits(model, "glmerMod")) {
    trans <- slot(model, "resp")$family$linkinv
    if (slot(model, "resp")$family$family  %in%  c("poisson", "quasipoisson")) {
      count=TRUE
    }
  }
  if (identical(trans, I) && inherits(model, "glmmTMB")) {
    trans <- model$modelInfo$family$linkinv
    if (model$modelInfo$family$family  %in%  c("poisson", "quasipoisson")) {
      count=TRUE
    }
  }

  vars <- unlist(strsplit(modelterm, "\\:"))
  mdf <- model.frame(model)
  if (any(!base::is.element(vars, names(mdf)[sapply(mdf,is.factor)]))) {
    stop(paste(vars, "must be factor(s)!"))
  }
  # option checking
  if (length(vars)==1) {
    atvar <- NULL
  }
  if (!is.null(permlist) && !any(c("NULL", "") %in% permlist)) {
    pairwise <- TRUE
    if (adj=="tukey") {
      stop("The p-value can't be adjusted by Tukey methd!")
    }
  }
  if (!is.null(atvar) && !any(c("NULL", "") == atvar)) {
    pairwise <- TRUE
  }
  if (adj != "none") {
    pairwise <- TRUE
  }
  # if (letterCI) {pairwise <- TRUE; adj <- "none"; if (is.null(level)) slevel <- level <- 0.166 else slevel <- level}
  # if (is.null(level)) slevel <- level <- 0.05 else slevel <- level
  if (letterCI) {
    pairwise <- TRUE
    slevel <- round(2*(1-pnorm(sqrt(2)*qnorm(1-level/2)/2)), 3)
    level <- pmlevel <- 0.05
  } else {
    pmlevel <- slevel <- level
  }

  if (!is.logical(meandecr)) {
    meandecr <- NULL
  }

  ctr.matrix <- Kmatrix(model, modelterm, covariate, data=data, prtnum=prtnum)
  KK <- ctr.matrix$K
  label <- ctr.matrix$fctnames
  rownames(label) <- rownames(KK)
  response <- ctr.matrix$response
  mp <- mymodelparm(model)
  n.table <- table(mdf[, vars, drop = FALSE])
  ndf <- data.frame(n.table)       ## To obtain info from model

  if (inherits(model, "glmmTMB")) {
    K <- KK
  } else {
    K <- KK[, mp$estimable, drop = FALSE]          # To match coef names
  }
  if (any(ndf$Freq == 0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]      # To delete any missing level in factor
    K <- K[rownames(K) %in% rnTrt,]
    label <- label[rownames(label) %in% rnTrt,]
  }
  pm <- K %*% mp$coef
  vcovm <- mp$vcov
  # ses <- sqrt(diag(K %*% tcrossprod(vcovm, K)))
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
  mt <- data.frame(pm, ses, label)
  LL <- UL <- NULL
  bkmt <- mt  # for back transformed

  mean.table <- round(xtabs(pm ~ ., mt[, c("pm", vars)], drop.unused.levels = TRUE), ndecimal)
  se.table <- round(xtabs(ses ~ ., mt[, c("ses", vars)], drop.unused.levels = TRUE), ndecimal+1)
  mean.table[!(n.table)] <- NA
  se.table[!(n.table)] <- NA
  if (length(vars) > 1) {
    varsnlevel <- numeric(0)
    for (i in vars) varsnlevel[i] <- nlevels(mdf[, i])
    tbvars <- names(sort(varsnlevel, decreasing = TRUE))
    mean.table <- ftable(mean.table, row.vars =tbvars[1], col.var=tbvars[-1])
    se.table <- ftable(se.table, row.vars =tbvars[1], col.var=tbvars[-1])
  }
  if (length(na.omit(unique(se.table))) == 1) {
    se.table <- min(se.table, na.rm=TRUE)
    names(se.table) <- "All means have the same SE"
  }

  nK <- nrow(K)          # To setup various matrix and row, col names
  rnK <- rownames(K)
  if (nK == 1) {
    SED.out <- NA
    LSD <- NA
  } else {
    K_indices <- t(utils::combn(nK, 2))
    K_indices <- K_indices[order(K_indices[ , 2]), ,drop=FALSE]
    Knum_pairs <- nrow(K_indices)
    rows <- rep(1:Knum_pairs, 2)  # The rows are the same for each pair of indices
    cols <- c(K_indices[, 1], K_indices[, 2])  # Indices of the pairs
    values <- c(rep(1, Knum_pairs), rep(-1, Knum_pairs))  # Corrected to 1 for the second species and -1 for the first
    # Create the sparse matrix CMatrix
    CMatrix <- Matrix::sparseMatrix(i = rows, j = cols, x = values, dims = c(Knum_pairs, nK))
    rK <- as.matrix(CMatrix %*% K)
    varn1 <- rnK[K_indices[, 2]]
    varn2 <- rnK[K_indices[, 1]]

    nKK <- nrow(KK)          # To setup various matrix and row, col names
    rnKK <- rownames(KK)
    KK_indices <- t(utils::combn(nKK, 2))
    KK_indices <- KK_indices[order(KK_indices[ , 2]), , drop=FALSE]
    KKvarn1 <- rnKK[KK_indices[, 2]]
    KKvarn2 <- rnKK[KK_indices[, 1]]

    KKvarndiff <- data.frame(matrix(unlist(strsplit(KKvarn1, "\\:")), byrow=T, nrow=length(KKvarn1)),
                             matrix(unlist(strsplit(KKvarn2, "\\:")), byrow=T, nrow=length(KKvarn2)))

    cm <- rK%*%mp$coef
    vcov.contr <- rK %*% tcrossprod(vcovm, rK)
    dses <- sqrt(diag(vcov.contr))
    if (adj == "bonferroni") {
      level <- level/length(dses) # The reson we need slevel
    }

    SED.out <- c(Max.SED = max(dses[!is.nan(dses)]), Min.SED = min(dses[!is.nan(dses)]), Aveg.SED = mean(dses[!is.nan(dses)]))
    dses.df <- data.frame(matrix(unlist(strsplit(varn1, "\\:")), byrow=T, nrow=length(varn1)),
                          matrix(unlist(strsplit(varn2, "\\:")), byrow=T, nrow=length(varn2)), dses)

    if (length(vars) > 1) { # all(length(vars) > 1, SED.out[1]!=SED.out[2])
      dses.m <- matrix(0, nrow=3, ncol=length(vars))
      colnames(dses.m) <- vars
      rownames(dses.m) <- c("Aveg.SED", "Min.SED", "Max.SED")
      for (i in 1:length(vars)) {
        varsindx <- as.character(dses.df[,i]) == as.character(dses.df[, length(vars)+i])
        if (any(varsindx)) { # in case of A/B
          dses.m[,i] <- summary.default(dses.df[varsindx,"dses"])[c(4,1,6)]
        } else {
          dses.m[,i] <- NA
          if (identical(atvar, vars[i])) {
            atvar <- NULL
          }
        }
      }
      attr(SED.out, "For the Same Level of Factor") <- dses.m
    }

    if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
      if (length(Df) == 0) {
        if (inherits(model, "lme")) {
          Df <- terms(model$fixDF)[modelterm]
          names(Df) <- NULL
          mt$Df <- round(Df, 2)
          pairDf <- Df
          Df_diff <- Df
        } else if (inherits(model, "lmerMod")) {
          #  Df <- median(df_term(model, modelterm), na.rm=TRUE)
          L_term <- get_contrasts_type1(model)[[modelterm]]
          Df <- mean(df_term(model, ctrmatrix=L_term))
          if (is.nan(Df)) {
            stop("You need provide Df for this model term!")
          }
          # terms_df <- round(df_term(model, modelterm), 2)
          terms_df <- round(df_term(model, ctrmatrix=K), 2)
          terms_df <- terms_df[!is.na(terms_df)]
          terms_df[terms_df <= 0] <- Df
          mt$Df <- terms_df
          Df_diff <- df_term(model, ctrmatrix=rK)
          Df_diff <- round(Df_diff[!is.na(Df_diff)], 2)
          Df_diff[Df_diff <= 0] <- Df
          pairDf <- round(mean(Df_diff), 2)
          # Df_diff[Df_diff <= 0] <- 1
        } else {
          Df <- mp$df
          mt$Df <- round(Df, 2)
          pairDf <- Df
          Df_diff <- Df
        }

        if (Df == 0) {
          stop("You need provide Df for this model term!")
        }
      } else {
        mt$Df <- Df
        Df_diff <- Df
        pairDf <- Df
      }

      bkmt <- mt  # update bkmt

      LSD <- round(qt(1 - level/2, df = Df) * SED.out, ndecimal+1)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      attr(LSD, "For the Same Level of Factor") <- NULL
      if (length(vars) > 1) {
        rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
        attr(LSD, "For the Same Level of Factor") <- round(qt(1 - level/2, df = Df) * dses.m, ndecimal+1)
      } # end of if LSD
      attr(LSD, "Significant level") <- pmlevel
      attr(LSD, "Degree of freedom") <- round(Df, 2)
    } else {
      LSD <- round(2 * SED.out[1:3], ndecimal+1)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      if (length(vars) > 1) {
        rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
        attr(LSD, "For the Same Level of Factor") <- round(2 * dses.m, ndecimal+1)
      }
      attr(LSD, "Note") <- "This is a approximate LSD (i.e. 2*SED) at 0.05 level."
    }

    if (pairwise) {
      p_valueMatrix <- NULL
      tvm <- t.p.valuem <- LSDm <- Diffm <- matrix(0, ncol = nK, nrow = nK)
      rownames(tvm) <- colnames(tvm) <- rownames(t.p.valuem) <- colnames(t.p.valuem) <- rownames(LSDm) <- colnames(LSDm) <- rnK
      t.v <- cm/dses
      if (all(is.null(permlist) || any(c("NULL", "") %in% permlist), adj=="tukey")) {
        if (inherits(model, "lmerMod")) {
          if (length(Df_diff) == 1) {
            Df_diff <- rep(Df_diff, length(t.v))
          }
          p.tukey <- sapply(1:length(t.v), function(m) ptukey(sqrt(2)*abs(t.v[m]), nK, Df_diff[m], lower.tail=FALSE))
        } else {
          p.tukey <- ptukey(sqrt(2)*abs(t.v), nK, Df, lower.tail=FALSE)
        }
      }
      tvm[upper.tri(tvm)] <- t.v
      if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
        if (inherits(model, "lmerMod")) {
          if (length(Df_diff) == 1) {
            Df_diff <- rep(Df_diff, length(t.v))
          }
          t.p.values <- sapply(1:length(t.v), function(m) 2 * pt(-abs(t.v[m]), Df_diff[m]))
        } else {
          t.p.values <- 2 * pt(-abs(t.v), Df)
        }
      } else {
        nsim <- length(permlist[[1]])
        tValue <- function(x, rK){
          cm <- rK %*% x$coef
          vcovm <- x$vcov
          vcov.contr <- rK %*% tcrossprod(vcovm, rK)
          ses <- sqrt(diag(vcov.contr))
          t.v <- cm/ses
          return(t.v)
        }

        if (.Platform$OS.type == "windows") {
          cl <- makeCluster(ncore)
          clusterExport(cl, c("tValue", "rK", "t.v"), envir = environment())
          t.TableL <- parLapplyLB(cl, permlist[[1]], function(x) round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6))
          stopCluster(cl)
        } else {
          t.TableL <- mclapply(permlist[[1]], function(x) {
            round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6)
          }, mc.cores=ncore)
        }

        t.TableL <- matrix(unlist(t.TableL), ncol = length(t.TableL))
        t.p.values <-  (rowSums(t.TableL)+1)/(nsim+1)

      } # end of if (is.null(permlist))

      if (is.null(atvar) || any(c("NULL", "") == atvar)) {
        if (adj == "tukey") {
          t.p.valuem[upper.tri(t.p.valuem)] <- p.tukey
        } else {
          t.p.valuem[upper.tri(t.p.valuem)] <- p.adjust(t.p.values, adj)
        }
        t.p.valuep <- t.p.valuem    # for plot
        t.p.valuem <- t(t.p.valuem) + tvm
        names(t.p.valuem) <- NULL
        diag(t.p.valuem) <- 1
        if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
          if (!inherits(model, "lmerMod")) {
            attr(t.p.valuem, "Degree of freedom") <- pairDf
          }
          attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, p-value (adjusted by '",
                                            adj, "' method) below the diagonal", sep="")
          LSDm.up <- qt(1-level/2, df = Df)*dses
          LSDm[upper.tri(LSDm)] <- LSDm.up
          Diffm[upper.tri(Diffm)] <- cm
          LSDm <- t(LSDm)+Diffm
          names(LSDm) <- NULL
          attr(LSDm,"Significant level") <- pmlevel
          attr(LSDm,"Degree of freedom") <- Df
          attr(LSDm,"Note") <- paste("LSDs matrix has mean differences (row-col) above the diagonal, LSDs (adjusted by '",
                                     adj, "' method) below the diagonal", sep="")
        } else {
          attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, and ", nsim, " times permutation p-value (adjusted by '", adj, "' method) below the diagonal", sep="")
        } # end of if (is.null(permlist))

        if (!is.null(meandecr) && is.logical(meandecr)) {
          groupRn <- rnK[order(mt$pm, decreasing = meandecr)]
        } else {
          groupRn <- NULL
        }
        t.p.valuemGrp <- t.p.valuem
        t.p.valuemGrp[upper.tri(t.p.valuemGrp)] <- t(t.p.valuemGrp)[upper.tri(t.p.valuemGrp)]
        if (!is.null(groupRn)) {
          t.p.valuemGrp <- t.p.valuemGrp[groupRn, groupRn]
        }
        if (nrow(t.p.valuep) > 2) {
          p_valueMatrix <- round(t(t.p.valuep), 4)
        }
        # if (all(nrow(t.p.valuep) > 2, pplot, plot, prtplt)) {
        #   mtitle <- plottitle
        #   if (is.null(plottitle) || any(c("NULL", "") == plottitle)) {
        #     mtitle <- paste("Level Plot of p-value (adjusted by '", adj, "' method)\n for Pairwise Comparison", sep="")
        #   }
        #   PMplot(t(t.p.valuep), level=pmlevel, legendx=0.69, mtitle=mtitle, newwd=newwd)
        # }
      } else {
        dses.df$tvalue <- t.v
        dses.df$pvalue <- t.p.values
        for (i in which(vars %in% atvar)) { # To find rows relating atvar
          KKvarndiff <- KKvarndiff[which(as.character(KKvarndiff[, i]) == as.character(KKvarndiff[, length(vars) + i])), ]
        }
        atvar.df <- f_loj_krc(KKvarndiff, dses.df, by.x=names(KKvarndiff), by.y=names(KKvarndiff))
        names(atvar.df)[1:length(vars)] <- vars
        for (i in vars) {       # To ensure factor vars  have the same level as original
          atvar.df[,i] <- factor(atvar.df[,i], levels=levels(mdf[, i][[1]]))
        }
        atvar_grp <- do.call("paste", c(atvar.df[, atvar, drop=FALSE], sep="_"))
        atvar.df$adj.pvalue <- unlist(lapply(split(atvar.df, atvar_grp)[unique(atvar_grp)], function(x) {
          t_value <- x$tvalue
          t_valueN <- length(na.omit(t_value))
          if (t_valueN < 2) {
            x$adj.pvalue <- x$pvalue
          } else {
            if (adj=="tukey") {
              x$adj.pvalue <- ptukey(sqrt(2)*abs(x$tvalue), t_valueN, Df, lower.tail=FALSE) # Df need to be updated
            } else {
              x$adj.pvalue <- p.adjust(x$pvalue, adj)
            }
          }
        }))

        rnK.df <- as.data.frame(matrix(unlist(strsplit(rnKK, "\\:")), byrow=T, nrow=nKK)) # To find the suitable names for pm
        colnames(rnK.df) <- vars
        for (i in vars) {       # To ensure factor vars  have the same level as original
          rnK.df[,i] <- factor(rnK.df[,i], levels=levels(mdf[, i][[1]]))
        }
        rnK.df  <- rnK.df [do.call(order, rnK.df[, atvar, drop=FALSE]),]
        resvar <- vars[!vars %in% atvar]   # The rest of vars rather than at var
        rnK.df <- rnK.df[, c(atvar, resvar)]   # To ensure the right matrix name later
        atvar.levels <- unique(do.call("paste", c(rnK.df[, atvar, drop=FALSE], sep=" : ")))
        resvar.levels <- unique(do.call("paste", c(rnK.df[, resvar, drop=FALSE], sep=" : ")))
        rcnplotm <- do.call("paste", c(rnK.df[, , drop=FALSE], sep=" : ")) # row col names of image plot

        # mean table arrange by atvar
        mt.atvar <- mt[do.call(order, mt[, atvar, drop=FALSE]),]
        if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
          bkmt <- mt.atvar[, c(atvar, resvar, "pm", "ses", "Df")]
        } else {
          bkmt <- mt.atvar[, c(atvar, resvar, "pm", "ses")]
        }
        listlength <- length(atvar.levels)
        pmlist <- pmlistTab <- pmlistLetter <- vector("list", listlength)
        indexlength <- nrow(atvar.df)/listlength
        nrow.pm <- length(resvar.levels)

        for (i in 1:listlength) {           # extract pvalues for each factor level
          atvar.pm <- matrix(0, nrow=nrow.pm, ncol=nrow.pm)
          atvar.pm[upper.tri(atvar.pm)] <- atvar.df$adj.pvalue[(indexlength*i-(indexlength-1)):(indexlength*i)]
          rownames(atvar.pm) <- colnames(atvar.pm) <- resvar.levels
          pmlist[[i]] <- t(atvar.pm)
          outtab <- round(as.table(t(atvar.pm)), 4)
          outtab[outtab < 0.0001] <- "0.0001"
          outtab[col(outtab) == row(outtab)] <- 1.0000
          outtab[upper.tri(outtab)] <-""
          Grpatvar.pm <- t(atvar.pm)+ atvar.pm
          if (all(!is.na(outtab))) {
            Grpatvar.pm[is.nan(Grpatvar.pm)] <- -1
            Grpatvar.letter <- multcompLetters(Grpatvar.pm, Letters=LETTERS, threshold=slevel)[resvar.levels]
            Grpatvar.letter[Grpatvar.pm[,1] == -1] <- ""
            outtab <- as.table(cbind(outtab, Group=Grpatvar.letter))
          } else {
            Grpatvar.letterM <- rep(NA, nrow.pm)
            names(Grpatvar.letterM) <- resvar.levels
            for (col_i in 1:nrow.pm) { # in case of missing levels
              atvar.rowname <- resvar.levels[!is.na(Grpatvar.pm[, col_i])]
              if (length(atvar.rowname) > 1) {
                break
              }
            }
            Grpatvar.pm[is.nan(Grpatvar.pm)] <- -1
            if (length(atvar.rowname) == 1) {
              Grpatvar.letterM[atvar.rowname] <- ""
            } else {
              Grpatvar.letterM[atvar.rowname] <- multcompLetters(Grpatvar.pm[atvar.rowname, atvar.rowname],
                                                                 Letters=LETTERS, threshold=slevel)[atvar.rowname]
            }
            Grpatvar.letterM[Grpatvar.pm[,1]==-1] <- ""
            outtab <- as.table(cbind(outtab, Group=Grpatvar.letterM))
            Grpatvar.letter <- Grpatvar.letterM[atvar.rowname]
          }
          pmlistLetter[[i]] <- Grpatvar.letter
          pmlistTab[[i]] <- outtab
        }
        if (nrow.pm > 2) {
          p_valueMatrix <- pmlist
        }
        # if (all(nrow.pm > 2, pplot, plot, prtplt)) {
        #   mtitle <- plottitle
        #   if (is.null(plottitle) || any(c("NULL", "") == plottitle)) {
        #     mtitle <- paste("Adjusted p-value (by '", adj, "' method)\n for Pairwise Comparison at Each Level of '",
        #                     paste(atvar, collapse =" and "), "'", sep="")
        #   }
        #   PMplot(pmlist, level=pmlevel, xylabel=rcnplotm, legendx=0.69, mtitle=mtitle, newwd=newwd)
        # }
      } # if (is.null(atvar))
    }# end of if (pairwise)

    meanTable <- mt
    if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
      meanTable$LL <- meanTable$pm - qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
      meanTable$UL <- meanTable$pm + qt(1 - slevel/2, df = meanTable$Df) * meanTable$ses
    } else {
      meanTable$LL <- meanTable$pm - 2 * meanTable$ses
      meanTable$UL <- meanTable$pm + 2 * meanTable$ses
      slevel <- 0.05
      meanTable$Df <- NA
    }

    rownames(meanTable) <- NULL
    meanTable <- meanTable[c(vars, "pm", "ses", "Df", "LL", "UL")]
    meanTable[c("pm", "LL", "UL")] <- round(meanTable[c("pm", "LL", "UL")], ndecimal)
    meanTable$ses <- round(meanTable$ses, ndecimal+1)
    colnames(meanTable) <- c(vars, "Mean", "SE", "Df", paste("LL(", (1 - slevel) * 100, "%)", sep = ""),
                             paste("UL(", (1 - slevel) * 100, "%)", sep = ""))
  }

  if (!is.null(trans)) {
    Mean <- Trt <- ciPlot <- NULL
    bkmt$Mean <- trans(bkmt$pm)-transOff
    if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)) {
      bkmt$Mean <- exp(bkmt$pm)-transOff
    }
    if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
      bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = bkmt$Df) * bkmt$ses)-transOff
      bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = bkmt$Df) * bkmt$ses)-transOff
    } else {
      bkmt$LL <- trans(bkmt$pm - 2 * bkmt$ses)-transOff
      bkmt$UL <- trans(bkmt$pm + 2 * bkmt$ses)-transOff
    }

    bkmt$pm <- bkmt$ses <- bkmt$Df <- NULL
    nc <- ncol(bkmt)
    bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], ndecimal)
    if (count) {
      bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 0)
      bkmt[, (nc - 2):nc][bkmt[, (nc - 2):nc] < 0] <- 0
    }

    rownames(bkmt) <- NULL
    bkmt$Trt <- NULL
    if (!(is.null(atvar) || any(c("NULL", "") == atvar))) {
      all_var_names <- colnames(meanTable)
      meanTable <- meanTable[c(atvar, resvar, setdiff(all_var_names, vars))]
      meanTable <- meanTable[do.call(order, meanTable[, atvar, drop=FALSE]), ]
      # meanTable <- with(meanTable, meanTable[order(eval(parse(text=c(atvar, resvar)[length(vars):1]))),])
    }

    # if (all(!identical(trans, function(x) x, ignore.environment=TRUE), !identical(trans, I, ignore.environment=TRUE))) {
    if (any(meanTable$Mean != bkmt$Mean)) {
      meanTable <- cbind(meanTable, round(bkmt[, (ncol(bkmt)-2):ncol(bkmt)], ndecimal))
      colnames(meanTable)[(ncol(meanTable)-2):ncol(meanTable)] <-  c("Bk_Mean", paste("Bk_LL(", (1 - slevel) * 100, "%)", sep = ""), paste("Bk_UL(", (1 - slevel) * 100, "%)", sep = ""))
    }

    if (pairwise) {
      if ((is.null(atvar) || any(c("NULL", "") == atvar))) {
        if (!is.null(meandecr) && is.logical(meandecr)) {
          meanTable <- meanTable[order(meanTable$Mean, decreasing = meandecr), ]
        }
        if (letterCI) {
          meanTable$LetterGrp <- ci_mcp(meanTable[, grepl("^LL", names(meanTable))], meanTable[, grepl("^UL", names(meanTable))])
          meanTable <- list(Table=meanTable,
                            Note=paste("Letter-based representation of pairwise comparisons based on (",
                                       (1 - slevel) * 100, "%) CI at ", level, " significant level.", sep=""))
        } else {
          t.p.valuemGrp[is.nan(t.p.valuemGrp)] <- -1
          meanTable$LetterGrp <- na.omit(multcompLetters(t.p.valuemGrp, Letters=LETTERS, threshold=slevel))
          meanTable$LetterGrp[is.nan(meanTable$SE)] <- ""
          meanTable <- list(Table=meanTable,
                            Note=paste("Letter-based representation of pairwise comparisons at significant level '", slevel, "'", sep=""))
        }
      } else {
        if (letterCI) {
          meanTable$LetterGrp <- unlist(lapply(split(meanTable, meanTable[, atvar]),
                                               function (x) {
                                                 ci_mcp(x[,grepl("^LL", names(meanTable))], x[,grepl("^UL", names(meanTable))])
                                               }))
          meanTable <- list(Table=meanTable,
                            Note=paste("Letter-based representation of pairwise comparisons based on (",
                                       (1 - slevel) * 100, "%) CI at ", level, " significant level for each ", atvar, " group.", sep=""))
        } else {
          # pmlistLetter <- pmlistLetter[unlist(lapply(pmlistLetter, function(x) !all(x=="")))]
          pmlistLetter <- unlist(pmlistLetter)[!is.na(unlist(pmlistLetter))]
          if (length(pmlistLetter)!=dim(meanTable)[1]) {
            pmlistLetter <- pmlistLetter[pmlistLetter!=""]
          }
          meanTable$LetterGrp <- pmlistLetter
          meanTable <- list(Table=meanTable,
                            Note=paste("Letter-based representation of pairwise comparisons at significant level '",
                                       slevel, "' at each level of ", paste(atvar, collapse =" and "), sep=""))
        }
      }
    }
  }

  if (pairwise) {
    if (is.null(atvar) || any(c("NULL", "") == atvar)) {
      if (all(is.null(permlist) || any(c("NULL", "") %in% permlist), any(c("none", "bonferroni") == adj))) {
        outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                           "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm,ndecimal+1),
                           "Pairwise p-value" = round(t.p.valuem, 4), mean_table=meanTable, p_valueMatrix=p_valueMatrix)
        class(outputlist) = "pdmlist"
        return(outputlist)
      } else {
        outputlist <- list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                           "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise p-value" = round(t.p.valuem, 4),
                           mean_table=meanTable, p_valueMatrix=p_valueMatrix)
        class(outputlist) = "pdmlist"
        return(outputlist)
      }
    } else {
      outputlist <- vector("list", listlength+7)
      outputlist[[1]] <- mean.table
      outputlist[[2]] <- se.table
      outputlist[[3]] <- SED.out
      outputlist[[4]] <- LSD
      outputlist[[5]] <- paste("For variable '", paste(resvar, collapse =" and "), "' at each level of '", paste(atvar, collapse =" and "), "'", sep="")
      for (i in 6: (listlength+5))  outputlist[[i]] <- pmlistTab[[i-5]]
      outputlist[[listlength+6]] <- meanTable
      outputlist[[listlength+7]] <- p_valueMatrix
      if (is.null(permlist) || any(c("NULL", "") %in% permlist)) {
        names(outputlist)<- c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                              "LSD", paste("Pairwise comparison p-value (adjusted by '", adj, "' method)", sep=""),
                              atvar.levels, "mean_table", "p_valueMatrix")[1:length(outputlist)]
      } else {
        names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
                                 "Approximated LSD"), paste("Pairwise '", nsim, "' times permuted p-value (adjusted by '", adj, "' method)", sep=""),
                               atvar.levels, "mean_table", "p_valueMatrix")[1:length(outputlist)]
      }
      class(outputlist) <- "pdmlist"
      return(outputlist)
    }
  } else {
    outputlist=list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                    "Standard Error of Differences" = SED.out, LSD = LSD, mean_table=meanTable)
    class(outputlist) <- "pdmlist"
    return(outputlist)
  }
  #  }
}
