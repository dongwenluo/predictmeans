#' Permutation Test of Linear Model
#' 
#' This function provides permutation t-tests for coefficients of (fixed)
#' effects and permutation F-tests for the terms in a linear model such as
#' \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' 
#' 
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, and \code{lmer}.
#' @param nperm The number of permutations. The default is 999.
#' @param type type of ANOVA test, "I", "II", "III", 1, 2, or 3. Roman numerals
#' are equivalent to the corresponding Arabic numerals.
#' @param test.statistic For type I ANOVA, F test is applied to all models,
#' while for type II and III ANOVA, F test is performed for \code{lm}, Chisq
#' test for \code{lme} and \code{gls} model, Chisq or F tests for \code{lmer}
#' model and Likelihood ratio, Wald or F tests for \code{glm} model.
#' @param exact A logical variable to indicate whether or not exact no. of
#' permutations will be used (applicable only to free the permutation case).
#' The default is FALSE.
#' @param data In some cases, you need to provide the data set used in model
#' fitting, especially when you have applied some variable trnasformation in
#' the model.
#' @param fo A model formula used in the \code{model}; \code{fo!=NULL} when the
#' formula is specified by function \code{formula}.
#' @param prt A logical variable to indicate whether or not to print output on
#' the screen. The default is TRUE.
#' @param ncore Number of core for parallel computing, the default value is 3.
#' @param seed Specify a random number generator seed, for reproducible
#' results.
#' @return The function produces permutation t-test table for coefficients of
#' (fixed) effects, permutation ANOVA table for model terms and a model
#' parameter list \code{permlist}, a list containing \code{nsim=4999} times
#' permutation refitted \code{model} parameters which are used in functions
#' \code{predictmeans} and \code{contrastmeans}.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' 
#' # # Not run for simplifying process of submiting pkg to CRAN
#' # library(predictmeans)
#' # Oats$nitro <- factor(Oats$nitro) 
#' # fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # # library(lme4)
#' # # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' # 
#' # # Permutation Test for model terms
#' # system.time(
#' #   permlme <- permmodels(model=fm, nperm=999)
#' # )  
#' # 
#' # # Permutation Test for multiple comparisons
#' # predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH", 
#' #              permlist=permlme, plot=FALSE)
#' # 
#' # # Permutation Test for specified contrasts
#' # cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#' #             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' # contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
#' 
permmodels <- function(model, nperm=999, type=c("I", "II", "III", 1, 2, 3),
                       test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
                       data=NULL, fo=NULL, prt=TRUE, ncore=3, seed) {
  
  options(scipen=6)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) {
    test.statistic <- "LR"
  }
  if (any(inherits(model, "gls"), inherits(model, "lme"))) {
    test.statistic <- "Chisq"
  }
  if (type %in% c("I", "1")) {
    test.statistic <- "F"
  }
  if (class(model)[1]=="aovlist") {
    stop("Plese use model 'lme' instead of 'aov'!")
  }
  if (inherits(model, "glmerMod")) {
    stop("This function is not applied to 'glmer' yet!")
  }
  
  if (any(inherits(model, "lm"), inherits(model, "aov"), inherits(model, "glm"))) {
    if (!is.null(data)) {
      mod_df <- data 
    } else {
      if (inherits(model, "glm")) {
        mod_df <- as.data.frame(model$data) 
      } else {
        mod_df <- as.data.frame(model.frame(model))
      }
    }
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]    
    if (grepl("[,]", yname)) {
      yname <- unlist(strsplit(yname, "[,] "))[2]
      yname <- gsub("\\)", "", unlist(strsplit(yname, " - "))) # ynames
    }
    diff_yname <- setdiff(colnames(mod_df), yname)
    if(!missing(seed)) {
      set.seed(seed)  
    }
    permy <- replicate(nperm, {
      mod_df[sample(1:nrow(mod_df)), yname, drop=FALSE]
    }, simplify = !(length(yname) > 1))
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        if (length(yname) > 1) {
          mod_df <- cbind(mod_df[, diff_yname], x) 
        } else {
          mod_df[, yname] <- x
        }
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) {
            aT <- anova(rfitmodel) 
          } else {
            aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          }
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
    } else {	
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) {
          clusterEvalQ(cl, library(car))
        }
        clusterExport(cl, c("mymodelparm", "mymodelparm.default", "model", 
                            "mod_df", "yname", "diff_yname", "type", 
                            "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          if (length(yname) > 1) {
            mod_df <- cbind(mod_df[, diff_yname], x) 
          } else { 
            mod_df[, yname] <- x
          }		
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1] == class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) { 
              aT <- anova(rfitmodel)
            } else {
              aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            }
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        
        stopCluster(cl)
      } else {  
        permmod <- mclapply(permy, function(x) {
          if (length(yname) > 1) {
            mod_df <- cbind(mod_df[, diff_yname], x) 
          } else {
            mod_df[, yname] <- x
          }
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1] == class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) {
              aT <- anova(rfitmodel) 
            } else {
              aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            }
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }    
  }
  
  if (any(inherits(model, "gls"), inherits(model, "lme"))) { 
    # fm$call$fixed
    # fm$call$random  
    model <- update(model, na.action=na.omit)
    mod_df <- as.data.frame(nlme::getData(model))
    X <- model.matrix(model)
    if (inherits(model, "gls")) {
      beta <- coef(model) 
    } else {
      beta <- nlme::fixef(model)
    }
    # check for columns dropped from model
    col_names <- names(beta)
    if (ncol(X) != length(col_names)) {
      X <- X[,col_names,drop=FALSE]
    }
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]
    y <- nlme::getResponse(model)
    X_y <- cbind(X, y)
    V_list <- build_Sigma_mats(model, sigma_scale =FALSE)
    ord_X_y <- get_order(A=V_list, B=X_y)
    X_y <- ord_X_y[[1]]
    mod_df <- mod_df[ord_X_y[[2]], ]
    xbeta <- as.vector(X_y[, col_names] %*% beta)
    errors <- X_y[, "y"] - xbeta
    V_matrix <- Matrix::bdiag(V_list)
    Ut <- t(chol(V_matrix))
    wt <- solve(Ut)
    
    #Weighting the residuals.
    wterrors <- wt%*%errors    
    
    if (!missing(seed)) {
      set.seed(seed)  
    }
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      as.data.frame(perm_y <- xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) {
            aT <- anova(rfitmodel) 
          } else {
            aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          }
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
    } else {
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) {
          clusterEvalQ(cl, library(car)) 
        } else {
          clusterEvalQ(cl, library(nlme))
        }
        clusterExport(cl, c("mymodelparm", "mymodelparm.lme", "mymodelparm.gls",
                            "mymodelparm.default", "model", "mod_df", "yname", 
                            "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) {
              aT <- anova(rfitmodel) 
            } else {
              aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            }
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      } else {  
        permmod <- mclapply(permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) {
              aT <- anova(rfitmodel) 
            } else {
              aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            }
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }
  }
  
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"))) { 
    lmer_var_names <- all.vars(model@call$formula)
    lmer_var_names <- lmer_var_names[lmer_var_names!="pi"]
    mod_df <- as.data.frame(nlme::getData(model))
    mod_df <- na.omit(mod_df[, lmer_var_names])	
    model <- update(model, data=mod_df)	
    theta <- getME(model, "theta")
    fixef <- fixef(model) 
    Lambda <- getME(model, "Lambda")
    Lambdac <- Matrix::tcrossprod(Lambda)
    V <- getME(model, "Z") %*% Lambdac %*% getME(model, "Zt")+diag(dim(getME(model, "Z"))[1])
    Ut <- t(chol(V))
    wt <- solve(Ut)
    xbeta <- as.vector(getME(model, "X") %*% fixef)
    if (inherits(model, "glmerMod")) {
      errors <- slot(model, "resp")$family$linkfun(getME(model, "y")) - xbeta 
    } else {
      errors <- getME(model, "y") - xbeta
    }
    
    #Weighting the residuals.
    wterrors <- wt%*%errors
    
    if (!missing(seed)) {
      set.seed(seed)  
    }
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      if (inherits(model, "glmerMod")) {
        perm_y <- slot(model, "resp")$family$linkinv(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) 
        if (slot(model, "resp")$family$family == "poisson") { 
          perm_y[perm_y<=0] <- 0
          perm_y <- round(perm_y, 0)
        }
        if (slot(model, "resp")$family$family == "binomial") {
          perm_y[perm_y<=0] <- 0
          perm_y[perm_y>=1] <- 1		
        }
        if (slot(model, "resp")$family$family == "gamma") {
          perm_y[perm_y<=0] <- 0		
        }
        as.data.frame(perm_y)
      } else {
        as.data.frame(perm_y <- xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))
      }
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) {
            aT <- anova(rfitmodel) 
          } else {
            aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          }
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    } else {
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        clusterEvalQ(cl, library(lme4))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lmerMod", "mymodelparm.default", 
                            "model", "type", "test.statistic"), envir = environment())       
        permmod <- parLapplyLB(cl, permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) {
              aT <- anova(rfitmodel) 
            } else {
              aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            }
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      } else {  
        permmod <- mclapply(permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) {
              aT <- anova(rfitmodel) 
            } else {
              aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            }
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }	
  }    
  permmod <- permmod[!(sapply(permmod,is.null))]
  nperm <- length(permmod)
  if (nperm==0) stop("permutation can't produce useful data!")
  if (!(inherits(model, "aov"))) {
    # tTable <- function(x) { # function to construct t-table for the model
    # if (class(model)[1]=="lmerMod"){
    # mp <- mymodelparm(x)
    # tTable <- cbind(mp$coef, sqrt(base::diag(mp$vcov)), mp$coef/sqrt(base::diag(mp$vcov)))
    # colnames(tTable) <- c("Estimate", "Std. Error", "t value")
    # return(round(tTable, 4))		
    # }else{
    # summ <- summary(x)
    # if (class(x)[1] %in% c("lme", "gls")) {
    # tTable <- summ$tTable
    # }else{
    # tTable <- coef(summ)
    # }
    # return(tTable)
    # }
    # }
    
    # model_tTable <- tTable(model)
    # Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    
    # if (nrow(model_tTable)==1) {
    # Tper.p <- (sum(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    # }else{
    # Tper.p <- (rowSums(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    # }
    
    model_tTable <- coef(summary(model))
    Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    if (nrow(model_tTable) == 1) {
      Tper.p <- (sum(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    } else {
      Tper.p <- (rowSums(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    }
    
    if ("(Intercept)" %in% rownames(model_tTable)) {
      Tper.p[1] <- NA
    }
    COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p), 4)
    if (prt) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEFFICENTS)
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "permutations.\n")
    }
  } else {
    COEFFICENTS <- NULL
  }
  
  if (type %in% c("I", "1")) {
    key_anova <- anova(model)
    Fvalue <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
    if (class(model)[1] %in% c("glm")) {
      Fvalue <- "Deviance"
    }
    if (nrow(key_anova) == 1) {
      Fper.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    } else {
      Fper.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    } 
    if ("(Intercept)" %in% rownames(key_anova)) {
      Fper.p[1] <- NA	
    }
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=Fper.p), 4)
  } else {
    key_anova <- car::Anova(model, type=type, test.statistic=test.statistic)
    if (class(model)[1] %in% c("lm", "aov")) {
      stats_value <- "F value"  
    } else {
      stats_value <-  test.statistic 
    }
    if (class(model)[1] %in% c("glm")) {
      stats_value <- switch(test.statistic, 
                            LR = "LR Chisq", 
                            Wald = "Chisq", 
                            F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
    }
    if (nrow(key_anova) == 1) {
      stats.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    } else {
      stats.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }  
    if ("(Intercept)" %in% rownames(key_anova)) {
      stats.p[1] <- NA
    }
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4) 
    attr(ANOVA, "heading") <- attr(key_anova, "heading")
  }
  
  ANOVA[is.na(ANOVA)] <-""
  if (prt) {
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) {
      cat(paste(attr(ANOVA, "heading"), "\n", sep=""))
    }
    print(ANOVA)
    cat(paste("\nNote: Perm_p_value of", test.statistic, "test is obtained using"), sQuote(nperm), "permutations.\n\n")
  }
  permlist <- vector("list", nperm)
  for (i in 1:nperm) {
    permlist[[i]] <- permmod[[i]][[1]]
  }
  return(invisible(list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA)))
} 
