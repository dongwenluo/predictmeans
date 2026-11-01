permmodels <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                       test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
					   data=NULL, fo=NULL, prt=TRUE, ncore=3, seed) {
  
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) test.statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test.statistic <- "Chisq"
  if (type %in% c("I", "1")) test.statistic <- "F"	  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "glmerMod")) stop("This function is not applied to 'glmer' yet!")
  
  if (any(inherits(model, "lm"), inherits(model, "aov"), inherits(model, "glm"))) {
    if (!is.null(data)) mod_df <- data else {
    if (inherits(model, "glm")) mod_df <- as.data.frame(model$data) else mod_df <- as.data.frame(model.frame(model))
	}
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]    
    if (grepl("[,]", yname)) {
      yname <- unlist(strsplit(yname, "[,] "))[2]
      yname <- gsub("\\)", "", unlist(strsplit(yname, " - "))) # ynames
    }
	diff_yname <- setdiff(colnames(mod_df), yname)
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      mod_df[sample(1:nrow(mod_df)), yname, drop=FALSE]
    }, simplify = !(length(yname) > 1))
		
	if (!is.null(fo)) {
	      permmod <- lapply(permy, function(x) {
        if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
        rfitmodel <- try(pm_refit(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
		  coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
    
    } else{	
	if (.Platform$OS.type=="windows") {
      cl <- makeCluster(ncore)
      clusterExport(cl, c("mymodelparm", "mymodelparm.default", "model", "mod_df", "yname", "diff_yname", "type", "test.statistic"), envir = environment()) 
      permmod <- parLapplyLB(cl, permy, function(x) {
		if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x		
        rfitmodel <- try(pm_refit(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
		  coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      stopCluster(cl)
    }else{  
      permmod <- mclapply(permy, function(x) {
        if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
        rfitmodel <- try(pm_refit(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
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
    model <- pm_refit(model, data=as.data.frame(pm_get_data(model)), na.action=na.omit, check=TRUE)
    mod_df <- as.data.frame(pm_get_data(model))
    X <- model.matrix(model)
    if (inherits(model, "gls")) beta <- coef(model) else beta <- nlme::fixef(model)
    # check for columns dropped from model
    col_names <- names(beta)
    if (ncol(X) != length(col_names)) X <- X[,col_names,drop=FALSE]    
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]
    y <- nlme::getResponse(model)
    X_y <- cbind(X, y)
    V_list <- build_Sigma_mats(model, sigma_scale =FALSE)
    ord_X_y <- get_order(A=V_list, B=X_y)
    X_y <- ord_X_y[[1]]
    mod_df <- mod_df[ord_X_y[[2]], ]
    xbeta <- as.vector(X_y[, col_names]%*%beta)
    errors <- X_y[, "y"] - xbeta
    V_matrix <- Matrix::bdiag(V_list)
    Ut <- t(chol(V_matrix))
    
    #Weighting the residuals.
    wterrors <- solve(Ut, errors)    
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      as.data.frame(perm_y <- xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))
    })
    
	if (!is.null(fo)) {
	        permmod <- lapply(permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(pm_refit(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
		  coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
	
	}else{
    if (.Platform$OS.type=="windows") {
      cl <- makeCluster(ncore)
      clusterEvalQ(cl, library(nlme))
      clusterExport(cl, c("mymodelparm", "mymodelparm.lme", "mymodelparm.gls","mymodelparm.default", "model", "mod_df", "yname", "type", "test.statistic"), envir = environment()) 
      permmod <- parLapplyLB(cl, permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(pm_refit(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
		  coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      stopCluster(cl)
    }else{  
      permmod <- mclapply(permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(pm_refit(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
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
	model <- pm_refit(model, data=mod_df, check=TRUE)	
    theta <- getME(model, "theta")
    fixef <- fixef(model) 
    Lambda <- getME(model, "Lambda")
    # V = Z Lambda Lambda' Z' + I, kept sparse; Ut is its lower Cholesky factor.
    ZL <- getME(model, "Z") %*% Lambda
    V <- Matrix::forceSymmetric(Matrix::tcrossprod(ZL) + Matrix::Diagonal(nrow(ZL)))
    Ut <- Matrix::t(Matrix::chol(V))
    xbeta <- as.vector(getME(model, "X")%*%fixef)
    if (inherits(model, "glmerMod")) {
      errors <- slot(model, "resp")$family$linkfun(getME(model, "y")) - xbeta 
	  }else errors <- getME(model, "y") - xbeta
    
    #Weighting the residuals.
    wterrors <- Matrix::solve(Ut, errors)
    
    if(!missing(seed)) set.seed(seed)  
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
		if (slot(model, "resp")$family$family == "Gamma") {
		  perm_y[perm_y<=0] <- 0		
		}
		as.data.frame(perm_y)
		}else as.data.frame(perm_y <- xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))
    })
    
	if (!is.null(fo)) {
	      permmod <- lapply(permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
		  coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
	
	}else{
    if (.Platform$OS.type=="windows") {
      cl <- makeCluster(ncore)	  
      clusterEvalQ(cl, library(lme4))
      clusterExport(cl, c("mymodelparm", "mymodelparm.lmerMod", "mymodelparm.default", "model", "type", "test.statistic"), envir = environment())       
      permmod <- parLapplyLB(cl, permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")){
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
		  coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      stopCluster(cl)
    }else{  
      permmod <- mclapply(permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
          mp <- mymodelparm(rfitmodel, df = 0)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
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
    if (nrow(model_tTable)==1) {
      Tper.p <- (sum(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    }else{
      Tper.p <- (rowSums(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    }

	if ("(Intercept)" %in% rownames(model_tTable)) Tper.p[1] <- NA
	COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p),4)
    if (prt) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEFFICENTS)
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "permutations.\n")
    }
  }else COEFFICENTS <- NULL
  
  if (type %in% c("I", "1")) {
    key_anova <- anova(model)
    Fvalue <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
    if (class(model)[1] %in% c("glm")) Fvalue <- "Deviance"
    if (nrow(key_anova)==1) {
      Fper.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    }else{
      Fper.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    } 
    if ("(Intercept)" %in% rownames(key_anova)) Fper.p[1] <- NA	
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=Fper.p),4)
  }else{
    key_anova <- pm_anova(model, type=type, test.statistic=test.statistic)
    if (class(model)[1] %in% c("lm", "aov")) stats_value <- "F value"  else stats_value <-  test.statistic 
    if (class(model)[1] %in% c("glm")) stats_value <- switch(test.statistic, LR = "LR Chisq", Wald = "Chisq", F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
    if (nrow(key_anova)==1) {
      stats.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }else{
      stats.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }  
	if ("(Intercept)" %in% rownames(key_anova)) stats.p[1] <- NA
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4) 
    attr(ANOVA, "heading") <- attr(key_anova, "heading")
  }
  
  ANOVA[is.na(ANOVA)] <-""
  if (prt) {
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep=""))
    print(ANOVA)
    cat(paste("\nNote: Perm_p_value of", test.statistic, "test is obtained using"), sQuote(nperm), "permutations.\n\n")
  }
  permlist <- vector("list", nperm)
  for (i in 1:nperm) {
    permlist[[i]] <- permmod[[i]][[1]]
  }
  return(invisible(list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA)))
} 
