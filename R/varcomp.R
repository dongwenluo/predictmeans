varcomp <- function(model, ci=TRUE, level=0.95) {
  # vcov_lmerMod(model, full = TRUE, ranpar = "var"), vcov_lmerMod(model, full = TRUE, ranpar = "sd")?
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"), inherits(model, "lmerModLmerTest"))) {
    varcomp <- as.data.frame(VarCorr(model), order = "lower.tri")
    vcov_index <- -1:-length(fixef(model))
	if (inherits(model, "glmerMod")) {
	  random_vcov <- as.matrix(vcov_glmerMod(model, full = TRUE, ranpar = "var")[vcov_index, vcov_index]) 
	  ci <- FALSE
	}else random_vcov <- as.matrix(vcov_lmerMod(model, full = TRUE, ranpar = "var")[vcov_index, vcov_index])
    if (inherits(model, "glmerMod")) random_terms <- gsub("cov_", "", rownames(random_vcov)) else random_terms <- gsub("cov_", "", colnames(random_vcov)) 
    varcomp$SE <- sqrt(diag(random_vcov))
	if (ci){
    varcomp <- cbind(varcomp, confint(model, level=level)[1:length(random_terms), ])
    kk <- abs(varcomp$sdcor^2 - varcomp$vcov) > 0.0000000000001
    varcomp[, 7] <- ifelse(!kk, ifelse(varcomp[, 7] < 0, 0, varcomp[, 7]^2), varcomp[, 7]*varcomp$vcov/varcomp$sdcor)
    varcomp[, 8] <- ifelse(!kk, varcomp[, 8]^2, varcomp[, 8]*varcomp$vcov/varcomp$sdcor)	
    varcomp <- round(varcomp[,c(4,6:8)], 4)
	}else{
	  varcomp <- round(varcomp[,c(4,6)], 4)	
	}	
    rownames(varcomp) <- random_terms
  }else if(inherits(model, "lme")){    
    pm_need("lmeInfo", "variance components of 'lme' models")
    vcov=unlist(lmeInfo::extract_varcomp(model))
    SE=sqrt(diag(lmeInfo::varcomp_vcov(model)))
    
    # cbind(Std.error, vcov)
    confint_list <- try(intervals(model, which="var-cov", level=level), silent = TRUE)
    if (!inherits(confint_list, "try-error")){
      confint_matrix <- as.matrix(rbind(do.call("rbind", confint_list[[1]]), confint_list$sigma))
      confint_names <- c(as.vector(sapply(confint_list[[1]], rownames)), "residual")
      rownames(confint_matrix) <- confint_names
      confint_matrixN <- confint_matrix
      confint_matrixN["residual", ] <- confint_matrix["residual",]^2 
      sd_detect <- grepl("^sd", confint_names, fixed = FALSE)
      if (any(sd_detect)) {
        sd_names <- sub("\\)$", "", sub("^sd\\(", "", confint_names[sd_detect]))
        confint_matrixN[sd_detect, ] <- confint_matrix[sd_detect,]^2
      }
      
      cor_detect <- grepl("^cor", confint_names, fixed = FALSE) 
	  est. <- NULL
      if (any(cor_detect)) {
        cor_names <- sub("\\)$", "", sub("^cor\\(", "", confint_names[cor_detect]))
        sd_value <- confint_matrix[sd_detect, "est."]
        names(sd_value) <- sd_names
        sd_prod <- sapply(strsplit(cor_names, ","), function(x){sdv <- sd_value[x]; sdv[1]*sdv[2]})
        # print(sd_prod)
        # print(confint_matrixN[cor_detect, ])
        confint_matrixN[cor_detect, ] <- confint_matrix[cor_detect,]*sd_prod
      }
      rownames(confint_matrixN) <- as.character(round(confint_matrixN[,"est."], 4))
      confint_matrixN <- confint_matrixN[as.character(round(vcov, 4)),]
      varcomp <- cbind(vcov, SE, confint_matrixN)
	  varcomp <- subset(varcomp, select=-c(est.))
	  rownames(varcomp) <- sub("^Tau.", "", rownames(varcomp))
    }else{
      varcomp <- cbind(vcov, SE)      
    }
  }else if(inherits(model, "glmmTMB")){
    theta_sd <- sqrt(diag(vcov(model, full=TRUE)))
    theta_sd <- theta_sd[grepl("theta", names(theta_sd))]
    theta <- getME(model, "theta")
    vcov <- exp(2*theta)
    SE=theta_sd*2*exp(2*theta)
    theta_UL <- theta+qnorm(p=(1-level)/2, lower.tail=FALSE)*theta_sd
    upper <- theta_sd*2*exp(2*theta_UL)
    theta_LL <- theta-qnorm(p=(1-level)/2, lower.tail=FALSE)*theta_sd
    lower <- theta_sd*2*exp(2*theta_LL)
    varcomp <- cbind(vcov, SE, lower, upper)   
  } else {
    stop("The model must be a lme, lmer, glmer or glmmTMB object!")
  }
  return(varcomp)  
}