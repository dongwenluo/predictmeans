varcomp <- function(model, ci=TRUE, level=0.95) {
  # vcov_lmerMod(model, full = TRUE, ranpar = "var"), vcov_lmerMod(model, full = TRUE, ranpar = "sd")?
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"), inherits(model, "lmerModLmerTest"))) {
    varcomp <- as.data.frame(VarCorr(model), order = "lower.tri")
    vcov_index <- -1:-length(fixef(model))
	if (inherits(model, "glmerMod")) {
	  random_vcov <- as.matrix(vcov_glmerMod(model, full = TRUE, ranpar = "var")[vcov_index, vcov_index, drop = FALSE])
	  ci <- FALSE
	}else random_vcov <- as.matrix(vcov_lmerMod(model, full = TRUE, ranpar = "var")[vcov_index, vcov_index, drop = FALSE])
    random_terms <- gsub("cov_", "", colnames(random_vcov))   # vcov_glmerMod() has no row names
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
    varcomp <- pm_tmb_varcomp(model, ci = ci, level = level)
  } else {
    stop("The model must be a lme, lmer, glmer or glmmTMB object!")
  }
  return(varcomp)  
}

## Variance components of a glmmTMB fit (conditional model, then the
## zero-inflation model with prefix "zi."), laid out as for lmer: variances
## and covariances of each term (lower triangle, by column), then the
## residual variance of a gaussian model. SEs by the delta method from
## vcov(full = TRUE): the components are recomputed by VarCorr() from
## perturbed parameters, so this works for every covariance structure.
## CIs are Wald intervals on the log scale for variances (glmmTMB's own
## scale, so always positive) and on the atanh scale for correlations
## (times the two SDs for a covariance).
pm_tmb_varcomp <- function(model, ci = TRUE, level = 0.95) {
  pf <- model$fit$parfull
  keep <- !names(pf) %in% c("b", "bzi")   # vcov(full = TRUE) has these, in this order
  pn <- names(pf)[keep]
  V <- vcov(model, full = TRUE)
  if (nrow(V) != length(pn)) stop("cannot match the parameters of this glmmTMB model to vcov(full = TRUE)")
  gaus <- stats::family(model)$family == "gaussian"
  use <- pn %in% c("theta", "thetazi", if (gaus) c("betad", "betadisp"))
  pos <- which(keep)[use]
  on.exit(model$obj$env$report(model$fit$parfull))   # leave the TMB object as it was

  components <- function(x) {
    m <- model
    m$fit$parfull[pos] <- x
    for (nm in unique(pn[use])) m$fit$par[names(m$fit$par) == nm] <- x[pn[use] == nm]
    vc <- VarCorr(m)
    est <- r <- sdp <- numeric(0)
    for (part in c("cond", "zi")) {
      for (j in seq_along(vc[[part]])) {
        v <- vc[[part]][[j]]
        is_diag <- identical(names(attr(v, "blockCode")), "diag")   # no covariances
        grp <- paste0(if (part == "zi") "zi.", names(vc[[part]])[j])
        cn <- rownames(v)
        v <- matrix(as.numeric(v), nrow(v))
        for (b in seq_len(ncol(v))) for (a in b:nrow(v)) {
          if (a != b && is_diag) next
          nm <- paste(c(grp, cn[unique(c(a, b))]), collapse = ".")
          est[nm] <- v[a, b]
          if (a != b) {
            sdp[nm] <- sqrt(v[a, a]*v[b, b])
            r[paste0("cor:", nm)] <- v[a, b]/sdp[nm]
          }
        }
      }
    }
    if (gaus) est["residual"] <- sigma(m)^2
    structure(c(est, r), sdprod = sdp)
  }

  x0 <- pf[pos]
  val <- components(x0)
  J <- numDeriv::jacobian(function(x) as.numeric(components(x)), x0)
  se <- sqrt(pmax(diag(J %*% V[use, use, drop = FALSE] %*% t(J)), 0))
  names(se) <- names(val)
  is_r <- startsWith(names(val), "cor:")
  est <- val[!is_r]
  out <- cbind(vcov = est, SE = se[!is_r])
  if (ci) {
    z <- stats::qnorm(1 - (1 - level)/2)
    lower <- ifelse(est > 0, est*exp(-z*se[!is_r]/est), 0)
    upper <- ifelse(est > 0, est*exp(z*se[!is_r]/est), NA)
    r <- val[is_r]
    nm <- sub("^cor:", "", names(r))
    sdprod <- attr(val, "sdprod")[nm]
    lower[nm] <- tanh(atanh(r) - z*se[is_r]/(1 - r^2))*sdprod
    upper[nm] <- tanh(atanh(r) + z*se[is_r]/(1 - r^2))*sdprod
    out <- cbind(out, lower, upper)
    colnames(out)[3:4] <- paste(format(100*c((1 - level)/2, 1 - (1 - level)/2), trim = TRUE,
                                       scientific = FALSE, digits = 3), "%")
  }
  as.data.frame(round(out, 4), check.names = FALSE)   # as for lmer
}
