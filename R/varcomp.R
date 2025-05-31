#' Calculate SE and CI of variance components for \code{lmer}, \code{glmer},
#' \code{lme} model
#'
#' This function calculates SE and CI of variance components for \code{lmer},
#' \code{glmer}, \code{lme}, \code{glmmTMB} model.
#'
#'
#' @param model Model object returned by \code{lmer}, \code{glmer}, \code{lme},
#' \code{glmmTMB}.
#' @param ci a logical value to indicates wheather or not to simulate a
#' confidence interval for \code{lmer} model, the default value is TRUE.
#' @param level level of confidence of CI, the default value is 0.95.
#' @return Variance components table.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' \dontrun{varcomp(fm)}
#' fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' varcomp(fm1)
#'
#' data(Orthodont, package="nlme")
#' mod <- lmer(distance ~ age + (age|Subject), data=Orthodont)
#' \dontrun{varcomp(mod)}
#' mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
#' varcomp(mod1)
#' @importFrom lmeInfo extract_varcomp varcomp_vcov
#' @importFrom lme4 VarCorr
#' @importFrom nlme intervals
#' @importFrom stats confint
#' @importFrom utils combn
#'
#' @export
varcomp <- function(model, ci=TRUE, level=0.95) {
  # vcov_lmerMod(model, full = TRUE, ranpar = "var"), vcov_lmerMod(model, full = TRUE, ranpar = "sd")?
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"), inherits(model, "lmerModLmerTest"))) {
    varcomp <- as.data.frame(VarCorr(model), order = "lower.tri")
    vcov_index <- -1:-length(fixef(model))
    if (inherits(model, "glmerMod")) {
      random_vcov <- as.matrix(vcov_glmerMod(model, full = TRUE, ranpar = "var")[vcov_index, vcov_index])
      ci <- FALSE
    } else {
      random_vcov <- as.matrix(vcov_lmerMod(model, full = TRUE, ranpar = "var")[vcov_index, vcov_index])
    }
    if (inherits(model, "glmerMod")) {
      random_terms <- gsub("cov_", "", rownames(random_vcov))
    } else {
      random_terms <- gsub("cov_", "", colnames(random_vcov))
    }
    varcomp$SE <- sqrt(diag(random_vcov))
    if (ci){
      varcomp <- cbind(varcomp, confint(model, level=level)[1:length(random_terms), ])
      kk <- abs(varcomp$sdcor^2 - varcomp$vcov) > 0.0000000000001
      varcomp[, 7] <- ifelse(!kk, ifelse(varcomp[, 7] < 0, 0, varcomp[, 7]^2), varcomp[, 7]*varcomp$vcov/varcomp$sdcor)
      varcomp[, 8] <- ifelse(!kk, varcomp[, 8]^2, varcomp[, 8]*varcomp$vcov/varcomp$sdcor)
      varcomp <- round(varcomp[,c(4,6:8)], 4)
    } else {
      varcomp <- round(varcomp[,c(4,6)], 4)
    }
    rownames(varcomp) <- random_terms
  } else if (inherits(model, "lme")) {
    vcov=unlist(extract_varcomp(model))
    SE=sqrt(diag(varcomp_vcov(model)))

    # cbind(Std.error, vcov)
    confint_list <- try(intervals(model, which="var-cov", level=level), silent = TRUE)
    if (!inherits(confint_list, "try-error")) {
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
    } else {
      varcomp <- cbind(vcov, SE)
    }
  } else {
  # else if(inherits(model, "glmmTMB")){
  # theta_sd <- sqrt(diag(vcov(model, full=TRUE)))
  # theta_sd <- theta_sd[grepl("theta", names(theta_sd))]
  # theta <- getME(model, "theta")
  # vc <- VarCorr(model)$cond
  # theta_id <- unlist(lapply(vc, function(x) {
  # n <- dim(x)[1]
  # v_ind <- rep(2, n*(n+1)/2)
  # v_ind[1:n] <- 1
  # return(v_ind)
  # }
  # ))
  # # vcov <- exp(2 * theta)
  # # vcov[theta_id==2] <- theta[theta_id==2]

  # vcov <- unlist(lapply(vc, function(x) c(diag(x), x[lower.tri(x, diag = FALSE)])))
  # SE=theta_sd*2*vcov
  # theta_UL <- theta+qnorm(p=(1-level)/2, lower.tail=FALSE)*theta_sd
  # upper <- theta_sd*2*exp(2*theta_UL)
  # upperN <- theta_sd*2*theta_UL
  # upper[theta_id==2] <- upperN[theta_id==2]
  # theta_LL <- theta-qnorm(p=(1-level)/2, lower.tail=FALSE)*theta_sd
  # lower <- theta_sd*2*exp(2*theta_LL)
  # lowerN <- theta_sd*2*theta_LL
  # lower[theta_id==2] <- lowerN[theta_id==2]
  # varcomp <- round(cbind(vcov, SE, lower, upper), 4)
  # }

    stop("The model must be a lme, lmer, glmer or glmmTMB object!")
  }
  return(varcomp)
}
