# Help page varcomp (roxygen2 source; man/varcomp.Rd is generated).
#' Calculate SE and CI of variance components for \code{lmer}, \code{glmer}, \code{lme}, \code{glmmTMB} model
#'
#' @description
#' This function calculates SE and CI of variance components for \code{lmer}, \code{glmer}, \code{lme}, \code{glmmTMB} model. For \code{lme} models it needs package \pkg{lmeInfo}.
#' @details
#' For \code{glmmTMB} models the table holds the variances and covariances of each random-effect term
#' (conditional model, then the zero-inflation model with prefix \code{zi.}) and, for gaussian models, the
#' residual variance. SEs are from the delta method. CIs are Wald intervals on the log scale for variances
#' and on the Fisher z (atanh) scale for correlations, times the two standard deviations for covariances.
#' @usage
#' varcomp(model, ci=TRUE, level=0.95)
#' @param model Model object returned by \code{lmer}, \code{glmer}, \code{lme}, \code{glmmTMB}.
#' @param ci a logical value to indicates wheather or not to simulate a confidence interval for \code{lmer} and \code{glmmTMB} models, the default value is TRUE.
#' @param level level of confidence of CI, the default value is 0.95.
#' @return
#' Variance components table.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
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
#' @name varcomp
NULL
