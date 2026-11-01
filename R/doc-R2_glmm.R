# Help page R2_glmm (roxygen2 source; man/R2_glmm.Rd is generated).
#' An adjusted coefficient of determination (R2) for generalized linear mixed models
#'
#' @description
#' This function produces adjusted R2 for generalized linear mixed models which was crafted following the guidance provided by Professor Hans-Peter Piepho.
#' @usage
#' R2_glmm(model, over_disp=FALSE)
#' @param model An object returned by \code{lmer}, \code{glmer} or \code{glmmTMB}.
#' @param over_disp A logical scalar to indicate whether \code{model} with over-dispersion or not. The default value is FALSE.
#' @return
#' Adjusted R2 in percentage for Total (fixed + random), Fiexd, Random and individual random term.
#' @references
#' Piepho HP. An adjusted coefficient of determination (R2 ) for generalized linear mixed models in one go. Biom J. 2023 Oct;65(7):e2200290. doi: 10.1002/bimj.202200290. Epub 2023 May 1. PMID: 37127864.
#' @examples
#'   library(predictmeans)
#'   Oats$nitro <- factor(Oats$nitro)
#'   (fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats))
#'   R2_glmm(fm)
#'   (gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
#'               data = cbpp, family = binomial))
#'   R2_glmm(gm)
#' @name R2_glmm
NULL
