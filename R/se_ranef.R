# This function calculate SE for random effects of lmer, glmer, glmmTMB models
# modified from function se.ranef at package 'arm'

se_ranef <- function (object, rand_term=NULL) {
  stopifnot(inherits(object, c("glmmTMB", "lmerModLmerTest", "merMod")))
  if (inherits(object, "glmmTMB")) se.bygroup <- ranef(object, condVar = TRUE)[["cond"]] else se.bygroup <- ranef(object, condVar = TRUE)
  rand_names <- names(se.bygroup)
  
  if (!is.null(rand_term)) {
    stopifnot(rand_term %in% names(se.bygroup))
    rand_names <- rand_term
  }
  
  for (m in rand_names) {
    if (inherits(object, "glmmTMB")) vars.m <- attr(se.bygroup[[m]], "condVar") else vars.m <- attr(se.bygroup[[m]], "postVar")
    K <- dim(vars.m)[1]
    J <- dim(vars.m)[3]
    names.full <- dimnames(se.bygroup[[m]])
    se.bygroup[[m]] <- array(NA, c(J, K))
    for (j in 1:J) {
      se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[ , , j])))
    }
    dimnames(se.bygroup[[m]]) <- list(names.full[[1]], names.full[[2]])
  }
  if (!is.null(rand_term)) se.bygroup <- se.bygroup[[m]]
  return(se.bygroup)
}
