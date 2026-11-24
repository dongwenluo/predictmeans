## permanova.lmer() was removed in 2.0: permmodels(model, type = ) tests the
## terms of lmer and glmer fits (Freedman-Lane, the default; permanova.lmer()
## permuted through permlmer() with a refit with and without each term, so
## its p-values differ).
#' @export
permanova.lmer <- function(...) pm_defunct_fun("permanova.lmer")
