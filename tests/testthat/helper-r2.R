## R2_model()'s values as R2_glmm() (removed in 2.0) printed them: percentages
## rounded to 2 decimals, so that the 1.1.1 reference values compare exactly.
r2_pct <- function(r) stats::setNames(paste0(round(100 * as.numeric(unclass(r)), 2), "%"), names(r))
