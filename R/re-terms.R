## Random-effect terms of a fitted lmer/glmer or glmmTMB model (conditional
## model), in the column order of Z. For each term: the grouping factor, the
## names of its components (e.g. "(Intercept)", "Days"), its number of levels,
## its covariance matrix (from VarCorr()) and its columns in Z. Within a term
## the columns run by level, then component, so component k of a term with nc
## components is every nc-th column.
pm_re_terms <- function(model) {
  if (inherits(model, "glmmTMB")) {
    vc <- VarCorr(model)$cond
    rs <- model$modelInfo$reStruc$condReStruc
    nl <- vapply(rs, function(r) as.numeric(r$blockReps), numeric(1))
    groups <- trimws(sub(".*\\|", "", names(rs)))
  } else {
    vc <- VarCorr(model)
    cnms <- getME(model, "cnms")
    nl <- diff(getME(model, "Gp")) / lengths(cnms)
    groups <- names(cnms)
  }
  start <- 0
  out <- vector("list", length(vc))
  for (j in seq_along(vc)) {
    v <- vc[[j]]
    v <- matrix(as.numeric(v), nrow(v), dimnames = list(rownames(v), colnames(v)))
    q <- nl[j] * nrow(v)
    out[[j]] <- list(group = groups[j], comps = rownames(v), nl = nl[j], vcov = v,
                     cols = start + seq_len(q))
    start <- start + q
  }
  out
}

## Lower-triangular L with L L' = v, for a covariance matrix v that may be
## singular (a zero variance, or a correlation of +-1): a column with no
## variance left is set to zero instead of failing as chol() does.
pm_chol_lower <- function(v) {
  k <- nrow(v)
  L <- matrix(0, k, k)
  tol <- 1e-12 * max(diag(v), 0)
  for (j in seq_len(k)) {
    p <- seq_len(j - 1)
    d <- v[j, j] - sum(L[j, p]^2)
    if (d > tol) {
      L[j, j] <- sqrt(d)
      for (i in seq_len(k - j) + j) {
        L[i, j] <- (v[i, j] - sum(L[i, p] * L[j, p])) / L[j, j]
      }
    }
  }
  L
}

## Relative covariance factor Lambda of the random effects of a glmmTMB fit
## (G = sigma^2 Lambda Lambda', as getME(<merMod>, "Lambda")), built from
## VarCorr(). glmmTMB's own theta can't be used: its correlation parameters
## are on another scale, and glmmTMB >= 1.1.12 no longer stores Lambdat.
pm_tmb_lambda <- function(model, terms = pm_re_terms(model)) {
  s <- sigma(model)
  blocks <- lapply(terms, function(tm) {
    L <- Matrix::Matrix(pm_chol_lower(tm$vcov) / s, sparse = TRUE)
    kronecker(Matrix::Diagonal(tm$nl), L)
  })
  Matrix::bdiag(blocks)
}
