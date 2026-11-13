## Denominator df by grouping level ("inner-outer", df_method = "inner-outer";
## the default for GLMMs). nlme's rule (Pinheiro and Bates 2000, section
## 4.2.2; nlme's getFixDF()): each fixed-effect coefficient is estimated in
## the stratum of the random-effect grouping within whose groups its column
## is constant, and gets that stratum's df; a contrast gets the smallest df
## of the coefficients it uses. It depends on the design only, so it works
## for GLMMs, where Kenward-Roger and Satterthwaite do not. In simulations
## (dev/sim-df/README.md) it held 4-8% for GLMMs where the residual df gave
## 8-15% for effects that vary between clusters.
##
## Computed from ranks, which gives nlme's values for nested grouping factors
## and extends the rule to crossed ones:
##   a column constant within the groups of factor g (g the coarsest such
##   factor: a site-level treatment belongs to the site stratum even though
##   it is also constant within blocks) gets  n_g - rank([columns constant
##   within g, groups of the factor g is nested in]);
##   a column that varies within the groups of every factor gets
##   N - rank([X, indicators of all the grouping factors]);
##   the intercept (constant everywhere) gets the largest stratum df, as nlme.
## A term gets the df of the stratum of its innermost column (nlme's
## fixDF$terms); means (rows that use the intercept) take the df of their
## term, so all the means of a term share one df whatever the coding, while
## contrasts take the smallest df of the coefficients they use.

## The fixed-effect design and the grouping factors of a model.
pm_io_parts <- function(model) {
  if (inherits(model, "merMod")) {
    return(list(X = as.matrix(lme4::getME(model, "X")), flist = as.list(lme4::getME(model, "flist"))))
  }
  if (inherits(model, "glmmTMB")) {
    X <- as.matrix(glmmTMB::getME(model, "X"))
    fl <- model$modelInfo$reTrms$cond$flist
    return(list(X = X, flist = as.list(fl)))
  }
  stop("df_method = \"inner-outer\" needs a model with random effects (lmer, glmer, glmmTMB, lme).", call. = FALSE)
}

## df of each coefficient (named vector).
## Value: df of each coefficient (named), with attributes "intercept" (which
## coefficient is the intercept) and "terms" (df of each term).
pm_inner_outer_df <- function(model) {
  if (inherits(model, "lme")) {                                # nlme's own
    d <- model$fixDF$X
    return(structure(d, intercept = names(d) == "(Intercept)", terms = model$fixDF$terms))
  }
  parts <- pm_io_parts(model)
  X <- parts$X
  fl <- lapply(parts$flist, function(f) factor(f))
  if (!length(fl)) stop("the model has no random effects.", call. = FALSE)
  N <- nrow(X)
  rk <- function(M) qr(M, tol = 1e-7)$rank
  const <- function(x, g) all(abs(x - stats::ave(x, g)) <= 1e-8 * max(1, abs(x)))
  nested_in <- function(g, f) all(tapply(as.integer(f), g, function(v) length(unique(v)) == 1L))
  nlev <- vapply(fl, nlevels, integer(1))
  # columns constant within each factor's groups
  C <- vapply(fl, function(g) vapply(seq_len(ncol(X)), function(j) const(X[, j], g), logical(1)),
              logical(ncol(X)))
  C <- matrix(C, ncol(X), length(fl), dimnames = list(colnames(X), names(fl)))
  stratum_df <- function(k) {
    g <- fl[[k]]
    # the finest factor that g is nested in (its groups are unions of g's)
    up <- which(vapply(seq_along(fl), function(i) i != k && nlev[i] < nlev[k] && nested_in(g, fl[[i]]), logical(1)))
    M <- X[, C[, k], drop = FALSE]
    if (length(up)) {
      P <- fl[[up[which.max(nlev[up])]]]
      M <- cbind(M, stats::model.matrix(~ P - 1))
    }
    nlev[k] - rk(M)
  }
  Zall <- do.call(cbind, lapply(fl, function(g) stats::model.matrix(~ g - 1)))
  res_df <- N - rk(cbind(X, Zall))
  sdf <- vapply(seq_along(fl), stratum_df, numeric(1))
  # stratum of each column: the coarsest factor within which it is constant
  # (its number of groups), or the residual stratum (Inf)
  strat <- vapply(seq_len(ncol(X)), function(j) {
    k <- which(C[j, ])
    if (!length(k)) Inf else k[which.min(nlev[k])]
  }, numeric(1))
  df <- ifelse(is.infinite(strat), res_df, sdf[pmin(strat, length(sdf))])
  fine <- ifelse(is.infinite(strat), Inf, nlev[pmin(strat, length(nlev))])
  int <- apply(X, 2, function(x) all(abs(x - x[1]) <= 1e-12))
  df[int] <- max(c(sdf, res_df))
  names(df) <- colnames(X)
  # term df: the stratum of the term's innermost column
  asg <- attr(X, "assign")
  tl <- tryCatch(attr(stats::terms(model, fixed.only = TRUE), "term.labels"), error = function(e) NULL)
  if (is.null(tl)) tl <- tryCatch(attr(stats::terms(model), "term.labels"), error = function(e) NULL)
  tdf <- if (!is.null(asg) && length(asg) == ncol(X) && length(tl) >= max(asg)) {
    vapply(seq_along(tl), function(a) {
      j <- which(asg == a)
      if (!length(j)) NA_real_ else df[j[which.max(fine[j])]]
    }, numeric(1))
  }
  if (!is.null(tdf)) names(tdf) <- tl
  structure(df, intercept = int, terms = tdf)
}
