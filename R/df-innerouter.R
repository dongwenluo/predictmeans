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

## The fixed-effect design and the grouping factors of a model; 'drop':
## grouping factors left out (the smooth blocks of a semireg() fit, which
## are not strata of the design).
pm_io_parts <- function(model, drop = character(0)) {
  if (inherits(model, "merMod")) {
    fl <- as.list(lme4::getME(model, "flist"))
    return(list(X = as.matrix(lme4::getME(model, "X")), flist = fl[!names(fl) %in% drop]))
  }
  if (inherits(model, "glmmTMB")) {
    X <- as.matrix(glmmTMB::getME(model, "X"))
    fl <- as.list(model$modelInfo$reTrms$cond$flist)
    return(list(X = X, flist = fl[!names(fl) %in% drop]))
  }
  if (inherits(model, "clmm")) {                    # ordinal::clmm: the fixed design with its intercept (R/ordinal.R)
    tt <- stats::delete.response(stats::terms(pm_nobars(stats::formula(model))))
    X <- stats::model.matrix(tt, model$model, contrasts.arg = model$contrasts)
    fl <- as.list(model$gfList)
    # the J - 1 cut points replace the one intercept: J - 2 more parameters, estimated within
    # (the residual stratum; = SAS GLIMMIX's df for terms within clusters, dev/sas/ordinal.sas)
    return(list(X = X, flist = fl[!names(fl) %in% drop], extra = length(model$alpha) - 1))
  }
  stop("df_method = \"inner-outer\" needs a model with random effects (lmer, glmer, glmmTMB, clmm, lme).", call. = FALSE)
}

## The inner-outer df of the rows of L (columns: the coefficients named in
## 'coef_names'): contrasts take the smallest df of the coefficients they
## use; means (rows that use the intercept) the df of their term ('terms').
## Used by pm_inference() and for semireg() fits.
pm_io_row_df <- function(io, coef_names, L, terms) {
  L <- rbind(L)
  d <- io[coef_names]
  int <- attr(io, "intercept")[match(coef_names, names(io))]
  tdf <- attr(io, "terms")
  tdf <- tdf[intersect(terms, names(tdf))]
  apply(L, 1, function(l) {
    nz <- abs(l) > 1e-10
    if (!any(nz)) return(NA_real_)
    if (any(nz & int) && length(tdf) && all(is.finite(tdf))) return(min(tdf))
    min(d[nz])
  })
}

## The cases of the residual stratum: when each row holds several cases
## (binomial trials, case weights of clmm; pm_freq_w()) and the counts are
## whole numbers, their sum, so aggregated counts and one row per case get
## the same df; rows with no cases are left out (as a fit to one row per case
## has none). Weights that are not whole numbers are not counts (sampling or
## precision weights): rows. Counting rows made within-cluster tests far too
## conservative for binomial data (12 rows of 25 trials: 1.1% of true null
## hypotheses rejected at 5%, trials 4.3%; dev/sim-df-io-trials.R), but it
## hid overdispersion the fit leaves out (row SD 0.3: rows 3.5%, trials
## 9.25%). options(predictmeans.inner_outer_count = "rows") counts rows.
## Returns list(n, keep (rows used), count: "rows", "trials" or
## "case weights").
pm_io_cases <- function(model, N) {
  how <- getOption("predictmeans.inner_outer_count", "cases")
  if (!identical(how, "cases") && !identical(how, "rows")) {
    stop("options(predictmeans.inner_outer_count = ) must be \"cases\" or \"rows\".", call. = FALSE)
  }
  rows <- list(n = N, keep = rep(TRUE, N), count = "rows")
  if (how == "rows") return(rows)
  w <- tryCatch(pm_freq_w(model), error = function(e) e)
  if (is.null(w)) return(rows)
  if (inherits(w, "error") || !is.numeric(w) || length(w) != N || !all(is.finite(w)) || any(w < 0)) {
    warning("predictmeans could not match the case counts of the fit to its rows; the inner-outer df count rows.",
            call. = FALSE)
    return(rows)
  }
  if (!pm_counts(w)) return(rows)
  list(n = sum(w), keep = w > 0, count = if (inherits(model, "clmm")) "case weights" else "trials")
}

## TRUE when weights are counts of cases: whole numbers, not all 1 (one rule
## for the inner-outer df and the dispersion test, pm_check_dispersion()).
pm_counts <- function(w) is.numeric(w) && !all(w == 1) && all(abs(w - round(w)) <= 1e-8)

## The cases counted by the inner-outer df of an inference (pm_inference()),
## for settings$df_cases: NA (other df, or df given by the user: 'user'),
## "rows", "trials" or "case weights", with attribute "df" the residual
## stratum's df when it counts cases and a coefficient uses it (the
## overdispersion advice compares it with the df of the result).
pm_df_cases <- function(inf, user = FALSE) {
  ct <- inf$df_cases
  if (is.null(ct) || isTRUE(user)) return(NA_character_)
  structure(ct, df = inf$df_trials)
}

## df of each coefficient (named vector).
## Value: df of each coefficient (named), with attributes "intercept" (which
## coefficient is the intercept) and "terms" (df of each term).
pm_inner_outer_df <- function(model, drop = character(0)) {
  if (inherits(model, "lme")) {                                # nlme's own
    d <- model$fixDF$X
    return(structure(d, intercept = names(d) == "(Intercept)", terms = model$fixDF$terms, count = "rows"))
  }
  parts <- pm_io_parts(model, drop)
  X <- parts$X
  cs <- pm_io_cases(model, nrow(X))
  if (!all(cs$keep)) {                                       # rows with no cases (binomial n = 0)
    a <- attr(X, "assign")
    X <- X[cs$keep, , drop = FALSE]
    attr(X, "assign") <- a
  }
  fl <- lapply(parts$flist, function(f) factor(f[cs$keep]))
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
  res_df <- cs$n - rk(cbind(X, Zall)) - if (is.null(parts$extra)) 0 else parts$extra
  sdf <- vapply(seq_along(fl), stratum_df, numeric(1))
  # stratum of each column: the coarsest factor within which it is constant
  # (its number of groups), or the residual stratum (Inf)
  strat <- vapply(seq_len(ncol(X)), function(j) {
    k <- which(C[j, ])
    if (!length(k)) Inf else k[which.min(nlev[k])]
  }, numeric(1))
  df <- ifelse(is.infinite(strat), res_df, sdf[pmin(strat, length(sdf))])
  int <- apply(X, 2, function(x) all(abs(x - x[1]) <= 1e-12))
  # cases counted only where a coefficient is in the residual stratum (with an observation-level random
  # effect none is: the effects within clusters get the rows' stratum)
  if (cs$count != "rows" && !any(is.infinite(strat) & !int)) {
    cs$count <- "rows"
    res_df <- N - rk(cbind(X, Zall)) - if (is.null(parts$extra)) 0 else parts$extra
    df[is.infinite(strat)] <- res_df
  }
  fine <- ifelse(is.infinite(strat), Inf, nlev[pmin(strat, length(nlev))])
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
  # a df below 1 would give NaN limits: not available, so pm_df_method() says so and uses its fallback
  if (any(df < 1, na.rm = TRUE) || any(tdf < 1, na.rm = TRUE)) {
    stop("some are below 1 (few rows or cases per stratum)", call. = FALSE)
  }
  structure(df, intercept = int, terms = tdf, count = cs$count, trials_df = if (cs$count != "rows") res_df)
}
