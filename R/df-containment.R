## Containment denominator df for lmer models (SAS PROC MIXED DDFM=CONTAIN).
##
## SAS rule: a random effect contains a fixed effect when it holds all the
## fixed effect's variables (and possibly more). The df of a fixed effect is
## the smallest rank contribution, to [X Z], of the random effects that
## contain it; when none does, the residual df N - rank([X Z]) is used.
## Rank contributions are sequential, as in SAS: the rank that a random
## effect adds to X and the random effects before it. Effects are ordered with
## fewer variables first (Block before Block*Variety, as one writes them in a
## RANDOM statement), then by their order in the formula (lme4 reorders them).
##
## Random effects follow SAS: in (1 + x | g) the intercept part is the effect
## g and the slope part is x*g, so (1 | Block/Variety) gives Block and
## Block*Variety.

## Random effects of an lmer model: list of (vars, Z column index).
pm_random_effects <- function(model) {
  bars <- lme4::findbars(stats::formula(model))
  cnms <- lme4::getME(model, "cnms")
  Gp <- lme4::getME(model, "Gp")
  out <- list()
  # lme4 reorders the terms (most levels first): match each to its term in
  # the formula by the grouping variables, to get the slope variables and
  # the formula position
  bar_g <- lapply(bars, function(b) all.vars(b[[3]]))
  used <- logical(length(bars))
  for (i in seq_along(cnms)) {
    gvars <- unlist(strsplit(names(cnms)[i], ":", fixed = TRUE))
    hit <- which(!used & vapply(bar_g, function(g) setequal(g, gvars), logical(1)))[1]
    if (!is.na(hit)) used[hit] <- TRUE
    lhs_vars <- if (!is.na(hit)) all.vars(bars[[hit]][[2]]) else character(0)
    pos <- if (is.na(hit)) i else hit
    nc <- length(cnms[[i]])
    cols <- (Gp[i] + 1):Gp[i + 1]
    comp <- (seq_along(cols) - 1) %% nc + 1
    is_int <- cnms[[i]][comp] == "(Intercept)"
    if (any(is_int)) {
      out[[length(out) + 1]] <- list(name = names(cnms)[i], vars = gvars, cols = cols[is_int], pos = pos)
    }
    if (any(!is_int)) {
      out[[length(out) + 1]] <- list(name = paste(c(lhs_vars, names(cnms)[i]), collapse = ":"),
                                     vars = union(lhs_vars, gvars), cols = cols[!is_int], pos = pos)
    }
  }
  out
}

## Containment df of every fixed term (named, "(Intercept)" included), with
## attribute "residual" (N - rank([X Z])).
pm_containment_df <- function(model) {
  stopifnot(inherits(model, "merMod"))
  X <- Matrix::Matrix(lme4::getME(model, "X"), sparse = TRUE)
  Z <- methods::as(lme4::getME(model, "Z"), "CsparseMatrix")
  # rank(M) = rank(M'M): a (p+q) x (p+q) problem instead of n x (p+q)
  rk <- function(M) qr(as.matrix(Matrix::crossprod(M)), tol = 1e-10)$rank
  r_full <- rk(cbind(X, Z))
  res_df <- nrow(X) - r_full
  re <- pm_random_effects(model)
  re <- re[order(vapply(re, function(e) length(e$vars), integer(1)),
                 vapply(re, function(e) e$pos, numeric(1)))]
  contrib <- numeric(length(re))
  r_prev <- rk(X)
  M <- X
  for (k in seq_along(re)) {
    M <- cbind(M, Z[, re[[k]]$cols, drop = FALSE])
    r_now <- rk(M)
    contrib[k] <- r_now - r_prev
    r_prev <- r_now
  }
  re_vars <- lapply(re, `[[`, "vars")

  tt <- stats::terms(model)
  labs <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  term_vars <- lapply(labs, function(l) rownames(fac)[fac[, l] > 0])
  # variable names as they appear in the random effects (e.g. factor(x) -> x)
  term_vars <- lapply(term_vars, function(v) unique(unlist(lapply(v, function(s) {
    av <- all.vars(str2lang(s)); if (length(av)) av else s
  }))))
  if (attr(tt, "intercept") == 1) {
    labs <- c("(Intercept)", labs)
    term_vars <- c(list(character(0)), term_vars)
  }
  df <- vapply(term_vars, function(v) {
    ok <- vapply(re_vars, function(rv) all(v %in% rv), logical(1))
    if (any(ok)) min(contrib[ok]) else res_df
  }, numeric(1))
  names(df) <- labs
  attr(df, "residual") <- res_df
  df
}
