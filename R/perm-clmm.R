## permmodels() for ordinal::clmm (B4 G3, dev/review-b4.md): randomisation
## tests of the terms applied to whole clusters, one cluster-level variable
## at a time. For a variable v (a main-effect term of the model, constant
## within the clusters of a random-effect factor g; the finest such g),
## v's values are relabelled among the clusters within strata: the
## clusters that share the other variables constant within g (fixed blocks,
## other cluster-level factors) and the coarser grouping factors g is nested
## in (farms of plots). Under the null hypothesis that v has no effect at all
## (main effect and interactions) every such allocation is equally likely,
## so the test is exact whatever the within-cluster effects and cluster
## sizes. The statistic is v's Wald chi-square (Type II or III).
##   - Allocations: the product over strata of multinomial counts. At most
##     nperm: all are used (exact p-values). For a factor whose levels have
##     equal counts in every stratum, relabelling its levels gives the same
##     statistic: one allocation per class is refitted (weight k!), which
##     removes ties broken by optimiser noise and divides the work by k!.
##     Otherwise nperm random allocations (seed, else the user's stream);
##     near-ties (relative 1e-6) count as at least as large.
##   - A refit that does not converge (pm_clmm_ok()) counts as at least as
##     extreme as the observed allocation (conservative), with a note.
## Returns the permmodels() result (tables, notes, advice).
pm_perm_clmm <- function(model, nperm, type, data, seed, ncore) {
  w <- pm_freq_w(model)
  if (!is.null(w) && any(w != 1)) {
    stop("permmodels() needs a clmm fit with one row per observation (no case weights).", call. = FALSE)
  }
  if (!pm_clmm_ok(model)) stop("The clmm fit has not converged (gradient not near 0); refit it first.", call. = FALSE)
  d <- if (is.null(data)) pm_clmm_rows(model) else as.data.frame(data)
  gf <- lapply(model$gfList, factor)
  if (nrow(d) != length(gf[[1]])) {
    stop("'data' must hold the rows the clmm was fitted to (", length(gf[[1]]), " rows; without those dropped for ",
         "missing values).", call. = FALSE)
  }
  ff <- pm_nobars(stats::formula(model))
  tl <- attr(stats::terms(ff), "term.labels")
  fvars <- all.vars(stats::delete.response(stats::terms(ff)))
  const <- function(x, g) all(tapply(as.character(x), g, function(v) length(unique(v)) == 1L))
  nested_in <- function(g, f) nlevels(f) < nlevels(g) && const(f, g)
  aT <- pm_anova(model, type = type, test_statistic = "Chisq")
  tests <- list()
  for (v in intersect(tl, fvars)) {                      # main-effect terms of one variable
    ks <- which(vapply(gf, function(g) const(d[[v]], g), logical(1)))
    if (!length(ks)) next
    k <- ks[which.max(vapply(gf[ks], nlevels, integer(1)))]       # the finest grouping that holds v
    g <- gf[[k]]
    first <- !duplicated(g)
    others <- setdiff(fvars[vapply(fvars, function(u) const(d[[u]], g), logical(1))], v)
    coarse <- gf[vapply(gf, function(f) nested_in(g, f), logical(1))]
    strata_df <- c(lapply(d[others], function(x) x[first]), lapply(coarse, function(f) f[first]))
    strat <- if (length(strata_df)) do.call(paste, c(lapply(strata_df, as.character), sep = "\r")) else
      rep("", sum(first))
    tests[[v]] <- list(v = v, g = names(gf)[k], unit = g, val = d[[v]][first], strat = strat,
                       others = c(others, names(coarse)))
  }
  if (!length(tests)) {
    stop("No main-effect term of this clmm is constant within the clusters of a random effect: there is no ",
         "treatment of whole clusters to permute (use bootmodels()).", call. = FALSE)
  }
  res <- lapply(tests, function(tt) pm_perm_clmm_one(model, d, tt, aT, nperm, type, seed, ncore))
  p <- vapply(res, `[[`, 0, "p")
  ANOVA <- cbind(round(aT, 4), Perm_p_value = NA_real_)
  ANOVA[names(p), "Perm_p_value"] <- round(p, 4)
  ctab <- pm_coef_table(model)
  COEF <- round(cbind(ctab, Perm_p_value = NA_real_), 4)
  for (v in names(p)) {                                  # one coefficient: its |z| test is the term's
    j <- which(attr(pm_basis(model, d[1, , drop = FALSE]), "assign") == match(v, tl))
    if (length(j) == 1) COEF[colnames(pm_basis(model, d[1, , drop = FALSE]))[j], "Perm_p_value"] <- round(p[[v]], 4)
  }
  attr(ANOVA, "heading") <- attr(aT, "heading")
  notes <- vapply(res, `[[`, "", "note")
  out <- pm_permmodels_result(list(permlist = NULL, COEFFICENTS = COEF, ANOVA = ANOVA),
                              "Note: Perm_p_value of a single-coefficient term is that of its term (below).",
                              paste(c("Note: Perm_p_value from randomisation of whole clusters:", notes),
                                    collapse = " "))
  attr(out, "resampling") <- "cluster randomisation"
  out
}

## One cluster-level variable: its allocations, refits and p-value.
pm_perm_clmm_one <- function(model, d, tt, aT, nperm, type, seed, ncore) {
  units <- seq_along(tt$val)
  by_s <- split(units, tt$strat)
  lev <- as.character(tt$val)
  count <- exp(sum(vapply(by_s, function(i) lfactorial(length(i)) - sum(lfactorial(table(lev[i]))), 0)))
  if (round(count) <= 1) {
    return(list(p = NA_real_, note = paste0(tt$v, ": no other allocation within ",
                                           paste(tt$others, collapse = ", "), " (not tested).")))
  }
  # a factor with equal counts of its levels in every stratum: the statistic does not change when its levels
  # are relabelled, so one allocation per relabelling class (levels in order of first appearance)
  isf <- is.factor(tt$val) || is.character(tt$val)
  lv <- unique(lev)
  sym <- isf && all(vapply(by_s, function(i) {
    tb <- table(factor(lev[i], levels = lv))
    all(tb == tb[1])
  }, logical(1)))
  wt <- if (sym) factorial(length(lv)) else 1
  enumerate <- count <= nperm
  alloc <- if (enumerate) {
    # all allocations: the product over strata of the distinct orderings of each stratum's values
    per <- lapply(by_s, function(i) pm_multiset_perms(match(lev[i], lv)))
    ix <- as.matrix(expand.grid(lapply(per, function(m) seq_len(ncol(m)))))
    A <- apply(ix, 1, function(r) {
      out <- integer(length(units))
      for (s in seq_along(by_s)) out[by_s[[s]]] <- per[[s]][, r[s]]
      out
    })
    A <- matrix(A, nrow = length(units))
    if (sym) A <- A[, apply(A, 2, function(a) identical(unique(a), seq_along(unique(a)))), drop = FALSE]
    A
  } else {
    draw <- function() replicate(nperm, {
      out <- match(lev, lv)
      for (i in by_s) out[i] <- out[i][sample.int(length(i))]
      out
    })
    if (is.null(seed)) draw() else pm_with_seed(seed, draw())
  }
  alloc <- matrix(alloc, nrow = length(units))
  pos <- match(tt$unit, unique(tt$unit))
  vals <- tt$val
  stat <- function(a) {
    dd <- d
    nv <- vals[match(lv[a], lev)]                         # a value with that label (keeps factor levels)
    dd[[tt$v]] <- nv[pos]
    f <- suppressWarnings(tryCatch(pm_refit(model, data = dd), error = function(e) NULL))
    if (!pm_clmm_ok(f)) return(NA_real_)
    tryCatch(suppressMessages(pm_anova(f, type = type, test_statistic = "Chisq"))[tt$v, "Chisq"],
             error = function(e) NA_real_)
  }
  sims <- unlist(pm_perm_map(lapply(seq_len(ncol(alloc)), function(j) alloc[, j]), stat, ncore,
                             packages = "ordinal"))
  obs <- aT[tt$v, "Chisq"]
  ge <- is.na(sims) | sims >= obs * (1 - 1e-6) - 1e-10       # failed refits: as extreme (conservative)
  nfail <- sum(is.na(sims))
  total <- if (enumerate) round(count) else ncol(alloc) + 1
  p <- if (enumerate) wt * sum(ge) / total else (sum(ge) + 1) / total
  floor_p <- (if (enumerate) wt else 1) / total
  note <- paste0(tt$v, ": relabelled among the levels of ", tt$g,
                 if (length(tt$others)) paste0(" within ", paste(tt$others, collapse = ", ")), "; ",
                 format(round(count), big.mark = ","), " allocations, ",
                 if (enumerate) paste0("all used (exact p-value", if (sym) paste0(", ", wt, " relabellings of ",
                                                                                  tt$v, "'s levels each"), ")")
                 else paste0(ncol(alloc), " random ones"),
                 "; smallest attainable p-value ", format(signif(floor_p, 3), scientific = FALSE),
                 if (nfail) paste0("; ", nfail, " refits did not converge and count as at least as extreme"),
                 ". The null hypothesis is no effect of ", tt$v, " at all (with its interactions).")
  list(p = min(1, p), note = note)
}
