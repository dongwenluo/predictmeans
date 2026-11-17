## Means over empty cells, and comparisons of such means on common cells.
##
## With weights = "equal", a mean that averages over cells of the interacting
## factors that have no data (missing by accident, not absent by design as the
## combinations of nested factors are) is the average of its observed cells:
## the rule of predictmeans 1.x (ASReml and Genstat: !PRESENT). Two such means
## average over different cells, so their difference also carries the effects
## of the other factors: with an empty cell the comparison rejected 8-78% of
## true null hypotheses in dev/sim-empty-cells.R (dev/review-defaults.md).
##
## So a comparison (a contrast whose coefficients sum to zero) that involves
## such a mean is made on the cells that all its means have: at each
## combination of the other factors observed for every one of them, the
## contrast of the cell means, averaged with the usual (equal) weights; 4-6%
## in the same simulation. Factors nested in a factor of the term (or a factor
## of the term is nested in them) cannot be matched and are averaged within
## each mean, as for the means. With no common cell the comparison is NA.

## The variables of each term of the model: a list of character vectors, one
## per term (e.g. A:B gives c("A", "B"), A:poly(x, 2) gives c("A", "x")).
pm_term_sets <- function(Terms) {
  fac_mat <- attr(Terms, "factors")
  if (!length(fac_mat)) return(list())
  row_vars <- lapply(rownames(fac_mat), function(r) {
    e <- tryCatch(str2lang(r), error = function(err) NULL)
    if (is.null(e)) r else all.vars(e)
  })
  lapply(seq_len(ncol(fac_mat)), function(t) unique(unlist(row_vars[fac_mat[, t] > 0])))
}

## Pairs of nested factors, from the data alone: b is nested in a when each
## level of b occurs with one level of a (and a has several levels). Their
## other combinations cannot exist. Returns a list of list(a, b, seen), 'seen'
## the combinations that occur.
pm_nested_pairs <- function(data, vars) {
  out <- list()
  if (length(vars) < 2) return(out)
  code <- lapply(data[vars], function(z) as.integer(factor(z)))
  nlev <- vapply(code, function(z) length(unique(z)), integer(1))
  for (a in vars) for (b in vars) {
    if (a == b || nlev[[a]] < 2) next
    if (length(unique(code[[a]] * (max(code[[b]]) + 1) + code[[b]])) != nlev[[b]]) next
    out[[length(out) + 1]] <- list(a = a, b = b,
                                   seen = unique(paste(as.character(data[[a]]), as.character(data[[b]]), sep = "\r")))
  }
  out
}

## TRUE for the rows of 'grid' that hold a combination of nested factors that
## cannot exist ('pairs': pm_nested_pairs()).
pm_absent_by_design <- function(grid, pairs) {
  out <- rep(FALSE, nrow(grid))
  for (p in pairs) {
    out <- out | !(paste(as.character(grid[[p$a]]), as.character(grid[[p$b]]), sep = "\r") %in% p$seen)
  }
  out
}

## The cells of the model and its empty cells.
##
## A cell is a combination of the levels of the factors of one interaction
## term (the term with the most factors of each group of terms: A:B:C, not
## also A:B). Factors that only have main effects, or that interact with
## other factors only, do not make cells together: in y ~ A*B + C the cells
## are those of A:B, and a combination of A and C without data is not an
## empty cell (the model is additive in C). A cell without data that could
## exist (not a combination of nested factors, 'nested') is an empty cell:
## with weights = "equal" the means average over the observed cells only.
##
## Returns NULL (no cell is without data) or list(present = function(grid) ->
## observed?, inter_vars (the factors that interact), factors (variables x
## terms), design = function(grid) -> exists by design?, unmatched (factors
## nested in or nesting the term's factors), missing (labels of the empty
## cells), missing_term (the term of each)). apply_rule: FALSE when the rows
## of K are not defined by a factor.
pm_empty_rule <- function(Terms, factor_names, preddf, facs, nested = list(), apply_rule = TRUE,
                          confounded = list()) {
  if (!apply_rule) return(NULL)
  # 'confounded': pairs of factors that share no term but are confounded in
  # the data (Kmatrix()): their combinations make cells as a term's do
  sets <- unique(c(lapply(pm_term_sets(Terms), function(s) intersect(factor_names, s)), confounded))
  sets <- sets[lengths(sets) > 1]
  sets <- sets[!vapply(seq_along(sets), function(i) any(vapply(sets[-i], function(o) all(sets[[i]] %in% o) &&
                                                                 length(o) > length(sets[[i]]), logical(1))), logical(1))]
  if (!length(sets)) return(NULL)
  cell_key <- function(df, vs) do.call(paste, c(lapply(vs, function(v) as.character(df[[v]])), sep = "\r"))
  seen <- lapply(sets, function(vs) unique(cell_key(preddf, vs)))
  inter_vars <- factor_names[factor_names %in% unlist(sets)]
  inter_pairs <- nested[vapply(nested, function(p) all(c(p$a, p$b) %in% inter_vars), logical(1))]
  missing <- missing_term <- character(0)
  any_unobserved <- FALSE
  for (i in seq_along(sets)) {
    vs <- sets[[i]]
    all_cells <- do.call(expand.grid, c(lapply(preddf[vs], function(z) levels(factor(z))), KEEP.OUT.ATTRS = FALSE,
                                        stringsAsFactors = FALSE))
    unobs <- !(cell_key(all_cells, vs) %in% seen[[i]])
    any_unobserved <- any_unobserved || any(unobs)
    acc <- all_cells[unobs & !pm_absent_by_design(all_cells, inter_pairs[vapply(inter_pairs, function(p)
      all(c(p$a, p$b) %in% vs), logical(1))]), , drop = FALSE]
    if (nrow(acc)) {
      missing <- c(missing, do.call(paste, c(lapply(acc, as.character), sep = ":")))
      missing_term <- c(missing_term, rep(paste(vs, collapse = ":"), nrow(acc)))
    }
  }
  if (!any_unobserved) return(NULL)
  fac_mat <- attr(Terms, "factors")
  if (length(fac_mat)) rownames(fac_mat) <- pm_term_vars(rownames(fac_mat))
  list(present = function(grid) Reduce(`&`, lapply(seq_along(sets), function(i) cell_key(grid, sets[[i]]) %in% seen[[i]])),
       inter_vars = inter_vars,
       factors = if (length(fac_mat)) fac_mat,
       design = function(grid) !pm_absent_by_design(grid, inter_pairs),
       unmatched = pm_unmatched_factors(preddf, setdiff(inter_vars, facs), intersect(facs, factor_names)),
       missing = missing, missing_term = missing_term, confounded = confounded)
}

## The means over empty cells of a reference grid, with a warning naming them,
## or NULL when no mean averages over an empty cell.
##   rg           the reference grid (pm_make_refgrid())
##   rule         pm_empty_rule()
##   basis        function(grid rows) -> basis matrix, as for K
##   coef_names   the columns of K
##   N            null space of the design (rank-deficient fits) or NULL
##   labels       row names of K
##   modelterm    for the warning
## Value: list(affected (logical, by row of K), over (by row of K: the cells
## a mean over empty cells averages, e.g. "N 0, 0.2, 0.4"; NA for the other
## means), missing (labels of the empty cells), compare = function(C,
## together = FALSE)).
pm_empty_cells <- function(rg, rule, basis, coef_names, N, labels, modelterm) {
  if (!length(rule$missing)) return(NULL)
  observed <- rule$present(rg$grid)
  design <- rule$design(rg$grid)
  unmatched <- rule$unmatched
  inter_vars <- rule$inter_vars
  key <- rg$key
  accidental <- design & !observed
  affected <- rep(FALSE, nrow(rg$target))
  affected[unique(key[accidental])] <- TRUE
  affected[rg$empty] <- FALSE
  names(affected) <- labels
  if (!any(affected)) return(NULL)

  # The matching key of a comparison: the dimensions of the grid other than
  # the term's, less the factors that cannot be matched: those nested in (or
  # nesting) the term's factors, and those whose every model term also
  # contains a factor in which the compared means differ (e.g. Farm in
  # Sample_date/Farm when dates are compared: a farm's effect is a different
  # parameter at each date, so matching farms would remove nothing). The
  # factors left are those with effects shared by the compared means.
  dn <- names(rg$idx)
  base_dims <- setdiff(dn, rg$tdims)
  dim_unmatched <- vapply(rg$dim_vars[base_dims], function(v) any(v %in% unmatched), logical(1))
  base_dims <- base_dims[!dim_unmatched]
  tvars <- names(rg$target)
  fac <- rule$factors                                  # variables x terms (0/1)
  shared_with <- function(v, D) {                      # some term with v and no factor of D
    if (is.null(fac) || !v %in% rownames(fac)) return(TRUE)
    tv <- fac[v, ] > 0
    !length(D) || any(colSums(fac[intersect(D, rownames(fac)), tv, drop = FALSE] > 0) == 0)
  }
  matching_dims <- function(S) {
    D <- tvars[vapply(tvars, function(v) length(unique(as.character(rg$target[S, v]))) > 1, logical(1))]
    base_dims[vapply(rg$dim_vars[base_dims], function(vs) all(vapply(vs, shared_with, logical(1), D = D)),
                     logical(1))]
  }
  mkey_of <- function(mdims) {
    mk <- if (length(mdims)) do.call(paste, c(rg$idx[mdims], sep = "\r")) else rep("", length(key))
    match(mk, unique(mk))
  }
  w0 <- rg$w0

  describe <- function(rows, mdims, unm = unmatched, vars = NULL) {
    # the factors whose levels describe the common cells
    desc_vars <- if (!is.null(vars)) vars else intersect(setdiff(inter_vars, unm), unlist(rg$dim_vars[mdims]))
    if (!length(desc_vars)) return(NA_character_)
    g <- unique(rg$grid[rows, desc_vars, drop = FALSE])
    short <- function(x) if (length(x) > 6) c(x[1:6], "...") else x
    if (length(desc_vars) == 1) {
      return(paste0(desc_vars, " ", paste(short(as.character(sort(g[[1]]))), collapse = ", ")))
    }
    levs <- lapply(g, function(z) unique(as.character(z)))
    if (nrow(g) == prod(lengths(levs))) {             # all combinations: list each factor's levels
      return(paste(vapply(desc_vars, function(v) paste0(v, " ", paste(short(levs[[v]]), collapse = ", ")),
                          character(1)), collapse = "; "))
    }
    paste0(paste(desc_vars, collapse = ":"), " ",
           paste(short(do.call(paste, c(lapply(g, as.character), sep = ":"))), collapse = ", "))
  }

  compare <- function(C, together = FALSE) {
    C <- rbind(C)
    nr <- nrow(C)
    sub <- na <- rep(FALSE, nr)
    over <- rep(NA_character_, nr)
    ii <- jj <- integer(0)
    xx <- numeric(0)
    S_all <- which(colSums(abs(C)) > 0)
    for (r in seq_len(nr)) {
      cr <- C[r, ]
      S <- if (together) S_all else which(abs(cr) > 0)
      if (!length(S) || !any(affected[S]) || abs(sum(cr)) > 1e-8 * max(abs(cr))) next
      sub[r] <- TRUE
      mdims <- matching_dims(S)
      mkey <- mkey_of(mdims)
      sets <- lapply(S, function(k) unique(mkey[observed & key == k]))
      M <- Reduce(intersect, sets)
      if (!length(M)) {
        na[r] <- TRUE
        over[r] <- "none"
        next
      }
      rows <- which(observed & key %in% S & mkey %in% M)
      # within each mean and common cell: its observed grid rows (factors
      # that cannot be matched); across common cells: their base weights
      tot <- stats::ave(w0[rows], key[rows], mkey[rows], FUN = sum)
      first <- !duplicated(mkey[rows])
      om <- stats::setNames(w0[rows][first], mkey[rows][first])
      om <- om / sum(om)
      v <- cr[key[rows]] * om[as.character(mkey[rows])] * w0[rows] / tot
      ii <- c(ii, rep(r, length(rows)))
      jj <- c(jj, rows)
      xx <- c(xx, v)
      over[r] <- describe(rows[cr[key[rows]] != 0 | together], mdims)
    }
    L <- matrix(NA_real_, nr, length(coef_names), dimnames = list(rownames(C), coef_names))
    estimable <- rep(TRUE, nr)
    ok <- sub & !na
    if (any(ok)) {
      Wc <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(nr, nrow(rg$grid)))
      Lf <- as.matrix(pm_refgrid_average(list(grid = rg$grid, W = Wc, empty = rep(FALSE, nr)), basis))
      if (!is.null(N)) {
        np <- Lf[, rownames(N), drop = FALSE] %*% N
        estimable <- apply(abs(np), 1, max) <= 1e-6 * pmax(1, apply(abs(Lf), 1, max))
      }
      L[ok, ] <- Lf[ok, coef_names, drop = FALSE]
    }
    estimable[na] <- FALSE
    list(L = L, sub = sub, na = na, estimable = estimable, over = over)
  }
  short <- function(x) paste(if (length(x) > 6) c(x[1:6], "...") else x, collapse = ", ")
  # the means by their factor levels (without covariate values in the grid)
  fcols <- names(rg$target)[vapply(rg$target, function(v) is.factor(v) || is.character(v), logical(1))]
  fac_lab <- if (length(fcols)) do.call(paste, c(lapply(rg$target[fcols], as.character), sep = ":")) else labels
  aff_lab <- unique(fac_lab[affected])
  one <- length(aff_lab) == 1
  # the empty cells by their term: "A:B a3:b2" (several terms: joined by "; ")
  mt <- if (length(rule$missing_term) == length(rule$missing)) rule$missing_term else
    rep(paste(inter_vars, collapse = ":"), length(rule$missing))
  miss_lab <- paste(vapply(unique(mt), function(tm) paste(tm, short(rule$missing[mt == tm])), character(1)),
                    collapse = "; ")
  # factors confounded in the data: say so first, and how to say it in the model
  conf <- if (length(rule$confounded)) paste0(
    paste(vapply(rule$confounded, function(p) paste0("'", p[1], "' and '", p[2], "'"), character(1)), collapse = "; "),
    if (length(rule$confounded) == 1) " are" else ": each pair is", " confounded in the data: the means average over ",
    "the combinations that exist (to say so: ",
    paste(vapply(rule$confounded, function(p) paste0(p[1], "/", p[2]), character(1)), collapse = ", "),
    " in the formula, leave one factor out, or weights = \"cells\"). ")
  warning(conf, "No data for ", miss_lab, ". ",
          if (one) "The mean of " else "The means of ", modelterm, " = ", short(aff_lab),
          if (one) " averages only its observed cells; its pairs" else " average only their observed cells; their pairs",
          " use the cells both means have (column 'over'). See vignette(\"missing-cells\").",
          call. = FALSE)
  # the cells each mean over empty cells averages (as 'over' of the pairs)
  over <- stats::setNames(rep(NA_character_, length(affected)), labels)
  # (with the factors that cannot be matched; only the factors with levels
  # missing, unless only their combinations are)
  avg_vars <- intersect(inter_vars, unlist(rg$dim_vars[setdiff(dn, rg$tdims)]))
  for (k in which(affected)) {
    rows <- which(observed & key == k)
    part <- avg_vars[vapply(avg_vars, function(v) length(unique(rg$grid[rows, v])) < length(unique(rg$grid[[v]])),
                            logical(1))]
    over[k] <- describe(rows, NULL, vars = if (length(part)) part else avg_vars)
  }
  list(affected = affected, over = over, missing = rule$missing, compare = compare)
}

## Factors of 'cands' that cannot be matched across the levels of the factors
## 'facs': nested in one of them, or one of them nested in it.
pm_unmatched_factors <- function(data, cands, facs) {
  nested_in <- function(inner, outer) {
    a <- as.character(data[[outer]])
    b <- as.character(data[[inner]])
    length(unique(a)) > 1 && length(unique(b)) > 1 && all(tapply(a, b, function(v) length(unique(v)) == 1L))
  }
  cands[vapply(cands, function(b) any(vapply(facs, function(a) nested_in(b, a) || nested_in(a, b), logical(1))),
               logical(1))]
}

## Replace the rows of L (comparisons in coefficient space, from C %*% K) that
## involve means over empty cells by their comparisons on common cells.
## 'cells' is Kmatrix()$empty_cells; C has one column per row of that K;
## 'cols' are the coefficient columns of L. Returns list(L, sub, na, over).
pm_common_rows <- function(cells, C, L, cols, together = FALSE) {
  if (is.null(cells)) return(list(L = L, sub = rep(FALSE, nrow(L)), na = rep(FALSE, nrow(L)),
                                  over = rep(NA_character_, nrow(L))))
  cc <- cells$compare(C, together = together)
  ok <- cc$sub & !cc$na
  if (any(ok)) L[ok, ] <- cc$L[ok, cols, drop = FALSE]
  na <- cc$na | (cc$sub & !cc$estimable)
  list(L = L, sub = cc$sub, na = na, over = cc$over)
}

## Pairs of factors that may be confounded in the data: they share no term
## of the model (and are not declared nested) and not all their combinations
## occur. Returns a list of character vectors c(a, b). Whether a pair is the
## reason a mean cannot be estimated is found by trying (Kmatrix()).
pm_confounded_pairs <- function(data, factor_names, Terms, nested = list()) {
  out <- list()
  if (length(factor_names) < 2) return(out)
  sets <- pm_term_sets(Terms)
  together <- function(a, b) any(vapply(sets, function(s) all(c(a, b) %in% s), logical(1))) ||
    any(vapply(nested, function(p) setequal(c(p$a, p$b), c(a, b)), logical(1)))
  for (i in seq_len(length(factor_names) - 1)) for (j in (i + 1):length(factor_names)) {
    a <- factor_names[i]
    b <- factor_names[j]
    if (together(a, b)) next
    fa <- factor(data[[a]])
    fb <- factor(data[[b]])
    if (length(unique(paste(as.integer(fa), as.integer(fb)))) < nlevels(fa) * nlevels(fb)) {
      # outer factor first (for the hint "write a/b"): the levels of the inner
      # one occur with fewer levels of the other
      per_b <- mean(tapply(as.integer(fa), fb, function(v) length(unique(v))), na.rm = TRUE)
      per_a <- mean(tapply(as.integer(fb), fa, function(v) length(unique(v))), na.rm = TRUE)
      out[[length(out) + 1]] <- if (per_a < per_b) c(b, a) else c(a, b)
    }
  }
  out
}
