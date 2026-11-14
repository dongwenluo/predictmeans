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

## The cells: the combinations of the factors that interact in the model
## (blocks or other main-effect-only factors would make every cell look
## single, hiding structurally empty combinations such as nested
## treatments). When some have no data (and not every cell is single), means
## average only over the observed ones. Returns NULL (no empty cells: means
## average over all cells) or list(present = function(grid) -> observed?,
## inter_vars, design = function(grid) -> exists by design?, unmatched
## (factors nested in or nesting the term's factors), missing (labels of the
## cells missing not by design)). apply_rule: FALSE when the term is a
## covariate.
pm_empty_rule <- function(Terms, factor_names, preddf, facs, apply_rule = TRUE) {
  fac_mat <- attr(Terms, "factors")
  ord <- attr(Terms, "order")
  inter_rows <- if (length(fac_mat) && any(ord > 1))
    rownames(fac_mat)[rowSums(fac_mat[, ord > 1, drop = FALSE] > 0) > 0] else character(0)
  inter_vars <- factor_names[factor_names %in% inter_rows |
                               vapply(factor_names, function(v) any(grepl(paste0("\\b", v, "\\b"), inter_rows)),
                                      logical(1))]
  if (length(inter_vars) < 2) inter_vars <- factor_names
  ndf <- data.frame(table(preddf[, inter_vars, drop = FALSE]))
  if (!(apply_rule && max(c(0, ord)) > 1 && length(inter_vars) > 1 && any(ndf$Freq == 0) &&
        !setequal(unique(ndf$Freq), c(0, 1)))) return(NULL)
  cell_key <- function(df) do.call(paste, c(lapply(inter_vars, function(v) as.character(df[[v]])), sep = "\r"))
  seen <- unique(cell_key(preddf))
  # combinations absent by design (a factor nested in another: each of its
  # levels occurs with one level of the other) are not missing cells
  missing_cells <- ndf[ndf$Freq == 0 & !pm_nested_absent(ndf[, inter_vars, drop = FALSE], preddf, inter_vars),
                       inter_vars, drop = FALSE]
  if (length(fac_mat)) rownames(fac_mat) <- pm_term_vars(rownames(fac_mat))
  list(present = function(grid) cell_key(grid) %in% seen,
       inter_vars = inter_vars,
       factors = if (length(fac_mat)) fac_mat,
       design = function(grid) !pm_nested_absent(grid[, inter_vars, drop = FALSE], preddf, inter_vars),
       unmatched = pm_unmatched_factors(preddf, setdiff(inter_vars, facs), intersect(facs, factor_names)),
       missing = do.call(paste, c(lapply(missing_cells, as.character), sep = ":")))
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
## Value: list(affected (logical, by row of K), missing (labels of the empty
## cells), compare = function(C, together = FALSE)).
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

  describe <- function(rows, mdims) {
    # the factors whose levels describe the common cells
    desc_vars <- intersect(setdiff(inter_vars, unmatched), unlist(rg$dim_vars[mdims]))
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
  warning("No data for ", paste(inter_vars, collapse = ":"), " ", short(rule$missing), ". ",
          if (one) "The mean of " else "The means of ", modelterm, " = ", short(aff_lab),
          if (one) " averages only its observed cells; comparisons with it" else
            " average only their observed cells; comparisons with them",
          " use only the cells that the compared means all have (column 'over' of the pairs).",
          call. = FALSE)
  list(affected = affected, missing = rule$missing, compare = compare)
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
