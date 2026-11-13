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

## The empty-cell bookkeeping of a reference grid, or NULL when no mean
## averages over an empty cell.
##   rg           the reference grid (pm_make_refgrid())
##   observed     grid rows whose cell (of 'inter_vars') has data
##   design       grid rows whose cell exists by design (not nested-absent)
##   unmatched    factors that cannot be matched across means (nesting)
##   inter_vars   the interacting factors (the cells)
##   basis        function(grid rows) -> basis matrix, as for K
##   coef_names   the columns of K
##   N            null space of the design (rank-deficient fits) or NULL
##   labels       row names of K
## Value: list(affected (logical, by row of K), missing_cells (labels),
## compare = function(C, together = FALSE)).
pm_empty_cells <- function(rg, observed, design, unmatched, inter_vars, basis, coef_names, N, labels) {
  key <- rg$key
  accidental <- design & !observed
  affected <- rep(FALSE, nrow(rg$target))
  affected[unique(key[accidental])] <- TRUE
  affected[rg$empty] <- FALSE
  names(affected) <- labels
  if (!any(affected)) return(NULL)

  # the matching key: every dimension of the grid except the term's and the
  # unmatched (nesting) factors
  dn <- names(rg$idx)
  nested_dims <- dn[vapply(rg$dim_vars[dn], function(v) any(v %in% unmatched), logical(1))]
  mdims <- setdiff(dn, c(rg$tdims, nested_dims))
  mk <- if (length(mdims)) do.call(paste, c(rg$idx[mdims], sep = "\r")) else rep("", length(key))
  mkey <- match(mk, unique(mk))
  w0 <- rg$w0
  # the factors whose levels describe the common cells
  desc_vars <- intersect(setdiff(inter_vars, unmatched), unlist(rg$dim_vars[mdims]))

  describe <- function(rows) {
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
      over[r] <- describe(rows[cr[key[rows]] != 0 | together])
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
  list(affected = affected, compare = compare)
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
