##------------------------------------------------------------------------
## Block-diagonal covariance matrices of nlme fits: lists of dense blocks,
## one per group (attribute "groups": the group of each observation, in the
## data's order). Used by build-variance-matrices.R and the permutation code.
##------------------------------------------------------------------------

## Split a matrix by a factor: its rows ("row"), columns ("col") or diagonal
## blocks ("both"); a vector gives diagonal matrices.
matrix_list <- function(x, fac, dim) {
  if (is.vector(x)) {
    if (dim != "both") stop("Object must be a matrix in order to subset by ", dim, ".")
    return(lapply(split(x, fac), function(v) diag(v, nrow = length(v))))
  }
  lapply(levels(fac), function(f) switch(dim,
                                         row = x[fac == f, , drop = FALSE],
                                         col = x[, fac == f, drop = FALSE],
                                         both = x[fac == f, fac == f, drop = FALSE]))
}

## A sparse n x n matrix with the blocks 'A' (a list, named as 'pos') at the
## observations 'pos' (a list of row indices per block).
pm_place_blocks <- function(A, pos, n) {
  nm <- if (is.null(names(A))) names(pos) else names(A)
  parts <- lapply(seq_along(A), function(k) {
    i <- pos[[nm[k]]]
    B <- as.matrix(A[[k]])
    data.frame(i = rep(i, times = length(i)), j = rep(i, each = length(i)), x = as.vector(B))
  })
  d <- do.call(rbind, parts)
  Matrix::sparseMatrix(i = d$i, j = d$j, x = d$x, dims = c(n, n))
}

## The connected components of the non-zero pattern of a symmetric matrix:
## an integer label per row (observations in one component are correlated).
pm_sparse_components <- function(S) {
  S <- methods::as(methods::as(S, "generalMatrix"), "TsparseMatrix")
  n <- nrow(S)
  lab <- seq_len(n)
  i <- S@i + 1L
  j <- S@j + 1L
  repeat {                                       # each row takes the smallest label among its neighbours
    m <- pmin(lab[i], lab[j])
    agg <- tapply(m, i, min)
    at <- as.integer(names(agg))
    new <- lab
    new[at] <- pmin(new[at], agg)
    new <- new[new]                              # labels are row numbers: jump to the label's label
    if (identical(new, lab)) break
    lab <- new
  }
  match(lab, unique(lab))
}

## Add the blocks 'small_mat' (one per level of 'indices') into 'big_mat'.
add_submatrices <- function(indices, small_mat, big_mat) {
  levs <- levels(indices)
  if (nlevels(indices) != length(small_mat)) stop("Levels of indices do not match entries of small_mat.")
  for (i in seq_along(levs)) {
    ind <- levs[i] == indices
    big_mat[ind, ind] <- big_mat[ind, ind] + small_mat[[i]]
  }
  big_mat
}

## Add the blocks of a finer grouping into those of a coarser one that nests
## it (crosswalk: data frame of the fine and the coarse group of each row).
add_bdiag <- function(small_mats, big_mats, crosswalk) {
  small_indices <- lapply(split(crosswalk[[1]], crosswalk[[2]]), droplevels)
  big_indices <- unique(crosswalk)
  big_indices <- big_indices[[2]][order(big_indices[[1]])]
  small_mats_list <- split(small_mats, big_indices)
  Map(add_submatrices, indices = small_indices, small_mat = small_mats_list, big_mat = big_mats)
}

## Rows of B reordered to follow the blocks of A (the block order); returns
## list(B reordered, the original row of each new row).
get_order <- function(A, B, block = attr(A, "groups")) {
  if (is.null(names(A))) names(A) <- seq_along(A)
  A_names <- names(A)
  C_indx <- rep(A_names, times = vapply(A, nrow, 1L))
  if (is.null(block)) block <- C_indx
  C <- B
  B_index <- rep(0, nrow(B))
  for (b in A_names) {
    ind <- block == b
    ind_C <- C_indx == b
    C[ind_C, ] <- B[ind, ]
    B_index[ind_C] <- which(ind)
  }
  list(C, B_index)
}
