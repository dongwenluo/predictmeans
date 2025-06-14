



##------------------------------------------------------------------------
## block-diagonal matrix addition, multiplication, and trace functions matrix_list
##------------------------------------------------------------------------

# turn matrix into a list of sub-matrices

sub_f <- function(x, fac, dim) {
  function(f) {
    switch(dim,
           row = x[fac == f, , drop = FALSE],
           col = x[, fac == f, drop = FALSE],
           both = x[fac == f, fac == f, drop = FALSE])
  }
}

matrix_list <- function(x, fac, dim) {
  if (is.vector(x)) {
    if (dim != "both") {
      stop(paste0("Object must be a matrix in order to subset by ", dim, "."))
    }
    x_list <- split(x, fac)
    lapply(x_list, function(x) {
      diag(x, nrow = length(x))
    })
  } else {
    lapply(levels(fac), sub_f(x, fac, dim))
  }
}

# turn block-diagonal into regular matrix

unblock <- function(A, block = attr(A, "groups")) {
  if (is.null(block)) {
    block <- factor(rep(names(A), times = sapply(A, function(x) {
      dim(x)[1]
    })))
  }
  n <- length(block)
  mat <- matrix(0, n, n)
  for (i in levels(block)) {
    index <- i == block
    mat[index, index] <- A[[i]]
  }
  return(mat)
}


# sum of two conformable block-diagonal matrices

sum_blockblock <- function(A, B) {
  mapply(function(a, b) {
    a + b
  },
  a = A,
  b = B,
  SIMPLIFY = FALSE)
}

# generic matrix minus block-diagonal

matrix_minus_block <- function(A, B, block = attr(B, "groups")) {
  if (is.null(block)) {
    block <- rep(names(B), times = sapply(B, function(x) {
      dim(x)[1]
    }))
  }

  mat <- A
  for (i in unique(block)) {
    index <- i == block
    mat[index, index] <- mat[index, index] - B[[i]]
  }
  return(mat)
}


# block-diagonal minus generic matrix

block_minus_matrix <- function(A, B, block = attr(A, "groups")) {
  if (is.null(block)) {
    block <- rep(names(A), times = sapply(A, function(x) {
      dim(x)[1]
    }))
  }
  mat <- -B
  for (i in unique(block)) {
    index <- i == block
    mat[index, index] <- mat[index, index] + A[[i]]
  }
  return(mat)
}

add_submatrices <- function(indices, small_mat, big_mat) {
  levs <- levels(indices)
  if (nlevels(indices) != length(small_mat)) {
    stop("Levels of indices do not match entries of small_mat.")
  }
  for (i in seq_along(levs)) {
    ind <- levs[i] == indices
    big_mat[ind, ind] <- big_mat[ind, ind] + small_mat[[i]]
  }
  big_mat
}

add_bdiag <- function(small_mats, big_mats, crosswalk) {
  small_indices <- lapply(split(crosswalk[[1]], crosswalk[[2]]), droplevels)
  big_indices <- unique(crosswalk)
  big_indices <- big_indices[[2]][order(big_indices[[1]])]
  small_mats_list <- split(small_mats, big_indices)
  Map(
    add_submatrices,
    indices = small_indices,
    small_mat = small_mats_list,
    big_mat = big_mats
  )
}

# sum of conformable diagonal matrix and block-diagonal matrix

add_diag <- function(d, M) {
  diag(M) <- diag(M) + d
  M
}

add_diag_bdiag <- function(diag_mats, big_mats) {
  Map(add_diag, d = diag_mats, M = big_mats)
}

get_order <- function(A, B, block = attr(A, "groups")) {
  if (is.null(names(A))) {
    names(A) <- seq_along(A)
  }
  A_names <- names(A)
  C_indx <- rep(A_names, times = sapply(A, function(x) {
    dim(x)[1]
  }))

  if (is.null(block)) {
    block <- C_indx
  }
  C <- B
  B_index <- rep(0, nrow(B))

  for (b in A_names) {
    ind <- block == b
    ind_C <- C_indx == b
    C[ind_C, ] <-  B[ind, ]
    B_index[ind_C] <- which(ind)
  }

  return(list(C, B_index))
}
