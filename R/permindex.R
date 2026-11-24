#' Permutation Index
#'
#' @description
#' Permutations allowed by the design: free, within blocks, of whole groups (e.g. whole plots), or of whole groups
#' within blocks. Give the result to \code{\link{permmodels}} or \code{\link{permlmer}} (argument \code{perms}) for
#' restricted (design-based) permutation tests; when its columns are every permutation the design allows, use
#' \code{exact = TRUE} in \code{permmodels()}.
#' @usage
#' permindex(data, block=NULL, group=NULL, nsim=4999, seed)
#' @param data The data the model was fitted to.
#' @param block Optional name (in "quotes") of the blocking factor in \code{data}.
#' @param group Optional name (in "quotes") of the factor whose levels are permuted as whole units.
#' @param nsim The number of permutations. The default is 4999.
#' @param seed Specify a random number generator seed, for reproducible results. The user's random number stream is
#'   left as it was.
#' @return
#' A matrix with one row per row of \code{data} and up to \code{nsim} columns, each a permutation of the row numbers
#' of \code{data} (row i of the permuted data is row \code{index[i]} of \code{data}). Duplicated permutations are
#' dropped, so there are fewer columns when the design allows fewer distinct permutations.
#' @details
#' Groups are permuted as whole units, so the groups that change places (within a block, if \code{block} is given)
#' must have the same number of rows. The data need not be sorted: before predictmeans 2.0 the indices referred to
#' the data sorted by \code{block} and \code{group}, so they were right only for data already in that order (the
#' permutations are the same as before for sorted data and the same \code{seed}).
#' @seealso \code{\link{permmodels}}, \code{\link{permlmer}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   block <- rep(1:3, each=12)
#'   group <- rep(rep(1:3, each=4), 3)
#'   data <- data.frame(block, group)
#'   cbind(data, permindex(data, block="block", group="group", nsim=5))  
#'                         # Permute group as a whole within each block first, 
#'                         # then permute obs within each group.
#'   cbind(data, permindex(data, block="block",  nsim=5)) 
#'                         # Permute obs within each block only.
#'   cbind(data, permindex(data, group="group", nsim=5)) 
#'                         # Permute groups as a whole block first, 
#'                         # then permute obs within each group.
#'   cbind(data, permindex(data, nsim=5))  # Free permutation.
#'
#'   ## A randomisation test of nitrogen in a randomised block design:
#'   ## permutations within blocks only
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lm(yield ~ Block + nitro, data = Oats)
#'   idx <- permindex(Oats, block = "Block", nsim = 199, seed = 1)
#'   permmodels(fm, perms = idx, ncore = 1)
#' @name permindex
NULL

#' @export
permindex <- function(data, block=NULL, group=NULL, nsim=4999, seed) {
  data <- as.data.frame(data)
  n <- nrow(data)
  # The permutations are drawn for the data sorted by block and group (as in
  # predictmeans 1.x) and then mapped back to the data's own row order.
  ord <- if (any(!is.null(block), !is.null(group))) {
    do.call(order, unname(as.list(data[, c(block, group), drop = FALSE])))
  } else seq_len(n)
  data <- data[ord, , drop = FALSE]
  rownames(data) <- NULL
  seed <- if (missing(seed)) NULL else seed
  draw <- function(expr, offset) if (is.null(seed)) expr else pm_with_seed(seed + offset, expr)
  # groups that change places must have the same number of rows
  same_size <- function(g) {
    sz <- table(as.character(g))
    if (length(unique(sz)) > 1) {
      stop("Whole groups can only change places when they have the same number of rows; the group sizes",
           if (!is.null(block)) " within a block", " are ", paste(unique(sz), collapse = ", "), ".", call. = FALSE)
    }
  }
  if (all(is.null(block), is.null(group))) {
    permIndex <- draw(replicate(nsim, sample.int(n)), 0)
  } else if (all(!is.null(block), is.null(group))) {
    Bsplit <- split(seq_len(n), as.character(data[, block]))    # split 1:n seq into nB blocks
    permIndex <- draw(replicate(nsim, unlist(lapply(Bsplit, function (x) x[sample.int(length(x))]), recursive = FALSE, use.names = FALSE)), 1)
  } else if (all(is.null(block), !is.null(group))) {
    same_size(data[, group])
    permIndex <- draw(replicate(nsim, {
      Gsplit <- split(seq_len(n), as.character(data[, group]))   # split 1:n seq into separate groups
      nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
      unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)
      }), 2)          # permute the order of the whole groups and repeat process nsim times
  } else {
    for (b in unique(as.character(data[, block]))) same_size(data[as.character(data[, block]) == b, group])
    ind <- seq_len(n)
    names(ind) <- as.character(data[, group])
    Bsplit <- split(ind, as.character(data[, block]))    # split 1:n seq into blocks
    permIndex <- draw(replicate(nsim, unlist(lapply(Bsplit, function (x) {    # within each block
        Gsplit <- split(x, names(x))  # split each block seq into groups
        nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
        unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)          # permute the order of groups
      }), recursive = FALSE, use.names = FALSE)), 3)
  }
  permIndex <- matrix(permIndex, nrow = n)
  # drop duplicated permutations (the separator keeps 1,12 apart from 11,2)
  permIndex <- permIndex[, !duplicated(apply(permIndex, 2, paste, collapse = ",")), drop = FALSE]
  # back to the data's own order: original row ord[i] takes original row ord[permIndex[i, ]]
  out <- permIndex
  out[ord, ] <- ord[permIndex]
  out
}
