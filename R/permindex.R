#' Permutation Index
#'
#' This function obtains permutation index for a dataset.
#'
#'
#' @param data Data object used in the \code{model} fitting.
#' @param block Name (in "quotes") for the blocking factor in the \code{data}.
#' @param group Name (in "quotes") for the group factor in the \code{data}.
#' @param nsim The number of permutations. The default is 4999.
#' @param seed Specify a random number generator seed, for reproducible
#' results.
#' @return A matrix has 'nsim' columns of permuted index.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
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
#' @export
permindex <- function(data, block=NULL, group=NULL, nsim=4999, seed) {
  if (any(!is.null(block), !is.null(group))) {
    data <- data[do.call(order, data[, c(block, group), drop=FALSE]),]
  }
  rownames(data) <- NULL
  n <- nrow(data)
  if (all(is.null(block), is.null(group))) {
    if (!missing(seed)) {
      set.seed(seed)
    }
    permIndex <- replicate(nsim, sample.int(n))
  } else if (all(!is.null(block), is.null(group))) {
    Bsplit <- split(seq_len(n), as.character(data[, block]))    # split 1:n seq into nB blocks
    if (!missing(seed)) {
      set.seed(seed+1)
    }
    permIndex <- replicate(nsim, unlist(lapply(Bsplit, function (x) x[sample.int(length(x))]), recursive = FALSE, use.names = FALSE))
  } else if (all(is.null(block), !is.null(group))) {   # data=data; group="group"; nsim=5
    if (!missing(seed)) {
      set.seed(seed+2)
    }
    permIndex <- replicate(nsim, {
      Gsplit <- split(seq_len(n), as.character(data[, group]))   # split 1:n seq into separate groups
      nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
      unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)
    })          # permute the order of the whole groups and repeat process nsim times
  } else {
    ind <- seq_len(n)
    names(ind) <- as.character(data[, group])
    Bsplit <- split(ind, as.character(data[, block]))    # split 1:n seq into blocks
    if (!missing(seed)) {
      set.seed(seed+3)
    }
    permIndex <- replicate(nsim, unlist(lapply(Bsplit, function (x) {    # within each block
      Gsplit <- split(x, names(x))  # split each block seq into groups
      nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
      unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)          # permute the order of groups
    }), recursive = FALSE, use.names = FALSE))
  }

  colnIndex <- apply(permIndex, 2, function(x) paste(x, collapse=""))
  colnames(permIndex) <- colnIndex
  permIndex <- permIndex[, unique(colnIndex)]
  colnames(permIndex) <- NULL
  return(permIndex)
}
