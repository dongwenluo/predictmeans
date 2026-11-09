#' Families of Contrasts among Means
#'
#' @description
#' Build a standard family of contrasts among the means of a model term, as a matrix for the \code{ctrmatrix} argument
#' of \code{\link{contrastmeans}}. The families are those of \pkg{emmeans}, with the same coefficients and names:
#' \describe{
#' \item{\code{"pairwise"}, \code{"revpairwise"}}{all pairs, earlier minus later level (or later minus earlier);}
#' \item{\code{"trt.vs.ctrl"}}{each level minus the reference level \code{ref} (Dunnett's comparisons; use
#'   \code{adj = "mvt"} in \code{contrastmeans} for the exact Dunnett adjustment);}
#' \item{\code{"consec"}}{each level minus the one before;}
#' \item{\code{"poly"}}{orthogonal polynomials (linear, quadratic, ...) with integer coefficients, for equally spaced
#'   levels;}
#' \item{\code{"eff"}}{each level minus the average of all levels;}
#' \item{\code{"del.eff"}}{each level minus the average of the other levels.}
#' }
#' With \code{by}, the family is built within each level (or combination of levels) of the \code{by} factor(s), as
#' \code{atvar} does for pairwise comparisons. Note that \code{contrastmeans} adjusts the p-values of all rows as one
#' family.
#' @usage
#' contrast_family(.x, family = c("pairwise", "revpairwise", "trt.vs.ctrl", "consec", "poly", "eff", "del.eff"),
#'                 ref = 1, by = NULL, max.degree = NULL)
#' @param .x A \code{\link{Kmatrix}} or \code{\link{predictmeans}} result for the model term, or a character vector of
#'   the level names (then \code{by} can't be used).
#' @param family The family of contrasts (see Description).
#' @param ref For \code{"trt.vs.ctrl"}: the reference level, by position or name (within each \code{by} group), the
#'   default is the first level.
#' @param by Optional name(s) of factor(s) of the term: build the family within each of their levels.
#' @param max.degree For \code{"poly"}: the highest degree, the default is the number of levels - 1 (at most 6).
#' @return
#' A matrix with one row per contrast (named as in \pkg{emmeans}, e.g. \code{"0.2 - 0"}, \code{"linear"}; with
#' \code{by} the group comes first, e.g. \code{"Victory: 0.2 - 0"}) and one column per level of the term. Levels whose
#' means can't be estimated are left out of the families.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   # Dunnett: each nitrogen level against no nitrogen, exact adjustment
#'   K <- Kmatrix(fm, "nitro")
#'   contrastmeans(fm, "nitro", ctrmatrix = contrast_family(K, "trt.vs.ctrl"), adj = "mvt")
#'   # linear and quadratic trends in nitrogen within each variety
#'   K2 <- Kmatrix(fm, "nitro:Variety")
#'   cf <- contrast_family(K2, "poly", by = "Variety", max.degree = 2)
#'   contrastmeans(fm, "nitro:Variety", ctrmatrix = cf)
#' @name contrast_family
NULL

contrast_family <- function(.x, family = c("pairwise", "revpairwise", "trt.vs.ctrl", "consec", "poly", "eff", "del.eff"),
                            ref = 1, by = NULL, max.degree = NULL) {
  family <- match.arg(family)
  x <- .x
  if (is.list(x) && !is.null(x$K)) {                 # Kmatrix() result
    levs <- rownames(x$K)
    fac <- x$fctnames
    estimable <- attr(x$K, "estimable")
    if (is.null(estimable)) estimable <- rep(TRUE, length(levs))
    estimable <- estimable & !apply(x$K, 1, function(r) any(is.na(r)))
  } else if (inherits(x, "predictmeans")) {
    vars <- unlist(strsplit(x$settings$modelterm, ":"))
    fac <- x$means[vars]
    levs <- do.call(paste, c(lapply(fac, as.character), sep = ":"))
    estimable <- x$means$estimable
  } else if (is.character(x)) {
    levs <- x
    fac <- NULL
    estimable <- rep(TRUE, length(levs))
  } else {
    stop("'.x' must be a Kmatrix() or predictmeans() result, or a character vector of level names.")
  }
  if (anyDuplicated(levs)) stop("The level names are not unique.")

  if (is.null(by)) {
    groups <- list(seq_along(levs))
    within <- levs
    gname <- ""
  } else {
    if (is.null(fac)) stop("'by' needs a Kmatrix() or predictmeans() result.")
    if (!all(by %in% names(fac))) {
      stop("'by' must name factor(s) of the term: ", paste(sQuote(names(fac)), collapse = ", "))
    }
    rest <- setdiff(names(fac), by)
    if (length(rest) == 0) stop("'by' can't hold all the factors of the term.")
    gkey <- do.call(paste, c(lapply(fac[by], as.character), sep = ":"))
    within <- do.call(paste, c(lapply(fac[rest], as.character), sep = ":"))
    gname <- unique(gkey)
    groups <- lapply(gname, function(g) which(gkey == g))
  }
  if (any(!estimable)) {
    message("Left out of the contrasts (means not estimable): ", paste(levs[!estimable], collapse = ", "))
  }

  rows <- list()
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    idx <- idx[estimable[idx]]
    if (length(idx) < 2) next
    C <- pm_family_coefs(family, within[idx], ref, max.degree)       # levels x contrasts
    M <- matrix(0, ncol(C), length(levs), dimnames = list(colnames(C), levs))
    M[, idx] <- t(C)
    if (!is.null(by)) rownames(M) <- paste(gname[k], rownames(M), sep = ": ")
    rows[[k]] <- M
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) stop("There are no contrasts: each group needs at least two estimable means.")
  do.call(rbind, rows)
}

## Coefficients of one family for levels 'levs' (levels x contrasts), as
## emmeans' *.emmc functions.
pm_family_coefs <- function(family, levs, ref = 1, max.degree = NULL) {
  k <- length(levs)
  pair <- function(i, j) { v <- numeric(k); v[i] <- 1; v[j] <- -1; v }
  if (family %in% c("pairwise", "revpairwise")) {
    ij <- which(upper.tri(diag(k)), arr.ind = TRUE)
    ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
    if (family == "revpairwise") ij <- ij[order(ij[, 2], ij[, 1]), c(2, 1), drop = FALSE]
    C <- vapply(seq_len(nrow(ij)), function(r) pair(ij[r, 1], ij[r, 2]), numeric(k))
    nm <- paste(levs[ij[, 1]], levs[ij[, 2]], sep = " - ")
  } else if (family == "trt.vs.ctrl") {
    r <- if (is.character(ref)) match(ref, levs) else as.integer(ref)
    if (length(r) != 1 || is.na(r) || r < 1 || r > k) {
      stop("'ref' must be one of the levels (by name or position): ", paste(levs, collapse = ", "))
    }
    others <- setdiff(seq_len(k), r)
    C <- vapply(others, function(i) pair(i, r), numeric(k))
    nm <- paste(levs[others], levs[r], sep = " - ")
  } else if (family == "consec") {
    C <- vapply(2:k, function(i) pair(i, i - 1), numeric(k))
    nm <- paste(levs[2:k], levs[1:(k - 1)], sep = " - ")
  } else if (family == "poly") {
    deg <- if (is.null(max.degree)) min(6, k - 1) else min(max.degree, k - 1)
    P <- stats::poly(seq_len(k), deg)
    C <- apply(P, 2, function(con) {     # integer coefficients, as emmeans
      con <- con/min(con[con > 0.01])
      z <- max(abs(con - round(con)))
      while (z > 0.05) {
        con <- con/z
        z <- max(abs(con - round(con)))
      }
      round(con)
    })
    nm <- c("linear", "quadratic", "cubic", "quartic", paste("degree", 5:20))[seq_len(deg)]
  } else {
    w <- if (family == "eff") 1/k else 1/(k - 1)
    C <- -w + diag(if (family == "eff") 1 else 1 + w, k)
    nm <- paste(levs, "effect")
  }
  C <- matrix(C, nrow = k)
  dimnames(C) <- list(levs, nm)
  C
}
