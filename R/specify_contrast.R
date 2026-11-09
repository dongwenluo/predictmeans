#' Contrasts Written with Level Names
#'
#' @description
#' Build a contrast matrix for \code{\link{contrastmeans}} by writing each contrast as an expression in the level names
#' of a model term, so there is no need to count columns, even when some treatment combinations are missing.
#' Each level name (in backticks when it is not a syntactic name, e.g. \code{`Wet:Bio`}) stands for that mean.
#' @details
#' The arguments start with a dot so that contrasts can be given any name, including \code{x}.
#' @usage
#' specify_contrast(.x, ..., .zero_sum = TRUE)
#' @param .x A result of \code{\link{Kmatrix}} or \code{\link{predictmeans}} for the term, or a character vector of its level
#'   names (as in the row names of \code{Kmatrix(model, modelterm)$K}).
#' @param \dots Named expressions, one per contrast, e.g. \code{Bio_vs_None = (`Wet:Bio` + `Dry:Bio`)/2 - (`Wet:None` + `Dry:None`)/2}.
#' @param .zero_sum If \code{TRUE} (the default), every contrast's weights must sum to zero. Use \code{FALSE} for other linear
#'   combinations of the means, such as the average of two means.
#' @return
#' A matrix with one row per contrast (named as given) and one column per level. Columns of levels that are not
#' mentioned are zero. Means that cannot be estimated (e.g. missing cells) cannot be used.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @seealso
#' \code{\link{contrastmeans}}
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   cm <- specify_contrast(Kmatrix(fm, "nitro:Variety"),
#'     N_effect_Victory = `0.6:Victory` - `0:Victory`,
#'     N_effect_vs_others = (`0.6:Victory` - `0:Victory`) -
#'                          ((`0.6:Golden Rain` - `0:Golden Rain`) + (`0.6:Marvellous` - `0:Marvellous`))/2)
#'   contrastmeans(fm, "nitro:Variety", ctrmatrix = cm)
#' @name specify_contrast
NULL

## Build a contrast matrix from expressions in the level names of a term.
##
## specify_contrast(Kmatrix(fm, "A:B"),
##   bio_vs_untreated = (`Wet:Bio` + `Dry:Bio`)/2 - (`Wet:None` + `Dry:None`)/2)
##
## Each level name stands for that mean (a unit vector over the levels), so
## any linear combination can be written without counting columns. The result
## can be passed to contrastmeans(ctrmatrix = ).
specify_contrast <- function(.x, ..., .zero_sum = TRUE) {
  # dot-prefixed arguments, so contrasts can be named anything (even "x")
  x <- .x
  zero_sum <- .zero_sum
  if (is.list(x) && !is.null(x$K)) {
    # a Kmatrix() result
    levs <- rownames(x$K)
    estimable <- attr(x$K, "estimable")
    if (!is.null(estimable)) estimable <- estimable & !apply(x$K, 1, function(r) any(is.na(r)))
  } else if (inherits(x, "predictmeans")) {
    vars <- unlist(strsplit(x$settings$modelterm, ":"))
    levs <- do.call(paste, c(lapply(x$means[vars], as.character), sep = ":"))
    estimable <- x$means$estimable
  } else if (is.matrix(x) && !is.null(rownames(x))) {
    levs <- rownames(x)
    estimable <- attr(x, "estimable")
  } else if (is.character(x)) {
    levs <- x
    estimable <- NULL
  } else {
    stop("'.x' must be a Kmatrix() or predictmeans() result, or a character vector of level names.")
  }
  if (anyDuplicated(levs)) stop("The level names are not unique.")

  exprs <- as.list(substitute(list(...)))[-1]
  if (length(exprs) == 0) stop("Please give at least one contrast, e.g. A_vs_B = `a` - `b`.")
  ctr_names <- names(exprs)
  if (is.null(ctr_names) || any(ctr_names == "")) {
    stop("Please name every contrast, e.g. specify_contrast(K, A_vs_B = `a` - `b`).")
  }

  # Each level name is bound to its unit vector; other names come from the caller.
  n <- length(levs)
  env <- new.env(parent = parent.frame())
  for (i in seq_len(n)) assign(levs[i], replace(numeric(n), i, 1), envir = env)
  M <- matrix(0, length(exprs), n, dimnames = list(ctr_names, levs))
  for (k in seq_along(exprs)) {
    v <- tryCatch(eval(exprs[[k]], env), error = function(e) {
      stop("In contrast ", sQuote(ctr_names[k]), ": ", conditionMessage(e),
           ". Level names are: ", paste(sprintf("`%s`", levs), collapse = ", "), call. = FALSE)
    })
    if (!is.numeric(v) || length(v) != n) {
      stop("Contrast ", sQuote(ctr_names[k]), " must be a combination of the level names, such as `a` - `b`.")
    }
    M[k, ] <- v
  }

  if (zero_sum) {
    bad <- abs(rowSums(M)) > sqrt(.Machine$double.eps) * pmax(1, rowSums(abs(M)))
    if (any(bad)) {
      stop("The weights of ", paste(sQuote(ctr_names[bad]), collapse = ", "),
           " do not sum to zero. Use .zero_sum = FALSE for other linear combinations of the means.")
    }
  }
  if (!is.null(estimable)) {
    used <- colSums(abs(M)) > 0
    if (any(used & !estimable)) {
      stop("These means cannot be estimated and cannot be used in a contrast: ",
           paste(levs[used & !estimable], collapse = ", "))
    }
  }
  M
}
