adiag <- function(..., pad = as.integer(0), do.dimnames = TRUE) {
  ## function from package 'magic'
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1], list(pad = pad)))
    return(do.call("Recall", c(
      list(args[[1]]), list(jj),
      list(pad = pad)
    )))
  }
  a <- args[[1]]
  b <- args[[2]]
  if (is.null(b)) {
    return(a)
  }
  if (is.null(dim(a)) && is.null(dim(b))) {
    dim(a) <- rep(1, 2)
    dim(b) <- rep(1, 2)
  }
  if (is.null(dim(a)) && length(a) == 1) {
    dim(a) <- rep(1, length(dim(b)))
  }
  if (is.null(dim(b)) && length(b) == 1) {
    dim(b) <- rep(1, length(dim(a)))
  }
  if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
    stop("a and b must have identical number of dimensions")
  }
  s <- array(pad, dim.a + dim.b)
  s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
  ind <- lapply(seq(dim.b), function(i) {
    seq_len(dim.b[[i]]) +
      dim.a[[i]]
  })
  out <- do.call("[<-", c(list(s), ind, list(b)))
  n.a <- dimnames(a)
  n.b <- dimnames(b)
  if (do.dimnames && !is.null(n.a) && !is.null(n.b)) {
    dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
    names(dimnames(out)) <- names(n.a)
  }
  return(out)
}

vec2mat2 <- function(x, sep = "-") {
  splits <- strsplit(x, sep)
  n.spl <- sapply(splits, length)
  if (any(n.spl != 2)) {
    stop(
      "Names must contain exactly one '", sep, "' each;  instead got ",
      paste(x, collapse = ", ")
    )
  }
  x2 <- t(as.matrix(as.data.frame(splits)))
  dimnames(x2) <- list(x, NULL)
  x2
}

multcompLetters <- function(x, compare = "<", threshold = 0.05,
                            Letters = c(letters, LETTERS, "."),
                            reversed = FALSE) {
  # Fast input validation and conversion
  if (!is.logical(x)) {
    x <- do.call(compare, list(x, threshold))
  }

  Lvls <- rownames(x)
  if (is.null(Lvls)) {
    stop("Matrix must have row names")
  }
  n <- nrow(x)

  # Quick return for no differences
  if (!any(x[lower.tri(x)])) {
    return(structure(rep(Letters[1], n), names = Lvls))
  }

  # Initialize LetMat more efficiently
  LetMat <- matrix(TRUE, n, 1L, dimnames = list(Lvls, NULL))

  # Get significant pairs more efficiently using which() with arr.ind = TRUE
  sig_indices <- which(x[lower.tri(x)])
  if (length(sig_indices) == 0) {
    return(structure(rep(Letters[1], n), names = Lvls))
  }

  lt_indices <- which(lower.tri(x), arr.ind = TRUE)
  sig_pairs <- lt_indices[sig_indices, , drop = FALSE]

  # Optimized absorb function
  absorb <- function(mat) {
    nc <- ncol(mat)
    if (nc <= 1) {
      return(mat)
    }

    # Pre-compute column sums for efficiency
    col_sums <- colSums(mat)

    i <- 1
    while (i < nc) {
      j <- i + 1
      while (j <= nc) {
        # Quick check using column sums
        if (col_sums[i] >= col_sums[j]) {
          if (all(mat[mat[, j], i])) {
            mat <- mat[, -j, drop = FALSE]
            col_sums <- col_sums[-j]
            nc <- nc - 1
            next
          }
        } else {
          if (all(mat[mat[, i], j])) {
            mat <- mat[, -i, drop = FALSE]
            col_sums <- col_sums[-i]
            nc <- nc - 1
            i <- i - 1
            break
          }
        }
        j <- j + 1
      }
      i <- i + 1
    }
    mat
  }

  # Process pairs more efficiently
  for (i in seq_len(nrow(sig_pairs))) {
    idx1 <- sig_pairs[i, 1]
    idx2 <- sig_pairs[i, 2]

    common_cols <- which(LetMat[idx1, ] & LetMat[idx2, ])

    if (length(common_cols)) {
      new_col <- LetMat[, common_cols, drop = FALSE]
      new_col[idx1, ] <- FALSE
      LetMat[idx2, common_cols] <- FALSE
      LetMat <- cbind(LetMat, new_col)

      if (ncol(LetMat) > 2) {
        LetMat <- absorb(LetMat)
      }
    }
  }

  # Optimized column sorting
  sort_cols <- function(mat) {
    if (ncol(mat) <= 1) {
      return(mat)
    }

    first_true <- apply(mat, 2, function(x) match(TRUE, x, nomatch = nrow(mat) + 1))
    ord <- order(first_true)
    mat <- mat[, ord, drop = FALSE]

    # Handle ties more efficiently
    rle_result <- rle(first_true[ord])
    if (any(rle_result$lengths > 1)) {
      pos <- 1
      for (i in seq_along(rle_result$lengths)) {
        len <- rle_result$lengths[i]
        if (len > 1) {
          start_row <- rle_result$values[i]
          if (start_row < nrow(mat)) {
            idx <- pos:(pos + len - 1)
            sub_mat <- mat[(start_row + 1):nrow(mat), idx, drop = FALSE]
            if (nrow(sub_mat) > 0) {
              mat[(start_row + 1):nrow(mat), idx] <- sort_cols(sub_mat)
            }
          }
        }
        pos <- pos + len
      }
    }
    mat
  }

  # Sort and potentially reverse columns
  LetMat <- sort_cols(LetMat)
  if (reversed) {
    LetMat <- LetMat[, rev(seq_len(ncol(LetMat))), drop = FALSE]
  }

  # Efficient letter assignment
  k_ltrs <- ncol(LetMat)
  Letters <- if (k_ltrs <= length(Letters)) {
    Letters[seq_len(k_ltrs)]
  } else {
    make_letters <- function(k, ltrs) {
      if (k <= length(ltrs)) {
        return(ltrs[seq_len(k)])
      }
      extra <- paste0(
        ltrs[length(ltrs)],
        c(ltrs[-length(ltrs)], ltrs[length(ltrs)])
      )
      c(
        ltrs[-length(ltrs)],
        make_letters(k - length(ltrs) + 1, extra)
      )
    }
    make_letters(k_ltrs, Letters)
  }

  # Pre-compute letter widths and blanks
  letter_widths <- nchar(Letters)
  blanks <- vapply(letter_widths, function(w) {
    paste(rep(" ", w), collapse = "")
  }, character(1))

  # Create result vector efficiently
  result <- character(n)
  names(result) <- Lvls

  # Use vectorized operations where possible
  for (i in seq_len(n)) {
    letter_pos <- which(LetMat[i, ])
    if (length(letter_pos)) {
      chars <- character(k_ltrs)
      chars[] <- blanks
      chars[letter_pos] <- Letters[letter_pos]
      result[i] <- paste(chars, collapse = "")
    } else {
      result[i] <- paste(rep(" ", sum(letter_widths)), collapse = "")
    }
  }

  result
}

######################## Functions from lmerTest begin #################

get_contrasts_type1 <- function(model) {
  terms <- terms(model)
  X <- model.matrix(model)
  p <- ncol(X)
  if (p == 0L) {
    return(list(matrix(numeric(0L), nrow = 0L))) # no fixef
  }
  if (p == 1L && attr(terms, "intercept")) { # intercept-only model
    return(list(matrix(numeric(0L), ncol = 1L)))
  }
  # Compute 'normalized' doolittle factorization of XtX:
  L <- if (p == 1L) {
    matrix(1L)
  } else {
    t(doolittle(crossprod(X))$L)
  }
  dimnames(L) <- list(colnames(X), colnames(X))
  # Determine which rows of L belong to which term:
  ind.list <- term2colX(terms, X)[attr(terms, "term.labels")]
  lapply(ind.list, function(rows) L[rows, , drop = FALSE])
}

term2colX <- function(terms, X) {
  # Compute map from terms to columns in X using the assign attribute of X.
  # Returns a list with one element for each term containing indices of columns
  #   in X belonging to that term.
  if (is.null(asgn <- attr(X, "assign"))) {
    stop("Invalid design matrix:",
      "design matrix 'X' should have a non-null 'assign' attribute",
      call. = FALSE
    )
  }
  term_names <- attr(terms, "term.labels")
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if (has_intercept) {
    c("(Intercept)", term_names)[asgn + 1]
  } else {
    term_names[asgn[asgn > 0]]
  }
  if (!length(col_terms) == ncol(X)) {
    # should never happen.
    stop("An error happended when mapping terms to columns of X")
  }
  # get names of terms (including aliased terms)
  nm <- union(unique(col_terms), term_names)
  res <- lapply(setNames(as.list(nm), nm), function(x) numeric(0L))
  map <- split(seq_along(col_terms), col_terms)
  res[names(map)] <- map
  res[nm] # order appropriately
}

doolittle <- function(x, eps = 1e-6) {
  if (!is.matrix(x) || ncol(x) != nrow(x) || !is.numeric(x)) {
    stop("argument 'x' should be a numeric square matrix")
  }
  stopifnot(ncol(x) > 1L)
  n <- nrow(x)
  L <- U <- matrix(0, nrow = n, ncol = n)
  diag(L) <- rep(1, n)
  for (i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for (j in 1:n) {
      U[i, j] <- x[i, j]
      if (im1 > 0) {
        for (k in 1:im1) {
          U[i, j] <- U[i, j] - L[i, k] * U[k, j]
        }
      }
    }
    if (ip1 <= n) {
      for (j in ip1:n) {
        L[j, i] <- x[j, i]
        if (im1 > 0) {
          for (k in 1:im1) {
            L[j, i] <- L[j, i] - L[j, k] * U[k, i]
          }
        }
        L[j, i] <- if (abs(U[i, i]) < eps) {
          0
        } else {
          L[j, i] / U[i, i]
        }
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list(L = L, U = U)
}

########################

as_lmerModLT <- function(model, devfun, tol = 1e-8) {
  is_reml <- getME(model, "is_REML")
  # Coerce 'lme4-model' to 'lmerModLmerTest':
  res <- as(model, "lmerModLmerTest")
  # Set relevant slots of the new model object:
  res@sigma <- sigma(model)
  res@vcov_beta <- as.matrix(vcov(model))
  varpar_opt <- unname(c(res@theta, res@sigma))
  # Compute Hessian:
  h <- numDeriv::hessian(
    func = devfun_vp, x = varpar_opt, devfun = devfun,
    reml = is_reml
  )
  # Eigen decompose the Hessian:
  eig_h <- eigen(h, symmetric = TRUE)
  evals <- eig_h$values
  neg <- evals < -tol
  pos <- evals > tol
  zero <- evals > -tol & evals < tol
  if (sum(neg) > 0) { # negative eigenvalues
    eval_chr <- ifelse(sum(neg) > 1, "eigenvalues", "eigenvalue")
    evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
    warning(sprintf(
      "Model failed to converge with %d negative %s: %s",
      sum(neg), eval_chr, evals_num
    ), call. = FALSE)
  }
  # Note: we warn about negative AND zero eigenvalues:
  if (sum(zero) > 0) { # some eigenvalues are zero
    eval_chr <- ifelse(sum(zero) > 1, "eigenvalues", "eigenvalue")
    evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
    warning(sprintf(
      "Model may not have converged with %d %s close to zero: %s",
      sum(zero), eval_chr, evals_num
    ))
  }
  # Compute vcov(varpar):
  pos <- eig_h$values > tol
  q <- sum(pos)
  # Using the Moore-Penrose generalized inverse for h:
  h_inv <- with(eig_h, {
    vectors[, pos, drop = FALSE] %*% diag(1 / values[pos], nrow = q) %*%
      t(vectors[, pos, drop = FALSE])
  })
  res@vcov_varpar <- 2 * h_inv # vcov(varpar)
  # Compute Jacobian of cov(beta) for each varpar and save in list:
  Jac <- numDeriv::jacobian(func = get_covbeta, x = varpar_opt, devfun = devfun)
  res@Jac_list <- lapply(seq_len(ncol(Jac)), function(i) {
    array(Jac[, i], dim = rep(length(res@beta), 2))
  }) # k-list of jacobian matrices
  res
}

##########
devfun_vp <- function(varpar, devfun, reml) {
  nvarpar <- length(varpar)
  sigma2 <- varpar[nvarpar]^2
  theta <- varpar[-nvarpar]
  df_envir <- environment(devfun)
  devfun(theta) # Evaluate deviance function at varpar
  n <- nrow(df_envir$pp$V)
  # Compute deviance for ML:
  dev <- df_envir$pp$ldL2() + (df_envir$resp$wrss() + df_envir$pp$sqrL(1)) / sigma2 +
    n * log(2 * pi * sigma2)
  if (!reml) {
    return(dev)
  }
  # Adjust if REML is used:
  RX <- df_envir$pp$RX() # X'V^{-1}X ~ crossprod(RX^{-1}) = cov(beta)^{-1} / sigma^2
  dev + 2 * c(determinant(RX)$modulus) - ncol(RX) * log(2 * pi * sigma2)
}

##########
get_covbeta <- function(varpar, devfun) {
  nvarpar <- length(varpar)
  sigma <- varpar[nvarpar] # residual std.dev.
  theta <- varpar[-nvarpar] # ranef var-par
  devfun(theta) # evaluate REML or ML deviance 'criterion'
  df_envir <- environment(devfun) # extract model environment
  sigma^2 * tcrossprod(df_envir$pp$RXi()) # vcov(beta)
}

######################## Functions from lmerTest end #################

#######################
# https://rpubs.com/bbolker/waldvar
# http://rstudio-pubs-static.s3.amazonaws.com/28864_dd1f084207d54f5ea67c9d1a9c845d01.html

#######################################################
# from package merDeriv vcov.lmerMod.R

vcov_lmerMod <- function(object, ...) {
  if (!is(object, "lmerMod")) {
    stop("vcov.lmerMod() only works for lmer() models.")
  }
  dotdotdot <- list(...)
  if ("full" %in% names(dotdotdot)) {
    full <- dotdotdot$full
    if (!is.logical(full)) {
      stop("Argument full should be a logical value: TRUE or FALSE.")
    }
  } else {
    full <- FALSE
  }

  if ("information" %in% names(dotdotdot)) {
    ## allows partial argument matching
    ## doesn't need further argument checking
    information <- match.arg(
      dotdotdot$information,
      c("expected", "observed")
    )
  } else {
    information <- "expected"
  }

  if ("ranpar" %in% names(dotdotdot)) {
    ranpar <- dotdotdot$ranpar
  } else {
    ranpar <- "var"
  }

  parts <- getME(object, "ALL")
  yXbe <- parts$y - tcrossprod(parts$X, t(parts$beta))
  uluti <- length(parts$theta)
  Zlam <- tcrossprod(parts$Z, parts$Lambdat)
  V <- (tcrossprod(Zlam, Zlam) + Matrix::Diagonal(parts$n, 1)) * (parts$sigma)^2
  M <- solve(chol(V))
  invV <- tcrossprod(M, M)
  LambdaInd <- parts$Lambda
  LambdaInd@x[] <- parts$Lind
  invVX <- crossprod(parts$X, invV)
  Pmid <- solve(crossprod(parts$X, t(invVX)))
  P <- invV - tcrossprod(crossprod(invVX, Pmid), t(invVX))
  fixvar <- solve(tcrossprod(crossprod(parts$X, invV), t(parts$X)))

  if (full == FALSE) {
    fixvar
  } else {
    fixhes <- tcrossprod(crossprod(parts$X, invV), t(parts$X))
    uluti <- length(parts$theta)
    devV <- vector("list", (uluti + 1))
    devLambda <- vector("list", uluti)
    score_varcov <- matrix(NA, nrow = length(parts$y), ncol = uluti)

    for (i in 1:uluti) {
      devLambda[[i]] <- Matrix::forceSymmetric(LambdaInd == i, uplo = "L")
      devV[[i]] <- tcrossprod(tcrossprod(parts$Z, t(devLambda[[i]])), parts$Z)
    }

    devV[[(uluti + 1)]] <- Matrix::Diagonal(nrow(parts$X), 1)
    ranhes <- matrix(NA, nrow = (uluti + 1), ncol = (uluti + 1))
    entries <- rbind(
      matrix(rep(1:(uluti + 1), each = 2), (uluti + 1), 2, byrow = TRUE),
      t(combn((uluti + 1), 2))
    )
    entries <- entries[order(entries[, 1], entries[, 2]), ]

    if (parts$devcomp$dims[["REML"]] == 0) {
      if (information == "expected") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <-
          apply(entries, 1, function(x) {
            as.numeric(0.5 * lav_matrix_trace(tcrossprod(
              tcrossprod(crossprod(invV, devV[[x[1]]]), invV),
              t(devV[[x[2]]])
            )))
          })
      } else if (information == "observed") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- unlist(apply(
          entries,
          1, function(x) {
            as.vector(-as.numeric((1 / 2) *
              lav_matrix_trace(tcrossprod(tcrossprod(crossprod(
                invV,
                devV[[x[1]]]
              ), invV), t(devV[[x[2]]])))) +
              tcrossprod((tcrossprod(
                (crossprod(yXbe, tcrossprod(tcrossprod(crossprod(
                  invV,
                  devV[[x[1]]]
                ), invV), t(devV[[x[2]]])))),
                invV
              )), t(yXbe)))
          }
        ))
      }
    }

    if (parts$devcomp$dims[["REML"]] > 0) {
      if (information == "expected") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(
          entries,
          1, function(x) {
            as.numeric((1 / 2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(
              P,
              devV[[x[1]]]
            ), P), t(devV[[x[2]]]))))
          }
        )
      } else if (information == "observed") {
        ranhes[lower.tri(ranhes, diag = TRUE)] <- apply(
          entries,
          1, function(x) {
            -as.numeric((1 / 2) * lav_matrix_trace(tcrossprod(tcrossprod(crossprod(
              P,
              devV[[x[1]]]
            ), P), t(devV[[x[2]]])))) + tcrossprod((tcrossprod((crossprod(
              yXbe,
              tcrossprod(tcrossprod(
                crossprod(invV, devV[[x[1]]]),
                P
              ), t(devV[[x[2]]]))
            )), invV)), t(yXbe))
          }
        )
      }
    }

    ranhes <- Matrix::forceSymmetric(ranhes, uplo = "L")
    if (information == "expected") {
      varcov_beta <- matrix(0, length(devV), length(parts$beta))
    } else if (information == "observed") {
      varcov_beta <- matrix(NA, length(devV), length(parts$beta))
      for (j in seq_along(devV)) {
        varcov_beta[j, ] <- as.vector(tcrossprod(
          crossprod(
            parts$X,
            (tcrossprod(crossprod(invV, devV[[j]]), invV))
          ),
          t(yXbe)
        ))
      }
    }

    if (ranpar == "var") {
      ranhes <- ranhes
      varcov_beta <- varcov_beta
    } else if (ranpar == "sd") {
      sdcormat <- as.data.frame(VarCorr(object, comp = "Std.Dev"),
        order = "lower.tri"
      )
      ## TODO: Not sure we really need to call which. The logicals should take care of the subsetting
      idx <- is.na(sdcormat$var2)
      sdcormat$sdcor2[which(idx)] <- sdcormat$sdcor[which(idx)] * 2
      sdcormat$sdcor2[which(!idx)] <- sdcormat$vcov[which(!idx)] / sdcormat$sdcor[which(!idx)]
      varcov_beta <- sweep(varcov_beta, MARGIN = 1, sdcormat$sdcor2, `*`)
      weight <- apply(entries, 1, function(x) {
        sdcormat$sdcor2[x[1]] *
          sdcormat$sdcor2[x[2]]
      })
      ranhes[lower.tri(ranhes, diag = TRUE)] <- weight *
        ranhes[lower.tri(ranhes, diag = TRUE)]
      ranhes <- Matrix::forceSymmetric(ranhes, uplo = "L")
    } else {
      stop("ranpar needs to be var or sd for lmerMod object.")
    }

    full_varcov <- solve(rbind(
      cbind(fixhes, t(varcov_beta)),
      cbind(varcov_beta, ranhes)
    ))

    colnames(full_varcov) <- c(
      names(parts$fixef),
      paste("cov", names(parts$theta), sep = "_"),
      "residual"
    )

    callingFun <- try(deparse(sys.call(-2)), silent = TRUE)
    if (length(callingFun) > 1) {
      callingFun <- paste(callingFun, collapse = "")
    }

    if (!inherits(callingFun, "try-error") && grepl("summary.merMod", callingFun)) {
      return(fixvar)
    } else {
      return(full_varcov)
    }
  }
}

vcov_glmerMod <- function(object, ...) {
  if (!(family(object)$family %in% c("binomial", "poisson"))) {
    stop("family has to be binomial or poisson")
  }

  dotdotdot <- list(...)
  if ("full" %in% names(dotdotdot)) {
    full <- dotdotdot$full
  } else {
    full <- FALSE
  }

  if ("ranpar" %in% names(dotdotdot)) {
    ranpar <- dotdotdot$ranpar
  } else {
    ranpar <- "var"
  }

  if (full == FALSE) {
    full_vcov <- vcov.merMod(object)
  } else {
    if (length(getME(object, "l_i")) > 1L) {
      stop("Multiple cluster variables detected. This type of model is currently not supported for full vcov.")
    }

    ## Hessian was based on deviance function, which is the
    ## -2*LogLik. That's why divided by -2
    if (object@devcomp$dims[["nAGQ"]] == 0L) {
      stop("For full vcov, nAGQ of at least 1 is required.")
    }
    full_vcov_noorder <- -solve(object@optinfo$derivs$Hessian / (-2))

    ## Block order in Hessian was theta, beta. Reorganize to
    ## put fixed parameter block first to match with score
    ## matrix order.

    ## count parameter numbers
    p <- nrow(full_vcov_noorder)
    pfix <- length(object@beta)
    pran <- length(object@theta)
    ## reorder four blocks
    full_vcov <- matrix(NA, nrow(full_vcov_noorder), ncol(full_vcov_noorder))
    full_vcov[1:pfix, 1:pfix] <- full_vcov_noorder[(pran + 1):p, (pran + 1):p]
    full_vcov[(pfix + 1):p, (pfix + 1):p] <- full_vcov_noorder[1:pran, 1:pran]
    full_vcov[(pfix + 1):p, 1:pfix] <- full_vcov_noorder[1:pran, (pran + 1):p]
    full_vcov[1:pfix, (pfix + 1):p] <- full_vcov_noorder[(pran + 1):p, 1:pran]


    ## reparameterize for sd and var for random variance/covariance parameters.
    if (ranpar == "theta") {
      full_vcov <- full_vcov
    }

    if (ranpar == "sd") {
      dd <- devfun2(object, useSc = FALSE, signames = TRUE)
      nvp <- length(attr(dd, "thopt"))
      pars <- attr(dd, "optimum")
      pars <- pars[!is.na(names(pars))]
      hh <- hessian(dd, pars) / (-2)

      full_vcov_noorder <- -solve(hh)
      full_vcov <- matrix(
        NA, nrow(full_vcov_noorder),
        ncol(full_vcov_noorder)
      )
      full_vcov[1:pfix, 1:pfix] <-
        full_vcov_noorder[(pran + 1):p, (pran + 1):p]
      full_vcov[(pfix + 1):p, (pfix + 1):p] <-
        full_vcov_noorder[1:pran, 1:pran]
      full_vcov[(pfix + 1):p, 1:pfix] <-
        full_vcov_noorder[1:pran, (pran + 1):p]
      full_vcov[1:pfix, (pfix + 1):p] <-
        full_vcov_noorder[(pran + 1):p, 1:pran]
    }

    if (ranpar == "var") {
      dd <- devfun2(object, useSc = FALSE, signames = TRUE)
      nvp <- length(attr(dd, "thopt"))
      pars <- attr(dd, "optimum")
      pars <- pars[!is.na(names(pars))]
      hh <- hessian(dd, pars) / (-2)

      sdcormat <- as.data.frame(VarCorr(object, comp = "Std.Dev"),
        order = "lower.tri"
      )
      idx <- is.na(sdcormat$var2)
      ## TODO: Couldn't we just replace this with 0.5 / sqrt(sdcormat$sdcor[which(dx)])?
      ## TODO: Also, I don't think we need which
      sdcormat$sdcor2[which(idx)] <- 0.5 * (sdcormat$sdcor[which(idx)])^(-0.5)
      ## TODO: this looks like - (x/y)^-1 which we could just write out as -y/x
      sdcormat$sdcor2[which(!idx)] <- (-1) * (sdcormat$vcov[which(!idx)] / sdcormat$sdcor[which(!idx)])^(-1)
      hh[((pran + 1):p), (1:pran)] <- sweep(as.matrix(hh[((pran + 1):p), (1:pran)]), MARGIN = 2, sdcormat$sdcor2, `*`)
      hh[(1:pran), ((pran + 1):p)] <- t(hh[((pran + 1):p), (1:pran)])
      ## ranhes reparameterization
      if (pran == 1) {
        entries <- matrix(1, 1, 1)
        weight <- (sdcormat$sdcor2)^2
      } else {
        entries <- rbind(matrix(rep(1:pran, each = 2), pran, 2, byrow = TRUE), t(combn(pran, 2)))
        entries <- entries[order(entries[, 1], entries[, 2]), ]
        weight <- apply(entries, 1, function(x) {
          sdcormat$sdcor2[x[1]] * sdcormat$sdcor2[x[2]]
        })
      }

      hh[1:pran, 1:pran][lower.tri(hh[1:pran, 1:pran], diag = TRUE)] <-
        weight * hh[1:pran, 1:pran][lower.tri(hh[1:pran, 1:pran], diag = TRUE)]
      if (pran > 1) {
        hh[1:pran, 1:pran] <- as.matrix(forceSymmetric(hh[1:pran, 1:pran], uplo = "L"))
      }

      full_vcov_noorder <- -solve(hh)
      full_vcov <- matrix(NA, nrow(full_vcov_noorder), ncol(full_vcov_noorder))
      full_vcov[1:pfix, 1:pfix] <- full_vcov_noorder[(pran + 1):p, (pran + 1):p]
      full_vcov[(pfix + 1):p, (pfix + 1):p] <- full_vcov_noorder[1:pran, 1:pran]
      full_vcov[(pfix + 1):p, 1:pfix] <- full_vcov_noorder[1:pran, (pran + 1):p]
      full_vcov[1:pfix, (pfix + 1):p] <- full_vcov_noorder[(pran + 1):p, 1:pran]
    }
    if (!(ranpar %in% c("sd", "theta", "var"))) {
      stop("ranpar needs to be sd, theta or var for glmerMod object.")
    }

    ## name the matrix
    parts <- getME(object, c("fixef", "theta"))
    colnames(full_vcov) <- c(names(parts$fixef), paste("cov", names(parts$theta), sep = "_"))
  }
  return(full_vcov)
}

lav_matrix_trace <- function(..., check = TRUE) {
  if (nargs() == 0L) {
    return(as.numeric(NA))
  }
  dots <- list(...)
  if (is.list(dots[[1]])) {
    mlist <- dots[[1]]
  } else {
    mlist <- dots
  }
  nMat <- length(mlist)
  if (nMat == 1L) {
    S <- mlist[[1]]
    if (check) {
      stopifnot(NROW(S) == NCOL(S))
    }
    out <- sum(S[lav_matrix_diag_idx(n = NROW(S))])
  } else if (nMat == 2L) {
    out <- sum(mlist[[1]] * t(mlist[[2]]))
  } else if (nMat == 3L) {
    A <- mlist[[1]]
    B <- mlist[[2]]
    C <- mlist[[3]]
    B2 <- B %*% C
    out <- sum(A * t(B2))
  } else {
    M1 <- mlist[[1]]
    M2 <- mlist[[2]]
    for (m in 3L:nMat) {
      M2 <- M2 %*% mlist[[m]]
    }
    out <- sum(M1 * t(M2))
  }
  out
}

lav_matrix_diag_idx <- function(n = 1L) {
  1L + (seq_len(n) - 1L) * (n + 1L)
}

########## R-function: ZOSull ##########
# For creation of O'Sullivan-type Z matrices.
# Last changed: 04 OCT 2021 by M.P.Wand.

ZOSull <- function(x, intKnots, range.x, drv = 0) {
  if (drv > 2) {
    stop("splines not smooth enough for more than 2 derivatives")
  }

  # Set defaults for `range.x' and `intKnots'
  if (missing(range.x)) {
    range.x <- c(1.05 * min(x) - 0.05 * max(x), 1.05 * max(x) - 0.05 * min(x))
  }

  if (missing(intKnots)) {
    numIntKnots <- min(length(unique(x)), 35)
    intKnots <- quantile(unique(x), seq(0, 1, length = (numIntKnots + 2))[-c(1, (numIntKnots + 2))])
  }
  numIntKnots <- length(intKnots)

  # Obtain the penalty matrix.
  allKnots <- c(rep(range.x[1], 4), intKnots, rep(range.x[2], 4))
  K <- length(intKnots)
  L <- 3 * (K + 8)
  xtilde <- (rep(allKnots, each = 3)[-c(1, (L - 1), L)] + rep(allKnots, each = 3)[-c(1, 2, L)]) / 2
  wts <- rep(diff(allKnots), each = 3) * rep(c(1, 4, 1) / 6, K + 7)
  Bdd <- splines::spline.des(allKnots, xtilde,
    derivs = rep(2, length(xtilde)),
    outer.ok = TRUE
  )$design
  Omega <- t(Bdd * wts) %*% Bdd

  # Use the spectral decomposition of Omega to obtain Z.
  eigOmega <- eigen(Omega)
  indsZ <- 1:(numIntKnots + 2)
  UZ <- eigOmega$vectors[, indsZ]
  LZ <- t(t(UZ) / sqrt(eigOmega$values[indsZ]))

  # Perform stability check.
  indsX <- (numIntKnots + 3):(numIntKnots + 4)
  UX <- eigOmega$vectors[, indsX]
  Lmat <- cbind(UX, LZ)
  stabCheck <- t(crossprod(Lmat, t(crossprod(Lmat, Omega))))
  if (sum(stabCheck^2) > 1.0001 * (numIntKnots + 2)) {
    print("WARNING: NUMERICAL INSTABILITY ARISING\\
              FROM SPECTRAL DECOMPOSITION")
  }

  # Obtain B and post-multiply by LZ matrix to get Z.
  B <- splines::spline.des(allKnots, x,
    derivs = rep(drv, length(x)),
    outer.ok = TRUE
  )$design

  Z <- crossprod(t(B), LZ)

  # Order the columns of Z with respect to the eigenvalues of "Omega":
  Z <- Z[, order(eigOmega$values[indsZ])]

  # Add the `range.x' and 'intKnots' as attributes
  # of the return object.
  attr(Z, "range.x") <- range.x
  attr(Z, "knots") <- intKnots

  # Return Z matrix with 2 attributes.
  return(Z)
}

############ R-function: Ztps ############
Ztps <- function(x, k, knots = NULL, range.x = NULL) {
  # Set up thin plate spline generalised
  # covariance function:
  tps.cov <- function(r, m = 2, d = 1) {
    r <- as.matrix(r)
    num.row <- nrow(r)
    num.col <- ncol(r)
    r <- as.vector(r)
    nzi <- seq_along(r)[r != 0]
    ans <- rep(0, length(r))
    if ((d + 1) %% 2 != 0) {
      ans[nzi] <- (abs(r[nzi]))^(2 * m - d) * log(abs(r[nzi])) # d is even
    } else {
      ans[nzi] <- (abs(r[nzi]))^(2 * m - d)
    }

    if (num.col > 1) {
      ans <- matrix(ans, num.row, num.col) # d is odd
    }
    return(ans)
  }

  # Set up function for matrix square-roots:
  matrix.sqrt <- function(A) {
    sva <- svd(A)
    if (min(sva$d) >= 0) {
      Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
    } else {
      stop("Matrix square root is not defined")
    }
    return(Asqrt)
  }

  if (is.null(knots)) {
    x1_grid <- seq(range(x[, 1])[1], range(x[, 1])[2], length = k)
    x2_grid <- seq(range(x[, 2])[1], range(x[, 2])[2], length = k)
    knots <- expand.grid(x1_grid, x2_grid)
    names(knots) <- colnames(x)
  }

  if (!is.null(range.x)) {
    inBdry <- HRW::pointsInPoly(knots, range.x)
    knots <- knots[inBdry, ]
  }

  # Obtain  matrix of inter-knot distances:
  numKnots <- nrow(knots)

  dist.mat <- matrix(0, numKnots, numKnots)
  dist.mat[lower.tri(dist.mat)] <- dist(as.matrix(knots))
  dist.mat <- dist.mat + t(dist.mat)

  Omega <- tps.cov(dist.mat, d = 2)

  # Obtain preliminary Z matrix of knot to data covariances:
  x.knot.diffs.1 <- outer(x[, 1], knots[, 1], "-")
  x.knot.diffs.2 <- outer(x[, 2], knots[, 2], "-")
  x.knot.dists <- sqrt(x.knot.diffs.1^2 + x.knot.diffs.2^2)

  prelim.Z <- tps.cov(x.knot.dists, m = 2, d = 2)

  # Transform to canonical form:
  sqrt.Omega <- matrix.sqrt(Omega)
  Z <- t(solve(sqrt.Omega, t(prelim.Z)))
  attr(Z, "knots") <- knots
  if (!is.null(range.x)) attr(Z, "range.x") <- range.x
  return(Z)
}
########################################################################
# ==========================================================================
#  https://www.r-bloggers.com/2021/08/r-dataframe-merge-while-keeping-orders-of-row-and-column/
# —————————————————————–
# Function : f_loj_krc
# —————————————————————–
# Left outer join while keeping orders of input rows and columns
# Meaning of input arguments are the same as those of merge()
# —————————————————————–
f_loj_krc <- function(x, y, by.x, by.y) {
  # save row id
  x.temp <- x
  x.temp$temp.id <- seq_len(nrow(x.temp))

  # each column names
  x.cn <- colnames(x)
  y.cn <- colnames(y)

  # replace column names of y with same names of x
  # to avoid duplicate fields
  for (i in seq_along(by.y)) {
    colnames(y)[which(y.cn == by.y[i])] <- by.x[i]
  }
  by.y <- by.x # since two fields are the same now

  # new column names of y
  y.cn <- colnames(y)

  # remove only joining key fields which are redundant
  # and keep only new informative fields
  y.cn.not.key <- setdiff(y.cn, by.y)

  # left outer join
  df <- merge(x = x.temp, y = y, by.x = by.x, by.y = by.y, all.x = TRUE)

  # recover the original rows and columns orders
  df <- df[order(df$temp.id), c(x.cn, y.cn.not.key)]
  rownames(df) <- NULL

  return(df)
}

########################################################
# To handle aovlist object refit by lmer function
aovlist_lmer <- function(object) {
  stopifnot(inherits(object, "aovlist"))
  mod_df <- model.frame(object)
  lmer_call <- match.call(aov, attr(object, "call"))

  trms <- terms(object)
  response <- as.character(attr(trms, "variables"))[[2]]

  # Find the Error terms
  trms_label <- attr(trms, "term.labels")
  err.idx <- grep("^Error\\(", trms_label)

  # Original Error term
  error_term <- trms_label[err.idx]
  error_term <- gsub("Error\\(|\\)", "", error_term)
  error_terms <- strsplit(error_term, "\\+\\s+")[[1]]

  # Reformat into lmer format
  lmer_random_parts <- paste0("(1|", error_terms, ")", collapse = " + ")
  lmer_fix_parts <- paste0(trms_label[-err.idx], collapse = " + ")

  lmer_call$formula <- as.formula(paste(response, " ~ ", lmer_fix_parts, "+", lmer_random_parts, sep = ""))
  lmer_call[[1]] <- as.name("lmer")
  assign(as.character(lmer_call$data), mod_df)
  lmer_object <- eval(lmer_call)
  return(lmer_object)
}

########################################################
reTrms_tmb <- function(model, ...) {
  form <- formula(model)
  cond_reTrms <- reformulas::mkReTrms(
    reformulas::findbars(form, ...),
    model.frame(model),
    reorder.terms = FALSE, calc.lambdat = TRUE
  )

  rmattr <- function(x, a = c("correlation", "blockCode", "stddev")) {
    for (aa in a) attr(x, aa) <- NULL
    x
  }

  mktheta <- function(model) {
    vc <- VarCorr(model)$cond
    get_chol <- function(v) {
      cc <- t(chol(rmattr(v)))
      cc[lower.tri(cc, diag = TRUE)]
    }
    theta <- lapply(vc, get_chol)
    return(unlist(theta) / sigma(model))
  }

  cond_reTrms$Lambdat@x <- unname(mktheta(model)[cond_reTrms$Lind])
  return(cond_reTrms)
}

########################################################
###################### for print
# Define print method for objects of class 'pdmlist'
#' @exportS3Method print pdmlist
print.pdmlist <- function(x, ...) {
  pos <- grep("predictmeansPlot|predictmeansciPlot|predictmeansBKPlot|predictmeansBarPlot|p_valueMatrix", names(x))
  if (length(pos) == 0) {
    x <- x[seq_along(x)]
  } else {
    x <- x[names(x)[-pos]]
  }
  NextMethod()
}

###################### for plot
# Define plot method for objects of class 'pdmlist'
#' @exportS3Method plot pdmlist
plot.pdmlist <- function(x, ...) {
  plotmt <- x$mean_table
  if (!is.data.frame(plotmt)) {
    plotmt <- plotmt[[1]]
  }

  if (!all(c("Mean", "SE") %in% names(plotmt))) {
    stop("Mean table must have columns named 'Mean' and 'SE'")
  }

  facts <- names(plotmt)[sapply(plotmt, is.factor)]
  if (length(facts) > 3) {
    stop("There is no plot for more than three-way interaction!")
  }
  names(facts) <- facts
  dots <- list(...)
  if ("plotord" %in% names(dots)) {
    plotord <- dots$plotord
    if (length(plotord) != length(facts)) {
      stop(paste("plotord must be a vector of length", length(facts)))
    }
    if (is.character(plotord)) {
      i <- match(plotord, facts)
      if (any(is.na(i))) {
        stop(paste(
          "The entries", paste0(plotord[is.na(i)], sep = ", "), "do not match any term in ",
          paste(facts, collapse = ", "), "\n"
        ))
      } else {
        plotord <- i
      }
    }
    facts <- facts[plotord]
    # dots$plotord <- NULL
  }

  facts <- lapply(facts, I)
  mpdots <- dots[intersect(names(dots), c("title", "xlab", "ylab", "scales", "basesz", "line"))]
  mpdots$plot_mt <- plotmt
  mpdots$x_var <- facts[1][[1]]
  mpdots$y_var <- "Mean"
  mpdots$col_var <- facts[2][[1]]
  mpdots$panel_var <- facts[3][[1]]
  if (!is.null(x$LSD)) {
    mpdots$bar_value <- x$LSD[3]
    mpdots$level <- attr(x$LSD, "Significant level")
  }

  meanPlot <- do.call(mean_plot, mpdots)

  # mean_plot(plotmt, facts[1][[1]], y_var="Mean", col_var = facts[2][[1]], panel_var = facts[3][[1]], ...)
  bpdots <- dots[intersect(names(dots), c("title", "xlab", "ylab", "scales", "basesz"))]
  bpdots$plot_mt <- plotmt
  bpdots$x_var <- facts[1][[1]]
  bpdots$y_var <- "Mean"
  bpdots$se_var <- "SE"
  bpdots$col_var <- facts[2][[1]]
  bpdots$panel_var <- facts[3][[1]]
  barPlot <- do.call(bar_plot, bpdots)

  ## barPlot <- bar_plot(plotmt, facts[1][[1]], y_var="Mean", se_var="SE", col_var = facts[2][[1]],
  ## panel_var = facts[3][[1]], ...)

  if ("mod_df" %in% names(dots)) {
    mod_df <- dots$mod_df
    dots$mod_df <- NULL
  } else {
    mod_df <- NULL
  }
  if ("resp_name" %in% names(dots)) {
    resp_name <- dots$resp_name
    dots$resp_name <- NULL
  } else {
    resp_name <- NULL
  }
  cidots <- dots[intersect(names(dots), c("jitterv", "basesz"))]
  cidots$plot_mt <- plotmt
  cidots$mod_df <- mod_df
  cidots$resp_name <- resp_name
  ciPlot <- do.call(ci_plot, cidots)

  # ciPlot <- ci_plot(plotmt, mod_df=mod_df, resp_name=resp_name, ...)
  if (!is.null(x$p_valueMatrix)) {
    pmPlot <- PMplot(x$p_valueMatrix, level = 0.05, legendx = 0.69)
  } else {
    pmPlot <- NULL
  }

  return(list(meanPlot = meanPlot, barPlot = barPlot, ciPlot = ciPlot, pmPlot = pmPlot))
}
