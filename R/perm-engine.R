## Shared parts of the permutation and bootstrap tests (permmodels(),
## permlmer(), permanova.lmer(), bootmodels(), and permutation p-values of
## pairs and contrasts in predictmeans() and contrastmeans()):
##   pm_perm_stream()  permutation indices, with a seed of their own;
##   pm_perm_map()     refits, in parallel if asked, failures kept as NULL;
##   pm_perm_drop()    drops failed refits and counts them;
##   pm_perm_p()       p-values from rounded statistics;
##   pm_perm_tstats()  |t| of contrasts from a list of refits (permlist);
##   pm_perm_fl_t()    Freedman-Lane p-values of contrasts (lm/aov), from the
##                     permutation object of permmodels(method = "freedman-lane").
## See dev/permutation-plan.md, section 8.

## Permutation indices. pm_perm_stream(n, seed) returns a function that gives
## the next m permutations of 1:n as an n x m integer matrix; successive calls
## continue the same stream, so drawing in chunks gives the same permutations
## as drawing all at once. With a seed (or a saved random number 'state') the
## stream has its own random number state and the user's stream is left as it
## was; without one it draws from the user's stream (as set.seed() before the
## call would).
pm_perm_stream <- function(n, seed = NULL, state = NULL) {
  genv <- globalenv()
  shared <- is.null(seed) && is.null(state)
  if (!is.null(seed)) state <- pm_with_seed(seed, get(".Random.seed", envir = genv, inherits = FALSE))
  if (shared && !exists(".Random.seed", envir = genv, inherits = FALSE)) set.seed(NULL)
  start <- if (shared) get(".Random.seed", envir = genv, inherits = FALSE) else state
  draw <- function(m) matrix(vapply(seq_len(m), function(i) sample.int(n), integer(n)), nrow = n)
  f <- function(m) {
    if (shared) return(draw(m))
    old <- get0(".Random.seed", envir = genv, inherits = FALSE)
    on.exit(if (is.null(old)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old, envir = genv))
    assign(".Random.seed", state, envir = genv)
    out <- draw(m)
    state <<- get(".Random.seed", envir = genv, inherits = FALSE)
    out
  }
  # the random number state before the first draw: pm_perm_stream(n, state =
  # attr(f, "start")) gives the same permutations again
  attr(f, "start") <- start
  f
}

## Permutations from a matrix given by the user (e.g. permindex()), or from
## pm_perm_stream(): a function giving the next m permutations. 'perms' must
## have one row per observation and permutations of 1:n as columns.
pm_perm_source <- function(n, seed = NULL, perms = NULL, state = NULL) {
  if (is.null(perms)) return(pm_perm_stream(n, seed, state))
  perms <- pm_perm_check_index(perms, n)
  pos <- 0L
  f <- function(m) {
    out <- perms[, pos + seq_len(m), drop = FALSE]
    pos <<- pos + m
    out
  }
  attr(f, "matrix") <- perms
  f
}

## Check a matrix of permutation indices: one row per observation used in
## the model, each column a permutation of 1:n.
pm_perm_check_index <- function(perms, n) {
  perms <- as.matrix(perms)
  if (nrow(perms) != n) {
    stop("'perms' must have one row per observation used in the model (", n, "), not ", nrow(perms), ".",
         call. = FALSE)
  }
  if (!is.numeric(perms) || any(perms != round(perms)) ||
      any(apply(perms, 2, function(p) anyDuplicated(p) > 0 || min(p) < 1 || max(p) > n))) {
    stop("each column of 'perms' must be a permutation of 1 to ", n, ".", call. = FALSE)
  }
  storage.mode(perms) <- "integer"
  perms
}

## Apply FUN to each element of X (usually one permuted or simulated data set
## each), on 'ncore' cores: forked processes (mclapply) on Unix, a socket
## cluster on Windows with 'packages' attached on each worker. An error in
## FUN gives NULL for that element, which pm_perm_drop() counts as a failure.
pm_perm_map <- function(X, FUN, ncore = 1, packages = character()) {
  safe <- function(x) tryCatch(FUN(x), error = function(e) NULL)
  ncore <- max(1L, as.integer(ncore))
  if (ncore == 1L || length(X) < 2L) return(lapply(X, safe))
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    for (p in packages) parallel::clusterCall(cl, library, p, character.only = TRUE)
    parallel::parLapplyLB(cl, X, safe)
  } else {
    parallel::mclapply(X, safe, mc.cores = ncore)
  }
}

## Drop failed refits (NULL, or a try-error from a forked process). Returns
## the remaining results and the number dropped.
pm_perm_drop <- function(res) {
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  list(res = res[!bad], failures = sum(bad))
}

## Permutation (or bootstrap) p-values. 'obs' holds k observed statistics and
## 'sims' the statistics from the permuted data sets (a k x B matrix, or a
## vector of length B when k = 1). Both are rounded to 'digits' decimals
## before they are compared, so numerical noise in statistics that are equal
## (an LRT of 0, the same F twice) does not decide the comparison. A
## p-value is (number of simulated >= observed + 1) / (number of simulated
## + 1), the observed data counting as one of the permutations, or without
## the +1 terms when 'exact' is TRUE (all permutations enumerated).
## Simulated statistics that are NA are left out; an NA observed statistic
## gives an NA p-value.
pm_perm_p <- function(obs, sims, exact = FALSE, digits = 6) {
  sims <- matrix(sims, nrow = length(obs))
  ok <- !is.na(sims)
  ge <- round(sims, digits) >= round(obs, digits)
  ge[!ok] <- FALSE
  p <- (rowSums(ge) + !exact) / (rowSums(ok) + !exact)
  p[is.na(obs)] <- NA
  p
}

## |t| of the contrasts in the rows of rK for each refit in a permlist (a list
## of lists with $coef and $vcov): a k x B matrix.
pm_perm_tstats <- function(refits, rK) {
  rK <- as.matrix(rK)
  matrix(vapply(refits, function(x) {
    vc <- rK %*% tcrossprod(x$vcov, rK)
    abs(as.vector(rK %*% x$coef) / sqrt(diag(vc)))
  }, numeric(nrow(rK))), nrow = nrow(rK))
}

## Number of refits that failed, as text for the notes under printed tables
## ("" when none failed).
pm_perm_failnote <- function(failures, what = "refits") {
  if (failures > 0) paste0(" (", failures, " ", what, " failed and were left out)") else ""
}

## Was a permlist given (predictmeans(), contrastmeans())? Old scripts pass
## "" or "NULL" for none.
pm_perm_given <- function(permlist) {
  if (is.null(permlist)) return(FALSE)
  if (is.character(permlist)) return(!all(permlist %in% c("NULL", "")))
  TRUE
}

## Number of permutations behind a permlist (refits, or a permutation object).
pm_perm_n <- function(permlist) {
  if (!is.null(permlist[["perm"]])) permlist[["perm"]]$nperm else length(permlist[[1]])
}

## Permutation object of permmodels(method = "freedman-lane"): what is needed
## to test any linear hypothesis later (pairs, contrasts) with the same
## permutations: the (weighted) design matrix X without aliased columns, the
## response z (offset removed, weighted), the random number state that
## regenerates the permutations, and names to check that a later call uses
## the same model.
pm_perm_fl_object <- function(X, z, nperm, exact, start, coef_names, idx = NULL) {
  structure(list(scheme = "freedman-lane", X = X, z = z, nperm = as.integer(nperm), exact = exact,
                 start = start, idx = idx, coef_names = coef_names),
            class = "pm_perm")
}

## Freedman-Lane p-values of the contrasts in the rows of L (columns: the
## model's non-aliased coefficients), one test per row: the model with that
## contrast set to 0 (X N, N a basis of the null space of the row) gives
## fitted values and residuals; the residuals are permuted, added back, and
## |t| of the row is computed again by least squares. The statistic is the
## t value of the fitted model (as predictmeans() shows it), so the
## p-value tests that comparison alone, whatever the other means are.
pm_perm_fl_t <- function(perm, L, model = NULL) {
  X <- perm$X
  z <- perm$z
  L <- matrix(as.numeric(as.matrix(L)), ncol = ncol(as.matrix(L)))
  if (ncol(L) != ncol(X) || (!is.null(model) && !identical(pm_perm_coef_names(model), perm$coef_names))) {
    stop("'permlist' was made for a different model.", call. = FALSE)
  }
  n <- nrow(X)
  qrX <- qr(X)
  df_res <- n - qrX$rank
  XtXi <- solve(crossprod(X))
  se1 <- sqrt(rowSums((L %*% XtXi) * L))
  abs_t <- function(j, Y) {
    B <- qr.coef(qrX, Y)
    s <- sqrt(colSums(qr.resid(qrX, Y)^2) / df_res)
    abs(as.vector(L[j, , drop = FALSE] %*% B)) / (se1[j] * s)
  }
  k <- nrow(L)
  red <- lapply(seq_len(k), function(j) {
    q <- qr(t(L[j, , drop = FALSE]))
    N <- qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
    f0 <- if (ncol(N)) qr.fitted(qr(X %*% N), z) else rep(0, n)
    list(f0 = f0, e0 = z - f0, obs = round(abs_t(j, matrix(z)), 6))
  })
  count <- valid <- numeric(k)
  perms <- pm_perm_source(n, perms = perm$idx, state = perm$start)
  chunk <- max(1L, min(perm$nperm, floor(2e6 / n)))
  done <- 0L
  while (done < perm$nperm) {
    m <- min(chunk, perm$nperm - done)
    idx <- perms(m)
    for (j in seq_len(k)) {
      tj <- round(abs_t(j, red[[j]]$f0 + matrix(red[[j]]$e0[idx], n)), 6)
      ok <- !is.na(tj)
      count[j] <- count[j] + sum(tj[ok] >= red[[j]]$obs)
      valid[j] <- valid[j] + sum(ok)
    }
    done <- done + m
  }
  p <- (count + !perm$exact) / (valid + !perm$exact)
  p[se1 == 0 | !is.finite(se1)] <- NA
  p
}

## Names of the non-aliased coefficients of an lm/aov model.
pm_perm_coef_names <- function(model) {
  b <- stats::coef(model)
  names(b)[!is.na(b)]
}

## With a permlist, permutation p-values compare the observed t statistics
## with those of the refits, which use the unadjusted covariance; the
## observed ones must use it too, so the Kenward-Roger adjustment
## (vcov_adjust) is not used. Returns the vcov_adjust to use.
pm_perm_vcov_adjust <- function(model, permlist, vcov_adjust, vcov_adjust_default) {
  if (!pm_perm_given(permlist) || !isTRUE(vcov_adjust)) return(vcov_adjust)
  if (!vcov_adjust_default && inherits(model, "lmerMod")) {
    warning("'vcov_adjust = TRUE' is not used with 'permlist': permutation p-values compare t statistics ",
            "from the unadjusted covariance, as the refits do.", call. = FALSE)
  }
  FALSE
}
