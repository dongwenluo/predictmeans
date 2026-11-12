## Shared parts of the permutation and bootstrap tests (permmodels(),
## permlmer(), permanova.lmer(), bootmodels(), and permutation p-values of
## pairs and contrasts in predictmeans() and contrastmeans()):
##   pm_perm_stream()  permutation indices, with a seed of their own;
##   pm_perm_map()     refits, in parallel if asked, failures kept as NULL;
##   pm_perm_drop()    drops failed refits and counts them;
##   pm_perm_p()       p-values from rounded statistics;
##   pm_perm_tstats()  |t| of contrasts from a list of refits (permlist).
## See dev/permutation-plan.md, section 8.

## Permutation indices. pm_perm_stream(n, seed) returns a function that gives
## the next m permutations of 1:n as an n x m integer matrix; successive calls
## continue the same stream, so drawing in chunks gives the same permutations
## as drawing all at once. With a seed the stream has its own random number
## state and the user's stream is left as it was; without one it draws from
## the user's stream (as set.seed() before the call would).
pm_perm_stream <- function(n, seed = NULL) {
  state <- if (is.null(seed)) NULL else
    pm_with_seed(seed, get(".Random.seed", envir = globalenv(), inherits = FALSE))
  draw <- function(m) matrix(vapply(seq_len(m), function(i) sample.int(n), integer(n)), nrow = n)
  function(m) {
    if (is.null(state)) return(draw(m))
    genv <- globalenv()
    old <- get0(".Random.seed", envir = genv, inherits = FALSE)
    on.exit(if (is.null(old)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old, envir = genv))
    assign(".Random.seed", state, envir = genv)
    out <- draw(m)
    state <<- get(".Random.seed", envir = genv, inherits = FALSE)
    out
  }
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
