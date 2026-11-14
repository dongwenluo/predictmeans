## Multiplicity adjustments, in one place for predictmeans() (the pairs table,
## the printed p-value matrices and the letters), contrastmeans(),
## pm_slopes() and the means of semireg fits.
##
## A family is the set of contrasts adjusted together: all the pairs, or with
## 'atvar' the pairs within each of its levels. pm_families() gives each
## family's contrasts and its size for each method; pm_adjust_p() the adjusted
## p-values and pm_adjust_crit() the simultaneous limits that go with them.
## Contrasts that cannot be estimated (NA) are left out of their family and
## stay NA.
##   tukey      Tukey-Kramer: |t| sqrt(2) against the studentized range of k
##              means (those taking part in the family's estimable pairs),
##              each pair with its own df (SAS ADJDFE=ROW)
##   sidak      1 - (1 - p)^m for m contrasts
##   scheffe    F = t^2/r on (r, df) df, r the rank of the family's contrasts
##              (as emmeans; k - 1 for the pairs among k independent means)
##   mvt        single-step multivariate t, P(max_j |T_j| >= |t_i|) with T
##              multivariate t with the contrasts' correlation R and df. This
##              is what SAS ADJUST=SIMULATE estimates by simulation; here it is
##              computed by mvtnorm::pmvt() (error below about 1e-4). Each
##              contrast uses its own df, rounded as emmeans does.
##   with permutation p-values ('perm_T': observed and permuted |t|), "tukey",
##              "mvt" and "westfall" use the family's permuted max |t|
##              (single step; "westfall" step down)
##   others     stats::p.adjust() (holm, hochberg, hommel, bonferroni, BH, BY)

pm_adjust_extra_methods <- c("sidak", "scheffe", "mvt")

## The families: 'fam' labels the family of each contrast (NA: in none),
## 'est' marks the estimable ones. Sizes: m (estimable contrasts), k (means
## in the estimable pairs, from 'ij', the two means of each pair; without
## 'ij', from the number of pairs k(k-1)/2) and r (rank of the rows of 'L',
## the estimable contrasts in coefficient space; or 'rank'; else k - 1).
pm_families <- function(fam, est, ij = NULL, L = NULL, rank = NULL) {
  lev <- unique(fam[!is.na(fam)])
  lapply(stats::setNames(lev, lev), function(f) {
    idx <- which(fam %in% f)
    ok <- idx[est[idx]]
    k <- if (!is.null(ij)) length(unique(c(ij[ok, 1], ij[ok, 2]))) else (1 + sqrt(1 + 8*length(idx)))/2
    r <- if (!is.null(rank)) rank else if (!is.null(L) && length(ok)) qr(as.matrix(L[ok, , drop = FALSE]))$rank else k - 1
    list(idx = idx, ok = ok, m = length(ok), k = k, r = max(r, 1))
  })
}

## Adjusted p-values of contrasts with t statistics 't', p-values 'p' and df
## 'df', within the families 'fams' (pm_families()). R: correlation of the
## contrasts ("mvt"); perm_T: list(obs, T) of permutation p-values.
pm_adjust_p <- function(adj, t, p, df, fams, R = NULL, perm_T = NULL) {
  t <- as.vector(t)
  p <- as.vector(p)
  out <- rep(NA_real_, length(p))
  df <- rep_len(as.numeric(df), length(p))
  for (f in fams) {
    ok <- f$ok
    if (!length(ok)) next
    out[ok] <- if (adj == "none" || f$m < 2) {
      p[ok]
    } else if (!is.null(perm_T) && adj %in% pm_perm_maxt_methods) {
      pm_perm_maxT(perm_T$obs[ok], perm_T$T[ok, , drop = FALSE], stepdown = adj == "westfall", p_raw = p[ok])
    } else switch(adj,
      tukey = stats::ptukey(sqrt(2)*abs(t[ok]), f$k, df[ok], lower.tail = FALSE),
      sidak = 1 - (1 - p[ok])^f$m,
      scheffe = stats::pf(t[ok]^2/f$r, f$r, df[ok], lower.tail = FALSE),
      mvt = pm_p_mvt(t[ok], df[ok], R[ok, ok, drop = FALSE]),
      stats::p.adjust(p[ok], adj))
  }
  out
}

pm_p_mvt <- function(t, df, R) {
  pm_need("mvtnorm", "adj = \"mvt\"")
  if (length(t) > 1000) stop("adj = \"mvt\" handles at most 1000 contrasts in a family.", call. = FALSE)
  R <- stats::cov2cor(as.matrix(R))
  k <- length(t)
  # pmvt() is a Monte Carlo integral: a fixed seed makes results repeat
  alg <- mvtnorm::GenzBretz(maxpts = 5e4, abseps = 1e-4)
  pm_with_seed(20240403, vapply(seq_len(k), function(i) {
    a <- abs(t[i])
    pr <- mvtnorm::pmvt(lower = rep(-a, k), upper = rep(a, k), df = pm_mvt_df(df[i]), corr = R, algorithm = alg)
    min(1, max(0, 1 - as.numeric(pr)))
  }, numeric(1)))
}

## Critical values c of simultaneous confidence limits (estimate +- c SE) that
## agree with the adjusted p-values: a contrast's adjusted p-value is below
## alpha exactly when its limits exclude 0 (as SAS AdjLower/AdjUpper and
## emmeans), with the family sizes of pm_families():
##   tukey       qtukey(1 - alpha, k, df)/sqrt(2), each pair's df
##   bonferroni  qt(1 - alpha/(2m), df)
##   sidak       qt(1 - (1 - (1 - alpha)^(1/m))/2, df)
##   scheffe     sqrt(r qf(1 - alpha, r, df))
##   mvt         the two-sided (1 - alpha) quantile of max |T| (mvtnorm::qmvt)
## With a permlist ('T': permuted |t|, one row per contrast): "tukey", "mvt"
## and "westfall" the (1 - alpha) quantile of the family's permuted max |t|
## (single-step; Westfall-Young's step-down has no simultaneous limits),
## "bonferroni" and "sidak" each contrast's permutation quantile at the
## per-contrast level. NULL for the methods with no simultaneous limits
## (holm, hochberg, hommel and the false discovery rate methods).
pm_adjust_crit <- function(adj, alpha, df, fams, R = NULL, T = NULL) {
  perm <- !is.null(T)
  have <- if (perm) c("tukey", "mvt", "westfall", "bonferroni", "sidak") else
    c("tukey", "bonferroni", "sidak", "scheffe", "mvt")
  if (!adj %in% have) return(NULL)
  n_all <- max(c(0, unlist(lapply(fams, `[[`, "idx"))), length(df))
  df <- rep_len(as.numeric(df), n_all)
  out <- rep(NA_real_, n_all)
  for (f in fams) {
    ok <- f$ok
    m <- f$m
    if (!m) next
    a1 <- switch(adj, bonferroni = alpha/m, sidak = 1 - (1 - alpha)^(1/m), alpha)
    out[ok] <- if (perm) {
      if (adj %in% c("bonferroni", "sidak")) pm_perm_crit(T[ok, , drop = FALSE], a1) else
        pm_perm_crit(T[ok, , drop = FALSE], alpha, family = rep("f", m))
    } else switch(adj,
      tukey = stats::qtukey(1 - alpha, max(f$k, 2), df[ok])/sqrt(2),
      bonferroni = , sidak = stats::qt(1 - a1/2, df[ok]),
      scheffe = sqrt(f$r*stats::qf(1 - alpha, f$r, df[ok])),
      mvt = pm_q_mvt(alpha, df[ok], R[ok, ok, drop = FALSE]))
  }
  out
}

## mvtnorm's integer df, rounded as emmeans does (0: normal)
pm_mvt_df <- function(d) {
  if (d > 0) d <- max(1, d)
  if (is.infinite(d) || d > 9999) d <- 0
  as.integer(floor(d + 0.25))
}

## Two-sided multivariate t quantiles for each contrast's df (rounded by
## pm_mvt_df()); mvtnorm::qmvt() with a fixed seed.
pm_q_mvt <- function(alpha, df, R) {
  pm_need("mvtnorm", "adj = \"mvt\"")
  if (length(df) > 1000) stop("adj = \"mvt\" handles at most 1000 contrasts in a family.", call. = FALSE)
  R <- stats::cov2cor(as.matrix(R))
  d <- vapply(df, pm_mvt_df, integer(1))
  alg <- mvtnorm::GenzBretz(maxpts = 5e4, abseps = 1e-4)
  q <- vapply(unique(d), function(dd) pm_with_seed(20240403,
    mvtnorm::qmvt(1 - alpha, tail = "both.tails", df = dd, corr = R, algorithm = alg)$quantile), numeric(1))
  q[match(d, unique(d))]
}

## Evaluate 'expr' with a fixed random seed and leave the user's random
## number stream as it was.
pm_with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

## Critical value for simultaneous bands over the rows of D (Ruppert, Wand
## and Carroll): the (1 - level) quantile of max_g |D_g theta| / se_g for
## theta ~ N(0, C), from 10000 draws.
pm_sim_crit <- function(D, C, level, nsim = 10000) {
  D <- as.matrix(D)
  C <- as.matrix(C)
  se <- sqrt(rowSums((D %*% C) * D))
  ok <- is.finite(se) & se > 0
  if (!any(ok)) return(stats::qnorm(1 - level/2))
  e <- eigen(C, symmetric = TRUE)
  A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow = length(e$values))
  DA <- D[ok, , drop = FALSE] %*% A
  M <- pm_with_seed(20240403, {
    Z <- matrix(stats::rnorm(ncol(A)*nsim), ncol(A))
    apply(abs(DA %*% Z)/se[ok], 2, max)
  })
  as.numeric(stats::quantile(M, 1 - level))
}
