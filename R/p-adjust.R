## Multiplicity adjustments beyond p.adjust() and Tukey, for one family of
## contrasts, as emmeans' "sidak", "scheffe" and "mvt":
##   sidak    1 - (1 - p)^m for m contrasts;
##   scheffe  F = t^2/r on (r, df) degrees of freedom, r the rank of the
##            family (k - 1 for the pairs among k means);
##   mvt      single-step multivariate t, P(max_j |T_j| >= |t_i|) with T
##            multivariate t with the contrasts' correlation R and df. This is
##            what SAS ADJUST=SIMULATE estimates by simulation; here it is
##            computed by mvtnorm::pmvt() (error below about 1e-4). Each
##            contrast uses its own df, rounded as emmeans does.
## Contrasts with NA t or p are left out of the family and stay NA.
pm_adjust_extra <- function(adj, t, p, df, rank = NULL, R = NULL) {
  out <- rep(NA_real_, length(t))
  ok <- !is.na(t) & !is.na(p)
  m <- sum(ok)
  if (m == 0) return(out)
  df <- rep_len(as.numeric(df), length(t))
  out[ok] <- switch(adj,
    sidak = 1 - (1 - p[ok])^m,
    scheffe = stats::pf(t[ok]^2/rank, rank, df[ok], lower.tail = FALSE),
    mvt = pm_p_mvt(t[ok], df[ok], R[ok, ok, drop = FALSE]),
    stop("unknown adjustment '", adj, "'", call. = FALSE))
  out
}

pm_adjust_extra_methods <- c("sidak", "scheffe", "mvt")

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
## agree with the adjusted p-values: a pair's adjusted p-value is below alpha
## exactly when its limits exclude 0 (as SAS AdjLower/AdjUpper and emmeans).
## 'fam' labels the family of each contrast (NA: in no family, limits NA);
## 'est' marks the estimable contrasts (the others stay NA and, as in the
## p-values, do not count in m for "bonferroni" and "sidak").
##   tukey       qtukey(1 - alpha, k, df)/sqrt(2), k means, each pair's df
##   bonferroni  qt(1 - alpha/(2m), df)
##   sidak       qt(1 - (1 - (1 - alpha)^(1/m))/2, df)
##   scheffe     sqrt(r qf(1 - alpha, r, df)), r = k - 1
##   mvt         the two-sided (1 - alpha) quantile of max |T| (mvtnorm::qmvt)
## With a permlist ('T': permuted |t|, one row per contrast): "tukey", "mvt"
## and "westfall" the (1 - alpha) quantile of the family's permuted max |t|
## (single-step; Westfall-Young's step-down has no simultaneous limits),
## "bonferroni" and "sidak" each contrast's permutation quantile at the
## per-contrast level. NULL for the methods with no simultaneous limits
## (holm, hochberg, hommel and the false discovery rate methods).
pm_adjust_crit <- function(adj, alpha, df, fam, est, R = NULL, T = NULL) {
  perm <- !is.null(T)
  have <- if (perm) c("tukey", "mvt", "westfall", "bonferroni", "sidak") else
    c("tukey", "bonferroni", "sidak", "scheffe", "mvt")
  if (!adj %in% have) return(NULL)
  n_all <- length(fam)
  df <- rep_len(as.numeric(df), n_all)
  out <- rep(NA_real_, n_all)
  for (f in unique(stats::na.omit(fam))) {
    idx <- which(fam %in% f)
    ok <- idx[est[idx]]
    if (!length(ok)) next
    k <- (1 + sqrt(1 + 8*length(idx)))/2        # means in the family, from its k(k-1)/2 pairs
    m <- length(ok)
    a1 <- switch(adj, bonferroni = alpha/m, sidak = 1 - (1 - alpha)^(1/m), alpha)
    out[ok] <- if (perm) {
      if (adj %in% c("bonferroni", "sidak")) pm_perm_crit(T[ok, , drop = FALSE], a1) else
        pm_perm_crit(T[ok, , drop = FALSE], alpha, family = rep("f", m))
    } else switch(adj,
      tukey = stats::qtukey(1 - alpha, k, df[ok])/sqrt(2),
      bonferroni = , sidak = stats::qt(1 - a1/2, df[ok]),
      scheffe = sqrt((k - 1)*stats::qf(1 - alpha, k - 1, df[ok])),
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
