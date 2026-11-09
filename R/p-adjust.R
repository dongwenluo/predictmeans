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
  fix_df <- function(d) {
    if (d > 0) d <- max(1, d)
    if (is.infinite(d) || d > 9999) d <- 0          # 0: normal
    as.integer(floor(d + 0.25))
  }
  R <- stats::cov2cor(as.matrix(R))
  k <- length(t)
  # pmvt() is a Monte Carlo integral: a fixed seed makes results repeat
  alg <- mvtnorm::GenzBretz(maxpts = 5e4, abseps = 1e-4)
  pm_with_seed(20240403, vapply(seq_len(k), function(i) {
    a <- abs(t[i])
    pr <- mvtnorm::pmvt(lower = rep(-a, k), upper = rep(a, k), df = fix_df(df[i]), corr = R, algorithm = alg)
    min(1, max(0, 1 - as.numeric(pr)))
  }, numeric(1)))
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
