## Follow-up 2: does whitening with the symmetric square root of V (instead
## of its Cholesky factor) make the residual-permutation schemes exact under
## the LMM global null? The Cholesky factor treats the first unit of each
## block differently from the others; the symmetric root treats the units of
## a block alike (Commenges 2003 on exchangeability-preserving transforms).
## Same design and seeds rule as fu_lmm_global (dev/sim-perm-followup.R),
## different cell names, so different data sets.
##   joint_chol / fl_chol  as in dev/sim-perm-followup.R (Cholesky)
##   joint_sym  / fl_sym   the same with V^(1/2) symmetric (eigen)
##   within, satt          references
## Run: Rscript dev/sim-perm-followup2.R [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
nperm <- if (length(args) >= 1) as.integer(args[1]) else 99
cores <- if (length(args) >= 3) as.integer(args[3]) else 1
source("dev/sim-perm-lib.R")

sym_root <- function(V) {
  e <- eigen(as.matrix(V), symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
}

lmm_global_sym <- function(d, f0, f1, term) {
  m1 <- lme4::lmer(f1, data = d, REML = TRUE, control = ctl)
  m0 <- lme4::lmer(f0, data = d, REML = TRUE, control = ctl)
  Fs <- function(y) stats::anova(refit(m1, y))[term, "F value"]
  obs <- stats::anova(m1)[term, "F value"]
  n <- nrow(d)
  parts <- function(m, sym) {
    ZL <- getME(m, "Z") %*% getME(m, "Lambda")
    V <- Matrix::tcrossprod(ZL) + Matrix::Diagonal(n)
    U <- if (sym) sym_root(V) else as.matrix(lower_chol(ZL, rep(1, n)))
    xb <- as.vector(getME(m, "X") %*% fixef(m))
    list(U = U, xb = xb, e = as.vector(solve(U, d$y - xb)))
  }
  P <- list(joint_chol = parts(m1, FALSE), fl_chol = parts(m0, FALSE),
            joint_sym = parts(m1, TRUE), fl_sym = parts(m0, TRUE))
  blocks <- split(seq_len(n), d$block)
  s <- matrix(NA_real_, nperm, 5, dimnames = list(NULL, c(names(P), "within")))
  for (k in seq_len(nperm)) {
    i <- sample.int(n)
    w <- seq_len(n); for (b in blocks) w[b] <- b[sample.int(length(b))]
    for (nm in names(P)) {
      p <- P[[nm]]
      ys <- if (startsWith(nm, "joint")) p$xb[i] + as.vector(p$U %*% p$e[i]) else p$xb + as.vector(p$U %*% p$e[i])
      s[k, nm] <- tryCatch(Fs(ys), error = function(e) NA)
    }
    s[k, "within"] <- tryCatch(Fs(d$y[w]), error = function(e) NA)
  }
  satt <- stats::anova(lmerTest::as_lmerModLmerTest(m1))[term, "Pr(>F)"]
  c(apply(s, 2, pval_perm, obs = obs), satt = satt)
}

description <- list(
  fu2_lmm_global_sym = "LMM global null as fx_global (8 blocks); Cholesky vs symmetric-root whitening; 2000 data sets",
  fu2_lmm_lev_sym = "LMM, extreme covariate related to A (as sim-freedman-lane lmer leverage, 8 blocks), A null; Cholesky vs symmetric root; 1000 data sets")
R_of <- c(fu2_lmm_global_sym = 2000, fu2_lmm_lev_sym = 1000)

one_set <- function(cell) {
  switch(cell,
    fu2_lmm_global_sym = {
      d <- lev_design(8); d$x <- rnorm(nrow(d))
      d$y <- rnorm(8)[d$block] + rnorm(nrow(d))
      lmm_global_sym(d, y ~ x + (1 | block), y ~ x + A + (1 | block), "A")
    },
    fu2_lmm_lev_sym = {
      d <- lev_design(8)
      d$y <- 3 * (d$x - mean(d$x)) + rnorm(8)[d$block] + rnorm(nrow(d))
      lmm_global_sym(d, y ~ x + (1 | block), y ~ x + A + (1 | block), "A")
    },
    stop("unknown cell ", cell))
}

cells <- if (length(args) >= 2 && args[2] != "all") strsplit(args[2], ",")[[1]] else names(description)
for (cell in cells) run_cells(cell, one_set, R_of[[cell]], nperm, cores, "dev/sim-perm-followup", description)
cat("SIM DONE\n")
