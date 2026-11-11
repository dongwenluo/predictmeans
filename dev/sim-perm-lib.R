## Shared functions for the permutation simulation studies
## (dev/sim-perm-edge.R; dev/sim-perm-ranef.R has its own copy so that a run in
## progress is never affected by edits here).
##
## Conventions (keep them, the research log in dev/sim-study/ relies on them):
## - each cell is one design x truth x model; data set i of a cell uses
##   set.seed(sum(utf8ToInt(cell)) * 1000 + i), so results do not depend on
##   the number of cores or the platform;
## - p-values use the +1 rule: (#{T* >= T} + 1) / (#valid + 1);
## - every p-value of every data set is saved, so any level can be studied later.

suppressPackageStartupMessages({ library(lme4); library(predictmeans) })
ctl <- lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
gctl <- glmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")

## ---- designs and distributions ---------------------------------------------
## split plot: nblock blocks x 4 whole plots (factor A) x nsub subplots (factor B)
sp_design <- function(nsub = 3, nblock = 6) {
  d <- expand.grid(sub = seq_len(nsub), wp = 1:4, block = seq_len(nblock))
  d$A <- factor(d$wp); d$B <- factor(d$sub)
  d$block <- factor(d$block); d$wpid <- interaction(d$block, d$wp, drop = TRUE)
  d
}
skew <- function(n) rexp(n) - 1                  # mean 0, sd 1
heavy <- function(n) rt(n, df = 3) / sqrt(3)       # sd 1
lev_design <- function(nblock = 8) {               # covariate with extreme values, larger for A3
  d <- expand.grid(rep = 1:2, A = factor(1:3), block = factor(seq_len(nblock)))
  d$x <- exp(rnorm(nrow(d)) + 0.8 * (d$A == "3"))
  d
}
assoc_design <- function(nblock = 8, reps = 1) {   # A and B associated within blocks
  cl <- rbind(cbind(1, c(1, 1, 1, 2)), cbind(2, c(2, 2, 2, 3)), cbind(3, c(3, 3, 3, 1)))
  d <- do.call(rbind, lapply(seq_len(nblock), function(b) data.frame(block = b, A = cl[, 1], B = cl[, 2])))
  d <- d[rep(seq_len(nrow(d)), each = reps), ]
  d$block <- factor(d$block); d$A <- factor(d$A); d$B <- factor(d$B)
  d
}

pval_perm <- function(stat, obs) (sum(stat >= obs - 1e-10, na.rm = TRUE) + 1) / (sum(!is.na(stat)) + 1)
## 50:50 mixture of chi-square(k-1) and chi-square(k) (k = 1: chi-square(0) and (1))
chibar <- function(lrt, k = 1) {
  lrt <- max(lrt, 0)
  0.5 * stats::pchisq(lrt, k, lower.tail = FALSE) + 0.5 * (if (k == 1) 0 else stats::pchisq(lrt, k - 1, lower.tail = FALSE))
}
lower_chol <- function(ZL, d) Matrix::t(Matrix::chol(Matrix::forceSymmetric(Matrix::tcrossprod(ZL) + Matrix::Diagonal(x = d))))

## ---- random effects, linear mixed models ----------------------------------
## permlmer() (Lee-Braun) and the chi-bar reference
lmm_test <- function(d, f0, f1, k = 1) {
  m0 <- lmer(f0, data = d, REML = TRUE, control = ctl)
  m1 <- lmer(f1, data = d, REML = TRUE, control = ctl)
  pp <- suppressMessages(permlmer(m0, m1, nperm = nperm, ncore = 1))$`Perm-p`[2]
  lrt <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  c(perm = pp, chibar = chibar(lrt, k))
}

## ---- fixed effects, linear mixed models -----------------------------------
## permmodels() form: F of `term` (sequential, term last) from refits of the
## full model. "joint" = the maintainer's scheme (full-model fitted values and
## whitened residuals permuted together); "fl" = Freedman-Lane (reduced model's
## fitted values kept, its whitened residuals permuted); "satt" = Satterthwaite F.
lmm_fixed_F <- function(d, f0, f1, term) {
  m1 <- lme4::lmer(f1, data = d, REML = TRUE, control = ctl)
  m0 <- lme4::lmer(f0, data = d, REML = TRUE, control = ctl)
  Fs <- function(m) stats::anova(m)[term, "F value"]
  obs <- Fs(m1)
  U1 <- lower_chol(getME(m1, "Z") %*% getME(m1, "Lambda"), rep(1, nrow(d)))
  xb1 <- as.vector(getME(m1, "X") %*% fixef(m1)); e1 <- as.vector(Matrix::solve(U1, d$y - xb1))
  U0 <- lower_chol(getME(m0, "Z") %*% getME(m0, "Lambda"), rep(1, nrow(d)))
  xb0 <- as.vector(getME(m0, "X") %*% fixef(m0)); e0 <- as.vector(Matrix::solve(U0, d$y - xb0))
  sj <- sf <- numeric(nperm)
  for (k in seq_len(nperm)) {
    i <- sample.int(nrow(d))
    sj[k] <- tryCatch(Fs(refit(m1, xb1[i] + as.vector(U1 %*% e1[i]))), error = function(e) NA)
    sf[k] <- tryCatch(Fs(refit(m1, xb0 + as.vector(U0 %*% e0[i]))), error = function(e) NA)
  }
  satt <- stats::anova(lmerTest::as_lmerModLmerTest(m1))[term, "Pr(>F)"]
  c(joint = pval_perm(sj, obs), fl = pval_perm(sf, obs), satt = satt)
}

## ---- GLMMs on the working scale -------------------------------------------
## linearise a glmer fit: working response z, weights w, fitted fixed part xb,
## lower Cholesky factor U of V = Z G Z' + diag(1/w), whitened residuals e
linearise <- function(g) {
  eta <- g@resp$eta - g@resp$offset
  fam <- stats::family(g)
  w <- g@resp$sqrtWrkWt()^2
  z <- eta + (g@resp$y - g@resp$mu) / fam$mu.eta(eta)
  xb <- as.vector(getME(g, "X") %*% fixef(g))
  U <- lower_chol(getME(g, "Z") %*% getME(g, "Lambda"), 1 / w)
  list(z = z, w = w, xb = xb, U = U, e = as.vector(Matrix::solve(U, z - xb)))
}

## random effects: Lee-Braun on the working scale (linearised at H0)
glmm_test <- function(d, f0, f1, family, k = 1) {
  g0 <- glmer(f0, data = d, family = family, control = gctl)
  g1 <- glmer(f1, data = d, family = family, control = gctl)
  lrt_glmm <- 2 * (as.numeric(logLik(g1)) - as.numeric(logLik(g0)))
  L <- linearise(g0)
  dz <- d; dz$z <- L$z; dz$w <- L$w
  l0 <- lme4::lmer(stats::update(f0, z ~ .), data = dz, weights = w, REML = TRUE, control = ctl)
  l1 <- lme4::lmer(stats::update(f1, z ~ .), data = dz, weights = w, REML = TRUE, control = ctl)
  lr <- function(y) max(2 * (as.numeric(logLik(refit(l1, y))) - as.numeric(logLik(refit(l0, y)))), 0)
  obs <- lr(L$z)
  s <- vapply(seq_len(nperm), function(j) {
    i <- sample.int(nrow(d))
    tryCatch(lr(L$xb + as.vector(L$U %*% L$e[i])), error = function(e) NA)
  }, numeric(1))
  c(work_perm = pval_perm(s, obs), work_chibar = chibar(obs, k), glmm_chibar = chibar(lrt_glmm, k))
}

## fixed effects: the maintainer's scheme ("joint", linearised at the full fit)
## and Freedman-Lane ("fl", linearised at the reduced fit); F of `term` from a
## weighted LMM of the working response; reference: Laplace LRT, chi-square
glmm_fixed <- function(d, f0, f1, family, term) {
  g0 <- glmer(f0, data = d, family = family, control = gctl)
  g1 <- glmer(f1, data = d, family = family, control = gctl)
  lrt <- 2 * (as.numeric(logLik(g1)) - as.numeric(logLik(g0)))
  wf1 <- stats::update(f1, z ~ .)
  Fs <- function(m) stats::anova(m)[term, "F value"]
  setup <- function(L) {
    dz <- d; dz$z <- L$z; dz$w <- L$w
    lme4::lmer(wf1, data = dz, weights = w, REML = TRUE, control = ctl)
  }
  L1 <- linearise(g1); L0 <- linearise(g0)
  t1 <- setup(L1); t0 <- setup(L0)
  obs_j <- Fs(t1); obs_f <- Fs(t0)
  sj <- sf <- numeric(nperm)
  for (k in seq_len(nperm)) {
    i <- sample.int(nrow(d))
    sj[k] <- tryCatch(Fs(refit(t1, L1$xb[i] + as.vector(L1$U %*% L1$e[i]))), error = function(e) NA)
    sf[k] <- tryCatch(Fs(refit(t0, L0$xb + as.vector(L0$U %*% L0$e[i]))), error = function(e) NA)
  }
  c(joint = pval_perm(sj, obs_j), fl = pval_perm(sf, obs_f),
    lrt_chisq = stats::pchisq(lrt, df = length(fixef(g1)) - length(fixef(g0)), lower.tail = FALSE))
}

## ---- runner ------------------------------------------------------------------
## Runs each cell (skipping cells whose result file exists unless overwrite),
## saves dev/<outdir>/<cell>.rds with all p-values, settings and versions, and
## prints one line per cell. Forks on Unix; uses a socket cluster on Windows.
run_cells <- function(cells, one_set, R, nperm, cores, outdir, description = list(), overwrite = FALSE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (cell in cells) {
    f <- file.path(outdir, paste0(cell, ".rds"))
    if (file.exists(f) && !overwrite) { cat("CELL", cell, "exists, skipped\n"); next }
    t0 <- Sys.time()
    base <- sum(utf8ToInt(cell)) * 1000
    task <- function(i) {
      set.seed(base + i)
      tryCatch(suppressWarnings(suppressMessages(one_set(cell))), error = function(e) conditionMessage(e))
    }
    res <- if (cores > 1 && .Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, ls(globalenv()), envir = globalenv())
      parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(lme4); library(predictmeans) }))
      r <- parallel::parLapply(cl, seq_len(R), task)
      parallel::stopCluster(cl)
      r
    } else parallel::mclapply(seq_len(R), task, mc.cores = cores)
    ok <- vapply(res, is.numeric, logical(1))
    p <- do.call(rbind, res[ok])
    errors <- table(unlist(res[!ok]))
    rej <- colMeans(p <= 0.05)
    out <- list(cell = cell, description = description[[cell]], R = R, R_ok = sum(ok), nperm = nperm,
                p = p, reject = rej, errors = errors,
                seed_rule = "set.seed(sum(utf8ToInt(cell)) * 1000 + i)",
                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), date = format(Sys.time()),
                cores = cores, platform = R.version$platform,
                versions = c(R = R.version.string, lme4 = as.character(utils::packageVersion("lme4")),
                             Matrix = as.character(utils::packageVersion("Matrix")),
                             predictmeans = as.character(utils::packageVersion("predictmeans"))))
    saveRDS(out, f)
    cat(sprintf("CELL %-20s ok %4d/%d  %s  (%.1f min)\n", cell, sum(ok), R,
                paste(sprintf("%s %.3f", names(rej), rej), collapse = "  "), out$minutes))
  }
  invisible(NULL)
}
