## Simulation check of permutation tests for random effects (and one fixed-
## effect case), as implemented in permlmer() / permanova.lmer(), and a pilot
## of the proposed GLMM extension ("working-scale" permutation).
##
## Cells (each: type I error at 5% under H0, or power under H1):
##   lmm_sp_norm   split plot, test whole-plot variance, normal errors
##   lmm_sp_skew   as above, skewed errors and skewed block effects
##   lmm_sp_unbal  as above, normal errors, 20% of observations missing
##   lmm_slope     repeated measures, test a random slope (H0: random intercept)
##   lmm_sp_power  split plot, whole-plot SD 0.5 (power)
##   lmm_fix_lev   fixed term A with a covariate with extreme values, LRT by
##                 permanova.lmer's scheme (permute null fitted values with the
##                 whitened residuals) and by Freedman-Lane (keep the fitted values)
##   glmm_bin1     binary GLMM split plot, test whole-plot variance (working scale)
##   glmm_bin10    binomial (10 trials) GLMM split plot, same test
##   glmm_pois     Poisson GLMM split plot, same test
##   glmm_pois_power  Poisson, whole-plot SD 0.4 (power)
##   glmm_fix_lev_bin10, glmm_fix_lev_pois  fixed factor A with a covariate with
##                 extreme values in a GLMM (binomial 10 trials; Poisson), A null
##   glmm_fix_lev_pois_power  as above, A effect (power)
##   glmm_fix_assoc_bin1  binary GLMM, A tested, B (large effect) associated with A
## Methods:
##   perm      permlmer() as in the package (lmer cells)
##   chibar    REML LRT with the 50:50 chi-square(0):chi-square(1) reference
##   work_perm GLMM: linearise the H0 glmer fit (working response z, weights w),
##             whiten the marginal working residuals with V0 = Z G Z' + diag(1/w),
##             permute, re-correlate, refit weighted LMMs to z*, REML LRT
##   glmm_chibar  GLMM Laplace LRT with the 50:50 reference
##   joint     GLMM fixed term: the maintainer's scheme on the working scale
##             (full-model fitted values and whitened residuals permuted together)
##   fl        GLMM fixed term: Freedman-Lane on the working scale (reduced model)
##   lrt_chisq GLMM fixed term: Laplace LRT with the chi-square reference
##
## Run: Rscript dev/sim-perm-ranef.R [R] [nperm] [cells] [cores]
##   e.g. Rscript dev/sim-perm-ranef.R 1000 99 lmm_sp_norm,glmm_bin1 4
## Each cell writes dev/sim-perm-ranef/<cell>.rds and prints one line.

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
nperm <- if (length(args) >= 2) as.integer(args[2]) else 99
cells <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else
  c("lmm_sp_norm", "lmm_sp_skew", "lmm_sp_unbal", "lmm_slope", "lmm_sp_power", "lmm_fix_lev",
    "glmm_bin1", "glmm_bin10", "glmm_pois", "glmm_pois_power")
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
alpha <- 0.05
dir.create("dev/sim-perm-ranef", showWarnings = FALSE)
suppressPackageStartupMessages({ library(lme4); library(predictmeans) })
ctl <- lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
gctl <- glmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")

## ---- designs ---------------------------------------------------------------
## split plot: 6 blocks x 4 whole plots (trt A) x 3 subplots (trt B)
sp_design <- function(nsub = 3) {
  d <- expand.grid(sub = seq_len(nsub), wp = 1:4, block = 1:6)
  d$A <- factor(d$wp); d$B <- factor(d$sub)
  d$block <- factor(d$block); d$wpid <- interaction(d$block, d$wp)
  d
}
skew <- function(n) rexp(n) - 1

pval_perm <- function(stat, obs) (sum(stat >= obs - 1e-10, na.rm = TRUE) + 1) / (sum(!is.na(stat)) + 1)
chibar <- function(lrt) 0.5 * stats::pchisq(max(lrt, 0), 1, lower.tail = FALSE)

## permlmer() p-value and the chi-bar reference for one data set
lmm_test <- function(d, f0, f1) {
  m0 <- lmer(f0, data = d, REML = TRUE, control = ctl)
  m1 <- lmer(f1, data = d, REML = TRUE, control = ctl)
  pp <- suppressMessages(permlmer(m0, m1, nperm = nperm, ncore = 1))$`Perm-p`[2]
  lrt <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  c(perm = pp, chibar = chibar(lrt))
}

## fixed term by LRT (ML criterion at REML estimates, as permlmer does when the
## fixed effects differ): "joint" = permanova.lmer's scheme (null fitted values
## and whitened residuals permuted together); "fl" = fitted values kept
lmm_fixed <- function(d, f0, f1) {
  m0 <- lmer(f0, data = d, REML = TRUE, control = ctl)
  m1 <- lmer(f1, data = d, REML = TRUE, control = ctl)
  ZL <- getME(m0, "Z") %*% getME(m0, "Lambda")
  V <- Matrix::forceSymmetric(Matrix::tcrossprod(ZL) + Matrix::Diagonal(nrow(d)))
  U <- Matrix::t(Matrix::chol(V))
  xb <- as.vector(getME(m0, "X") %*% fixef(m0))
  w <- as.vector(Matrix::solve(U, d$y - xb))
  lr <- function(y) 2 * (as.numeric(logLik(refit(m1, y), REML = FALSE)) - as.numeric(logLik(refit(m0, y), REML = FALSE)))
  obs <- lr(d$y)
  sj <- sf <- numeric(nperm)
  for (k in seq_len(nperm)) {
    i <- sample.int(nrow(d))
    sj[k] <- tryCatch(lr(xb[i] + as.vector(U %*% w[i])), error = function(e) NA)
    sf[k] <- tryCatch(lr(xb + as.vector(U %*% w[i])), error = function(e) NA)
  }
  mf1 <- lmer(f1, data = d, REML = FALSE, control = ctl); mf0 <- lmer(f0, data = d, REML = FALSE, control = ctl)
  c(joint = pval_perm(sj, obs), fl = pval_perm(sf, obs),
    lrt_chisq = stats::pchisq(2 * (as.numeric(logLik(mf1)) - as.numeric(logLik(mf0))),
                              df = length(fixef(mf1)) - length(fixef(mf0)), lower.tail = FALSE))
}

## working-scale permutation test of a variance component in a GLMM
glmm_test <- function(d, f0, f1, family) {
  g0 <- glmer(f0, data = d, family = family, control = gctl)
  g1 <- glmer(f1, data = d, family = family, control = gctl)
  lrt_glmm <- 2 * (as.numeric(logLik(g1)) - as.numeric(logLik(g0)))
  ## linearise at the H0 fit
  eta <- g0@resp$eta - g0@resp$offset
  mu <- g0@resp$mu
  fam <- stats::family(g0)
  w <- g0@resp$sqrtWrkWt()^2                       # prior weights x mu.eta^2 / variance
  z <- eta + (g0@resp$y - mu) / fam$mu.eta(eta)
  xb <- as.vector(getME(g0, "X") %*% fixef(g0))
  ZL <- getME(g0, "Z") %*% getME(g0, "Lambda")      # G = Lambda Lambda' (no dispersion)
  V <- Matrix::forceSymmetric(Matrix::tcrossprod(ZL) + Matrix::Diagonal(x = 1 / w))
  U <- Matrix::t(Matrix::chol(V))
  e <- as.vector(Matrix::solve(U, z - xb))
  ## working LMMs (weights w), fitted once to z, then refitted
  dz <- d; dz$z <- z; dz$w <- w
  wf0 <- stats::update(f0, z ~ .); wf1 <- stats::update(f1, z ~ .)
  l0 <- lmer(wf0, data = dz, weights = w, REML = TRUE, control = ctl)
  l1 <- lmer(wf1, data = dz, weights = w, REML = TRUE, control = ctl)
  lr <- function(y) max(2 * (as.numeric(logLik(refit(l1, y))) - as.numeric(logLik(refit(l0, y)))), 0)
  obs <- lr(z)
  s <- vapply(seq_len(nperm), function(k) {
    i <- sample.int(nrow(d))
    tryCatch(lr(xb + as.vector(U %*% e[i])), error = function(e) NA)
  }, numeric(1))
  c(work_perm = pval_perm(s, obs), work_chibar = chibar(obs), glmm_chibar = chibar(lrt_glmm))
}

## linearise a glmer fit: working response z, working weights w, fitted
## fixed part xb, lower Cholesky factor U of V = Z G Z' + diag(1/w), and the
## whitened marginal working residuals e
linearise <- function(g) {
  eta <- g@resp$eta - g@resp$offset
  fam <- stats::family(g)
  w <- g@resp$sqrtWrkWt()^2
  z <- eta + (g@resp$y - g@resp$mu) / fam$mu.eta(eta)
  xb <- as.vector(getME(g, "X") %*% fixef(g))
  ZL <- getME(g, "Z") %*% getME(g, "Lambda")
  U <- Matrix::t(Matrix::chol(Matrix::forceSymmetric(Matrix::tcrossprod(ZL) + Matrix::Diagonal(x = 1 / w))))
  list(z = z, w = w, xb = xb, U = U, e = as.vector(Matrix::solve(U, z - xb)))
}

## fixed term of a GLMM on the working scale. "joint" = the maintainer's
## scheme (as permmodels() for lmer): linearise the full model, permute its
## fitted values and whitened residuals together. "fl" = Freedman-Lane:
## linearise the reduced model, keep its fitted values, permute its whitened
## residuals. Statistic: F of the term from a weighted LMM of the working
## response (term last, sequential), observed and permuted alike.
glmm_fixed <- function(d, f0, f1, family, term) {
  g0 <- glmer(f0, data = d, family = family, control = gctl)
  g1 <- glmer(f1, data = d, family = family, control = gctl)
  lrt <- 2 * (as.numeric(logLik(g1)) - as.numeric(logLik(g0)))
  wf1 <- stats::update(f1, z ~ .)
  Fstat <- function(m) stats::anova(m)[term, "F value"]
  setup <- function(L) {
    dz <- d; dz$z <- L$z; dz$w <- L$w
    lme4::lmer(wf1, data = dz, weights = w, REML = TRUE, control = ctl)
  }
  L1 <- linearise(g1); L0 <- linearise(g0)
  t1 <- setup(L1); t0 <- setup(L0)
  obs_j <- Fstat(t1); obs_f <- Fstat(t0)
  sj <- sf <- numeric(nperm)
  for (k in seq_len(nperm)) {
    i <- sample.int(nrow(d))
    sj[k] <- tryCatch(Fstat(refit(t1, L1$xb[i] + as.vector(L1$U %*% L1$e[i]))), error = function(e) NA)
    sf[k] <- tryCatch(Fstat(refit(t0, L0$xb + as.vector(L0$U %*% L0$e[i]))), error = function(e) NA)
  }
  c(joint = pval_perm(sj, obs_j), fl = pval_perm(sf, obs_f),
    lrt_chisq = stats::pchisq(lrt, df = length(fixef(g1)) - length(fixef(g0)), lower.tail = FALSE))
}

## ---- cells -----------------------------------------------------------------
one_set <- function(cell) {
  switch(cell,
    lmm_sp_norm = , lmm_sp_skew = , lmm_sp_unbal = , lmm_sp_power = {
      d <- sp_design()
      eb <- if (cell == "lmm_sp_skew") skew(6) else rnorm(6)
      ee <- if (cell == "lmm_sp_skew") skew(nrow(d)) else rnorm(nrow(d))
      sw <- if (cell == "lmm_sp_power") 0.5 else 0
      d$y <- 0.5 * as.numeric(d$A) + 0.3 * as.numeric(d$B) + eb[d$block] + sw * rnorm(24)[d$wpid] + ee
      if (cell == "lmm_sp_unbal") d <- d[sort(sample.int(nrow(d), round(0.8 * nrow(d)))), ]
      lmm_test(d, y ~ A * B + (1 | block), y ~ A * B + (1 | block) + (1 | wpid))
    },
    lmm_slope = {
      d <- expand.grid(time = 0:5, subj = factor(1:12))
      d$y <- 1 + 0.5 * d$time + rnorm(12)[d$subj] + rnorm(nrow(d))
      lmm_test(d, y ~ time + (1 | subj), y ~ time + (1 | subj) + (0 + time | subj))
    },
    lmm_fix_lev = {
      d <- expand.grid(rep = 1:2, A = factor(1:3), block = factor(1:8))
      d$x <- exp(rnorm(nrow(d)) + 0.8 * (d$A == "3"))
      d$y <- 3 * (d$x - mean(d$x)) + rnorm(8)[d$block] + rnorm(nrow(d))
      lmm_fixed(d, y ~ x + (1 | block), y ~ x + A + (1 | block))
    },
    glmm_fix_lev_bin10 = , glmm_fix_lev_pois = , glmm_fix_lev_pois_power = {
      d <- expand.grid(rep = 1:2, A = factor(1:3), block = factor(1:8))
      d$x <- exp(rnorm(nrow(d)) + 0.8 * (d$A == "3"))
      xc <- d$x - mean(d$x)
      bl <- 0.5 * rnorm(8)[d$block]
      if (cell == "glmm_fix_lev_bin10") {
        d$n <- 10; d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.3 * xc + bl))
        glmm_fixed(d, cbind(r, n - r) ~ x + (1 | block), cbind(r, n - r) ~ x + A + (1 | block), binomial, "A")
      } else {
        eff <- if (cell == "glmm_fix_lev_pois_power") c(0, 0.25, 0.5)[d$A] else 0
        d$y <- rpois(nrow(d), exp(1 + 0.15 * xc + eff + bl))
        glmm_fixed(d, y ~ x + (1 | block), y ~ x + A + (1 | block), poisson, "A")
      }
    },
    glmm_fix_assoc_bin1 = {
      cl <- rbind(cbind(1, c(1, 1, 1, 2)), cbind(2, c(2, 2, 2, 3)), cbind(3, c(3, 3, 3, 1)))
      d <- do.call(rbind, lapply(1:8, function(b) data.frame(block = b, A = cl[, 1], B = cl[, 2])))
      d <- d[rep(seq_len(nrow(d)), each = 2), ]
      d$block <- factor(d$block); d$A <- factor(d$A); d$B <- factor(d$B)
      d$y <- rbinom(nrow(d), 1, plogis(c(-1, 0, 1)[d$B] + 0.7 * rnorm(8)[d$block]))
      glmm_fixed(d, y ~ B + (1 | block), y ~ B + A + (1 | block), binomial, "A")
    },
    glmm_bin1 = , glmm_bin10 = , glmm_pois = , glmm_pois_power = {
      nsub <- if (cell == "glmm_bin1") 5 else 3
      d <- sp_design(nsub)
      sw <- if (cell == "glmm_pois_power") 0.4 else 0
      eta <- -0.3 + 0.3 * (as.numeric(d$A) - 2.5) + 0.2 * (as.numeric(d$B) - 2) + 0.7 * rnorm(6)[d$block] +
        sw * rnorm(24)[d$wpid]
      if (cell %in% c("glmm_pois", "glmm_pois_power")) {
        d$y <- rpois(nrow(d), exp(eta + 1))
        glmm_test(d, y ~ A + B + (1 | block), y ~ A + B + (1 | block) + (1 | wpid), poisson)
      } else {
        n <- if (cell == "glmm_bin1") 1 else 10
        d$n <- n; d$r <- rbinom(nrow(d), n, plogis(eta))
        glmm_test(d, cbind(r, n - r) ~ A + B + (1 | block), cbind(r, n - r) ~ A + B + (1 | block) + (1 | wpid), binomial)
      }
    })
}

for (cell in cells) {
  t0 <- Sys.time()
  base <- sum(utf8ToInt(cell)) * 1000
  res <- parallel::mclapply(seq_len(R), function(i) {
    set.seed(base + i)
    tryCatch(suppressWarnings(suppressMessages(one_set(cell))), error = function(e) NULL)
  }, mc.cores = cores)
  ok <- !vapply(res, is.null, logical(1))
  p <- do.call(rbind, res[ok])
  out <- list(cell = cell, R = R, R_ok = sum(ok), nperm = nperm, p = p, reject = colMeans(p <= alpha),
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, file.path("dev/sim-perm-ranef", paste0(cell, ".rds")))
  cat(sprintf("CELL %-15s ok %4d/%d  %s  (%.1f min)\n", cell, sum(ok), R,
              paste(sprintf("%s %.3f", names(out$reject), out$reject), collapse = "  "), out$minutes))
}
cat("SIM DONE\n")
