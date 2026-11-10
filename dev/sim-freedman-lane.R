## Simulation study: permutation schemes for testing one term when other terms
## have effects. Compares
##   raw   - the scheme of permmodels(): permute the response (lm), or for mixed
##           models permute the fitted values and decorrelated residuals of the
##           full model together (all effects are scrambled);
##   FL    - Freedman and Lane (1983): permute the residuals of the reduced model
##           (without the tested term), add them to its fitted values and refit
##           (for mixed models, residuals decorrelated with the reduced model's
##           covariance, as in permmodels(), then correlated again);
##   F     - the parametric F test (lm) or Satterthwaite F test (lmer).
## Type I error (tested term has no effect) and power, at the 5% level.
##
## Run: Rscript dev/sim-freedman-lane.R [R_lm] [R_lmer] [part]   (defaults 2000, 300, "all")
## part = "lmer-leverage" runs only the mixed-model design with a covariate
## with extreme values (dev/sim-freedman-lane-lmer-leverage.rds).
## Writes dev/sim-freedman-lane-results.rds and prints the tables.

args <- commandArgs(TRUE)
R_lm <- if (length(args) >= 1) as.integer(args[1]) else 2000
R_lmer <- if (length(args) >= 2) as.integer(args[2]) else 300
part <- if (length(args) >= 3) args[3] else "all"
nperm <- 199
alpha <- 0.05
set.seed(20261110)

## ---- errors ---------------------------------------------------------------
errs <- list(
  normal = function(n) rnorm(n),
  skewed = function(n) (rexp(n) - 1),                 # mean 0, sd 1
  heavy  = function(n) rt(n, df = 3) / sqrt(3))        # sd 1

## ---- lm designs -----------------------------------------------------------
## A: tested factor (3 levels); B: nuisance factor with a large effect.
make_design <- function(type) {
  if (type == "balanced") {
    d <- expand.grid(A = factor(1:3), B = factor(1:3), r = 1:4)
  } else if (type == "unbalanced") {
    # A and B associated: most units of A level i are in B level i
    counts <- matrix(c(8, 2, 2, 2, 8, 2, 2, 2, 8), 3)
    d <- do.call(rbind, lapply(1:3, function(i) do.call(rbind, lapply(1:3, function(j)
      data.frame(A = rep(i, counts[i, j]), B = rep(j, counts[i, j]))))))
    d$A <- factor(d$A); d$B <- factor(d$B)
  } else if (type == "covariate") {
    # covariate correlated with A
    d <- data.frame(A = factor(rep(1:3, each = 12)))
    d$x <- as.numeric(d$A) + rnorm(36, sd = 0.6)
  } else {
    # leverage: skewed covariate with extreme values, larger in A level 3
    # (the case where raw permutation is known to fail; Anderson & Legendre 1999)
    d <- data.frame(A = factor(rep(1:3, each = 8)))
    d$x <- exp(rnorm(24, sd = 1) + 0.8 * (d$A == "3"))
  }
  d
}
nuisance <- function(d, type) {
  if (type %in% c("covariate", "leverage")) 3 * (d$x - mean(d$x)) else c(0, 3, 6)[as.integer(d$B)]
}
red_formula <- function(type) if (type %in% c("covariate", "leverage")) ~ x else ~ B
full_formula <- function(type) if (type %in% c("covariate", "leverage")) ~ x + A else ~ B + A

## F statistics for many responses at once (columns of Y), from QR projections
f_stats <- function(X0, X1, Y) {
  q0 <- qr(X0); q1 <- qr(X1)
  r0 <- colSums(qr.resid(q0, Y)^2)
  r1 <- colSums(qr.resid(q1, Y)^2)
  df1 <- q1$rank - q0$rank; df2 <- nrow(X1) - q1$rank
  ((r0 - r1) / df1) / (r1 / df2)
}

lm_one <- function(d, y, type) {
  X0 <- model.matrix(red_formula(type), d)
  X1 <- model.matrix(full_formula(type), d)
  n <- length(y)
  idx <- replicate(nperm, sample.int(n))
  f_obs <- f_stats(X0, X1, matrix(y))
  df1 <- qr(X1)$rank - qr(X0)$rank; df2 <- n - qr(X1)$rank
  # raw: permute the response
  f_raw <- f_stats(X0, X1, matrix(y[idx], n))
  # Freedman-Lane: permute residuals of the reduced model
  q0 <- qr(X0); fit0 <- qr.fitted(q0, y); e0 <- y - fit0
  f_fl <- f_stats(X0, X1, fit0 + matrix(e0[idx], n))
  c(F = pf(f_obs, df1, df2, lower.tail = FALSE),
    raw = (sum(f_raw >= f_obs - 1e-12) + 1) / (nperm + 1),
    FL = (sum(f_fl >= f_obs - 1e-12) + 1) / (nperm + 1))
}

lm_results <- list()
if (part == "all") for (type in c("balanced", "unbalanced", "covariate", "leverage")) for (en in names(errs)) for (eff in c(0, 1)) {
  p <- t(replicate(R_lm, {
    d <- make_design(type)
    y <- nuisance(d, type) + eff * c(0, 0.5, 1)[as.integer(d$A)] + errs[[en]](nrow(d))
    lm_one(d, y, type)
  }))
  lm_results[[paste(type, en, eff)]] <- data.frame(model = "lm", design = type, errors = en,
    effect = if (eff == 0) "null (type I)" else "A effect (power)", t(colMeans(p <= alpha)))
  cat(sprintf("lm %-10s %-6s effect=%d: F %.3f raw %.3f FL %.3f\n", type, en, eff,
              mean(p[, "F"] <= alpha), mean(p[, "raw"] <= alpha), mean(p[, "FL"] <= alpha)))
}

## ---- lmer design ------------------------------------------------------------
## 8 blocks (random), A tested (3 levels), B nuisance (3 levels, large effect),
## A and B associated within blocks; y ~ B + A + (1 | block).
suppressPackageStartupMessages(library(lme4))
lmer_design <- function() {
  cells <- rbind(cbind(1, c(1, 1, 1, 2)), cbind(2, c(2, 2, 2, 3)), cbind(3, c(3, 3, 3, 1)))
  d <- do.call(rbind, lapply(1:8, function(b) data.frame(block = b, A = cells[, 1], B = cells[, 2])))
  d$block <- factor(d$block); d$A <- factor(d$A); d$B <- factor(d$B)
  d
}
chol_lower <- function(m) {          # lower Cholesky factor of V = Z L L' Z' + I
  ZL <- getME(m, "Z") %*% getME(m, "Lambda")
  V <- Matrix::forceSymmetric(Matrix::tcrossprod(ZL) + Matrix::Diagonal(nrow(ZL)))
  Matrix::t(Matrix::chol(V))
}
lmer_one <- function(d, full_f = y ~ B + A + (1 | block), red_f = y ~ B + (1 | block)) {
  ctl <- lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
  full <- lmer(full_f, data = d, REML = TRUE, control = ctl)
  red <- lmer(red_f, data = d, REML = TRUE, control = ctl)
  f_of <- function(m) anova(m)["A", "F value"]          # sequential: A after B
  f_obs <- f_of(full)
  p_F <- suppressMessages(lmerTest::as_lmerModLmerTest(full))
  p_F <- anova(p_F)["A", "Pr(>F)"]
  n <- nrow(d)
  # raw (permmodels): permute fitted values and decorrelated residuals together
  xb <- as.vector(getME(full, "X") %*% fixef(full)); U <- chol_lower(full)
  w <- as.vector(Matrix::solve(U, d$y - xb))
  # Freedman-Lane: reduced model's fitted values stay; its decorrelated residuals are permuted
  xb0 <- as.vector(getME(red, "X") %*% fixef(red)); U0 <- chol_lower(red)
  w0 <- as.vector(Matrix::solve(U0, d$y - xb0))
  f_raw <- f_fl <- numeric(nperm)
  for (k in seq_len(nperm)) {
    i <- sample.int(n)
    f_raw[k] <- tryCatch(f_of(refit(full, xb[i] + as.vector(U %*% w[i]))), error = function(e) NA)
    f_fl[k] <- tryCatch(f_of(refit(full, xb0 + as.vector(U0 %*% w0[i]))), error = function(e) NA)
  }
  c(F = p_F,
    raw = (sum(f_raw >= f_obs - 1e-10, na.rm = TRUE) + 1) / (sum(!is.na(f_raw)) + 1),
    FL = (sum(f_fl >= f_obs - 1e-10, na.rm = TRUE) + 1) / (sum(!is.na(f_fl)) + 1))
}
nperm <- 99
lmer_results <- list()
if (part == "all") for (en in c("normal", "skewed")) for (eff in if (en == "normal") c(0, 1) else 0) {
  p <- t(replicate(R_lmer, {
    d <- lmer_design()
    d$y <- c(0, 3, 6)[as.integer(d$B)] + eff * c(0, 0.6, 1.2)[as.integer(d$A)] +
      rnorm(8, sd = 1)[as.integer(d$block)] + errs[[en]](nrow(d))
    suppressWarnings(suppressMessages(lmer_one(d)))
  }))
  lmer_results[[paste(en, eff)]] <- data.frame(model = "lmer", design = "blocks, A~B associated", errors = en,
    effect = if (eff == 0) "null (type I)" else "A effect (power)", t(colMeans(p <= alpha)))
  cat(sprintf("lmer %-6s effect=%d: F %.3f raw %.3f FL %.3f\n", en, eff,
              mean(p[, "F"] <= alpha), mean(p[, "raw"] <= alpha), mean(p[, "FL"] <= alpha)))
}

## lmer with a covariate with extreme values: 8 blocks, each with two units of
## each A level; x log-normal, larger for A level 3, with a large slope;
## y ~ x + A + (1 | block), A tested.
lmer_lev_design <- function() {
  d <- expand.grid(rep = 1:2, A = factor(1:3), block = factor(1:8))
  d$x <- exp(rnorm(nrow(d), sd = 1) + 0.8 * (d$A == "3"))
  d
}
if (part == "lmer-leverage") {
  lev_results <- list()
  for (en in names(errs)) for (eff in if (en == "normal") c(0, 1) else 0) {
    p <- t(replicate(R_lmer, {
      d <- lmer_lev_design()
      d$y <- 3 * (d$x - mean(d$x)) + eff * c(0, 0.5, 1)[as.integer(d$A)] +
        rnorm(8, sd = 1)[as.integer(d$block)] + errs[[en]](nrow(d))
      suppressWarnings(suppressMessages(lmer_one(d, y ~ x + A + (1 | block), y ~ x + (1 | block))))
    }))
    lev_results[[paste(en, eff)]] <- data.frame(model = "lmer", design = "blocks, leverage covariate", errors = en,
      effect = if (eff == 0) "null (type I)" else "A effect (power)", t(colMeans(p <= alpha)), R = R_lmer)
    cat(sprintf("lmer leverage %-6s effect=%d: F %.3f raw %.3f FL %.3f\n", en, eff,
                mean(p[, "F"] <= alpha), mean(p[, "raw"] <= alpha), mean(p[, "FL"] <= alpha)))
  }
  lev <- do.call(rbind, lev_results); rownames(lev) <- NULL
  saveRDS(lev, "dev/sim-freedman-lane-lmer-leverage.rds")
  print(lev, digits = 3)
  quit(save = "no")
}

res <- rbind(do.call(rbind, lm_results), do.call(rbind, lmer_results))
rownames(res) <- NULL
res$R <- ifelse(res$model == "lm", R_lm, R_lmer)
saveRDS(res, "dev/sim-freedman-lane-results.rds")
print(res, digits = 3)
