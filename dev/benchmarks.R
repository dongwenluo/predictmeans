# Timings for Kenward-Roger and Satterthwaite df and predictmeans().
# Usage: Rscript --no-init-file dev/benchmarks.R [library]
# Compares with pbkrtest and lmerTest when they are installed.

lib <- commandArgs(TRUE)[1]
suppressPackageStartupMessages({
  if (is.na(lib)) library(predictmeans) else library(predictmeans, lib.loc = lib)
})
tm <- function(expr, reps = 3) {
  e <- substitute(expr); pf <- parent.frame()
  min(vapply(seq_len(reps), function(i) system.time(eval(e, pf))[["elapsed"]], numeric(1)))
}
sim <- function(n_block, n_rep = 2) {
  set.seed(1)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), Block = factor(seq_len(n_block)), rep = seq_len(n_rep))
  d$WP <- interaction(d$Block, d$A)
  d$y <- as.numeric(d$A) + as.numeric(d$B) + rnorm(n_block)[d$Block] + rnorm(nlevels(d$WP))[d$WP] + rnorm(nrow(d))
  d
}
has <- function(p) requireNamespace(p, quietly = TRUE)
rows <- list()
for (nb in c(10, 50, 250)) {
  d <- sim(nb)
  m <- lme4::lmer(y ~ A*B + (1|Block) + (1|WP), data = d)
  L <- Kmatrix(m, "A:B")$K
  r <- data.frame(n = nrow(d),
    kr_vcov = tm(predictmeans:::pm_vcov_kr(m)),
    kr_df = tm(df_term(m, ctrmatrix = L, type = "Kenward-Roger")),
    satt_df = tm(df_term(m, ctrmatrix = L, type = "Satterthwaite")),
    predictmeans = tm(predictmeans(m, "A:B", adj = "BH", plot = FALSE, prtnum = FALSE)),
    containment = tm(predictmeans(m, "A:B", df_method = "containment", plot = FALSE, prtnum = FALSE)),
    pbkrtest_vcov = if (has("pbkrtest") && nrow(d) <= 3000) tm(pbkrtest::vcovAdj(m), 1) else NA,
    lmerTest_satt = if (has("lmerTest")) tm({mt <- lmerTest::as_lmerModLmerTest(m); for (i in seq_len(nrow(L))) lmerTest::contest(mt, L[i, ], ddf = "Satterthwaite")}, 1) else NA)
  rows[[length(rows) + 1]] <- r
  print(r, row.names = FALSE)
}
res <- do.call(rbind, rows)
cat("\nSeconds (minimum of 3 runs; pbkrtest/lmerTest 1 run):\n"); print(res, row.names = FALSE)

## Permutation refits (type III tables, every refit) and means over an empty cell
oats <- as.data.frame(nlme::Oats)
oats$nitro <- factor(oats$nitro)
oats$Variety <- factor(oats$Variety)
oe <- droplevels(oats[!(oats$Variety == "Victory" & oats$nitro == "0.6"), ])   # empty cell
mtc <- mtcars; mtc$cyl <- factor(mtc$cyl); mtc$am <- factor(mtc$am)
ml <- lm(mpg ~ wt + cyl * am, data = mtc)
mm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = oats)
me <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = oe)
mg <- glm(round(mpg) ~ wt + cyl, family = poisson, data = mtc)
q <- function(x) suppressWarnings(suppressMessages(x))
perm <- data.frame(
  lm_raw_III = tm(q(permmodels(ml, nperm = 999, type = "III", seed = 1, method = "raw", prt = FALSE))),
  lm_fl_III = tm(q(permmodels(ml, nperm = 999, type = "III", seed = 1, method = "freedman-lane", prt = FALSE))),
  glm_III = tm(q(permmodels(mg, nperm = 199, type = "III", seed = 1, prt = FALSE))),
  lmer_refit_III = tm(q(permmodels(mm, nperm = 49, type = "III", seed = 1, prt = FALSE)), 1),
  lmer_norefit_III = tm(q(permmodels(mm, nperm = 999, type = "III", seed = 1, refit = FALSE, prt = FALSE))))
cat("\nPermutation tests, seconds:\n"); print(perm, row.names = FALSE)
pl <- q(permmodels(ml, nperm = 999, type = "II", seed = 1, method = "freedman-lane", prt = FALSE))
cells <- data.frame(
  lmer_full = tm(q(predictmeans(mm, "nitro:Variety", plot = FALSE, prtnum = FALSE))),
  lmer_empty_cell = tm(q(predictmeans(me, "nitro:Variety", plot = FALSE, prtnum = FALSE))),
  lmer_empty_main = tm(q(predictmeans(me, "nitro", plot = FALSE, prtnum = FALSE))),
  lm_permlist = tm(q(predictmeans(ml, "cyl:am", permlist = pl, plot = FALSE, prtnum = FALSE))))
cat("\npredictmeans(), seconds:\n"); print(cells, row.names = FALSE)
