## Plan step 13: when is refit = FALSE (Freedman-Lane with V fixed) safe?
## Split plot, test of the whole-plot factor A (null), B with effects, with
## 6 (study 1b: V fixed 8.6%), 12, 24 and 48 blocks (whole-plot error df 15,
## 33, 69, 141). Methods: Freedman-Lane with V fixed (fl), with refits (flr),
## Kenward-Roger F (kr). Uses dev/sim-perm-lib.R. Results: dev/sim-perm-refit/.
##
## Run: Rscript dev/sim-perm-refit.R [R] [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500
nperm <- if (length(args) >= 2) as.integer(args[2]) else 199
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")
pA <- function(tab, term = "A") as.numeric(tab[term, "Perm_p_value"])

description <- list(
  rf_sp_b12 = "split plot, 12 blocks x 4 whole plots (A) x 3 (B), A null, B (0,1,2); test of A",
  rf_sp_b24 = "as rf_sp_b12 with 24 blocks",
  rf_sp_b48 = "as rf_sp_b12 with 48 blocks")

one_set <- function(cell) {
  nb <- switch(cell, rf_sp_b12 = 12, rf_sp_b24 = 24, rf_sp_b48 = 48, stop("unknown cell ", cell))
  d <- sp_design(nblock = nb)
  d <- d[order(d$block, d$wp, d$sub), ]
  d$y <- c(0, 1, 2)[d$B] + rnorm(nb)[d$block] + rnorm(nlevels(d$wpid))[d$wpid] + rnorm(nrow(d))
  m <- lmer(y ~ A * B + (1 | block) + (1 | block:A), data = d, control = ctl)
  fl <- permmodels(m, nperm = nperm, type = "II", prt = FALSE, method = "freedman-lane", refit = FALSE)$ANOVA
  flr <- predictmeans:::pm_perm_fl_refit(m, nperm, "II", FALSE, NULL, FALSE, NULL, 1, only = "A")$ANOVA
  kr <- anova(lmerTest::as_lmerModLmerTest(m), type = "II", ddf = "Kenward-Roger")["A", "Pr(>F)"]
  c(fl = pA(fl), flr = pA(flr), kr = kr)
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-perm-refit", description)
cat("SIM DONE\n")
