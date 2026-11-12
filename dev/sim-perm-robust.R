## Plan step 12: heteroscedasticity-robust (HC3) Freedman-Lane permutations,
## permmodels(robust = TRUE), for lm with error variances the model ignores.
## Methods: usual F test (fs), Freedman-Lane with the usual F (fl), robust
## Freedman-Lane (rob), robust HC3 Wald F against F(q, residual df) (hc3).
## Uses dev/sim-perm-lib.R. Results: dev/sim-perm-robust/.
##
## Run: Rscript dev/sim-perm-robust.R [R] [nperm] [cells|all] [cores]

args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
nperm <- if (length(args) >= 2) as.integer(args[2]) else 199
cores <- if (length(args) >= 4) as.integer(args[4]) else 1
source("dev/sim-perm-lib.R")

description <- list(
  rb_bf = "lm y ~ x + A; A groups of 6/12/18 with error SD 3/1/1 (Behrens-Fisher-like), A null; test A",
  rb_bf_small = "as rb_bf with groups of 3/6/9",
  rb_bf_large = "as rb_bf with groups of 12/24/36",
  rb_eq = "as rb_bf with equal SDs (control)",
  rb_hetx = "lm y ~ x + A, 36 units, error SD |x| + 0.2, slope of x null; test x",
  rb_hetx_large = "as rb_hetx with 108 units")

one_set <- function(cell) {
  if (startsWith(cell, "rb_hetx")) {
    n <- if (cell == "rb_hetx_large") 108 else 36
    d <- data.frame(A = factor(rep(1:3, length.out = n)), x = rnorm(n))
    d$y <- c(0, 1, 2)[d$A] + (abs(d$x) + 0.2) * rnorm(n)
    term <- "x"
  } else {
    k <- switch(cell, rb_bf_small = 3, rb_bf_large = 12, 6)
    d <- data.frame(A = factor(rep(1:3, k * (1:3))))
    sds <- if (cell == "rb_eq") c(1, 1, 1) else c(3, 1, 1)
    d$x <- rnorm(nrow(d))
    d$y <- d$x + sds[d$A] * rnorm(nrow(d))
    term <- "A"
  }
  fm <- lm(y ~ x + A, data = d)
  fs <- drop1(fm, test = "F")[term, "Pr(>F)"]
  fl <- as.numeric(permmodels(fm, nperm = nperm, type = "II", prt = FALSE)$ANOVA[term, "Perm_p_value"])
  rb <- permmodels(fm, nperm = nperm, type = "II", prt = FALSE, robust = TRUE)$ANOVA
  hc3 <- pf(rb[term, "F value"], rb[term, "Df"], df.residual(fm), lower.tail = FALSE)
  c(fs = fs, fl = fl, rob = rb[term, "Perm_p_value"], hc3 = hc3)
}

cells <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else names(description)
run_cells(cells, one_set, R, nperm, cores, "dev/sim-perm-robust", description)
cat("SIM DONE\n")
