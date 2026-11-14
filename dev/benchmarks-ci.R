# Timing check for pull requests: key timings of this version against the
# base branch, both run on the same machine (so the runner's speed cancels).
# Fails (exit status 1) when a timing is more than 1.5 times the base's and
# at least 0.1 s longer.
# Usage: Rscript --no-init-file dev/benchmarks-ci.R <lib-base> <lib-new>
# (run by .github/workflows/benchmarks.yaml; each library is timed in its own
# R process, since one session cannot load two versions of a package).

args <- commandArgs(TRUE)

if (identical(args[1], "--time")) {
  lib <- args[2]
  suppressPackageStartupMessages(library(predictmeans, lib.loc = lib))
  tm <- function(expr, reps = 3) {
    e <- substitute(expr)
    pf <- parent.frame()
    r <- tryCatch(suppressWarnings(suppressMessages(eval(e, pf))), error = function(e) e)
    if (inherits(r, "error")) return(NA_real_)     # not in this version
    min(vapply(seq_len(reps), function(i) system.time(suppressWarnings(suppressMessages(eval(e, pf))))[["elapsed"]],
               numeric(1)))
  }
  set.seed(1)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), Block = factor(1:50), rep = 1:2)
  d$WP <- interaction(d$Block, d$A)
  d$y <- as.numeric(d$A) + as.numeric(d$B) + rnorm(50)[d$Block] + rnorm(nlevels(d$WP))[d$WP] + rnorm(nrow(d))
  m <- lme4::lmer(y ~ A*B + (1|Block) + (1|WP), data = d)                       # n = 1200
  de <- droplevels(d[!(d$A == "3" & d$B == "4"), ])
  me <- lme4::lmer(y ~ A*B + (1|Block) + (1|WP), data = de)                     # an empty cell
  oats <- as.data.frame(nlme::Oats); oats$nitro <- factor(oats$nitro)
  mo <- lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats)
  mtc <- mtcars; mtc$cyl <- factor(mtc$cyl); mtc$am <- factor(mtc$am)
  ml <- lm(mpg ~ wt + cyl*am, data = mtc)
  big <- expand.grid(S = factor(1:6), W = factor(1:4), B = factor(1:40))
  big$y <- rnorm(nrow(big)) + as.numeric(big$W)
  big$y[sample(nrow(big), 20)] <- NA
  ma <- aov(y ~ W*S + Error(B/W), data = big)                                  # n = 960, 20 missing
  cm <- rbind(c(-1, 1, 0, 0), c(0, -1, 1, 0), c(0, 0, -1, 1))
  q <- function(...) predictmeans(..., plot = FALSE, prtnum = FALSE)
  res <- c(
    predictmeans_lmer_1200 = tm(q(m, "A:B", adj = "tukey")),
    predictmeans_empty_cell = tm(q(me, "A:B", atvar = "A", adj = "holm")),
    predictmeans_satterthwaite = tm(q(m, "A:B", df_method = "Satterthwaite")),
    contrastmeans_lmer = tm(contrastmeans(mo, "nitro", ctrmatrix = cm)),
    predictmeans_aov_missing = tm(q(ma, "W:S", pairwise = TRUE)),
    permmodels_lm_III = tm(permmodels(ml, nperm = 499, type = "III", seed = 1, method = "raw", prt = FALSE), 1),
    permmodels_lmer_refit = tm(permmodels(mo, nperm = 19, seed = 1, prt = FALSE), 1),
    residplot_lmer = tm(residplot(m, plot = FALSE)),
    varcomp_lmer = tm(varcomp(m)))
  saveRDS(res, args[3])
  quit(save = "no")
}

lib_base <- args[1]
lib_new <- args[2]
# this script, to run it again for each library
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
run <- function(lib) {
  f <- tempfile(fileext = ".rds")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--no-init-file", shQuote(self), "--time", shQuote(lib), shQuote(f)))
  if (status != 0 || !file.exists(f)) stop("timing failed for ", lib, call. = FALSE)
  readRDS(f)
}
b <- run(lib_base)
n <- run(lib_new)
tab <- data.frame(base = b[names(n)], new = n, row.names = names(n))
tab$ratio <- round(tab$new / tab$base, 2)
tab$slower <- !is.na(tab$ratio) & tab$ratio > 1.5 & (tab$new - tab$base) >= 0.1
cat("Seconds (minimum of 3 runs; permutation tests 1 run):\n")
print(tab)
if (any(tab$slower)) {
  cat("\nSlower than the base branch:", paste(rownames(tab)[tab$slower], collapse = ", "), "\n")
  quit(save = "no", status = 1)
}
cat("\nNo timing is clearly slower than the base branch.\n")
