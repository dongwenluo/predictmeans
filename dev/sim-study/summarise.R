## Collect the results of the permutation simulation studies into one table.
## Run from the package root: Rscript dev/sim-study/summarise.R
## Writes dev/sim-study/results.csv (one row per cell x method) and
## dev/sim-study/results.md (the same as a markdown table).
##
## For each cell and method: rejection rate at 5% (and 1%), its Monte Carlo
## standard error, the number of usable data sets and failures, and a flag:
##   "ok"      within the 95% Monte Carlo range of the nominal level,
##   "liberal" above it, "conservative" below it (null cells only),
##   "power"   for power cells (name ends in "_power").

dirs <- c("dev/sim-perm-ranef", "dev/sim-perm-edge")
files <- unlist(lapply(dirs, function(d) list.files(d, pattern = "\\.rds$", full.names = TRUE)))
rows <- lapply(files, function(f) {
  x <- readRDS(f)
  p <- x$p
  is_power <- grepl("_power$", x$cell)
  n_ok <- if (!is.null(x$R_ok)) x$R_ok else nrow(p)
  do.call(rbind, lapply(colnames(p), function(m) {
    r5 <- mean(p[, m] <= 0.05, na.rm = TRUE)
    r1 <- mean(p[, m] <= 0.01, na.rm = TRUE)
    lo <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / n_ok); hi <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / n_ok)
    flag <- if (is_power) "power" else if (r5 > hi) "liberal" else if (r5 < lo) "conservative" else "ok"
    data.frame(study = basename(dirname(f)), cell = x$cell, method = m,
               reject_05 = round(r5, 3), mc_se = round(sqrt(r5 * (1 - r5) / n_ok), 4),
               reject_01 = round(r1, 3), flag = flag, n_ok = n_ok, R = x$R,
               failures = x$R - n_ok, nperm = x$nperm, minutes = round(x$minutes, 1),
               description = if (is.null(x$description)) "" else x$description)
  }))
})
res <- do.call(rbind, rows)
res <- res[order(res$study, res$cell, res$method), ]
dir.create("dev/sim-study", showWarnings = FALSE)
utils::write.csv(res, "dev/sim-study/results.csv", row.names = FALSE)

md <- c("| study | cell | method | reject 5% | MC SE | reject 1% | flag | n | failures |",
        "|---|---|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %.3f | %.4f | %.3f | %s | %d | %d |", res$study, res$cell, res$method,
                res$reject_05, res$mc_se, res$reject_01, res$flag, res$n_ok, res$failures))
writeLines(md, "dev/sim-study/results.md")
print(res[, c("cell", "method", "reject_05", "reject_01", "flag", "n_ok", "failures")], row.names = FALSE)
