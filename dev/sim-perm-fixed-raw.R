## Fix-up for dev/sim-perm-fixed.R. Its first runs of the LMM cells
## (fx2_lmm_sp, fx2_lmm_lev, fx2_lmm_unbal) called permmodels() for the "raw"
## column without method = "raw"; with the 2.0 default (method = "auto") that
## was Freedman-Lane with refits for all terms, not the raw scheme. This
## script regenerates the same data sets (same seeds) and computes the raw
## scheme's p-value only, then stores it as "raw" in the cell's file; the old
## column is kept as "flr_all" (Freedman-Lane with refits, a second
## independent set of permutations).
## Run: Rscript dev/sim-perm-fixed-raw.R <cells> <cores>
args <- commandArgs(TRUE)
fix_cells <- strsplit(args[1], ",")[[1]]
cores <- as.integer(args[2])
raw_only <- TRUE
src <- readLines("dev/sim-perm-fixed.R")
src[grep("^cells <- ", src)] <- "cells <- character()"
src <- src[!grepl("^run_cells\\(|^cat\\(\"SIM DONE", src)]
commandArgs <- function(...) c("500", "199", "none", "1")
eval(parse(text = src))
for (cell in fix_cells) {
  f <- file.path("dev/sim-perm-fixed", paste0(cell, ".rds"))
  x <- readRDS(f)
  if ("flr_all" %in% colnames(x$p)) { cat("CELL", cell, "already fixed\n"); next }
  base <- sum(utf8ToInt(cell)) * 1000
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(x$R), function(i) {
    set.seed(base + i)
    tryCatch(suppressWarnings(suppressMessages(one_set(cell)))[["raw"]], error = function(e) NA_real_)
  }, mc.cores = cores)
  raw <- unlist(res)
  ok <- vapply(seq_len(x$R), function(i) TRUE, TRUE)   # the cell's data sets, in order
  p <- x$p
  colnames(p)[colnames(p) == "raw"] <- "flr_all"
  # rows of p are the data sets that succeeded in the first run (all, for these cells)
  stopifnot(nrow(p) == x$R)
  x$p <- cbind(raw = raw, p)
  x$reject <- colMeans(x$p <= 0.05, na.rm = TRUE)
  x$note <- paste("raw recomputed by dev/sim-perm-fixed-raw.R on", format(Sys.time()),
                  "; flr_all = the first run's mislabelled raw column (Freedman-Lane with refits)")
  saveRDS(x, f)
  cat(sprintf("CELL %-20s raw %.3f  (%.1f min)\n", cell, mean(raw <= 0.05, na.rm = TRUE),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
