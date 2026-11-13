## Tables of the df simulations (dev/sim-df.R, dev/sim-df-glmm.R): percent of
## true null hypotheses rejected at the 5% level, one table per cell, rows =
## contrasts or F tests, columns = df methods. With 1000 data sets a true 5%
## lies in 3.6-6.4% (99% interval). Writes <dir>/results.md.
##
## Run: Rscript dev/sim-df-summary.R dev/sim-df   (or dev/sim-df-glmm)

dir <- commandArgs(TRUE)[1]
files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
out <- c(paste0("# Results: ", basename(dir)), "",
         "Percent of true null hypotheses rejected at the 5% level (two-sided t tests of single contrasts,",
         "F tests of terms). With R data sets a true 5% lies in 5 +- 2.58 sqrt(0.05 x 0.95 / R) x 100 (99%).",
         "Bold: outside that range.", "")
for (f in files) {
  r <- readRDS(f)
  p <- r$p[, !startsWith(colnames(r$p), "df_"), drop = FALSE]
  rej <- colMeans(p <= 0.05, na.rm = TRUE) * 100
  nm <- do.call(rbind, strsplit(names(rej), ".", fixed = TRUE))
  tab <- tapply(rej, list(factor(nm[, 1], unique(nm[, 1])), factor(nm[, 2], unique(nm[, 2]))), identity)
  half <- 2.58 * sqrt(0.05 * 0.95 / r$R_ok) * 100
  fmt <- function(x) ifelse(is.na(x), "", ifelse(abs(x - 5) > half, sprintf("**%.1f**", x), sprintf("%.1f", x)))
  out <- c(out, paste0("## ", r$cell), "", r$description, "",
           sprintf("%d data sets (%d failed); a true 5%% lies in %.1f-%.1f%%.", r$R_ok, r$R - r$R_ok, 5 - half, 5 + half), "",
           paste0("| | ", paste(colnames(tab), collapse = " | "), " |"),
           paste0("|---|", paste(rep("---", ncol(tab)), collapse = "|"), "|"))
  for (i in seq_len(nrow(tab))) out <- c(out, paste0("| ", rownames(tab)[i], " | ", paste(fmt(tab[i, ]), collapse = " | "), " |"))
  dfc <- r$p[, startsWith(colnames(r$p), "df_"), drop = FALSE]
  if (ncol(dfc)) out <- c(out, "", paste("Mean df:", paste(sprintf("%s %.1f", sub("^df_", "", colnames(dfc)), colMeans(dfc)),
                                                         collapse = ", ")))
  if (length(r$errors)) out <- c(out, "", paste("Errors:", paste(names(r$errors), r$errors, collapse = "; ")))
  out <- c(out, "")
}
writeLines(out, file.path(dir, "results.md"))
cat(out, sep = "\n")
