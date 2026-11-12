## Summary of study 1b (pairs, dev/sim-perm-pairs/, and term tests,
## dev/sim-perm-fixed/). Run from the package root:
##   Rscript dev/sim-study/summarise-pairs.R
## Writes dev/sim-study/results-1b.md.
##
## Pairs: for each cell and method, the rejection rate at 5% of each true
## null pair (mean and largest over those pairs), the familywise error over
## them (Holm, max-|t| or Tukey: the smallest adjusted p-value of the true
## null pairs <= 0.05), and the power for the pairs whose means differ.
## Term tests: rejection rate at 5% of each method.

null_pairs <- list(
  pr_rcbd_global = c("1v2", "1v3", "1v4", "2v3", "2v4", "3v4"),
  pr_rcbd_partial = c("1v2", "1v3", "2v3"), pr_rcbd_skew_partial = c("1v2", "1v3", "2v3"),
  pr_sp_A_partial = c("1v2", "1v3", "2v3"), pr_sp_B_partial = "1v2",
  pr_lev_global = c("1v2", "1v3", "2v3"),
  pr_lme_varident = c("1v2", "1v3", "1v4", "2v3", "2v4", "3v4"),
  pr_gls_ar1 = c("1v2", "1v3", "2v3"), pr_rcbd_power = character())

out <- c("# Study 1b results", "",
         "Rejection rates at the 5% level. Monte Carlo range for a true rate of 5%: 0.036-0.064 (1000 data sets),",
         "0.031-0.069 (500).", "")

files <- list.files("dev/sim-perm-pairs", pattern = "\\.rds$", full.names = TRUE)
if (length(files)) {
  out <- c(out, "## Pairwise comparisons", "",
           "| cell | method | per pair (mean) | per pair (max) | FWER Holm | FWER max-t / Tukey | power (mean) |",
           "|---|---|---|---|---|---|---|")
  for (f in files) {
    x <- readRDS(f)
    p <- x$p
    nul <- null_pairs[[x$cell]]
    for (m in c("joint", "flpair", "flterm", "param")) {
      cols <- grep(paste0("^", m, "\\.[0-9]+v[0-9]+$"), colnames(p), value = TRUE)
      pairs <- sub(".*\\.", "", cols)
      r <- colMeans(p[, cols, drop = FALSE] <= 0.05, na.rm = TRUE)
      rn <- r[pairs %in% nul]
      rp <- r[!pairs %in% nul]
      fam <- function(col) if (col %in% colnames(p)) sprintf("%.3f", mean(p[, col] <= 0.05, na.rm = TRUE)) else ""
      mx <- if (m == "param") fam("param.tukey") else fam(paste0(m, ".maxT"))
      out <- c(out, sprintf("| %s | %s | %s | %s | %s | %s | %s |", x$cell, m,
                            if (length(rn)) sprintf("%.3f", mean(rn)) else "",
                            if (length(rn)) sprintf("%.3f", max(rn)) else "",
                            if (length(rn)) fam(paste0(m, ".holm")) else "", if (length(rn)) mx else "",
                            if (length(rp)) sprintf("%.3f", mean(rp)) else ""))
    }
  }
  out <- c(out, "")
}

files <- list.files("dev/sim-perm-fixed", pattern = "\\.rds$", full.names = TRUE)
if (length(files)) {
  out <- c(out, "## Tests of a fixed term", "", "| cell | n | method | reject 5% | reject 1% |", "|---|---|---|---|---|")
  for (f in files) {
    x <- readRDS(f)
    for (m in colnames(x$p)) {
      out <- c(out, sprintf("| %s | %d | %s | %.3f | %.3f |", x$cell, x$R_ok, m,
                            mean(x$p[, m] <= 0.05, na.rm = TRUE), mean(x$p[, m] <= 0.01, na.rm = TRUE)))
    }
  }
}
writeLines(out, "dev/sim-study/results-1b.md")
cat(out, sep = "\n")
