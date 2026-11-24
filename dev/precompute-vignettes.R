## Knit the vignettes whose code is too slow for CRAN's check time limit
## (CRAN pre-test 2026-11-24: overall check time 21 min > 10 min). The source
## is vignettes/<name>.Rmd.orig (left out of the package by .Rbuildignore);
## this script writes vignettes/<name>.Rmd with the results and figures in
## it, so the package's vignette build runs no model code. Their code is
## still tested: run this script after any change to the package that can
## change their output, and commit the .Rmd and the figure files.
## Usage, from the package root: Rscript dev/precompute-vignettes.R
pkgload::load_all(".", quiet = TRUE)   # the development version, not an installed one
precomputed <- c("resampling", "smooth")
old <- setwd("vignettes"); on.exit(setwd(old))
for (v in precomputed) {
  unlink(Sys.glob(paste0(v, "-[0-9a-z]*.png")))
  knitr::opts_chunk$set(fig.path = paste0(v, "-"), dpi = 96)
  knitr::knit(paste0(v, ".Rmd.orig"), paste0(v, ".Rmd"), quiet = TRUE, envir = new.env())
  knitr::opts_chunk$restore()
  cat(v, ":", length(Sys.glob(paste0(v, "-*.png"))), "figures\n")
}
