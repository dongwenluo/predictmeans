# Regenerate tests/testthat/_data/reference-v1.1.1.rds from predictmeans 1.1.1
# (commit 7f398a5), so the tests can check that later versions reproduce its
# numbers. From the package root:
#
#   mkdir -p /tmp/pm-1.1.1 /tmp/lib-1.1.1
#   git archive 7f398a5 | tar -x -C /tmp/pm-1.1.1
#   R CMD INSTALL --library=/tmp/lib-1.1.1 /tmp/pm-1.1.1
#   Rscript --no-init-file dev/make-reference.R /tmp/lib-1.1.1

lib <- commandArgs(TRUE)[1]
suppressPackageStartupMessages({
  library(predictmeans, lib.loc = lib)
  library(HRW)
})
stopifnot(packageVersion("predictmeans", lib.loc = lib) == "1.1.1")
source("tests/testthat/helper-reference-cases.R")
grDevices::pdf(NULL)
saveRDS(reference_cases(), "tests/testthat/_data/reference-v1.1.1.rds", version = 2)
