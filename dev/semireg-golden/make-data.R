## The data sets of the golden semi-parametric examples (M. P. Wand's lme
## scripts in original/), from the files in inst/extdata/semireg/, saved as
## package data (data/*.rda). tests/testthat/test-semireg-golden.R checks that
## the two agree.
##   Rscript dev/semireg-golden/make-data.R
rd <- function(f, ...) utils::read.table(file.path("inst/extdata/semireg", f), ...)
growthIndiana <- rd("growthIndiana.txt", header = TRUE)
SBMDfem <- rd("SBMDfem.txt", header = TRUE)
SBMDfem$ethnicity <- factor(SBMDfem$ethnicity)
ozoneCalif <- rd("ozoneCalif.txt", header = TRUE)
scallop <- rd("scallop.txt", header = TRUE)
scallopKnots <- stats::setNames(rd("scallopKnots.txt"), c("longitude", "latitude"))
scallopBdry <- stats::setNames(rd("scallopBdry.txt"), c("longitude", "latitude"))
copperKnots <- stats::setNames(rd("copperKnots.txt"), c("xcoord", "ycoord"))
wageTUnion <- rd("wageTUnion.txt", header = TRUE)
files <- character(0)
for (nm in c("growthIndiana", "SBMDfem", "ozoneCalif", "scallop", "scallopKnots", "scallopBdry", "copperKnots",
             "wageTUnion")) {
  files[nm] <- file.path("data", paste0(nm, ".rda"))
  save(list = nm, file = files[nm], compress = "bzip2", version = 2)
}
tools::resaveRdaFiles(files, compress = "auto")    # these files only, not the package's other data
