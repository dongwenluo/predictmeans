# Reverse dependencies of predictmeans on CRAN, checked with the CRAN version
# and with this version (the working directory: the package sources).
# Usage: Rscript dev/revdep-check.R   (needs CRAN; run by
# .github/workflows/revdep-check.yaml). Results in revdep-results/.

options(repos = c(CRAN = "https://cloud.r-project.org"))
out <- "revdep-results"
dir.create(out, showWarnings = FALSE)
db <- utils::available.packages()
rev <- tools::package_dependencies("predictmeans", db = db, reverse = TRUE,
                                   which = c("Depends", "Imports", "LinkingTo", "Suggests"))[[1]]
cat("Reverse dependencies:", if (length(rev)) rev else "none", "\n")
writeLines(rev, file.path(out, "revdeps.txt"))
if (!length(rev)) quit(save = "no")

# the reverse dependencies' own dependencies (with their Suggests), and the
# two versions of predictmeans in their own libraries
utils::install.packages(rev, dependencies = TRUE)
lib_cran <- normalizePath(file.path(out, "lib-cran"), mustWork = FALSE)
lib_new <- normalizePath(file.path(out, "lib-new"), mustWork = FALSE)
dir.create(lib_cran, showWarnings = FALSE)
dir.create(lib_new, showWarnings = FALSE)
utils::install.packages("predictmeans", lib = lib_cran)
stopifnot(system2("R", c("CMD", "INSTALL", "--no-multiarch", "-l", shQuote(lib_new), ".")) == 0)

src <- utils::download.packages(rev, destdir = out, type = "source")[, 2]
status <- function(dir) {
  log <- file.path(dir, paste0(sub("_.*", "", basename(dir)), ".Rcheck"), "00check.log")
  if (!file.exists(log)) return("no check log")
  st <- grep("^Status:", readLines(log), value = TRUE)
  if (length(st)) sub("^Status: ", "", st[length(st)]) else "unfinished"
}
res <- do.call(rbind, lapply(src, function(tarball) {
  pkg <- sub("_.*", "", basename(tarball))
  one <- function(lib, tag) {
    dir <- file.path(out, paste0(pkg, "_", tag))
    dir.create(dir, showWarnings = FALSE)
    libs <- paste(c(lib, .libPaths()), collapse = .Platform$path.sep)
    system2("R", c("CMD", "check", "--no-manual", "-o", shQuote(dir), shQuote(tarball)),
            env = c(paste0("R_LIBS=", libs), "_R_CHECK_FORCE_SUGGESTS_=false"),
            stdout = file.path(dir, "check-output.txt"), stderr = file.path(dir, "check-output.txt"))
    status(dir)
  }
  data.frame(package = pkg, cran = one(lib_cran, "cran"), new = one(lib_new, "new"))
}))
res$changed <- res$cran != res$new
print(res, row.names = FALSE)
utils::write.csv(res, file.path(out, "summary.csv"), row.names = FALSE)
if (any(res$changed)) cat("\nStatus changed for:", res$package[res$changed],
                          "- see revdep-results/<package>_new/*.Rcheck/00check.log\n")
