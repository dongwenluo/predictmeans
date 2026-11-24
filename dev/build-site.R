# Build the pkgdown site locally without the maintainer's private notes.
# pkgdown renders every top-level .md file as a page; AGENTS.md, CLAUDE.md and
# HANDOVER.md stay in the package root (Positron reads AGENTS.md as project
# memory) but must not be on the site. The GitHub workflow
# (.github/workflows/pkgdown.yaml) removes them the same way.
#
# Usage, from the package root:  Rscript dev/build-site.R
# (or source("dev/build-site.R") in a fresh R session)

private <- c("AGENTS", "CLAUDE", "HANDOVER")
pkgdown::build_site(preview = FALSE, new_process = FALSE, install = TRUE)
gone <- unlist(lapply(private, function(p) Sys.glob(file.path("docs", paste0(p, ".*")))))
unlink(gone)
if (length(gone)) message("Removed from the site: ", paste(basename(gone), collapse = ", "))
left <- unlist(lapply(private, function(p) Sys.glob(file.path("docs", paste0(p, ".*")))))
stopifnot(length(left) == 0)
