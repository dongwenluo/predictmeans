# predictmeans: start here

Package version 2.0.0, ready for CRAN. Before any work, read `CLAUDE.md`
(the maintainer's rules) and `HANDOVER.md` (state of the package, what was
done, the open next steps in order). Look up details on demand in row O16 of
`dev/predictmeans-2.0-plan.md` and in `dev/NEWS-2.0-development.md`; don't read
them in full.

## Working rules (short form)

- One step per message; run only the related tests while working.
- Every change: a short user-facing NEWS.md entry (2.0.0 release note) and an
  O16 update.
- Full test suite, R CMD check and one git commit at the end of a block.
- Argument names are snake_case; result columns keep broom names.
- Help pages come from roxygen in `R/*.R`: edit there, then `devtools::document()`.
- AGENTS.md, CLAUDE.md, HANDOVER.md stay in the root but are kept out of the
  package (`.Rbuildignore`) and the site: build it with `dev/build-site.R`
  (pkgdown publishes every top-level `.md`; the GitHub workflow removes them too).

- `vignettes/resampling.Rmd` and `vignettes/smooth.Rmd` are generated: edit the
  `.Rmd.orig` files, then `Rscript dev/precompute-vignettes.R` (CRAN check time).
- A new slow test (over 0.5 s) starts with `skip_on_cran()`.

## Testing traps

- Run `options(warnPartialMatchDollar = FALSE)` first: if an earlier session
  left it on, tests warn from inside glmmTMB and mgcv.
- Keep the global environment clear of test-like objects (`b`, `d`, ...):
  tests run in the live session can find them. A fresh R session is safest.
- `Sys.setenv(NOT_CRAN = "true")`; one file with
  `testthat::test_file(path, package = "predictmeans", load_package = "source")`;
  full suite with `devtools::test(reporter = "silent")`, printing totals only.
- Reinstall the dev build before rendering articles or the site
  (`devtools::install(dependencies = FALSE, upgrade = FALSE)`): eda4treeR pulls
  CRAN predictmeans 1.1.1 into the library.
