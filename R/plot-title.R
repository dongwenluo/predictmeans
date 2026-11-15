## Titles of the package's plots: one helper breaks long ones into lines.
##
## A title is broken at spaces so that no line is longer than 'width'
## characters; line breaks already in it are kept (also one at the end, used
## as space under the title). The default width is what fits a plot about 7
## inches wide at the title's font size (base size 'basesz' x 'rel': 1.2 for
## titles, 1 for subtitles), divided by the number of panels side by side;
## options(predictmeans.title_width = ) sets another number of characters.
pm_wrap_title <- function(x, basesz = 12, rel = 1.2, panels = 1, width = NULL) {
  if (is.null(x) || !is.character(x) || !length(x)) return(x)
  if (is.null(width)) width <- getOption("predictmeans.title_width", NULL)
  if (is.null(width)) width <- max(20, floor(900 / (basesz * rel) / panels))
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    tail_nl <- regmatches(s, regexpr("\n*$", s))
    parts <- strsplit(sub("\n*$", "", s), "\n", fixed = TRUE)[[1]]
    paste0(paste(vapply(parts, function(p) paste(strwrap(p, width + 1), collapse = "\n"), character(1)),
                 collapse = "\n"), tail_nl)
  }, character(1), USE.NAMES = FALSE)
}
