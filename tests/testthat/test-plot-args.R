## Plot settings: plot_args, the same names as arguments (1.x), and the stages.

pa_fit <- function() lm(yield ~ Block + nitro * Variety, data = transform(as.data.frame(nlme::Oats), nitro = factor(nitro)))

test_that("plot_args and the 1.x argument names give the same plots", {
  fm <- pa_fit()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  a <- predictmeans(fm, "nitro", which.plots = "mean", prtnum = FALSE,
                    plot_args = list(plottitle = "Yield by nitrogen", plotxlab = "N", basesz = 14, lsd_bar = FALSE))
  b <- predictmeans(fm, "nitro", which.plots = "mean", prtnum = FALSE,
                    plottitle = "Yield by nitrogen", plotxlab = "N", basesz = 14, lsd_bar = FALSE)
  ga <- a$predictmeansPlot$meanPlot; gb <- b$predictmeansPlot$meanPlot
  expect_identical(ga$labels$title, "Yield by nitrogen\n")
  expect_identical(ga$labels$x, "N")
  expect_equal(ggplot2::layer_data(ga, 1), ggplot2::layer_data(gb, 1))
  expect_identical(ga$labels, gb$labels)
  # abbreviated names, as before
  cc <- predictmeans(fm, "nitro", which.plots = "mean", prtnum = FALSE, plott = "Yield by nitrogen", prtplt = FALSE)
  expect_identical(cc$predictmeansPlot$meanPlot$labels$title, "Yield by nitrogen\n")
  # plot_args wins over the same name given as an argument
  d <- predictmeans(fm, "nitro", which.plots = "mean", prtnum = FALSE, plottitle = "one", plot_args = list(plottitle = "two"))
  expect_identical(d$predictmeansPlot$meanPlot$labels$title, "two\n")
  # an order for the plot asks for the mean plot
  e <- predictmeans(fm, "nitro:Variety", prtnum = FALSE, plot_args = list(plotord = c(2, 1), prtplt = FALSE))
  expect_s3_class(e$predictmeansPlot$meanPlot, "ggplot")
})

test_that("unknown names are errors; the 1.x plot choices still work with a warning", {
  fm <- pa_fit()
  expect_error(predictmeans(fm, "nitro", prtnum = FALSE, colour = "red"), "unused argument\\(s\\) in predictmeans\\(\\): colour")
  expect_error(predictmeans(fm, "nitro", prtnum = FALSE, plot_args = list(colour = "red")), "'plot_args' can hold")
  expect_error(predictmeans(fm, "nitro", prtnum = FALSE, plot_args = "big"), "named list")
  expect_error(predictmeans(fm, "nitro", prtnum = FALSE, plot_args = list(14)), "named list")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_warning(r <- predictmeans(fm, "nitro", barplot = TRUE, prtplt = FALSE, prtnum = FALSE),
                 'barplot is deprecated: use which.plots = c\\("mean", "bar", "pval", "back"\\)')
  expect_s3_class(r$predictmeansBarPlot, "ggplot")
  expect_identical(names(formals(predictmeans))[1:8],
                   c("model", "modelterm", "data", "pairwise", "atvar", "adj", "Df", "level"))
  expect_false(any(c("mplot", "barplot", "pplot", "bkplot", "plottitle", "newwd") %in% names(formals(predictmeans))))
  expect_identical(names(pm_plot_args()), c(names(pm_plot_defaults), "legacy"))
})

test_that("the stages are functions that return named lists; a later stage replaces a value", {
  fm <- pa_fit()
  pa <- pm_plot_args()
  op <- options(scipen = 6); on.exit(options(op))             # (as predictmeans() sets for its tables)
  opt <- pm_setup(list(model = fm, modelterm = "nitro:Variety", data = NULL, pairwise = FALSE, atvar = "Variety",
                       adj = "tukey", Df = NULL, level = NULL, covariate = NULL, meandecr = NULL, letterCI = FALSE,
                       trans = I, trans_default = TRUE, transOff = 0, count = FALSE, which.plots = c("mean", "pval", "back"),
                       plot = FALSE, prtnum = FALSE, permlist = NULL, ncore = 1, ndecimal = 4, cov_reduce = mean,
                       exclude = NULL, vcov_adjust = TRUE, vcov_adjust_default = TRUE, weights = "equal",
                       df_method = "auto", vcov = NULL, equivalence = NULL, bias_adjust = FALSE, sigma = NULL,
                       average = "link"), pa)
  expect_true(opt$pairwise)                                   # 'adj' and 'atvar' ask for comparisons
  expect_identical(opt$level, 0.05)
  est <- pm_estimate(opt)
  cmp <- pm_compare(opt, est)
  tab <- pm_tables(opt, est, cmp)
  expect_true(all(c("K", "mt", "inf", "mean.table") %in% names(est)))
  expect_true(all(c("pair_tab", "LSD", "SED.out", "cmp_int", "atvar", "mt", "Df") %in% names(cmp)))
  expect_identical(names(tab), c("meanTable", "ci_level", "bkmt", "nc"))
  expect_null(est$mt$Df)
  expect_false(is.null(cmp$mt$Df))                            # the df of the means are added by pm_compare()
  v <- pm_flat(opt, est, cmp, tab)
  expect_false(anyDuplicated(names(v)) > 0)
  expect_identical(v$mt, cmp$mt)
  expect_identical(v$bkmt, tab$bkmt)
  out <- pm_assemble(v, list(meanPlot = NULL, ciPlot = NULL, predictmeansBarPlot = NULL))
  ref <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", prtnum = FALSE)
  expect_equal(out, ref)
  # a NULL from a later stage replaces an earlier value
  expect_null(pm_flat(list(a = 1, b = 2), list(a = NULL))$a)
  expect_identical(names(pm_flat(list(a = 1, b = 2), list(a = NULL))), c("a", "b"))
})
