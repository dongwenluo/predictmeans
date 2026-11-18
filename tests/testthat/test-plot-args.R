## Plot settings: plot_args (the 1.x names as arguments stop), and the stages.

pa_fit <- function() lm(yield ~ Block + nitro * Variety, data = transform(as.data.frame(nlme::Oats), nitro = factor(nitro)))

test_that("plot_args sets the plots; the 1.x argument names stop with their 2.0 spelling", {
  fm <- pa_fit()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  a <- predictmeans(fm, "nitro", plot = "mean",
                    plot_args = list(plottitle = "Yield by nitrogen", plotxlab = "N", basesz = 14, lsd_bar = FALSE))
  ga <- a$predictmeansPlot$meanPlot
  expect_identical(ga$labels$title, "Yield by nitrogen\n")
  expect_identical(ga$labels$x, "N")
  # abbreviated names in plot_args, as before
  cc <- predictmeans(fm, "nitro", plot = "mean", plot_args = list(plott = "Yield by nitrogen"))
  expect_identical(cc$predictmeansPlot$meanPlot$labels$title, "Yield by nitrogen\n")
  # 1.x: settings as arguments of their own, also abbreviated
  expect_error(predictmeans(fm, "nitro", plottitle = "one"), "'plottitle' was removed in predictmeans 2.0: use plot_args = list(plottitle", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plott = "one"), "'plottitle' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", bar = TRUE), "'barplot' was removed in predictmeans 2.0: use plot = c(", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", which.plots = "bar"), "'which.plots' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", prtplt = FALSE), "'prtplt' was removed", fixed = TRUE)
  # an order for the plot asks for the mean plot
  e <- predictmeans(fm, "nitro:Variety", plot = TRUE, plot_args = list(plotord = c(2, 1)))
  expect_s3_class(e$predictmeansPlot$meanPlot, "ggplot")
  # the plots wanted by name
  r <- predictmeans(fm, "nitro", plot = c("mean", "bar"))
  expect_s3_class(r$predictmeansBarPlot, "ggplot")
})

test_that("unknown names are errors", {
  fm <- pa_fit()
  expect_error(predictmeans(fm, "nitro", colour = "red"), "predictmeans(): unused argument: colour", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plot_args = list(colour = "red")), "'plot_args' can hold")
  expect_error(predictmeans(fm, "nitro", plot_args = list(newwd = TRUE)), "'plot_args' can hold")
  expect_error(predictmeans(fm, "nitro", plot_args = "big"), "named list")
  expect_error(predictmeans(fm, "nitro", plot_args = list(14)), "named list")
  expect_identical(names(formals(predictmeans))[1:7],
                   c("model", "modelterm", "data", "pairwise", "atvar", "adj", "conf.level"))
  expect_false(any(c("Df", "level", "which.plots", "mplot", "barplot", "pplot", "bkplot", "plottitle", "newwd", "prtnum",
                     "ndecimal", "transOff", "count") %in% names(formals(predictmeans))))
  expect_identical(names(pm_plot_args()), names(pm_plot_defaults))
})

test_that("the stages are functions that return named lists; a later stage replaces a value", {
  fm <- pa_fit()
  pa <- pm_plot_args()
  op <- options(scipen = 6); on.exit(options(op))             # (as predictmeans() sets for its tables)
  opt <- pm_setup(list(model = fm, modelterm = "nitro:Variety", data = NULL, pairwise = FALSE, atvar = "Variety",
                       adj = "tukey", Df = NULL, level = 0.05, covariate = NULL, meandecr = NULL, letterCI = FALSE,
                       trans = I, trans_default = TRUE, which.plots = c("mean", "pval", "back"),
                       plot = FALSE, permlist = NULL, ncore = 1, cov_reduce = mean,
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
  ref <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  expect_equal(out, ref)
  # a NULL from a later stage replaces an earlier value
  expect_null(pm_flat(list(a = 1, b = 2), list(a = NULL))$a)
  expect_identical(names(pm_flat(list(a = 1, b = 2), list(a = NULL))), c("a", "b"))
})
