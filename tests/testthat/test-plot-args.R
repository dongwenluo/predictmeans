## Plots: predictmeans() draws none (2.0); the 1.x and earlier-2.0 plot
## arguments stop with the autoplot() call that replaces them.

pa_fit <- function() lm(yield ~ Block + nitro * Variety, data = transform(as.data.frame(nlme::Oats), nitro = factor(nitro)))

test_that("plot and plot_args stop with their autoplot() spelling; plot = FALSE is accepted", {
  fm <- pa_fit()
  expect_error(predictmeans(fm, "nitro", plot = TRUE), "'plot' was removed in predictmeans 2.0: use autoplot(x)",
               fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plot = "bar"), "\"bar\" (bar)", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plot_args = list(basesz = 14)), "'plot_args' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plottitle = "one"), "+ labs(title, x, y)", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plott = "one"), "'plottitle' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", bar = TRUE), "'barplot' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", which.plots = "bar"), "'which.plots' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", prtplt = FALSE), "'prtplt' was removed", fixed = TRUE)
  expect_error(predictmeans(fm, "nitro", plotord = c(2, 1)), "x, group, panel", fixed = TRUE)
  expect_equal(predictmeans(fm, "nitro", plot = FALSE), predictmeans(fm, "nitro"))
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  expect_false(any(c("predictmeansPlot", "predictmeansBarPlot") %in% names(pm)))
  expect_true("p_valueMatrix" %in% names(pm))
  # nothing is drawn
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev))
  rec <- grDevices::recordPlot()
  invisible(predictmeans(fm, "nitro"))
  expect_identical(length(grDevices::recordPlot()[[1]]), length(rec[[1]]))
})

test_that("unknown names are errors", {
  fm <- pa_fit()
  expect_error(predictmeans(fm, "nitro", colour = "red"), "predictmeans(): unused argument: colour", fixed = TRUE)
  expect_identical(names(formals(predictmeans))[1:7],
                   c("model", "modelterm", "data", "pairwise", "atvar", "adj", "conf_level"))
  expect_false(any(c("Df", "level", "plot", "plot_args", "which.plots", "mplot", "barplot", "pplot", "bkplot",
                     "plottitle", "newwd", "prtnum", "ndecimal", "transOff", "count") %in% names(formals(predictmeans))))
})

test_that("the stages are functions that return named lists; a later stage replaces a value", {
  fm <- pa_fit()
  op <- options(scipen = 6); on.exit(options(op))             # (as predictmeans() sets for its tables)
  opt <- pm_setup(list(model = fm, modelterm = "nitro:Variety", data = NULL, pairwise = NULL, atvar = "Variety",
                       adj = "tukey", Df = NULL, level = 0.05, covariate = NULL, meandecr = NULL, letter_ci = FALSE,
                       trans = I, trans_default = TRUE, permlist = NULL, ncore = 1, cov_reduce = mean,
                       exclude = NULL, vcov_adjust = TRUE, vcov_adjust_default = TRUE, weights = "equal",
                       df_method = "auto", vcov = NULL, equivalence = NULL, bias_adjust = FALSE, sigma = NULL,
                       average = "link"))
  expect_true(opt$pairwise && opt$compare)                    # 'adj' and 'atvar' ask for comparisons
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
  out <- pm_assemble(v)
  ref <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  expect_equal(out, ref)
  # a NULL from a later stage replaces an earlier value
  expect_null(pm_flat(list(a = 1, b = 2), list(a = NULL))$a)
  expect_identical(names(pm_flat(list(a = 1, b = 2), list(a = NULL))), c("a", "b"))
})
