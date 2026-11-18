## Long plot titles are broken into lines (pm_wrap_title(), R/plot-title.R).

test_that("pm_wrap_title() breaks long titles at spaces and keeps the line breaks given", {
  long <- "Predicted means for \"nitrogen_fertiliser_rate\" by \"oat_variety_grown\" with Aveg.LSD (5%) Bar"
  w <- pm_wrap_title(long, 12)
  lines <- strsplit(w, "\n")[[1]]
  expect_gt(length(lines), 1)
  expect_true(all(nchar(lines) <= floor(900/(12*1.2))))
  expect_identical(gsub("\n", " ", w), long)                    # nothing lost
  expect_identical(pm_wrap_title("short"), "short")
  expect_identical(pm_wrap_title("a\nb\n"), "a\nb\n")           # given breaks, also the one at the end
  expect_null(pm_wrap_title(NULL))
  expect_identical(pm_wrap_title(""), "")
  expect_gt(nchar(pm_wrap_title(long, 12, width = 200)), 0)
  expect_false(grepl("\n", pm_wrap_title(long, 12, width = 200)))
  # larger fonts and panels side by side: shorter lines
  expect_gt(length(strsplit(pm_wrap_title(long, 20), "\n")[[1]]), length(lines))
  expect_gt(length(strsplit(pm_wrap_title(long, 12, panels = 2), "\n")[[1]]), length(lines))
  op <- options(predictmeans.title_width = 30); on.exit(options(op))
  expect_true(all(nchar(strsplit(pm_wrap_title(long), "\n")[[1]]) <= 30))
})

test_that("the plots of the package wrap their titles", {
  O <- as.data.frame(nlme::Oats)
  names(O)[names(O) == "Variety"] <- "oat_variety_grown_in_the_field_trial"
  O$nitrogen_fertiliser_rate_applied <- factor(O$nitro)
  fm <- lm(yield ~ Block + nitrogen_fertiliser_rate_applied * oat_variety_grown_in_the_field_trial, data = O)
  term <- "nitrogen_fertiliser_rate_applied:oat_variety_grown_in_the_field_trial"
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  pm <- predictmeans(fm, term, adj = "tukey", plot = c("mean", "pval"))
  multi <- function(g) grepl("\n.", g$labels$title)
  expect_true(multi(pm$predictmeansPlot$meanPlot))
  expect_true(multi(ggplot2::autoplot(pm)))
  expect_true(multi(ggplot2::autoplot(pm, type = "pvalue")))
  expect_true(multi(PMplot(pm$p_valueMatrix, mtitle = paste(rep("a long title for the level plot", 4), collapse = " "))))
  cm <- covariatemeans(lm(yield ~ nitro * oat_variety_grown_in_the_field_trial, data = O), "oat_variety_grown_in_the_field_trial",
                       covariate = "nitro",
                       plot_args = list(mtitle = paste(rep("a long title for the fitted lines", 4), collapse = " ")))
  expect_true(multi(cm$plt))
})
