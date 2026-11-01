test_that("results match predictmeans 1.1.1 for everything 1.1.1 supported", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  ref <- readRDS(test_path("_data", "reference-v1.1.1.rds"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  new <- reference_cases()
  expect_setequal(names(new), names(ref))
  # Cases where 1.1.1 itself stopped with an error have nothing to reproduce;
  # their fixes are tested in test-fixes.R.
  failed_in_ref <- vapply(ref, function(r) is.character(r) && length(r) == 1 && grepl("^ERROR", r), logical(1))
  expect_identical(names(ref)[failed_in_ref], c("r2_lmer_od", "r2_tmb"))
  # 1.1.1 refitted this lme model on nlme::Oats instead of the local copy
  # (bug B12), so its frozen result is wrong; see test-fixes.R.
  wrong_in_ref <- "perm_lme"
  for (nm in setdiff(names(ref)[!failed_in_ref], wrong_in_ref)) {
    expect_equal(clean_reference(new[[nm]]), clean_reference(ref[[nm]]),
                 tolerance = 1e-6, label = nm)
  }
})
