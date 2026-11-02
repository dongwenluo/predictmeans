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
  # Results from iterative mixed-model fits (lme4, nlme, glmmTMB, permutation
  # refits) depend on those packages' optimizers, which change between
  # versions (e.g. lme4 2.x moves a glmer SE by 0.5%); they are compared to
  # 1%. Everything else (lm, glm, K matrices) must match to 1e-6.
  mixed <- grepl("lmer|lme_|_lme|glmer|clin|tmb|miss_lmer|atp|semi|perm|r2_|varcomp|contr_", names(ref)) &
    !grepl("_K($|_)", names(ref)) & !grepl("^perm_lm", names(ref))
  names(mixed) <- names(ref)
  for (nm in setdiff(names(ref)[!failed_in_ref], wrong_in_ref)) {
    expect_equal(clean_reference(new[[nm]]), clean_reference(ref[[nm]]),
                 tolerance = if (mixed[[nm]]) 1e-2 else 1e-6, label = nm)
  }
})
