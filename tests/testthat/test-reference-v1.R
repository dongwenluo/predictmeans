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
  # versions (e.g. lme4 2.x moves a glmer SE by 0.5%, which moves a
  # confidence limit near zero by 1.4% of the vector); they are compared to
  # 2%. Everything else (lm, glm, K matrices) must match to 1e-6.
  mixed <- grepl("lmer|lme_|_lme|glmer|clin|tmb|miss_lmer|atp|semi|perm|r2_|varcomp|contr_", names(ref)) &
    !grepl("_K($|_)", names(ref)) & !grepl("^perm_lm", names(ref))
  names(mixed) <- names(ref)
  # Since 2.0 permmodels() whitens mixed-model residuals with an
  # order-invariant square root instead of the Cholesky factor (NEWS), so the
  # permutation p-values and refits of these cases differ from 1.1.1; their
  # other parts (estimates, SEs, ANOVA tables) are still compared.
  rewhitened <- c("perm_lmer", "perm_lmer_III")
  strip_perm <- function(x) {
    x <- x[setdiff(names(x), "permlist")]
    lapply(x, function(tab) tab[, setdiff(colnames(tab), "Perm_p_value"), drop = FALSE])
  }
  # Since 2.0 contrastmeans() tables also have confidence limits; the 1.1.1
  # columns are compared.
  strip_ci <- function(x) {
    if (is.list(x) && !is.null(x[["The t tests of the specified contrasts"]])) {
      tab <- x[["The t tests of the specified contrasts"]]
      x[["The t tests of the specified contrasts"]] <-
        structure(tab[, setdiff(colnames(tab), c("lower.CL", "upper.CL")), drop = FALSE], Note = attr(tab, "Note"))
    }
    x
  }
  # Since 2.0 Type III tests do not depend on the contrasts (sum-to-zero
  # coding, as SAS Type 3); 1.1.1 tested each main effect at the first level
  # of the other factors under treatment contrasts (test-anova.R). The
  # interaction is the same either way.
  type3 <- "perm_lmer_III"
  # Since 2.0 varcomp() has one layout for every class and Satterthwaite's
  # intervals (1.1.1: profile intervals for lmer, nlme's for lme, a matrix for
  # lme); the estimates and SEs are compared.
  vc_cases <- c("varcomp_lmer", "varcomp_lme")
  vc_core <- function(x) { x <- as.matrix(x); unname(cbind(as.numeric(x[, 1]), as.numeric(x[, 2]))) }
  for (nm in setdiff(names(ref)[!failed_in_ref], wrong_in_ref)) {
    a <- strip_ci(new[[nm]]); b <- ref[[nm]]
    if (nm %in% rewhitened) { a <- strip_perm(a); b <- strip_perm(b) }
    if (nm %in% vc_cases) { a <- vc_core(a); b <- vc_core(b) }
    if (nm %in% type3) {
      a$ANOVA <- a$ANOVA[grepl(":", rownames(a$ANOVA)), , drop = FALSE]
      b$ANOVA <- b$ANOVA[grepl(":", rownames(b$ANOVA)), , drop = FALSE]
    }
    expect_equal(clean_reference(a), clean_reference(b),
                 tolerance = if (mixed[[nm]]) 2e-2 else 1e-6, label = nm)
  }
})
