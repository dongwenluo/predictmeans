## Satterthwaite df for gls and lme models (R/satt-nlme.R); reference values
## from the SAS run of dev/sas/pm-reference.sas (DDFM=SATTERTHWAITE)

rep_data <- function() {
  y <- c(20.1, 24.5, 28.2, 18.4, 22.0, 25.9, 22.3, 26.1, 30.4, 19.0, 23.4, 27.1, 21.5, 25.0, 29.0, 17.8, 21.2, 24.8,
         19.5, 21.0, 22.4, 21.0, 22.8, 23.9, 18.2, 19.5, 20.8, 20.4, 21.9, 23.1, 22.1, 23.5, 24.6, 19.0, 20.1, 21.5)
  d <- data.frame(Subj = factor(rep(1:12, each = 3)), Group = factor(rep(c("TrtA", "TrtB"), each = 18)),
                  Time = factor(rep(1:3, 12)), Y = y)
  set.seed(1)
  d[sample(nrow(d)), ]                       # rows not sorted by subject: tests the ordering
}

test_that("gls with an unstructured covariance: Satterthwaite df equal SAS", {
  r <- rep_data()
  g <- nlme::gls(Y ~ Group * Time, data = r, correlation = nlme::corSymm(form = ~ as.integer(Time) | Subj),
                 weights = nlme::varIdent(form = ~ 1 | Time))
  pm <- predictmeans(g, "Group:Time", df_method = "Satterthwaite", plot = FALSE, prtnum = FALSE)
  expect_equal(pm$means$df, rep(10, 6), tolerance = 1e-3)
  expect_identical(pm$settings$df_method, "Satterthwaite")
  # pairs: SAS 10 or 10.1 / 10.2 (printed to one decimal)
  expect_true(all(pm$pairs$df > 9.99 & pm$pairs$df < 10.25))
  expect_equal(sort(round(pm$pairs$df, 1)), sort(c(rep(10, 9), 10.1, 10.1, rep(10.2, 4))))
  jt <- pm_joint(g, "Time", df_method = "Satterthwaite")
  expect_equal(c(jt$df2, jt$F.ratio), c(10, 1795.59), tolerance = 1e-4)
  expect_identical(attr(jt, "df_method"), "Satterthwaite")
  ct <- contrastmeans(g, "Time", ctrmatrix = rbind(c(-1, 1, 0)), df_method = "Satterthwaite")[[1]]
  expect_equal(unname(ct[, "df"]), 10, tolerance = 1e-3)
  expect_error(pm_satt_parts_nlme(update(g, method = "ML")), "REML")
})

test_that("Kenward-Roger for gls and lme: SAS DDFM=KR F tests, and lme equals lmer", {
  r <- rep_data()
  g <- nlme::gls(Y ~ Group * Time, data = r, correlation = nlme::corSymm(form = ~ as.integer(Time) | Subj),
                 weights = nlme::varIdent(form = ~ 1 | Time))
  for (dm in c("Kenward-Roger", "KR-SAS")) {
    pk <- predictmeans(g, "Group:Time", df_method = dm, plot = FALSE, prtnum = FALSE)
    expect_equal(pk$means$df, rep(10, 6), tolerance = 1e-3)
    expect_equal(sort(round(pk$pairs$df, 1)), sort(c(rep(10, 9), 10.1, 10.1, rep(10.2, 4))))
    jt <- rbind(pm_joint(g, "Group", df_method = dm), pm_joint(g, "Time", df_method = dm),
                pm_joint(g, "Group:Time", ctrmatrix = t(kronecker(contr.helmert(2), contr.helmert(3))), df_method = dm))
    expect_equal(jt$F.ratio, c(5.55, 1616.03, 395.30), tolerance = 1e-3)     # SAS DDFM=KR
    expect_equal(jt$df2, c(10, 9, 9), tolerance = 1e-3)
  }
  expect_equal(pk$settings$df_method, "KR-SAS")
  # lme with random effects only: the same as lmer's Kenward-Roger (pbkrtest)
  d <- as.data.frame(nlme::Oats); d$nitro <- factor(d$nitro)
  d <- droplevels(subset(d, !(Block == "II" & Variety == "Victory" & nitro == "0.6")))
  l <- nlme::lme(yield ~ Variety * nitro, random = ~ 1 | Block/Variety, data = d)
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d)
  expect_equal(unclass(pm_vcov_kr_nlme(l))[, ], unclass(pm_vcov_kr(m))[, ], tolerance = 1e-5, ignore_attr = TRUE)
  pl <- predictmeans(l, "Variety", df_method = "KR-SAS", plot = FALSE, prtnum = FALSE)
  pr <- predictmeans(m, "Variety", df_method = "KR-SAS", plot = FALSE, prtnum = FALSE)
  expect_equal(pl$pairs[, c("std.error", "df")], pr$pairs[, c("std.error", "df")], tolerance = 1e-5)
  expect_equal(pl$pairs$df, c(9.79, 9.97, 9.97), tolerance = 2e-3)                 # SAS DDFM=KR
})

test_that("lme: Satterthwaite df equal lmer's (and SAS's) for the split plot with a plot missing", {
  d <- as.data.frame(nlme::Oats); d$nitro <- factor(d$nitro)
  d <- droplevels(subset(d, !(Block == "II" & Variety == "Victory" & nitro == "0.6")))
  l <- nlme::lme(yield ~ Variety * nitro, random = ~ 1 | Block/Variety, data = d)
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d)
  pl <- predictmeans(l, "Variety:nitro", df_method = "Satterthwaite", plot = FALSE, prtnum = FALSE)
  pr <- predictmeans(m, "Variety:nitro", df_method = "Satterthwaite", plot = FALSE, prtnum = FALSE)
  expect_equal(pl$pairs$df, pr$pairs$df, tolerance = 1e-3)
  expect_equal(range(pl$pairs$df), c(31.3, 44.7), tolerance = 2e-3)        # SAS 31.3 ... 44.7
  # with a within-subject correlation and unequal variances: the rebuilt covariance is checked inside
  r <- rep_data()
  la <- nlme::lme(Y ~ Group * Time, random = ~ 1 | Subj, data = r,
                  correlation = nlme::corAR1(form = ~ as.integer(Time) | Subj),
                  weights = nlme::varIdent(form = ~ 1 | Time))
  P <- pm_satt_parts_nlme(la)
  expect_s3_class(P, "pm_satt")
  df <- pm_satt_ddf(diag(6)[2, , drop = FALSE], P)
  expect_true(is.finite(df) && df > 0)
})
