## df_method = "between-within" (SAS DDFM=BETWITHIN); reference values from
## the SAS PROC MIXED run of dev/sas/pm-reference.sas (dev/sas/sas-results.csv)

oats_f <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}
rep_data <- function() {
  y <- c(20.1, 24.5, 28.2, 18.4, 22.0, 25.9, 22.3, 26.1, 30.4, 19.0, 23.4, 27.1, 21.5, 25.0, 29.0, 17.8, 21.2, 24.8,
         19.5, 21.0, 22.4, 21.0, 22.8, 23.9, 18.2, 19.5, 20.8, 20.4, 21.9, 23.1, 22.1, 23.5, 24.6, 19.0, 20.1, 21.5)
  data.frame(Subj = factor(rep(1:12, each = 3)), Group = factor(rep(c("TrtA", "TrtB"), each = 18)),
             Time = factor(rep(1:3, 12)), Y = y)
}

test_that("pm_betwithin_df() equals SAS DDFM=BETWITHIN (split plot, one plot missing, unstructured)", {
  d <- oats_f()
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d)
  expect_equal(c(pm_betwithin_df(m)), c("(Intercept)" = 5, Variety = 10, nitro = 15, "Variety:nitro" = 30))
  d71 <- droplevels(subset(d, !(Block == "II" & Variety == "Victory" & nitro == "0.6")))
  m71 <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d71)
  expect_equal(unname(pm_betwithin_df(m71)[-1]), c(10, 15, 29))
  l71 <- nlme::lme(yield ~ Variety * nitro, random = ~ 1 | Block/Variety, data = d71)
  expect_equal(unname(pm_betwithin_df(l71)[-1]), c(10, 15, 29))
  r <- rep_data()
  g <- nlme::gls(Y ~ Group * Time, data = r, correlation = nlme::corSymm(form = ~ as.integer(Time) | Subj),
                 weights = nlme::varIdent(form = ~ 1 | Time))
  expect_equal(as.vector(pm_betwithin_df(g)), rep(10, 4))           # SAS: unstructured, all between-subject
  # other covariance structures: within-subject df for Time (classical repeated measures)
  ga <- nlme::gls(Y ~ Group * Time, data = r, correlation = nlme::corAR1(form = ~ as.integer(Time) | Subj))
  expect_equal(as.vector(pm_betwithin_df(ga)), c(10, 10, 20, 20))
  expect_equal(as.vector(pm_betwithin_df(lme4::lmer(Y ~ Group * Time + (1 | Subj), data = r))), c(10, 10, 20, 20))
})

test_that("predictmeans(), contrastmeans() and pm_joint() use the between-within df and the model's covariance", {
  d <- oats_f()
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d)
  pm <- predictmeans(m, "nitro", df_method = "between-within", plot = FALSE, prtnum = FALSE)
  expect_equal(unique(pm$means$df), 15)
  expect_equal(unique(pm$pairs$df), 15)
  expect_equal(pm$means$std.error, rep(7.1747, 4), tolerance = 1e-4)          # SAS (model-based SE)
  expect_equal(pm$pairs$std.error, rep(4.4358, 6), tolerance = 1e-4)
  expect_identical(pm$settings$df_method, "between-within")
  expect_false(pm$settings$vcov_adjust)
  ct <- contrastmeans(m, "Variety:nitro", ctrmatrix = rbind(c(-1, 1, rep(0, 10))), df_method = "between-within")[[1]]
  expect_equal(unname(ct[, "df"]), 30)
  jt <- pm_joint(m, "Variety", df_method = "between-within")
  expect_equal(c(jt$df1, jt$df2), c(2, 10))
  expect_equal(jt$F.ratio, 1.49, tolerance = 0.01)                             # SAS Type 3 F
  expect_identical(attr(jt, "df_method"), "between-within")
  # repeated measures with gls: SAS gives 10 for every effect
  r <- rep_data()
  g <- nlme::gls(Y ~ Group * Time, data = r, correlation = nlme::corSymm(form = ~ as.integer(Time) | Subj),
                 weights = nlme::varIdent(form = ~ 1 | Time))
  pg <- predictmeans(g, "Group:Time", df_method = "between-within", plot = FALSE, prtnum = FALSE)
  expect_equal(unique(pg$means$df), 10)
  expect_equal(unique(pg$pairs$df), 10)
  # crossed random effects have no subject; other classes ignore it
  d$Plot <- factor(rep(1:12, 6))
  mc <- lme4::lmer(yield ~ nitro + (1 | Block) + (1 | Plot), data = d)
  expect_error(predictmeans(mc, "nitro", df_method = "between-within", plot = FALSE, prtnum = FALSE), "subject")
  expect_warning(predictmeans(lm(yield ~ nitro, data = d), "nitro", df_method = "between-within",
                              plot = FALSE, prtnum = FALSE), "lmer, lme and gls")
})
