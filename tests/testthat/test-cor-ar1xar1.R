## corAR1xAR1() (R/cor-ar1xar1.R, dev/review-ar1ar1.md). References: nlme's
## corAR1 (one field row) and the Genstat run of Slate Hall
## (dev/genstat/genstat-results-2026-11-26-ar1ar1-slatehall.txt; data in
## _data/ar1xar1/).

slatehall <- function() {
  d <- read.csv(test_path("_data", "variogram", "slatehall.csv"))
  d$variety <- factor(d$variety)
  d$rep <- factor(d$replicates)
  d
}
genstat <- function(what, lab) read.csv(test_path("_data", "ar1xar1", paste0(what, lab, ".csv")))

test_that("one field row: corAR1xAR1 equals nlme's corAR1", {
  d <- slatehall()
  d1 <- d[d$fieldrow == 1, ]
  a <- nlme::gls(yield ~ fieldcolumn, data = d1, correlation = nlme::corAR1(0.6, form = ~ fieldcolumn, fixed = TRUE))
  b <- nlme::gls(yield ~ fieldcolumn, data = d1,
                 correlation = corAR1xAR1(c(0.3, 0.6), form = ~ fieldrow + fieldcolumn, fixed = TRUE))
  expect_equal(as.numeric(logLik(b)), as.numeric(logLik(a)), tolerance = 1e-10)
  expect_equal(coef(b), coef(a), tolerance = 1e-10)
  expect_length(coef(b$modelStruct$corStruct), 0)                    # fixed: nothing to estimate
})

test_that("Slate Hall AR1 x AR1 (A) and with a nugget (B) equal Genstat", {
  skip_on_cran()
  d <- slatehall()
  fA <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
  fB <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn, nugget = TRUE))
  pA <- coef(fA$modelStruct$corStruct, unconstrained = FALSE)
  pB <- coef(fB$modelStruct$corStruct, unconstrained = FALSE)
  # variance parameters to Genstat's printed digits (phi to 1e-3 for A, 1e-2 for B: the nugget ridge)
  expect_equal(unname(pA), c(0.4586, 0.6838), tolerance = 1e-3)
  expect_equal(fA$sigma^2, 3.876, tolerance = 1e-3)
  expect_equal(unname(pB[1:2]), c(0.6827, 0.8438), tolerance = 1e-2)
  expect_equal(fB$sigma^2 * pB[["nugget"]], 0.486, tolerance = 2e-3)                    # units
  expect_equal(fB$sigma^2 * (1 - pB[["nugget"]]) / (fB$sigma^2 * pB[["nugget"]]), 9.420, tolerance = 2e-3)
  expect_equal(-2 * (as.numeric(logLik(fA)) - as.numeric(logLik(fB))), 249.35 - 242.35, tolerance = 1e-3)
  # the 25 means and all 300 SEDs (unadjusted, as Genstat's VPREDICT)
  for (lab in c("A", "B")) {
    f <- if (lab == "A") fA else fB
    r <- predictmeans(f, "variety", vcov_adjust = FALSE)
    m <- genstat("means", lab)
    expect_equal(r$means$estimate[match(m$variety, r$means$variety)], m$mean, tolerance = 1e-5)
    s <- genstat("sed", lab)
    key <- paste(pmin(s$variety1, s$variety2), pmax(s$variety1, s$variety2))
    a <- as.integer(sub(":.*", "", r$pairs$level1))
    b <- as.integer(sub(":.*", "", r$pairs$level2))
    expect_equal(r$pairs$std.error[match(key, paste(pmin(a, b), pmax(a, b)))], s$sed, tolerance = 1e-5)
  }
  # Wald statistic and d.d.f. of the variety test (Genstat: 313.04, 80.0; 245.39, 75.7 with KR-average)
  jA <- pm_joint(fA, "variety", df_method = "KR-average")
  jB <- pm_joint(fB, "variety", df_method = "KR-average")
  expect_equal(jA$df2, 80.0, tolerance = 1e-3)
  expect_equal(jB$df2, 75.7, tolerance = 1e-3)
  expect_equal(as.numeric(anova(fA, type = "marginal")["variety", "F-value"]) * 24, 313.04, tolerance = 1e-4)
})

test_that("lme: replicates in one group keep the AR1 x AR1 over the field; at the boundary it is the gls fit", {
  skip_on_cran()
  d <- slatehall()
  d$all <- factor(1)
  fA <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
  fC <- nlme::lme(yield ~ variety, random = list(all = nlme::pdIdent(~ rep - 1)), data = d,
                  correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn | all))
  # Genstat's replicate component is negative (-0.024); lme keeps it >= 0, so it reaches 0: model A
  expect_lt(as.numeric(nlme::VarCorr(fC)[1, "Variance"]), 1e-4)
  expect_equal(as.numeric(logLik(fC)), as.numeric(logLik(fA)), tolerance = 1e-6)
  expect_equal(unname(coef(fC$modelStruct$corStruct, FALSE)), unname(coef(fA$modelStruct$corStruct, FALSE)),
               tolerance = 1e-3)
})

test_that("missing plots, groups, checks", {
  d <- slatehall()
  dm <- d[-c(5, 20, 77, 130), ]
  cs <- nlme::Initialize(corAR1xAR1(c(0.5, 0.7), form = ~ fieldrow + fieldcolumn), dm)
  R <- nlme::corMatrix(cs)
  rc <- as.matrix(dm[, c("fieldrow", "fieldcolumn")])
  expect_equal(R, 0.5^abs(outer(rc[, 1], rc[, 1], "-")) * 0.7^abs(outer(rc[, 2], rc[, 2], "-")),
               ignore_attr = TRUE)
  cn <- nlme::Initialize(corAR1xAR1(c(0.5, 0.7, 0.2), form = ~ fieldrow + fieldcolumn, nugget = TRUE), dm)
  Rn <- nlme::corMatrix(cn)
  expect_equal(unname(diag(Rn)), rep(1, nrow(dm)))
  expect_equal(Rn[1, 2], 0.8 * R[1, 2])
  # groups: one AR1 x AR1 per replicate
  cg <- nlme::Initialize(corAR1xAR1(c(0.5, 0.7), form = ~ fieldrow + fieldcolumn | rep), d)
  expect_length(nlme::corMatrix(cg), 6L)
  expect_error(corAR1xAR1(c(1, 0.5)), "\\(-1, 1\\)")
  d$rf <- factor(d$fieldrow)
  expect_error(nlme::Initialize(corAR1xAR1(form = ~ rf + fieldcolumn), d), "integer positions")
})

test_that("residvariogram(): the AR1 x AR1 fit leaves no spatial pattern; a model without one is pointed to it", {
  skip_on_cran()
  d <- slatehall()
  fit <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
  v <- residvariogram(fit, row = "fieldrow", col = "fieldcolumn", data = d, nsim = 19, plot = FALSE)
  expect_false(any(grepl("^variogram_", advice(v)$id)))
  v0 <- residvariogram(lm(yield ~ variety, data = d), row = "fieldrow", col = "fieldcolumn", data = d, plot = FALSE)
  expect_match(advice(v0)$text[advice(v0)$id == "variogram_spatial"], "corAR1xAR1\\(form = ~ row \\+ col\\).*t2\\(row, col\\)")
})
