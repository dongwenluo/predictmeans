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

test_that("Genstat tier (i): phi and the variance ratio fixed at Genstat's estimates give its means", {
  d <- slatehall()
  f <- nlme::gls(yield ~ variety, data = d,
                 correlation = corAR1xAR1(c(0.4586, 0.6838), form = ~ fieldrow + fieldcolumn, fixed = TRUE))
  m <- genstat("means", "A")
  b <- coef(f)
  est <- c(b[1], b[1] + b[-1])
  expect_equal(unname(est[m$variety]), m$mean, tolerance = 5e-5)     # (Genstat's phi printed to 4 decimals)
})

test_that("review: positions, transformations, boundary df, na.exclude, registration", {
  skip_on_cran()
  d <- slatehall()
  # two trials stacked without '| trial': a clear stop, not a wrong model
  d2 <- rbind(transform(d, trial = "a"), transform(d, trial = "b", yield = yield + rnorm(150)))
  expect_error(nlme::gls(yield ~ variety, data = d2,
                         correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn, nugget = TRUE)), "\\| trial")
  f2 <- nlme::gls(yield ~ variety, data = d2, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn | trial))
  expect_s3_class(f2, "gls")
  expect_error(nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ I(3 * fieldrow) + fieldcolumn)),
               "no transformations")
  # lme with replicates at the boundary: Kenward-Roger leaves that variance out, so the df are the gls ones
  fA <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn | rep))
  fL <- nlme::lme(yield ~ variety, random = ~ 1 | rep, data = d,
                  correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
  rA <- predictmeans(fA, "variety")
  rL <- suppressMessages(predictmeans(fL, "variety"))
  expect_equal(rL$settings$df_method, "Kenward-Roger")
  expect_equal(rL$pairs$df, rA$pairs$df, tolerance = 0.02)
  expect_equal(rL$means$df, rA$means$df, tolerance = 0.02)
  # a nugget at the boundary (data without one): Kenward-Roger still works
  set.seed(3)
  sim <- expand.grid(fieldrow = 1:10, fieldcolumn = 1:15)
  sim$variety <- d$variety[order(d$fieldrow, d$fieldcolumn)]
  Rr <- 0.5^abs(outer(1:10, 1:10, "-")); Rc <- 0.7^abs(outer(1:15, 1:15, "-"))
  sim$yield <- as.vector(t(chol(kronecker(Rc, Rr))) %*% rnorm(150))
  fn <- nlme::gls(yield ~ variety, data = sim, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn, nugget = TRUE))
  rn <- predictmeans(fn, "variety")
  expect_equal(rn$settings$df_method, "Kenward-Roger")
  expect_true(all(rn$pairs$df < 150))
  # na.exclude: as the fit without the row
  dn <- d
  dn$yield[5] <- NA
  fe <- nlme::gls(yield ~ variety, data = dn, na.action = na.exclude,
                  correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
  fo <- nlme::gls(yield ~ variety, data = d[-5, ], correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
  expect_equal(predictmeans(fe, "variety")$means$df, predictmeans(fo, "variety")$means$df, tolerance = 1e-6)
  # the coef<- method is registered (the installed package uses it, not nlme's default)
  expect_false(is.null(utils::getS3method("coef<-", "corAR1xAR1", optional = TRUE, envir = asNamespace("nlme"))))
})

test_that("a simulated field with known phi is recovered", {
  skip_on_cran()
  set.seed(11)
  g <- expand.grid(row = 1:20, col = 1:20)
  Rr <- 0.6^abs(outer(1:20, 1:20, "-")); Rc <- 0.4^abs(outer(1:20, 1:20, "-"))
  L <- t(chol(kronecker(Rc, Rr)))
  est <- t(replicate(4, {
    g$y <- as.vector(L %*% rnorm(400))
    f <- nlme::gls(y ~ 1, data = g, correlation = corAR1xAR1(form = ~ row + col))
    coef(f$modelStruct$corStruct, unconstrained = FALSE)
  }))
  expect_equal(unname(colMeans(est)), c(0.6, 0.4), tolerance = 0.1)
})

test_that("normalized residuals: right with grouped correlations on unsorted data, and independent of row order", {
  skip_on_cran()
  d <- slatehall()
  d <- d[order(d$fieldrow, d$fieldcolumn), ]                  # not sorted by replicate
  d$t <- ave(seq_len(nrow(d)), d$rep, FUN = seq_along)
  ds <- d[order(d$rep, d$t), ]
  k <- match(rownames(d), rownames(ds))
  sr <- function(m) as.vector(predictmeans:::pm_std_resid(m))
  # nlme's own normalized residuals are right only for data sorted by the groups
  f <- nlme::gls(yield ~ variety, data = d, correlation = nlme::corAR1(form = ~ t | rep))
  fs <- nlme::gls(yield ~ variety, data = ds, correlation = nlme::corAR1(form = ~ t | rep))
  expect_equal(unname(sr(f)), as.vector(residuals(fs, type = "normalized"))[k], tolerance = 1e-10)
  l <- nlme::lme(yield ~ variety, random = ~ 1 | rep, data = d, correlation = nlme::corAR1(form = ~ t))
  ls <- nlme::lme(yield ~ variety, random = ~ 1 | rep, data = ds, correlation = nlme::corAR1(form = ~ t))
  expect_equal(unname(sr(l)), as.vector(residuals(ls, type = "normalized"))[k], tolerance = 1e-10)
  # corAR1xAR1: a plot's residual does not depend on the order of the rows
  set.seed(1)
  sh <- d[sample(nrow(d)), ]
  a <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn | rep))
  b <- nlme::gls(yield ~ variety, data = sh, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn | rep))
  key <- function(x) paste(x$fieldrow, x$fieldcolumn)
  expect_equal(unname(sr(b)[match(key(d), key(sh))]), unname(sr(a)), tolerance = 1e-6)
})
