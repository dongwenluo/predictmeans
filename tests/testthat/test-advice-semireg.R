## Advice rules of semireg() and semipred() (R/advice.R)

ids_of <- function(x) pm_advice_get(x)$id

sr_data <- function() {
  set.seed(3318)
  d <- data.frame(x = runif(120, 0, 10), g = factor(rep(1:6, 20)))
  d$y <- sin(d$x) + rnorm(120, 0, 0.3)
  d$z <- d$x + rnorm(120)
  d
}

test_that("semireg(): a smooth estimated as a straight line is advised, a curve is not", {
  d <- sr_data()
  curved <- suppressMessages(semireg(y ~ s(x) + (1 | g), data = d))
  expect_length(ids_of(curved), 0)
  straight <- suppressMessages(semireg(z ~ s(x), data = d))
  expect_identical(ids_of(straight), "smooth_linear")
  expect_true(any(grepl("^Advice:", capture.output(print(straight)))))
})

test_that("semipred(): extrapolation and pointwise regions", {
  skip_on_cran()
  d <- sr_data()
  fit <- suppressMessages(semireg(y ~ s(x) + (1 | g), data = d))
  expect_length(ids_of(semipred(fit, "x", "x")), 0)
  out <- semipred(fit, "x", "x", at = seq(-2, 12, length = 30))
  expect_identical(ids_of(out), "extrapolation")
  expect_true(any(grepl("^Advice:", capture.output(print(out)))))
  expect_identical(ids_of(semipred(fit, "x", "x", deriv = 1)), "pointwise_regions")
  expect_length(ids_of(semipred(fit, "x", "x", deriv = 1, band = "simultaneous")), 0)
})

test_that("semireg means of aggregated binomial data: df_cases recorded, trials advice, no singular-fit noise", {
  skip_on_cran()
  data(cbpp, package = "lme4", envir = environment())
  cbpp$x <- as.numeric(cbpp$period) + seq(0, 0.5, length.out = nrow(cbpp))
  f <- suppressMessages(semireg(cbind(incidence, size - incidence) ~ period + s(x) + (1 | herd), data = cbpp,
                                family = binomial))
  p <- suppressMessages(predictmeans(f, "period"))
  expect_identical(as.vector(p$settings$df_cases), "trials")
  expect_false("singular_fit" %in% advice(p)$id)
})

test_that("summary() of semireg means works (the means carry 'estimable')", {
  set.seed(1)
  d <- data.frame(x = runif(60), g = factor(rep(1:3, 20)))
  d$y <- sin(3 * d$x) + as.numeric(d$g) / 2 + rnorm(60, 0, 0.3)
  f <- suppressMessages(semireg(y ~ g + s(x), data = d))
  p <- predictmeans(f, "g")
  expect_true(all(p$means$estimable))
  expect_output(print(summary(p)), "Predicted means")
})

test_that("semireg means: pairs get the Kenward-Roger SEs and df of the mixed model; vcov_adjust = FALSE does not", {
  skip_on_cran()
  skip_if_not_installed("pbkrtest")
  set.seed(3)
  d <- expand.grid(g = factor(1:4), rep = 1:10)
  d$blk <- factor(d$rep)
  d$x <- runif(nrow(d))
  d$y <- as.numeric(d$g) / 3 + 2 * sin(6 * d$x) + rnorm(10)[d$blk] + rnorm(nrow(d), 0, 0.4)
  f <- suppressMessages(semireg(y ~ g + s(x) + (1 | blk), data = d))
  vc <- as.data.frame(lme4::VarCorr(f$semer))
  expect_gt(max(vc$sdcor[grepl("^s", vc$grp)]), 0.01)               # a fitted curve, not a straight line
  # outside reference: pbkrtest's adjusted covariance of the mixed model (expected information)
  expect_equal(as.matrix(pm_vcov_kr(f$semer, "expected")), as.matrix(pbkrtest::vcovAdj(methods::as(f$semer, "lmerMod"))),
               tolerance = 1e-6, ignore_attr = TRUE)
  p <- predictmeans(f, "g")
  pn <- predictmeans(f$semer, "g")                         # the mixed model: KR by default
  k <- match(paste(p$pairs$level1, p$pairs$level2), paste(pn$pairs$level1, pn$pairs$level2))
  expect_equal(p$pairs$std.error, pn$pairs$std.error[k], tolerance = 1e-8)
  expect_equal(p$pairs$df, pn$pairs$df[k], tolerance = 1e-8)
  expect_true(p$settings$vcov_adjust)
  expect_match(p$settings$df_method, "^Kenward-Roger \\(pairs")
  expect_output(print(summary(p)), "Kenward-Roger adjusted SEs \\(pairs of fixed effects\\)\\)")
  p0 <- predictmeans(f, "g", vcov_adjust = FALSE)
  expect_false(p0$settings$vcov_adjust)
  expect_true(all(p0$pairs$std.error < p$pairs$std.error))
  expect_equal(p0$means$std.error, p$means$std.error)        # the means keep Satterthwaite
  # a smooth for each level of g: pairs within g are fixed-effect contrasts (KR), the others are not
  d$h <- factor(d$rep %% 2)
  fb <- suppressMessages(semireg(y ~ g * h + s(x, by = g) + (1 | blk), data = d))
  pb <- predictmeans(fb, "g:h", atvar = "g")
  pbn <- suppressWarnings(predictmeans(fb$semer, "g:h", atvar = "g"))
  w <- which(!is.na(pb$pairs$by))
  kb <- match(paste(pb$pairs$level1[w], pb$pairs$level2[w]), paste(pbn$pairs$level1, pbn$pairs$level2))
  expect_false(anyNA(kb))
  expect_equal(pb$pairs$std.error[w], pbn$pairs$std.error[kb], tolerance = 1e-8)
  expect_match(pb$settings$df_method, "of fixed effects only")
  expect_error(contrastmeans(f, "g", ctrmatrix = rbind(c(1, -1, 0, 0))), "does not handle semireg")
  # cells without data: NA with a message (were NaN)
  d$h2 <- factor(rep(1:2, length.out = nrow(d)))                 # g and h2 confounded: half the cells empty
  fe <- suppressMessages(semireg(y ~ g * h2 + s(x) + (1 | blk), data = d))
  expect_message(pe <- predictmeans(fe, "g:h2"), "cells without data")
  expect_false(any(is.nan(pe$means$estimate)))
  expect_true(all(is.na(pe$means$estimate[!pe$means$estimable])))
})
