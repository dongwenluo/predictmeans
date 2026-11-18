## Comparison intervals (R/comparison.R): intervals whose overlap answers the
## pairwise tests; letterCI reads its letters from them.

cmp_oats <- function() transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
# do the comparison intervals of each pair overlap?
cmp_no_overlap <- function(r) {
  m <- r$means; p <- r$pairs
  if (!is.null(p$by)) p <- p[!is.na(p$by), ]
  key <- do.call(paste, c(lapply(m[unlist(strsplit(r$settings$modelterm, ":"))], as.character), sep = ":"))
  i <- match(p$level1, key); j <- match(p$level2, key)
  list(no = m$comparison.low[i] > m$comparison.high[j] | m$comparison.low[j] > m$comparison.high[i], p = p)
}

test_that("independent means with equal SEs: the 83.4% interval of the literature", {
  O <- cmp_oats()
  r <- predictmeans(lm(yield ~ nitro, O), "nitro", pairwise = TRUE)
  cm <- r$settings$comparison
  df <- 68
  expect_equal(cm$mult, qt(0.975, df)/sqrt(2), tolerance = 1e-10)
  expect_true(cm$exact)
  # z version: 83.42% (Austin & Hux 2002; Maghsoodloo & Huang 2010); here with t on 68 df
  expect_equal(1 - 2*pnorm(-qnorm(0.975)/sqrt(2)), 0.8342, tolerance = 1e-4)
  expect_equal(cm$conf.level, 1 - 2*pt(-qt(0.975, df)/sqrt(2), df), tolerance = 1e-10)
  expect_equal(r$means$comparison.high - r$means$estimate, cm$mult*r$means$std.error, tolerance = 1e-10)
  # half an LSD on each side of a mean: the same interval
  expect_equal(r$means$comparison.high - r$means$estimate, rep(r$pairs$lsd[1]/2, 4), tolerance = 1e-10)
})

test_that("independent means with unequal SEs: c = t sqrt(1 + K^2)/(1 + K) for each pair", {
  set.seed(8)
  d <- data.frame(g = factor(rep(c("a", "b"), c(4, 36))), y = rnorm(40))
  r <- predictmeans(lm(y ~ g, d), "g", pairwise = TRUE)
  K <- r$means$std.error[1]/r$means$std.error[2]                # 3
  expect_equal(K, 3, tolerance = 1e-10)
  expect_equal(r$settings$comparison$mult, qt(0.975, 38)*sqrt(1 + K^2)/(1 + K), tolerance = 1e-10)
  # Maghsoodloo & Huang (2010): 87.87% at K = 3 (normal quantiles)
  expect_equal(1 - 2*pnorm(-qnorm(0.975)*sqrt(1 + 9)/4), 0.87874, tolerance = 1e-5)
})

test_that("correlated means (random blocks): narrower than 83.4%, and overlap still equals the test", {
  O <- cmp_oats()
  fm <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = O)
  r <- predictmeans(fm, "nitro", pairwise = TRUE)
  cm <- r$settings$comparison
  expect_true(cm$exact)
  expect_lt(cm$conf.level, 0.6)                                  # the SED (5.1) is smaller than one SE (7.3)
  expect_equal(cm$mult, r$pairs$lsd[1]/(2*r$means$std.error[1]), tolerance = 1e-8)
  o <- cmp_no_overlap(r)
  expect_identical(o$no, o$p$p.value < 0.05)
  # 95% confidence limits overlap for a pair that differs with p < 0.001
  k <- which(o$p$level1 == "0" & o$p$level2 == "0.2")
  expect_lt(o$p$p.value[k], 0.001)
  expect_gt(r$means$conf.high[1], r$means$conf.low[2])
  # with an adjustment that has a critical value the intervals go with the adjusted tests
  rt <- predictmeans(fm, "nitro", adj = "tukey")
  ot <- cmp_no_overlap(rt)
  expect_identical(ot$no, ot$p$p.adj < 0.05)
  expect_true(rt$settings$comparison$adjusted)
  expect_gt(rt$settings$comparison$mult, cm$mult)
  # step-wise methods: the intervals of the unadjusted tests
  rh <- predictmeans(fm, "nitro", adj = "holm")
  expect_false(rh$settings$comparison$adjusted)
  expect_equal(rh$settings$comparison$mult, cm$mult)
})

test_that("unequal SEDs: one average multiplier, and the pairs judged differently are counted", {
  O <- cmp_oats()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = O)
  r <- predictmeans(fm, "nitro:Variety", pairwise = TRUE)
  cm <- r$settings$comparison
  expect_false(cm$exact)
  o <- cmp_no_overlap(r)
  expect_equal(cm$pairs, 66)
  expect_equal(cm$disagree, sum(o$no != (o$p$p.value < 0.05)))
  expect_gt(cm$disagree, 0)
  expect_match(capture.output(print(summary(r)))[grep("^Comparison intervals", capture.output(print(summary(r))))],
               "pairs are judged differently by their own test")
  # within the levels of 'atvar' the SEDs are equal: exact, from those pairs only
  ra <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  ca <- ra$settings$comparison
  expect_true(ca$exact)
  expect_equal(ca$pairs, 18)
  oa <- cmp_no_overlap(ra)
  expect_identical(oa$no, oa$p$p.adj < 0.05)
})

test_that("letterCI: letters from the comparison intervals shown in the mean table, at the test level", {
  O <- cmp_oats()
  fm <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = O)
  r <- predictmeans(fm, "nitro", letterCI = TRUE)
  tab <- r$mean_table$Table
  cm <- r$settings$comparison
  expect_equal(r$settings$level, 0.05)                          # (1.x: 0.166, the level of the intervals)
  lab <- paste0("LL(", round(100*cm$conf.level, 1), "%)")
  expect_true(lab %in% names(tab))
  expect_equal(tab[[lab]], round(r$means$comparison.low, 4))
  expect_match(r$mean_table$Note, "comparison intervals")
  # balanced: the same letters as the tests
  rp <- predictmeans(fm, "nitro", pairwise = TRUE)
  expect_identical(trimws(tab$LetterGrp), trimws(rp$mean_table$Table$LetterGrp))
  # the stable table keeps the 95% confidence limits
  expect_equal(r$means$conf.low, rp$means$conf.low)
  # within 'atvar' groups
  fs <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = O)
  ra <- predictmeans(fs, "nitro:Variety", atvar = "Variety", letterCI = TRUE)
  rb <- predictmeans(fs, "nitro:Variety", atvar = "Variety", pairwise = TRUE)
  expect_identical(trimws(ra$mean_table$Table$LetterGrp), trimws(rb$mean_table$Table$LetterGrp))
  # back-transformed limits follow
  O$ly <- log(O$yield)
  rl <- predictmeans(lm(ly ~ Block + nitro, O), "nitro", letterCI = TRUE, trans = exp)
  tl <- rl$mean_table$Table
  expect_equal(tl[[grep("^Bk_LL", names(tl))]], round(exp(rl$means$comparison.low), 4))
})

test_that("autoplot(interval = 'comparison') draws them; no pairs, no intervals", {
  O <- cmp_oats()
  fm <- lm(yield ~ Block + nitro * Variety, data = O)
  r <- predictmeans(fm, "nitro", adj = "tukey")
  g <- ggplot2::autoplot(r, interval = "comparison")
  eb <- ggplot2::layer_data(g, 3)
  expect_equal(eb$ymin, r$means$comparison.low)
  expect_match(g$labels$subtitle, "5% level \\(tukey\\) when their bars do not overlap")
  # half-LSD bars go with the adjusted p-values too
  gl <- ggplot2::layer_data(ggplot2::autoplot(r, interval = "lsd"), 3)
  expect_equal(gl$ymax - gl$ymin, rep(r$LSD[["Aveg.LSD"]], 4), tolerance = 1e-4)
  # a result without them (one mean, or the pairs not tabulated)
  none <- r
  none$means$comparison.low <- none$means$comparison.high <- NULL
  expect_error(ggplot2::autoplot(none, interval = "comparison"), "no comparison intervals")
  op <- options(predictmeans.max_means = 2)
  big <- suppressMessages(predictmeans(fm, "nitro"))
  options(op)
  expect_null(big$settings$comparison)
})
