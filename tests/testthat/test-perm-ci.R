## Confidence intervals of pairs and contrasts (t, and permutation-t with a permlist)

test_that("pm_perm_crit(): the k-th largest permuted |t|, per row or per family", {
  T <- rbind(1:19, (1:19) / 2)
  expect_equal(pm_perm_crit(T, 0.1), c(18, 9))            # k = floor(0.1 * 20) = 2
  expect_equal(pm_perm_crit(T, 0.1, family = c("a", "a")), c(18, 18))
  expect_equal(pm_perm_crit(T, 0.01), c(Inf, Inf))        # too few permutations
})

test_that("contrastmeans() limits: t limits, and permutation-t limits that agree with the test", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro + Variety, data = Oats)
  cn <- rbind(c(-1, 1, 0, 0), c(0, -1, 1, 0), c(0, 0, -1, 1))
  tab <- contrastmeans(fm, "nitro", ctrmatrix = cn)[[1]]
  expect_equal(unname(tab[, "lower.CL"]), round(unname(tab[, "Estimate"] - qt(0.975, tab[, "df"]) * tab[, "Std. Error"]), 4),
               tolerance = 1e-3)
  skip_on_cran()
  # refits: the permutation p-value is <= level exactly when 0 is outside the limits
  fl <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = Oats)
  r <- permmodels(fl, nperm = 99, seed = 4, prt = FALSE, ncore = 1, method = "raw")
  pt <- contrastmeans(fl, "nitro", ctrmatrix = cn, permlist = r, level = 0.1)[[1]]
  outside <- pt[, "lower.CL"] > 0 | pt[, "upper.CL"] < 0
  expect_identical(unname(outside), unname(pt[, "Permuted Pr(>|t|)"] <= 0.1))
  # simultaneous limits (max-|t|) are wider
  ps <- contrastmeans(fl, "nitro", ctrmatrix = cn, permlist = r, adj = "westfall", level = 0.1)[[1]]
  expect_true(all(ps[, "upper.CL"] - ps[, "lower.CL"] >= pt[, "upper.CL"] - pt[, "lower.CL"] - 1e-8))
})

test_that("predictmeans() pairs: conf.low, conf.high and lsd", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro + Variety, data = Oats)
  pm <- predictmeans(fm, "nitro", plot = FALSE, prtnum = FALSE)
  pr <- pm$pairs
  expect_equal(pr$conf.high - pr$estimate, qt(0.975, pr$df) * pr$std.error)
  expect_equal(pr$lsd, qt(0.975, pr$df) * pr$std.error)
  r <- permmodels(fm, nperm = 199, seed = 1, prt = FALSE)
  pp <- predictmeans(fm, "nitro", permlist = r, plot = FALSE, prtnum = FALSE, ncore = 1)$pairs
  expect_true(all(is.finite(pp$lsd)))
  expect_equal(pp$conf.high - pp$estimate, pp$lsd)
  # with permutation Tukey: simultaneous limits with one critical value for
  # the family (adj.conf.*); conf.low, conf.high and lsd do not depend on adj
  pk <- predictmeans(fm, "nitro", permlist = r, adj = "tukey", plot = FALSE, prtnum = FALSE, ncore = 1)$pairs
  crit <- (pk$adj.conf.high - pk$estimate) / pk$std.error
  expect_equal(max(crit) - min(crit), 0, tolerance = 1e-8)
  expect_equal(pk$estimate - pk$adj.conf.low, pk$adj.conf.high - pk$estimate)
  expect_true(all(pk$adj.conf.high - pk$estimate >= pp$lsd - 1e-8))
  expect_equal(pk[c("conf.low", "conf.high", "lsd")], pp[c("conf.low", "conf.high", "lsd")])
  expect_null(pp$adj.conf.low)
})
