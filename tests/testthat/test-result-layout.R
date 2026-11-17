## contrastmeans(), pm_slopes() and predictmeans() on semireg fits: tables in
## the layout of predictmeans()$pairs (estimate, std.error, df, statistic,
## p.value, p.adj, conf.low, conf.high, adj.conf.*), the 1.x parts unchanged.

test_that("contrastmeans(): the contrasts table, the 1.x table in its place, unadjusted and adjusted p-values", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fm <- lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  cm <- rbind(a = c(-1, 0, 1, rep(0, 9)), b = c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
  ct <- contrastmeans(fm, "nitro:Variety", ctrmatrix = cm, joint = TRUE, adj = "bonferroni")
  expect_s3_class(ct, "contrastmeans")
  expect_named(ct, c("The t tests of the specified contrasts", "K", "Joint F test of the contrasts", "contrasts",
                     "joint", "settings"))
  old <- ct[[1]]
  tab <- ct$contrasts
  expect_named(tab, c("contrast", "estimate", "std.error", "df", "statistic", "p.value", "p.adj", "conf.low",
                      "conf.high", "adj.conf.low", "adj.conf.high", "p.method"))
  expect_equal(tab$contrast, c("a", "b"))
  expect_equal(tab$estimate, unname(old[, "Estimate"]), tolerance = 1e-4)
  expect_equal(tab$p.adj, unname(old[, "Pr(>|t|)"]), tolerance = 1e-3)
  expect_equal(tab$p.value, 2*pt(-abs(tab$statistic), tab$df))
  expect_equal(tab$p.adj, pmin(1, 2*tab$p.value))                     # bonferroni, 2 contrasts
  expect_equal(tab$conf.low, unname(old[, "lower.CL"]), tolerance = 1e-4)
  expect_equal(tab$adj.conf.high, unname(old[, "adj.upper.CL"]), tolerance = 1e-4)
  expect_equal(tab$p.adj < 0.05, tab$adj.conf.high < 0 | tab$adj.conf.low > 0)
  expect_equal(ct$joint$F.ratio, ct[[3]]$F.ratio, tolerance = 1e-4)
  expect_identical(as.data.frame(ct), tab)
  expect_identical(as.data.frame(ct, what = "joint"), ct$joint)
  if (requireNamespace("generics", quietly = TRUE)) expect_identical(generics::tidy(ct), tab)
  expect_equal(ct$settings$df_method, "Kenward-Roger")
  expect_equal(ct$settings$adjust, "bonferroni")
  expect_output(print(ct), "Joint F test: F\\(2, 30.2\\)")
  expect_error(as.data.frame(contrastmeans(fm, "nitro:Variety", ctrmatrix = cm), what = "joint"), "joint = TRUE")
})

test_that("contrastmeans(): odds ratios for a logit link when the weights sum to zero; permutation p-values", {
  data(Drug, package = "predictmeans")
  g <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
  t2 <- contrastmeans(g, "Drug", ctrmatrix = rbind(c(1, -1, 0, 0), c(1, 1, 0, -2)))$contrasts
  expect_equal(t2$odds.ratio, exp(t2$estimate))
  expect_equal(t2$odds.ratio.conf.low, exp(t2$conf.low))
  expect_null(contrastmeans(g, "Drug", ctrmatrix = rbind(c(1, 0, 0, 0)))$contrasts$odds.ratio)   # not a contrast
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fl <- lm(yield ~ Block + nitro + Variety, data = Oats)
  r <- permmodels(fl, nperm = 29, seed = 1, prt = FALSE, ncore = 1)
  cp <- contrastmeans(fl, "nitro", ctrmatrix = rbind(c(-1, 1, 0, 0), c(-1, 0, 0, 1)), permlist = r, adj = "holm")
  expect_equal(unique(cp$contrasts$p.method), "permutation")
  expect_true(all(is.na(cp$contrasts$df)))
  expect_equal(cp$contrasts$p.adj, unname(cp[[1]][, "Permuted Pr(>|t|)"]), tolerance = 1e-3)
  expect_true(all(cp$contrasts$p.adj >= cp$contrasts$p.value))
  expect_equal(cp$settings$df_method, "permutation")
})

test_that("pm_slopes(): estimate column; pairs with limits that agree with the Tukey p-values", {
  Oats <- as.data.frame(nlme::Oats)
  fm <- lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  s <- pm_slopes(fm, "nitro", "Variety", adj = "tukey")
  expect_named(s$slopes, c("Variety", "estimate", "std.error", "df", "statistic", "p.value", "conf.low", "conf.high"))
  p <- s$pairs
  expect_named(p, c("level1", "level2", "estimate", "std.error", "df", "statistic", "p.value", "p.adj", "conf.low",
                    "conf.high", "adj.conf.low", "adj.conf.high"))
  expect_equal(p$conf.high - p$estimate, qt(0.975, p$df)*p$std.error)
  expect_equal(p$adj.conf.high - p$estimate, qtukey(0.95, 3, p$df)/sqrt(2)*p$std.error)
  expect_identical(as.data.frame(s), s$slopes)
  expect_identical(as.data.frame(s, what = "pairs"), p)
  expect_equal(s$settings$covariate, "nitro")
  expect_output(print(s), "Kenward-Roger df")
})

test_that("predictmeans() on a semireg fit: pairs with conf.*, adj.conf.* and ratios, as predictmeans()", {
  set.seed(3)
  d <- expand.grid(x = seq(0, 1, length = 15), A = factor(c("a1", "a2", "a3")), blk = factor(1:6))
  d$y <- rpois(nrow(d), exp(1 + sin(2*d$x) + 0.3*(d$A == "a2") + rnorm(6, 0, 0.2)[d$blk]))
  f <- semireg(y ~ A + s(x, bs = "os", k = 6) + (1 | blk), family = poisson, data = d, prt = FALSE)
  pm <- suppressMessages(predictmeans(f, "A", adj = "tukey"))
  p <- pm$pairs
  expect_true(all(c("conf.low", "conf.high", "lsd", "adj.conf.low", "adj.conf.high", "ratio", "ratio.conf.low",
                    "ratio.conf.high") %in% names(p)))
  expect_equal(p$conf.high - p$estimate, qt(0.975, p$df)*p$std.error)
  expect_equal(p$p.adj < 0.05, p$adj.conf.low > 0 | p$adj.conf.high < 0)
  expect_equal(p$ratio, exp(p$estimate))
})
