# Regression tests for bugs found while splitting predictmeans() into stages.

oats_lmer <- function() {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  suppressMessages(lmerTest::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats))
}

test_that("permlist works with atvar", {
  m <- oats_lmer()
  perm <- suppressWarnings(permmodels(m, nperm = 9, prt = FALSE, ncore = 1, seed = 1))
  r <- predictmeans(m, "nitro:Variety", permlist = perm, atvar = "Variety",
                    plot = FALSE, prtnum = FALSE)
  expect_equal(nrow(r$means), 12)
  expect_true(all(is.na(r$means$df)))
  expect_true(all(r$pairs$p.method == "permutation"))
  # p-values within each variety are adjusted and present
  within <- r$pairs[!is.na(r$pairs$by), ]
  expect_equal(nrow(within), 3 * 6)
  expect_false(anyNA(within$p.adj))
})

test_that("a user-supplied Df is used for every pair of an lmer model", {
  m <- oats_lmer()
  r <- predictmeans(m, "nitro", Df = 20, adj = "BH", plot = FALSE, prtnum = FALSE)
  pv <- r$`Pairwise p-value`
  tv <- pv[upper.tri(pv)]
  expect_equal(t(pv)[upper.tri(pv)], round(p.adjust(2 * pt(-abs(tv), 20), "BH"), 4),
               tolerance = 1e-4)
  expect_true(all(r$means$df == 20))
  r_t <- predictmeans(m, "nitro", Df = 20, adj = "tukey", plot = FALSE, prtnum = FALSE)
  pt_ <- r_t$`Pairwise p-value`
  expect_false(anyNA(t(pt_)[upper.tri(pt_)]))
  expect_equal(t(pt_)[upper.tri(pt_)],
               round(ptukey(sqrt(2) * abs(tv), 4, 20, lower.tail = FALSE), 4), tolerance = 1e-4)
})
