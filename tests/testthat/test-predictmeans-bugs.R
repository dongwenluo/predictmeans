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
