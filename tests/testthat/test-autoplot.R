oats_fac <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}

test_that("autoplot() draws the means from the means table", {
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE, prtnum = FALSE)
  p <- ggplot2::autoplot(pm)
  expect_s3_class(p, "ggplot")
  pts <- ggplot2::layer_data(p, 2)
  expect_equal(sort(pts$y), sort(pm$means$estimate))
  bars <- ggplot2::layer_data(p, 3)
  expect_equal(sort(bars$ymax), sort(pm$means$conf.high))
  pl <- ggplot2::autoplot(pm, interval = "lsd")
  bl <- ggplot2::layer_data(pl, 3)
  within <- pm$pairs[!is.na(pm$pairs$by), ]
  expect_equal(unique(round(bl$ymax - bl$ymin, 8)), round(mean(within$lsd), 8))
  expect_length(ggplot2::autoplot(pm, interval = "none")$layers, 2)
})

test_that("autoplot(type = 'pvalue') shows one tile per comparison", {
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro", adj = "BH", plot = FALSE, prtnum = FALSE)
  tiles <- ggplot2::layer_data(ggplot2::autoplot(pm, type = "pvalue"), 1)
  expect_equal(nrow(tiles), nrow(pm$pairs))
  pa <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "BH", plot = FALSE, prtnum = FALSE)
  pp <- ggplot2::autoplot(pa, type = "pvalue")
  expect_equal(nrow(ggplot2::layer_data(pp, 1)), 3 * 6)
  expect_s3_class(pp$facet, "FacetWrap")
})

test_that("autoplot(type = 'response') uses the back-transformed means", {
  data("Drug", package = "predictmeans", envir = environment())
  fm <- glm(cbind(R, N - R) ~ Drug + log10Dose, data = Drug, family = binomial)
  pm <- predictmeans(fm, "Drug", plot = FALSE, prtnum = FALSE)
  pts <- ggplot2::layer_data(ggplot2::autoplot(pm, type = "response"), 2)
  expect_equal(sort(pts$y), sort(pm$means$response))
  fl <- lm(yield ~ nitro, data = oats_fac())
  expect_error(ggplot2::autoplot(predictmeans(fl, "nitro", plot = FALSE, prtnum = FALSE), type = "response"),
               "not back-transformed")
})

test_that("summary() describes the result", {
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats_fac())
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE, prtnum = FALSE)
  out <- capture.output(print(summary(pm)))
  expect_true(any(grepl("Predicted means for 'nitro:Variety'.*Kenward-Roger", out)))
  expect_true(any(grepl("18 pairwise comparisons within each level of Variety", out)))
  expect_true(any(grepl("^SED: min", out)))
  expect_true(any(grepl("^LSD \\(5%\\): min", out)))
  expect_true(any(grepl("adjusted by 'tukey'", out)))
  expect_true(any(grepl("Means sharing a letter", out)))
  perm <- suppressWarnings(permmodels(fm, nperm = 9, prt = FALSE, ncore = 1, seed = 1))
  pp <- predictmeans(fm, "nitro", permlist = perm, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_true(any(grepl("permutation \\(9 permutations\\)", capture.output(print(summary(pp))))))
})
