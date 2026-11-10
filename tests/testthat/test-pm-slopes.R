oats_slopes <- function() {
  d <- as.data.frame(nlme::Oats)
  set.seed(7)
  d[-sample(nrow(d), 9), ]
}

test_that("pm_slopes() equals emmeans::emtrends() for lmer (Kenward-Roger), with Tukey pairs", {
  skip_if_not_installed("emmeans")
  d <- oats_slopes()
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = d)
  s <- pm_slopes(fm, "nitro", "Variety", adj = "tukey")
  et <- emmeans::emtrends(fm, ~Variety, var = "nitro", lmer.df = "kenward-roger")
  e <- summary(et)
  i <- match(as.character(s$slopes$Variety), as.character(e$Variety))
  expect_equal(s$slopes$slope, e$nitro.trend[i], tolerance = 1e-6)
  expect_equal(s$slopes$std.error, e$SE[i], tolerance = 1e-6)
  expect_equal(s$slopes$df, e$df[i], tolerance = 1e-6)
  ep <- summary(pairs(et))
  j <- match(paste(s$pairs$level1, "-", s$pairs$level2), as.character(ep$contrast))
  expect_equal(s$pairs$estimate, ep$estimate[j], tolerance = 1e-6)
  expect_equal(s$pairs$p.adj, ep$p.value[j], tolerance = 1e-6)
})

test_that("pm_slopes() with poly(): the slope at 'at' (central difference, exact for a quadratic)", {
  skip_if_not_installed("emmeans")
  d <- oats_slopes()
  fp <- lm(yield ~ poly(nitro, 2)*Variety + Block, data = d)
  s <- pm_slopes(fp, "nitro", "Variety", at = 0.2)
  e <- summary(emmeans::emtrends(fp, ~Variety, var = "nitro", at = list(nitro = 0.2), delta.var = 1e-7))
  i <- match(as.character(s$slopes$Variety), as.character(e$Variety))
  expect_equal(s$slopes$slope, e$nitro.trend[i], tolerance = 1e-5)
  expect_equal(s$slopes$std.error, e$SE[i], tolerance = 1e-5)
  expect_equal(unique(s$slopes$df), df.residual(fp))
})

test_that("pm_slopes() for lme uses nlme's df of the covariate terms", {
  d <- oats_slopes()
  fl <- nlme::lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = d)
  s <- pm_slopes(fl, "nitro", "Variety")
  tdf <- terms(fl$fixDF)
  expect_equal(unique(s$slopes$df), min(tdf[c("nitro", "nitro:Variety")]))
  expect_equal(nrow(s$pairs), 3)
  expect_error(pm_slopes(fl, "Variety", "Variety"), "numeric variable")
})
