test_that("contrast_family() gives emmeans' coefficient matrices", {
  skip_if_not_installed("emmeans")
  e <- getNamespace("emmeans")
  for (lv in list(letters[1:4], paste0("x", 1:6))) {
    for (f in c("pairwise", "revpairwise", "trt.vs.ctrl", "consec", "poly", "eff", "del.eff")) {
      expect_equal(contrast_family(lv, f), t(as.matrix(get(paste0(f, ".emmc"), e)(lv))), info = f)
    }
  }
  expect_equal(rownames(contrast_family(letters[1:3], "trt.vs.ctrl", ref = "b")), c("a - b", "c - b"))
  expect_error(contrast_family(letters[1:3], "trt.vs.ctrl", ref = "z"), "'ref' must be one of")
})

test_that("Dunnett and trends within 'by' groups through contrastmeans() match emmeans", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("mvtnorm")
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  set.seed(7)
  d <- d[-sample(nrow(d), 9), ]
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = d)
  r <- contrastmeans(fm, "nitro", ctrmatrix = contrast_family(Kmatrix(fm, "nitro"), "trt.vs.ctrl"), adj = "mvt")[[1]]
  em <- summary(emmeans::contrast(emmeans::emmeans(fm, ~nitro, lmer.df = "kenward-roger"), "trt.vs.ctrl", adjust = "mvt"))
  expect_equal(unname(r[, "Estimate"]), round(em$estimate, 4))
  expect_equal(unname(r[, "Pr(>|t|)"]), round(em$p.value, 4), tolerance = 1e-3)

  pm <- predictmeans(fm, "nitro:Variety", plot = FALSE)
  cf <- contrast_family(pm, "poly", by = "Variety", max.degree = 2)
  expect_equal(rownames(cf)[1:2], c("Golden Rain: linear", "Golden Rain: quadratic"))
  r2 <- contrastmeans(fm, "nitro:Variety", ctrmatrix = cf)[[1]]
  e2 <- summary(emmeans::contrast(emmeans::emmeans(fm, ~nitro | Variety, lmer.df = "kenward-roger"),
                                  "poly", max.degree = 2))
  expect_equal(unname(r2[, "Estimate"]), round(e2$estimate, 4))
  expect_equal(unname(r2[, "Std. Error"]), round(e2$SE, 4))
  expect_equal(unname(r2[, "df"]), round(e2$df, 4))
})
