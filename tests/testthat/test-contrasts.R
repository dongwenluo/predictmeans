oats_fac <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}

test_that("specify_contrast() builds the matrix from level names", {
  K <- Kmatrix(lm(yield ~ nitro*Variety, data = oats_fac()), "nitro:Variety")
  cm <- specify_contrast(K, a = `0.6:Victory` - `0:Victory`,
                         b = (`0.2:Marvellous` + `0.4:Marvellous`)/2 - `0:Marvellous`)
  expect_identical(dim(cm), c(2L, 12L))
  expect_identical(rownames(cm), c("a", "b"))
  expect_equal(unname(cm["a", c("0.6:Victory", "0:Victory")]), c(1, -1))
  expect_equal(sum(abs(cm["a", ])), 2)
  expect_equal(unname(cm["b", c("0.2:Marvellous", "0.4:Marvellous", "0:Marvellous")]), c(0.5, 0.5, -1))
  expect_error(specify_contrast(K, avg = (`0:Victory` + `0.2:Victory`)/2), "do not sum to zero")
  expect_equal(sum(specify_contrast(K, avg = (`0:Victory` + `0.2:Victory`)/2, .zero_sum = FALSE)), 1)
  expect_identical(rownames(specify_contrast(K, x = `0.6:Victory` - `0:Victory`)), "x")
  expect_error(specify_contrast(K, y = `0.8:Victory` - `0:Victory`), "Level names are")
  expect_error(specify_contrast(K, `0.6:Victory` - `0:Victory`), "name every contrast")
  pm <- predictmeans(lm(yield ~ nitro*Variety, data = oats_fac()), "nitro:Variety", plot = FALSE, prtnum = FALSE)
  expect_equal(specify_contrast(pm, a = `0.6:Victory` - `0:Victory`), cm["a", , drop = FALSE])
})

test_that("contrastmeans() with specify_contrast() matches emmeans", {
  skip_if_not_installed("emmeans")
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats_fac())
  cm <- specify_contrast(Kmatrix(fm, "nitro:Variety"),
                         N_Victory = `0.6:Victory` - `0:Victory`,
                         N_interaction = (`0.6:Victory` - `0:Victory`) - (`0.6:Marvellous` - `0:Marvellous`))
  res <- contrastmeans(fm, "nitro:Variety", ctrmatrix = cm)[[1]]
  expect_identical(rownames(res), c("N_Victory", "N_interaction"))
  em <- emmeans::emmeans(fm, ~nitro:Variety, lmer.df = "kenward-roger")
  lev <- paste(em@grid$nitro, em@grid$Variety, sep = ":")
  ec <- summary(emmeans::contrast(em, list(N_Victory = cm[1, lev], N_interaction = cm[2, lev])))
  expect_equal(unname(res[, "Estimate"]), round(ec$estimate, 4))
  expect_equal(unname(res[, "df"]), round(ec$df, 4))
  expect_equal(unname(res[, "Pr(>|t|)"]), round(ec$p.value, 4))
})

test_that("contrastmeans() works with a missing cell (B18)", {
  skip_if_not_installed("emmeans")
  Om <- oats_fac()
  Om <- Om[!(Om$nitro == "0.6" & Om$Variety == "Victory"), ]
  fm <- lm(yield ~ nitro*Variety, data = Om)
  K <- suppressWarnings(Kmatrix(fm, "nitro:Variety"))
  cm <- specify_contrast(K, GR_vs_M = `0:Golden Rain` - `0:Marvellous`)
  res <- suppressWarnings(contrastmeans(fm, "nitro:Variety", ctrmatrix = cm))[[1]]
  em <- emmeans::emmeans(fm, ~nitro:Variety)
  lev <- paste(em@grid$nitro, em@grid$Variety, sep = ":")
  ec <- summary(emmeans::contrast(em, list(GR_vs_M = cm[1, lev])))
  expect_equal(unname(res[, "Estimate"]), round(ec$estimate, 4))
  expect_equal(unname(res[, "Std. Error"]), round(ec$SE, 4))
  expect_error(specify_contrast(K, bad = `0.6:Victory` - `0:Victory`), "cannot be estimated")
  bad <- matrix(0, 1, 12, dimnames = list(NULL, rownames(K$K)))
  bad[1, c("0.6:Victory", "0:Victory")] <- c(1, -1)
  expect_error(suppressWarnings(contrastmeans(fm, "nitro:Variety", ctrmatrix = bad)), "cannot be estimated")
})
