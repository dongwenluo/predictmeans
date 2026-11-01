oats <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}

test_that("legacy plot arguments give a deprecation warning, not an error (B1)", {
  fm <- lm(yield ~ nitro*Variety, data = oats())
  expect_warning(predictmeans(fm, "nitro", mplot = TRUE, barplot = FALSE, plot = FALSE, prtnum = FALSE),
                 "mplot and barplot are deprecated")
  expect_warning(predictmeans(fm, "nitro", mplot = TRUE, barplot = FALSE, pplot = FALSE, bkplot = FALSE,
                              plot = FALSE, prtnum = FALSE),
                 "mplot, barplot, pplot and bkplot are deprecated")
})

test_that("trans = NULL means no back-transformation (B2)", {
  fm <- lm(yield ~ nitro*Variety, data = oats())
  a <- predictmeans(fm, "nitro", trans = NULL, pairwise = TRUE, plot = FALSE, prtnum = FALSE)
  b <- predictmeans(fm, "nitro", pairwise = TRUE, plot = FALSE, prtnum = FALSE)
  expect_equal(a$mean_table, b$mean_table)
})

test_that("contrastmeans() and predictmeans() give the same permutation p-value (B3)", {
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats())
  perm <- suppressWarnings(permmodels(fm, nperm = 49, prt = FALSE, ncore = 1, seed = 11))
  pm <- predictmeans(fm, "nitro", permlist = perm, adj = "none", plot = FALSE, prtnum = FALSE, ncore = 1)
  cm <- contrastmeans(fm, "nitro", ctrmatrix = rbind(c(1, -1, 0, 0)), permlist = perm)
  # lower triangle of 'Pairwise p-value' holds the p-values; row 0.2, column 0
  expect_equal(unname(cm[[1]][1, "Permuted Pr(>|t|)"]), unname(round(pm$`Pairwise p-value`["0.2", "0"], 4)))
})

test_that("functions don't change options() for the rest of the session", {
  fm <- lm(yield ~ nitro*Variety, data = oats())
  op <- options(scipen = 0)
  on.exit(options(op))
  invisible(predictmeans(fm, "nitro", plot = FALSE, prtnum = FALSE))
  invisible(contrastmeans(fm, "nitro", ctrmatrix = rbind(c(1, -1, 0, 0))))
  invisible(permmodels(fm, nperm = 5, prt = FALSE, ncore = 1, seed = 1))
  expect_identical(getOption("scipen"), 0)
})

test_that("contrastmeans() without ctrmatrix needs an interactive session", {
  skip_if(interactive())
  fm <- lm(yield ~ nitro*Variety, data = oats())
  expect_error(contrastmeans(fm, "nitro"), "interactive session")
})

test_that("CookD() works for lmer models and agrees with lme", {
  data("ATP", package = "predictmeans", envir = environment())
  fa <- suppressMessages(lmer(ATP ~ A*B + (1|heart), data = ATP))
  fl <- lme(ATP ~ A*B, random = ~1|heart, data = ATP)
  c1 <- CookD(fa, group = "heart", plot = FALSE)
  c2 <- CookD(fl, group = "heart", plot = FALSE)
  expect_equal(c1, c2[names(c1)], tolerance = 1e-5)
})

test_that("R2_glmm(): glmmTMB data found outside the global environment; clear over_disp error", {
  data("Clinical", package = "predictmeans", envir = environment())
  local_fit <- function(dat) glmmTMB(cbind(Favorable, Total - Favorable) ~ Treatment + (1|Clinic),
                                     family = binomial, data = dat)
  r <- suppressMessages(capture.output(val <- R2_glmm(local_fit(Clinical))))
  expect_identical(unname(val), c("75.71%", "4.77%", "70.94%", "70.94%"))
  fm <- lmer(yield ~ nitro + (1|Block), data = oats())
  expect_error(R2_glmm(fm, over_disp = TRUE), "observation-level random effect")
})

test_that("R2_glmm() over_disp leaves out an observation-level random effect", {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
  d$y <- rpois(120, exp(0.5 + 0.3*d$x + rnorm(20, sd = 0.4)[d$g]))
  d$obs <- factor(seq_len(120))
  fm <- suppressMessages(glmer(y ~ x + (1|obs) + (1|g), family = poisson, data = d))
  r <- capture.output(od <- R2_glmm(fm, over_disp = TRUE))
  r <- capture.output(full <- R2_glmm(fm))
  # 1.1.1 values for this model
  expect_identical(unname(od), c("11.3%", "9.64%", "1.66%"))
  expect_identical(unname(full), c("26.19%", "9.64%", "16.55%", "14.89%", "1.66%"))
})

test_that("permmodels() and CookD() refit the model on its own data (B12)", {
  # A local data set whose name also exists on the search path (nlme::Oats)
  Oats <- oats()
  fm <- lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = Oats)
  p <- permmodels(fm, nperm = 9, prt = FALSE, ncore = 1, seed = 3)
  expect_equal(p$ANOVA$numDF, anova(fm)$numDF)
  expect_identical(rownames(p$COEFFICENTS), names(nlme::fixef(fm)))
  # A local data set with a name found nowhere else
  local_dat <- Oats
  fm2 <- lme(yield ~ nitro + Variety, random = ~1|Block, data = local_dat)
  expect_length(CookD(fm2, group = "Block", plot = FALSE), 6)
  fg <- gls(yield ~ nitro + Variety, data = local_dat)
  expect_length(permmodels(fg, nperm = 5, prt = FALSE, ncore = 1, seed = 1)$permlist, 5)
})
