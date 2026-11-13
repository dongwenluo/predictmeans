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
  perm <- suppressWarnings(permmodels(fm, method = "raw", nperm = 49, prt = FALSE, ncore = 1, seed = 11))
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
  invisible(permmodels(fm, method = "raw", nperm = 5, prt = FALSE, ncore = 1, seed = 1))
  expect_equal(getOption("scipen"), 0)
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
  local_fit <- function(dat) glmmTMB::glmmTMB(cbind(Favorable, Total - Favorable) ~ Treatment + (1|Clinic),
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
  p <- permmodels(fm, method = "raw", nperm = 9, prt = FALSE, ncore = 1, seed = 3)
  expect_equal(p$ANOVA$numDF, anova(fm)$numDF)
  expect_identical(rownames(p$COEFFICENTS), names(nlme::fixef(fm)))
  # A local data set with a name found nowhere else
  local_dat <- Oats
  fm2 <- lme(yield ~ nitro + Variety, random = ~1|Block, data = local_dat)
  expect_length(CookD(fm2, group = "Block", plot = FALSE), 6)
  fg <- gls(yield ~ nitro + Variety, data = local_dat)
  expect_length(permmodels(fg, method = "raw", nperm = 5, prt = FALSE, ncore = 1, seed = 1)$permlist, 5)
})

test_that("predictmeans() works for lmer models with a missing cell (B15)", {
  Om <- oats()
  Om <- Om[!(Om$nitro == "0.6" & Om$Variety == "Victory"), ]
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Om)
  pm <- suppressWarnings(predictmeans(fm, "nitro:Variety", plot = FALSE, prtnum = FALSE))
  expect_equal(nrow(pm$mean_table), 11)
  K <- suppressWarnings(Kmatrix(fm, "nitro:Variety")$K)
  K <- K[rownames(K) != "0.6:Victory", ]
  expect_equal(pm$mean_table$Df, round(unname(df_term(fm, ctrmatrix = K)), 2))
})

test_that("predictmeans() works for lmer models fitted by ML, with Satterthwaite df (B16)", {
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats(), REML = FALSE)
  expect_message(pm <- predictmeans(fm, "nitro", plot = FALSE, prtnum = FALSE), "Satterthwaite")
  K <- Kmatrix(fm, "nitro")$K
  expect_equal(pm$mean_table$Df, round(unname(df_term(fm, ctrmatrix = K, type = "Satterthwaite")), 2))
})

test_that("vcov_adjust = TRUE gives Kenward-Roger adjusted SEs, as emmeans", {
  skip_if_not_installed("emmeans")
  set.seed(7)
  Ou <- oats()
  Ou <- Ou[-sample(nrow(Ou), 9), ]
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou)
  em <- summary(emmeans::emmeans(fm, ~Variety, lmer.df = "kenward-roger"))
  adj <- predictmeans(fm, "Variety", vcov_adjust = TRUE, df_method = "KR-expected", plot = FALSE, prtnum = FALSE)$mean_table
  unadj <- predictmeans(fm, "Variety", vcov_adjust = FALSE, df_method = "KR-expected", plot = FALSE, prtnum = FALSE)$mean_table
  expect_equal(adj$SE, round(em$SE, 5))
  expect_equal(adj$Df, round(em$df, 2))
  expect_equal(adj$Mean, unadj$Mean)
  expect_true(all(adj$SE > unadj$SE))
  cm <- contrastmeans(fm, "Variety", ctrmatrix = rbind(c(1, -1, 0)), vcov_adjust = TRUE, df_method = "KR-expected")
  ec <- summary(emmeans::contrast(emmeans::emmeans(fm, ~Variety, lmer.df = "kenward-roger"),
                                  list(c(1, -1, 0))))
  expect_equal(unname(cm[[1]][1, "Std. Error"]), round(ec$SE, 4))
})

test_that("predictmeans() passes cov_reduce and exclude to Kmatrix()", {
  skip_if_not_installed("mgcv")
  set.seed(9)
  d <- expand.grid(A = factor(c("a1", "a2")), blk = factor(1:5), rep = 1:4)
  d$x <- runif(nrow(d))
  d$y <- as.numeric(d$A) + sin(4*d$x) + rnorm(5)[d$blk] + rnorm(nrow(d), sd = 0.3)
  g <- mgcv::gam(y ~ A + s(x) + s(blk, bs = "re"), data = d)
  K <- Kmatrix(g, "A", cov_reduce = "average")$K
  pm <- predictmeans(g, "A", cov_reduce = "average", plot = FALSE, prtnum = FALSE)
  expect_equal(pm$mean_table$Mean, round(as.vector(K %*% coef(g)), 4))
})

test_that("Kenward-Roger adjusted SEs are the default for REML lmer fits", {
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats())
  a <- predictmeans(fm, "Variety", plot = FALSE, prtnum = FALSE)
  b <- predictmeans(fm, "Variety", vcov_adjust = TRUE, plot = FALSE, prtnum = FALSE)
  expect_equal(a$means, b$means)
  expect_true(a$settings$vcov_adjust)
  ml <- lmer(yield ~ nitro + (1|Block), data = oats(), REML = FALSE)
  expect_no_warning(suppressMessages(predictmeans(ml, "nitro", plot = FALSE, prtnum = FALSE)))
})

test_that("permmodels() gives the same tables for lme4 and lmerTest fits", {
  d <- oats()
  f1 <- lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = d)
  f2 <- lmerTest::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = d)
  p1 <- suppressWarnings(permmodels(f1, method = "raw", nperm = 9, prt = FALSE, ncore = 1, seed = 1))
  p2 <- suppressWarnings(permmodels(f2, method = "raw", nperm = 9, prt = FALSE, ncore = 1, seed = 1))
  expect_equal(p1$ANOVA, p2$ANOVA)
  expect_equal(p1$COEFFICENTS, p2$COEFFICENTS)
  expect_true("Pr(>F)" %in% names(p1$ANOVA))
})

test_that("residplot(): glmmTMB with aliased columns; random part in the column order of Z", {
  skip_if_not_installed("glmmTMB")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  set.seed(1)
  da <- data.frame(g = factor(rep(1:6, each = 8)), A = factor(rep(1:4, 12)))
  da$y <- rpois(48, exp(1 + rnorm(6, sd = 0.3)[da$g]))
  da$B <- da$A      # aliased with A: glmmTMB drops its columns
  ma <- suppressMessages(glmmTMB::glmmTMB(y ~ A + B + (1 | g), family = poisson, data = da))
  expect_no_error(residplot(ma))
  # correlated term: Z b equals the random part of the linear predictor
  d <- lme4::sleepstudy
  d$y <- rpois(nrow(d), exp(d$Reaction / 100))
  gm <- suppressWarnings(lme4::glmer(y ~ Days + (Days | Subject), family = poisson, data = d))
  Z <- lme4::getME(gm, "Z")
  zb <- as.vector(Z %*% lme4::getME(gm, "b"))
  expect_equal(zb, unname(predict(gm) - predict(gm, re.form = NA)), tolerance = 1e-8)
  expect_no_error(residplot(gm))
  tm <- glmmTMB::glmmTMB(y ~ Days + (Days | Subject), family = poisson, data = d)
  expect_no_error(residplot(tm))
})

test_that("ci_mcp() accepts LL = UL; model.matrix() of an lme model takes 'data'", {
  # [3, 3] lies inside [2, 4]; [10, 12] overlaps neither
  expect_identical(unname(trimws(ci_mcp(c(3, 2, 10), c(3, 4, 12)))), c("A", "A", "B"))
  # the lme methods are internal to predictmeans
  fm <- nlme::lme(distance ~ age + Sex, random = ~ 1 | Subject, data = nlme::Orthodont)
  X <- predictmeans:::model.matrix.lme(fm)
  expect_identical(dim(X), c(nrow(nlme::Orthodont), 3L))
  rows <- c(1:4, 65:68)   # both sexes
  expect_equal(predictmeans:::model.matrix.lme(fm, data = nlme::Orthodont[rows, ]), X[rows, ], ignore_attr = TRUE)
})
