## One place for the inference choices (R/inference.R): predictmeans(),
## contrastmeans(), pm_joint() and pm_slopes() give the same df, covariance
## and label for the same contrast.

inf_models <- function() {
  oats <- as.data.frame(nlme::Oats)
  oats$nitro <- factor(oats$nitro)
  oats$Variety <- factor(oats$Variety, ordered = FALSE)
  oats$Block <- factor(oats$Block, ordered = FALSE)
  miss <- droplevels(oats[!(oats$Block == "I" & oats$Variety == "Golden Rain"), ])[-c(5, 17, 30), ]
  dr <- data.frame(Subj = factor(rep(1:12, each = 3)), Group = factor(rep(c("A", "B"), each = 18)),
                   Time = factor(rep(1:3, 12)))
  set.seed(1)
  dr$Y <- rnorm(36, 20 + as.integer(dr$Time) + (dr$Group == "B"), 1) + rnorm(12)[dr$Subj]
  list(
    lmer = list(m = lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = miss), term = "Variety"),
    lme = list(m = nlme::lme(yield ~ Variety * nitro, random = ~ 1 | Block/Variety, data = miss), term = "Variety"),
    gls = list(m = nlme::gls(Y ~ Group * Time, data = dr, correlation = nlme::corSymm(form = ~ as.integer(Time) | Subj),
                             weights = nlme::varIdent(form = ~ 1 | Time)), term = "Time"),
    lm = list(m = lm(yield ~ Variety * nitro + Block, data = miss), term = "Variety")
  )
}

test_that("pairs df agree across predictmeans(), contrastmeans() and pm_joint() for every df method", {
  ms <- inf_models()
  for (nm in names(ms)) {
    m <- ms[[nm]]$m
    term <- ms[[nm]]$term
    dfms <- if (nm == "lm") "auto" else c("auto", "Kenward-Roger", "Satterthwaite", "between-within", "KR-expected",
                                         if (nm == "lmer") "containment")
    for (d in dfms) {
      pm <- suppressWarnings(suppressMessages(predictmeans(m, term, df_method = d, plot = FALSE, prtnum = FALSE)))
      pr <- pm$pairs
      lev <- levels(factor(pm$means[[term]]))
      C <- t(vapply(seq_len(nrow(pr)), function(i) (lev == pr$level1[i]) - (lev == pr$level2[i]), numeric(length(lev))))
      colnames(C) <- lev
      cm <- suppressWarnings(contrastmeans(m, term, ctrmatrix = C, df_method = d))[[1]]
      expect_equal(unname(cm[, "df"]), round(pr$df, 4), tolerance = 1e-4, info = paste(nm, d))
      expect_equal(unname(cm[, "Std. Error"]), round(pr$std.error, 4), tolerance = 1e-4, info = paste(nm, d))
      jt <- suppressWarnings(pm_joint(m, term, ctrmatrix = C[1, , drop = FALSE], df_method = d))
      expect_equal(jt$df2, pr$df[1], tolerance = 1e-6, info = paste(nm, d))
      expect_identical(attr(jt, "df_method"), pm$settings$df_method, info = paste(nm, d))
    }
  }
})

test_that("settings$vcov_adjust says whether the Kenward-Roger covariance was used", {
  ms <- inf_models()
  s <- function(m, ...) suppressWarnings(suppressMessages(
    predictmeans(m, "Variety", plot = FALSE, prtnum = FALSE, ...)))$settings$vcov_adjust
  expect_true(s(ms$lmer$m))
  expect_false(s(ms$lmer$m, df_method = "Satterthwaite"))
  expect_false(s(ms$lmer$m, df_method = "containment"))
  expect_true(s(ms$lmer$m, df_method = "containment", vcov_adjust = TRUE))
  expect_false(s(ms$lmer$m, vcov_adjust = FALSE))
  expect_false(s(ms$lm$m))
  expect_true(s(ms$lme$m))                                   # Kenward-Roger by default for REML lme/gls
  expect_false(s(ms$lme$m, df_method = "containment"))       # nlme's df, as 1.x
  expect_true(suppressMessages(predictmeans(ms$gls$m, "Time", plot = FALSE, prtnum = FALSE))$settings$vcov_adjust)
})

test_that("the legacy p-value matrix uses each pair's df (gls with Kenward-Roger)", {
  m <- inf_models()$gls$m
  pm <- predictmeans(m, "Group:Time", df_method = "Kenward-Roger", pairwise = TRUE, plot = FALSE, prtnum = FALSE)
  P <- pm$`Pairwise p-value`
  expect_gt(length(unique(round(pm$pairs$df, 6))), 1)
  expect_equal(unname(t(P)[upper.tri(P)]), round(pm$pairs$p.value, 4), tolerance = 1e-4)   # p-values below the diagonal
})

test_that("slopes: one df per term from the terms with the covariate, for every one-df method", {
  oats <- as.data.frame(nlme::Oats)
  oats$Variety <- factor(oats$Variety, ordered = FALSE)
  oats$Block <- factor(oats$Block, ordered = FALSE)
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = oats)
  cdf <- pm_containment_df(m)
  s <- pm_slopes(m, "nitro", "Variety", df_method = "containment")
  expect_equal(unique(s$slopes$df), min(cdf[c("nitro", "Variety:nitro")]))
  b <- pm_slopes(m, "nitro", "Variety", df_method = "between-within")
  expect_equal(unique(b$slopes$df), min(pm_betwithin_df(m)[c("nitro", "Variety:nitro")]))
})

test_that("covariatemeans() follows the same defaults: Kenward-Roger df per point, adjusted SEs", {
  oats <- as.data.frame(nlme::Oats)
  oats$Variety <- factor(oats$Variety, ordered = FALSE)
  oats$Block <- factor(oats$Block, ordered = FALSE)
  oats <- oats[-c(3, 20, 41), ]
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = oats)
  a <- covariatemeans(m, "Variety", covariate = "nitro", covariateV = c(0, 0.3, 0.6))$pltdf
  u <- covariatemeans(m, "Variety", covariate = "nitro", covariateV = c(0, 0.3, 0.6), vcov_adjust = FALSE)$pltdf
  inf <- pm_inference(m)
  K <- Kmatrix(m, "Variety", "nitro", covariateV = c(0, 0.3, 0.6))$K
  expect_equal(a$ses, unname(sqrt(diag(K %*% inf$mp$vcov %*% t(K)))), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(a$ses, u$ses)))
  crit <- (a$UL - a$Mean) / a$ses
  expect_equal(crit, stats::qt(0.975, inf$df(K, "Variety")), tolerance = 1e-6)
})
