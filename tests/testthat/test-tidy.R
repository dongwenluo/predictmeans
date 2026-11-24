# The stable result tables of predictmeans(), checked against emmeans.

oats_fac <- function() {
  d <- as.data.frame(nlme::Oats)
  # non-numeric labels, so emmeans doesn't prefix them in contrast names
  d$nitro <- factor(d$nitro, labels = c("n0", "n2", "n4", "n6"))
  d
}
key_em <- function(g, vars) do.call(paste, c(lapply(g[vars], as.character), sep = ":"))

test_that("means match emmeans: estimates, SEs, KR df and confidence limits", {
  skip_if_not_installed("emmeans")
  set.seed(7)
  Ou <- oats_fac()
  Ou <- Ou[-sample(nrow(Ou), 9), ]
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou)
  pm <- predictmeans(fm, "nitro:Variety", vcov_adjust = TRUE, df_method = "KR-expected", plot = FALSE)
  em <- summary(emmeans::emmeans(fm, ~nitro:Variety, lmer.df = "kenward-roger"))
  m <- tidy.predictmeans(pm)
  idx <- match(key_em(m, c("nitro", "Variety")), key_em(em, c("nitro", "Variety")))
  expect_equal(m$estimate, em$emmean[idx], tolerance = 1e-8)
  expect_equal(m$std.error, em$SE[idx], tolerance = 1e-8)
  expect_equal(m$df, em$df[idx], tolerance = 1e-6)
  expect_equal(m$conf.low, em$lower.CL[idx], tolerance = 1e-6)
  expect_equal(m$conf.high, em$upper.CL[idx], tolerance = 1e-6)
  expect_true(all(m$estimable))
  expect_equal(sum(m$n), nrow(Ou))
  # same numbers as the legacy table, before rounding
  mt <- pm$mean_table
  expect_equal(round(m$estimate, 4), mt$Mean[match(key_em(m, c("nitro", "Variety")), key_em(mt, c("nitro", "Variety")))])
})

test_that("pairs match emmeans: raw, Tukey and BH p-values", {
  skip_if_not_installed("emmeans")
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  em_means <- emmeans::emmeans(fm, ~nitro)
  for (a in c("none", "tukey", "BH")) {
    pm <- predictmeans(fm, "nitro", adj = a, plot = FALSE)
    pr <- tidy.predictmeans(pm, "pairs")
    em <- summary(pairs(em_means, adjust = if (a == "none") "none" else a))
    idx <- match(paste(pr$level1, "-", pr$level2), as.character(em$contrast))
    expect_false(anyNA(idx))
    expect_equal(pr$estimate, em$estimate[idx], tolerance = 1e-10)
    expect_equal(pr$std.error, em$SE[idx], tolerance = 1e-10)
    expect_equal(pr$p.adj, em$p.value[idx], tolerance = 1e-6, label = a)
  }
})

test_that("pairs within atvar levels match emmeans 'by' comparisons", {
  skip_if_not_installed("emmeans")
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  for (a in c("tukey", "holm")) {
    pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = a, plot = FALSE)
    pr <- tidy.predictmeans(pm, "pairs")
    pr <- pr[!is.na(pr$by), ]
    expect_equal(nrow(pr), 3 * 6)
    em <- summary(pairs(emmeans::emmeans(fm, ~nitro | Variety), adjust = a))
    lv1 <- sub(":.*", "", pr$level1)
    lv2 <- sub(":.*", "", pr$level2)
    idx <- match(paste(pr$by, lv1, "-", lv2), paste(em$Variety, em$contrast))
    expect_equal(pr$p.adj, em$p.value[idx], tolerance = 1e-6, label = a)
  }
})

test_that("back-transformed means are on the response scale", {
  skip_if_not_installed("emmeans")
  data("Drug", package = "predictmeans", envir = environment())
  fm <- glm(cbind(R, N - R) ~ Drug + log10Dose, data = Drug, family = binomial)
  m <- tidy.predictmeans(predictmeans(fm, "Drug", plot = FALSE))
  em <- summary(emmeans::emmeans(fm, ~Drug, type = "response"))
  expect_equal(m$response, em$prob[match(m$Drug, em$Drug)], tolerance = 1e-8)
  expect_true(all(m$response.conf.low < m$response & m$response < m$response.conf.high))
})

test_that("permutation p-values in pairs equal the legacy p-value matrix", {
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats_fac())
  perm <- suppressWarnings(permmodels(fm, method = "raw", nperm = 19, ncore = 1, seed = 2))
  pm <- predictmeans(fm, "nitro", permlist = perm, adj = "BH", plot = FALSE, ncore = 1)
  pr <- tidy.predictmeans(pm, "pairs")
  expect_true(all(pr$p.method == "permutation"))
  legacy <- pm$`Pairwise p-value`
  expect_equal(round(pr$p.adj, 4), unname(legacy[cbind(pr$level2, pr$level1)]))
  expect_identical(pm$settings$nperm, 19L)
})

test_that("printing is unchanged and the stable slots have fixed names", {
  fm <- lm(yield ~ nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "BH", plot = FALSE)
  expect_s3_class(pm, c("predictmeans", "pdmlist"))
  out <- capture.output(print(pm))
  expect_false(any(grepl("^\\$means|^\\$pairs|^\\$settings", out)))
  expect_identical(names(tidy.predictmeans(pm))[1:2], c("nitro", "Variety"))
  expect_true(all(c("estimate", "std.error", "df", "conf.low", "conf.high", "group", "n", "estimable") %in%
                    names(as.data.frame(pm))))
  expect_true(all(c("level1", "level2", "by", "estimate", "std.error", "df", "statistic", "p.value",
                    "p.adj", "lsd", "p.method") %in% names(as.data.frame(pm, what = "pairs"))))
  skip_if_not_installed("generics")
  expect_identical(generics::tidy(pm), pm$means)
})

test_that("coef() and vcov() of the means match emmeans and give the pairs' SEs", {
  skip_if_not_installed("emmeans")
  set.seed(7)
  Ou <- oats_fac()
  Ou <- Ou[-sample(nrow(Ou), 9), ]
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Ou)
  pm <- predictmeans(fm, "nitro:Variety", pairwise = TRUE, df_method = "KR-expected", plot = FALSE)
  em <- emmeans::emmeans(fm, ~nitro:Variety, lmer.df = "kenward-roger")
  b <- coef(pm)
  V <- vcov(pm)
  expect_equal(names(b), rownames(V))
  idx <- match(names(b), key_em(summary(em), c("nitro", "Variety")))
  expect_equal(unname(b), summary(em)$emmean[idx], tolerance = 1e-8)
  expect_equal(unname(V), unname(vcov(em)[idx, idx]), tolerance = 1e-8)
  # each pair from coef()/vcov()
  L <- outer(pm$pairs$level1, names(b), "==") - outer(pm$pairs$level2, names(b), "==")
  expect_equal(drop(L %*% b), pm$pairs$estimate, tolerance = 1e-10)
  expect_equal(sqrt(diag(L %*% V %*% t(L))), pm$pairs$std.error, tolerance = 1e-10)
})

test_that("marginaleffects::hypotheses() works on predictmeans results", {
  skip_if_not_installed("marginaleffects")
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro", pairwise = TRUE, plot = FALSE)
  # newer marginaleffects only accepts classes listed in this option
  op <- options(marginaleffects_model_classes = c(getOption("marginaleffects_model_classes"), "predictmeans"))
  on.exit(options(op), add = TRUE)
  h <- marginaleffects::hypotheses(pm, hypothesis = "n2 - n0 = 0")
  expect_equal(as.numeric(h$estimate), -pm$pairs$estimate[1], tolerance = 1e-8)
  expect_equal(as.numeric(h$std.error), pm$pairs$std.error[1], tolerance = 1e-6)
})

test_that("with a permlist the limits of the means are the model's t limits at conf_level (not +- 2 SE)", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro, Oats)
  pl <- permmodels(fm, nperm = 19, seed = 1)
  p0 <- predictmeans(fm, "nitro")
  for (cl in c(0.95, 0.99)) {
    pp <- predictmeans(fm, "nitro", permlist = pl, conf_level = cl)
    crit <- qt(1 - (1 - cl)/2, df.residual(fm))
    expect_equal(pp$means$conf.low, pp$means$estimate - crit * pp$means$std.error, label = cl)
    expect_true(all(is.na(pp$means$df)))
    expect_equal(pp$mean_table$Table[[paste0("LL(", cl * 100, "%)")]], round(pp$means$conf.low, 4))
  }
  expect_equal(predictmeans(fm, "nitro", permlist = pl)$means$conf.low, p0$means$conf.low)
  # back-transformed limits too
  fg <- glm(round(yield/10) ~ Block + nitro, poisson, Oats)
  pg <- predictmeans(fg, "nitro", permlist = permmodels(fg, nperm = 19, seed = 1), conf_level = 0.9)
  expect_equal(pg$means$response.conf.low, exp(pg$means$estimate - qnorm(0.95) * pg$means$std.error))
  expect_snapshot(print(summary(predictmeans(fm, "nitro", permlist = pl))))
})
