oats <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d$x <- as.numeric(as.character(d$nitro))
  d
}

test_that("pm_data(): fitted rows and the response on both scales for every model class", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  d <- oats()
  fits <- list(
    lm = lm(log(yield + 1) ~ nitro * Variety, data = d, subset = Block != "I"),
    aov = aov(yield ~ nitro * Variety + Error(Block / Variety), data = d),
    glm = glm(cbind(round(yield), 200 - round(yield)) ~ nitro, binomial, data = d),
    gls = nlme::gls(sqrt(yield) ~ nitro, data = d),
    lme = nlme::lme(yield ~ nitro, random = ~ 1 | Block, data = d),
    lmer = lme4::lmer(log(yield) ~ nitro + (1 | Block), data = d),
    glmer = lme4::glmer(round(yield) ~ nitro + (1 | Block), poisson, data = d),
    tmb = glmmTMB::glmmTMB(round(yield) ~ nitro + (1 | Block), family = glmmTMB::nbinom2, data = d),
    gam = mgcv::gam(yield ~ nitro + s(x, k = 4), data = d),
    semireg = semireg(yield ~ s(x, k = 4) + Variety + (1 | Block), data = d))
  for (n in names(fits)) {
    pd <- pm_data(fits[[n]])
    r <- attr(pd, "response")
    expect_identical(r$name, "yield", label = n)
    expect_equal(nrow(pd), if (n == "lm") 60L else 72L, label = n)
    expect_true(all(c(".y", ".y_orig") %in% names(pd)), label = n)
    expect_equal(pd$.y_orig, if (n == "glm") round(pd$yield) / 200 else if (n %in% c("glmer", "tmb")) round(pd$yield)
                 else pd$yield, label = n)
  }
  pd <- pm_data(fits$lm)
  expect_equal(pd$.y, log(pd$yield + 1))
  expect_identical(attr(pd, "response")[c("expr", "scale", "back")],
                   list(expr = "log(yield + 1)", scale = "transformed", back = "formula"))
  expect_identical(attr(pm_data(fits$glmer), "response")$scale, "link")
})

test_that("pm_data(): the original response from the inverse when the data hold only the transformed one", {
  d <- oats()
  d$ly <- log(d$yield)
  fm <- lm(ly ~ nitro, data = d)
  pd <- pm_data(fm, trans = exp)
  expect_equal(pd$.y_orig, d$yield)
  expect_identical(attr(pd, "response")$back, "trans")
  expect_identical(attr(pm_data(fm), "response")$back, NA_character_)
})

test_that("axis labels and titles say what the means are, and 'back-transformed' only when they are", {
  d <- oats()
  d$ly <- log(d$yield)
  lab <- function(fit, ...) {
    s <- suppressMessages(suppressWarnings(predictmeans(fit, "nitro", ...)))$settings
    c(pm_scale_label(s, "model"), pm_scale_label(s, "response"), pm_scale_label(s, "response", "title"))
  }
  expect_identical(lab(lm(yield ~ nitro, d)), c("yield", "yield", "Means of yield for 'nitro'"))
  expect_identical(lab(lm(log(yield + 1) ~ nitro, d)),
                   c("log(yield + 1)", "Median of yield",
                     "Medians of yield for 'nitro' (back-transformed from log(yield + 1))"))
  expect_identical(lab(lm(log(yield) ~ nitro, d), bias_adjust = TRUE)[2:3],
                   c("Mean of yield", "Means of yield for 'nitro' (back-transformed from log(yield), bias-adjusted)"))
  expect_identical(lab(lm(log(yield) ~ nitro, d), trans = I), rep(c("log(yield)", "Means of log(yield) for 'nitro'"),
                                                                  c(2, 1)))
  expect_identical(lab(lm(ly ~ nitro, d), trans = exp)[1:2], c("ly", "trans(ly)"))
  expect_identical(lab(glm(cbind(round(yield), 200 - round(yield)) ~ nitro, binomial, data = d)),
                   c("logit(probability of yield)", "Probability of yield",
                     "Probabilities of yield for 'nitro' (inverse logit link)"))
  expect_identical(lab(lme4::glmer(round(yield) ~ nitro + (1 | Block), poisson, data = d))[2],
                   "Mean of yield, typical cluster")
  pm <- suppressMessages(predictmeans(lm(log(yield + 1) ~ nitro, d), "nitro"))
  g <- ggplot2::autoplot(pm, type = "ci")
  expect_identical(g$labels$x, "Median of yield")
  expect_false(grepl("Back", ggplot2::autoplot(predictmeans(lm(yield ~ nitro, d), "nitro"))$labels$title))
})
