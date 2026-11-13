## residplot(): four ggplot2 plots for every model class predictmeans handles.

rp <- function(m, ...) residplot(m, plot = FALSE, ...)
oats_f <- function() {
  O <- as.data.frame(nlme::Oats)
  O$nitro <- factor(O$nitro)
  O$Variety <- factor(O$Variety, ordered = FALSE)
  O$Block <- factor(O$Block, ordered = FALSE)
  O
}

test_that("every model class gives four ggplots with the right panels", {
  O <- oats_f()
  ms <- list(lm = lm(yield ~ Block + nitro * Variety, O),
             aovlist = aov(yield ~ nitro * Variety + Error(Block/Variety), O),
             glm = glm(breaks ~ wool * tension, family = poisson, data = warpbreaks),
             gls = nlme::gls(yield ~ nitro * Variety, correlation = nlme::corCompSymm(form = ~ 1 | Block), data = O),
             lme = nlme::lme(yield ~ nitro * Variety, random = ~ 1 | Block/Variety, data = O),
             lmer = lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = O),
             glmer = lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial,
                                 data = lme4::cbpp),
             tmb = glmmTMB::glmmTMB(breaks ~ wool * tension + (1 | tension:wool), family = poisson, data = warpbreaks))
  third <- c(lm = "Scale-location", aovlist = "Scale-location", glm = "Scale-location",
             gls = "Autocorrelation of normalized residuals", lme = "Normal Q-Q of random effects: (Intercept) | Block",
             lmer = "Normal Q-Q of random effects: (Intercept) | Variety:Block",
             glmer = "Normal Q-Q of random effects: (Intercept) | herd",
             tmb = "Normal Q-Q of random effects: (Intercept) | tension:wool")
  for (nm in names(ms)) {
    p <- suppressWarnings(rp(ms[[nm]]))
    expect_s3_class(p, "pm_residplot")
    expect_length(p, 4)
    expect_true(all(vapply(p, ggplot2::is.ggplot, logical(1))), label = nm)
    expect_identical(p[[3]]$labels$title, unname(third[nm]), label = nm)
    expect_identical(p[[4]]$labels$title, if (nm %in% c("lm", "glm")) "Cook's distance" else "Fitted vs observed",
                     label = nm)
    expect_no_error(ggplot2::ggplot_build(p[[1]]))
  }
  # axis labels: residuals on y, standard normal quantiles on x
  p <- rp(ms$lmer)
  expect_identical(p[[2]]$labels$x, "Theoretical quantiles (standard normal)")
  expect_identical(p[[2]]$labels$y, "Standardized residuals")
  expect_identical(p[[1]]$labels$x, "Fitted values")
  # lmer residuals are scaled by sigma
  expect_equal(attr(p, "data")$resid, unname(residuals(ms$lmer)) / sigma(ms$lmer), tolerance = 1e-10)
  expect_identical(rp(ms$glm)[[1]]$labels$y, "Randomized quantile residuals")
})

test_that("randomized quantile residuals: standard normal for a right Poisson model, reproducible", {
  set.seed(5)
  d <- data.frame(x = runif(2000))
  d$y <- rpois(2000, exp(1 + d$x))
  g <- glm(y ~ x, family = poisson, data = d)
  r <- attr(rp(g), "data")$resid
  expect_equal(c(mean(r), sd(r)), c(0, 1), tolerance = 0.06)
  expect_gt(stats::shapiro.test(r[1:500])$p.value, 0.01)
  expect_identical(attr(rp(g), "data")$resid, r)          # fixed seed
  s <- .Random.seed; rp(g); expect_identical(.Random.seed, s)   # the user's random stream is untouched
})

test_that("group colours, labels of extreme points, binary calibration, semireg, plotly", {
  O <- oats_f()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = O)
  p <- rp(fm, group = "Variety", id = 2)
  expect_identical(levels(attr(p, "data")$group), levels(O$Variety))
  lab_layer <- Filter(function(l) inherits(l$geom, "GeomText"), p[[1]]$layers)
  expect_equal(nrow(lab_layer[[1]]$data), 2)
  expect_length(Filter(function(l) inherits(l$geom, "GeomText"), rp(fm, id = FALSE)[[1]]$layers), 0)
  expect_error(rp(fm, group = "nothing"), "'group'")
  set.seed(1)
  b <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
  b$y <- rbinom(200, 1, plogis(0.5 * b$x))
  gb <- lme4::glmer(y ~ x + (1 | g), family = binomial, data = b)
  expect_identical(suppressMessages(rp(gb))[[4]]$labels$title, "Calibration (tenths of fitted values)")
  skip_if_not_installed("plotly")
  expect_s3_class(residplot(fm, interactive = TRUE, plot = FALSE), "plotly")
})
