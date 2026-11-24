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
    # (long titles are broken into lines: compare the words)
    expect_identical(gsub("\n", " ", p[[3]]$labels$title), unname(third[nm]), label = nm)
    expect_identical(p[[4]]$labels$title, if (nm %in% c("lm", "glm")) "Cook's distance" else "Fitted vs observed",
                     label = nm)
    expect_no_error(ggplot2::ggplot_build(p[[1]]))
  }
  # axis labels: residuals on y, standard normal quantiles on x
  p <- rp(ms$lmer)
  expect_identical(p[[2]]$labels$x, "Normal quantiles")
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

test_that("fitted vs observed: response from pm_data(), axis labelled with the response", {
  O <- oats_f()
  m1 <- lme4::lmer(log(yield + 1) ~ nitro + (1 | Block), data = O)
  p <- rp(m1)
  expect_identical(p[[4]]$labels$x, "log(yield + 1)")
  expect_equal(attr(p, "data")$observed, log(O$yield + 1))
  m2 <- nlme::lme(sqrt(yield) ~ nitro, random = ~ 1 | Block, data = O)
  expect_identical(rp(m2)[[4]]$labels$x, "sqrt(yield)")
  expect_identical(rp(aov(yield ~ nitro + Error(Block), O))[[4]]$labels$x, "yield")
  cb <- lme4::cbpp
  g <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cb)
  expect_identical(rp(g)[[4]]$labels$x, "Observed proportion of incidence")
  # a two-level factor response: 0/1, so the calibration plot (was "response not available")
  cb$y01 <- factor(cb$incidence > 0)
  t1 <- glmmTMB::glmmTMB(y01 ~ period + (1 | herd), family = binomial, data = cb)
  p <- suppressWarnings(rp(t1))
  expect_identical(attr(p, "data")$observed, as.numeric(cb$y01 == "TRUE"))
  expect_identical(gsub("\n", " ", p[[4]]$labels$title), "Calibration (tenths of fitted values)")
})

test_that("points keep the original row numbers when rows with missing values are dropped", {
  O <- oats_f()
  O$yield[c(3, 10)] <- NA
  O$nitro[20] <- NA
  ms <- list(lm = lm(yield ~ nitro + Block, O),
             lme = nlme::lme(yield ~ nitro, random = ~ 1 | Block, data = O, na.action = na.omit),
             gls = nlme::gls(yield ~ nitro, data = O, na.action = na.omit),
             lmer = lme4::lmer(yield ~ nitro + (1 | Block), data = O),
             tmb_pois = glmmTMB::glmmTMB(round(yield) ~ nitro + (1 | Block), family = poisson, data = O),
             gam = mgcv::gam(yield ~ nitro + s(Block, bs = "re"), data = O))
  keep <- setdiff(as.character(1:72), c("3", "10", "20"))
  for (nm in names(ms)) {
    d <- attr(suppressWarnings(rp(ms[[nm]])), "data")
    expect_identical(d$label, keep, label = nm)
  }
  # the largest residual is labelled with its own row of the data
  O$yield[50] <- 400
  d <- attr(rp(lme4::lmer(yield ~ nitro + (1 | Block), data = O)), "data")
  expect_identical(d$label[which.max(abs(d$resid))], "50")
})

test_that("gam fits show fitted vs observed, not Cook's distance", {
  O <- oats_f()
  g <- mgcv::gam(log(yield) ~ nitro + s(Block, bs = "re"), data = O)
  p <- rp(g)
  expect_identical(p[[4]]$labels$title, "Fitted vs observed")
  expect_identical(p[[4]]$labels$x, "log(yield)")
  expect_equal(attr(p, "data")$observed, log(O$yield))
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
  expect_identical(gsub("\n", " ", suppressMessages(rp(gb))[[4]]$labels$title), "Calibration (tenths of fitted values)")
  skip_if_not_installed("plotly")
  expect_s3_class(residplot(fm, interactive = TRUE, plot = FALSE), "plotly")
})

test_that("residplot(): many points as open transparent circles, bins on request, no legends; Q-Q line through all quartiles", {
  set.seed(8124)
  big <- data.frame(g = factor(sample(1:200, 6000, TRUE)), x = runif(6000), f = factor(sample(letters[1:3], 6000, TRUE)))
  big$y <- 2 + 0.5 * big$x + rnorm(200)[big$g] + rnorm(6000)
  fb <- lme4::lmer(y ~ x + (1 | g), data = big)
  geoms <- function(q) vapply(q$layers, function(l) class(l$geom)[1], "")
  p <- rp(fb)
  # default: every point, open circles, transparent
  lay <- Filter(function(l) inherits(l$geom, "GeomPoint"), p[[1]]$layers)[[1]]
  expect_equal(nrow(lay$data), 6000)
  expect_identical(lay$aes_params$shape, 1)
  expect_lt(lay$aes_params$alpha, 0.5)
  # small data: filled points, nearly opaque
  O <- oats_f()
  ls <- Filter(function(l) inherits(l$geom, "GeomPoint"), rp(lm(yield ~ Block + nitro * Variety, O))[[1]]$layers)[[1]]
  expect_identical(ls$aes_params$shape, 16)
  # no legends, also with group colours
  pg <- rp(fb, group = big$f)
  expect_true(all(vapply(pg, function(z) identical(z$theme$legend.position, "none"), logical(1))))
  # bins only when asked, without a colour bar, outstanding points on top
  pb <- rp(fb, dense = "bins")
  expect_true(any(c("GeomTile", "GeomBin2d") %in% geoms(pb[[1]])))
  fs <- Filter(function(s) "fill" %in% s$aesthetics, pb[[1]]$scales$scales)[[1]]
  expect_identical(fs$guide, "none")
  # Q-Q of many residuals: thinned in the centre, every tail point kept, line from all quartiles
  qd <- Filter(function(l) inherits(l$geom, "GeomPoint"), p[[2]]$layers)[[1]]$data
  expect_lt(nrow(qd), 6000)
  r <- attr(p, "data")$resid
  expect_true(all(sort(r)[c(1:100, 5901:6000)] %in% qd$resid))
  ab <- Filter(function(l) inherits(l$geom, "GeomAbline"), p[[2]]$layers)[[1]]$data
  expect_equal(unname(ab$slope), unname(diff(stats::quantile(r, c(0.25, 0.75))) / diff(stats::qnorm(c(0.25, 0.75)))))
  expect_error(rp(fb, dense = "hex"), "should be one of")
})

