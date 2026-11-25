# Cumulative link models (polr, clm): latent means, pairs and class probabilities against emmeans
# (decisions D1-D6, dev/review-classes.md).

ord_fits <- function() {
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  list(polr = MASS::polr(rating ~ temp * contact, data = wine, Hess = TRUE),
       clm = ordinal::clm(rating ~ temp * contact, data = wine),
       clm_probit_eq = ordinal::clm(rating ~ temp * contact, data = wine, link = "probit",
                                    threshold = "equidistant"))
}

test_that("latent means, SEs and pairs equal emmeans (mode = 'latent', z tests)", {
  skip_if_not_installed("emmeans")
  for (m in ord_fits()) {
    r <- predictmeans(m, "temp:contact", pairwise = TRUE)
    e <- as.data.frame(emmeans::emmeans(m, ~ temp * contact))
    e <- e[order(e$temp, e$contact), ]
    expect_equal(r$means$estimate, e$emmean, tolerance = 1e-6)
    expect_equal(r$means$std.error, e$SE, tolerance = 1e-6)
    expect_true(all(is.infinite(r$means$df)))
    ep <- as.data.frame(pairs(emmeans::emmeans(m, ~ temp * contact), adjust = "none", reverse = TRUE))
    key <- function(a, b) paste(pmin(a, b), pmax(a, b))
    lv <- function(x) gsub(" ", ":", x)
    ep_l <- do.call(rbind, strsplit(as.character(ep$contrast), " - "))
    i <- match(key(r$pairs$level1, r$pairs$level2), key(lv(ep_l[, 1]), lv(ep_l[, 2])))
    expect_false(anyNA(i))
    expect_equal(abs(r$pairs$estimate), abs(ep$estimate[i]), tolerance = 1e-6)
    expect_equal(r$pairs$p.value, ep$p.value[i], tolerance = 1e-6)
  }
})

test_that("class probabilities and the mean class equal emmeans for the cells; averaged on the latent scale", {
  skip_if_not_installed("emmeans")
  for (m in ord_fits()) {
    r <- predictmeans(m, "temp:contact", mean_type = "prob")
    ep <- as.data.frame(emmeans::emmeans(m, ~ rating | temp * contact, mode = "prob"))
    ep <- ep[order(ep$temp, ep$contact, ep$rating), ]
    expect_equal(r$probs$estimate, ep$prob, tolerance = 1e-6)
    expect_equal(r$probs$std.error, ep$SE, tolerance = 1e-6)
    expect_true(all(r$probs$conf.low > 0 & r$probs$conf.high < 1))
    expect_equal(as.vector(tapply(r$probs$estimate, paste(r$probs$temp, r$probs$contact), sum)), rep(1, 4))
    r <- predictmeans(m, "temp:contact", mean_type = "mean_class")
    em <- as.data.frame(emmeans::emmeans(m, ~ temp * contact, mode = "mean.class"))
    em <- em[order(em$temp, em$contact), ]
    expect_equal(r$probs$estimate, em$mean.class, tolerance = 1e-6)
    expect_equal(r$probs$std.error, em$SE, tolerance = 1e-6)
  }
  # over contact: probabilities at the latent mean (D2), by hand
  m <- ord_fits()$polr
  r <- predictmeans(m, "temp", mean_type = "prob")
  lat <- r$means$estimate                                  # x'b - mean(theta)
  th <- m$zeta
  P <- t(sapply(lat + mean(th), function(eta) diff(c(0, plogis(th - eta), 1))))
  expect_equal(r$probs$estimate, as.vector(t(P)), tolerance = 1e-10)
})

test_that("the scale is labelled; mean_type checks; nominal effects stop", {
  m <- ord_fits()$clm
  r <- predictmeans(m, "temp")
  expect_output(print(summary(r)), "latent mean of rating \\(logit\\); differences are log cumulative odds ratios")
  expect_output(print(summary(predictmeans(m, "temp", mean_type = "prob"))), "Class probabilities")
  expect_error(predictmeans(lm(breaks ~ tension, warpbreaks), "tension", mean_type = "prob"), "cumulative link")
  data(wine, package = "ordinal", envir = environment())
  mn <- ordinal::clm(rating ~ temp, nominal = ~ contact, data = wine)
  expect_error(predictmeans(mn, "temp"), "location effects only")
})
