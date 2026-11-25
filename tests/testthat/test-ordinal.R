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
  skip_on_cran()
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
    sgn <- ifelse(r$pairs$level1 == lv(ep_l[i, 1]), 1, -1)               # emmeans may list the pair the other way
    expect_equal(r$pairs$estimate, sgn * ep$estimate[i], tolerance = 1e-6)
    expect_equal(r$pairs$p.value, ep$p.value[i], tolerance = 1e-6)
  }
})

test_that("class probabilities and the mean class equal emmeans for the cells; averaged on the latent scale", {
  skip_on_cran()
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
    expect_true(all(r$probs$conf.low >= 1 & r$probs$conf.high <= 5))
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
  skip_on_cran()
  m <- ord_fits()$clm
  r <- predictmeans(m, "temp")
  expect_output(print(summary(r)), "latent mean of rating \\(logit\\); differences are log odds ratios of a higher class")
  expect_output(print(summary(predictmeans(m, "temp", mean_type = "prob"))), "Class probabilities")
  expect_error(predictmeans(lm(breaks ~ tension, warpbreaks), "tension", mean_type = "prob"), "cumulative link")
  data(wine, package = "ordinal", envir = environment())
  mn <- ordinal::clm(rating ~ temp, nominal = ~ contact, data = wine)
  expect_error(predictmeans(mn, "temp"), "location effects only")
})

test_that("a positive difference of latent means: higher classes more likely (log odds of P(Y > j))", {
  data(wine, package = "ordinal", envir = environment())
  m <- MASS::polr(rating ~ temp, data = wine, Hess = TRUE)
  r <- predictmeans(m, "temp", pairwise = TRUE)
  P <- predict(m, newdata = data.frame(temp = c("cold", "warm")), type = "probs")
  lo <- function(p) log((1 - cumsum(p)[1:4]) / cumsum(p)[1:4])        # log odds of P(Y > j), j = 1..4
  d <- r$means$estimate[r$means$temp == "warm"] - r$means$estimate[r$means$temp == "cold"]
  expect_equal(unname(lo(P[2, ]) - lo(P[1, ])), rep(d, 4), tolerance = 1e-8)
  expect_output(print(summary(r)), "log odds ratios of a higher class")
})

test_that("atvar keeps each class row with its mean; aliased clm coefficients; links; stops", {
  skip_on_cran()
  data(wine, package = "ordinal", envir = environment())
  m <- MASS::polr(rating ~ temp * contact, data = wine, Hess = TRUE)
  a <- predictmeans(m, "temp:contact", atvar = "contact", mean_type = "mean_class")$probs
  b <- predictmeans(m, "temp:contact", mean_type = "mean_class")$probs
  ka <- paste(a$temp, a$contact); kb <- paste(b$temp, b$contact)
  expect_equal(a$estimate[match(kb, ka)], b$estimate)
  # an empty cell: the interaction coefficient aliased (NA), as for polr
  w <- wine[!(wine$temp == "warm" & wine$contact == "yes"), ]
  mc <- ordinal::clm(rating ~ temp * contact, data = w)
  mp <- suppressWarnings(MASS::polr(rating ~ temp * contact, data = w, Hess = TRUE))
  rc <- suppressWarnings(predictmeans(mc, "temp"))
  rp <- suppressWarnings(predictmeans(mp, "temp"))
  expect_equal(rc$means$estimate, rp$means$estimate, tolerance = 1e-4)
  # loglog: predict.polr is the reference (emmeans differs for polr's loglog)
  ml <- MASS::polr(rating ~ temp * contact, data = wine, Hess = TRUE, method = "loglog")
  r <- predictmeans(ml, "temp:contact", mean_type = "prob")$probs
  P <- predict(ml, newdata = unique(r[c("temp", "contact")]), type = "probs")
  expect_equal(r$estimate, as.vector(t(P)), tolerance = 1e-6)
  # stops
  expect_error(predictmeans(MASS::polr(rating ~ temp, data = wine), "temp"), "Hess = TRUE")
  expect_error(predictmeans(m, "temp", trans = exp), "'trans' cannot be used")
  expect_error(predictmeans(m, "temp", bias_adjust = TRUE), "'bias_adjust' cannot be used")
  ma <- suppressWarnings(ordinal::clm(rating ~ temp, data = wine, link = "Aranda-Ordaz"))
  expect_error(predictmeans(ma, "temp"), "lambda")
})

test_that("advice: probabilities at averaged latent means; mean class spacing; plots without data", {
  data(wine, package = "ordinal", envir = environment())
  m <- ordinal::clm(rating ~ temp + contact, data = wine)
  expect_true("ordinal_latent_average" %in% advice(predictmeans(m, "temp", mean_type = "prob"))$id)
  expect_false("ordinal_latent_average" %in% advice(predictmeans(m, "temp:contact", mean_type = "prob"))$id)
  expect_true("ordinal_mean_class" %in% advice(predictmeans(m, "temp:contact", mean_type = "mean_class"))$id)
  expect_message(ggplot2::autoplot(predictmeans(m, "temp"), data = TRUE), "not shown with latent means")
})

test_that("clmm: latent means, pairs and class probabilities for a typical subject equal emmeans", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  data(wine, package = "ordinal", envir = environment())
  m <- ordinal::clmm(rating ~ temp * contact + (1 | judge), data = wine)
  r <- predictmeans(m, "temp:contact", pairwise = TRUE)
  e <- as.data.frame(emmeans::emmeans(m, ~ temp * contact))
  e <- e[order(e$temp, e$contact), ]
  expect_equal(r$means$estimate, e$emmean, tolerance = 1e-6)
  expect_equal(r$means$std.error, e$SE, tolerance = 1e-6)
  ep <- as.data.frame(pairs(emmeans::emmeans(m, ~ temp * contact), adjust = "none"))
  expect_equal(sort(abs(r$pairs$estimate)), sort(abs(ep$estimate)), tolerance = 1e-6)
  expect_equal(sort(r$pairs$p.value), sort(ep$p.value), tolerance = 1e-6)
  p <- predictmeans(m, "temp:contact", mean_type = "prob")
  eq <- as.data.frame(emmeans::emmeans(m, ~ rating | temp * contact, mode = "prob"))
  eq <- eq[order(eq$temp, eq$contact, eq$rating), ]
  expect_equal(p$probs$estimate, eq$prob, tolerance = 1e-6)
  expect_equal(p$probs$std.error, eq$SE, tolerance = 1e-6)
  expect_true("ordinal_typical_subject" %in% advice(p)$id)
  expect_false("ordinal_latent_average" %in% advice(p)$id)              # judge is a random term, not averaged over
  expect_output(print(summary(p)), "random effects at 0 \\(a typical subject\\)")
  expect_error(predictmeans(m, "temp", bias_adjust = TRUE), "'bias_adjust' cannot be used")
  expect_error(se_ranef(m), "not covered")
})
