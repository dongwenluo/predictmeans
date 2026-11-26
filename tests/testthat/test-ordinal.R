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

test_that("the scale is labelled; mean_type checks; nominal effects give class probabilities", {
  skip_on_cran()
  m <- ord_fits()$clm
  r <- predictmeans(m, "temp")
  expect_output(print(summary(r)), "latent mean of rating \\(logit\\); differences are log odds ratios of a higher class")
  expect_output(print(summary(predictmeans(m, "temp", mean_type = "prob"))), "Class probabilities")
  expect_error(predictmeans(lm(breaks ~ tension, warpbreaks), "tension", mean_type = "prob"), "models of a class response")
  data(wine, package = "ordinal", envir = environment())
  mn <- ordinal::clm(rating ~ temp, nominal = ~ contact, data = wine)
  expect_identical(predictmeans(mn, "temp")$settings$response$classes$what, "prob")    # (R/categorical.R)
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
  expect_output(print(summary(r)), "p-values: z tests")
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

test_that("clmm: latent means and class probabilities equal emmeans; inner-outer df for the pairs", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  data(wine, package = "ordinal", envir = environment())
  m <- ordinal::clmm(rating ~ temp * contact + (1 | judge), data = wine)
  r <- predictmeans(m, "temp:contact", pairwise = TRUE)
  e <- as.data.frame(emmeans::emmeans(m, ~ temp * contact))
  e <- e[order(e$temp, e$contact), ]
  expect_equal(r$means$estimate, e$emmean, tolerance = 1e-6)
  expect_equal(r$means$std.error, e$SE, tolerance = 1e-6)
  # temp and contact vary within judges: the df of the within-judge stratum, N - rank(X, Z) less the
  # 3 extra cut points: 57, as SAS GLIMMIX (dev/sas/ordinal.sas)
  expect_equal(r$settings$df_method, "inner-outer")
  expect_equal(unique(r$pairs$df), 72 - (4 + 9 - 1) - 3)
  expect_equal(r$pairs$p.value, 2 * pt(-abs(r$pairs$estimate / r$pairs$std.error), 57), tolerance = 1e-10)
  ep <- as.data.frame(pairs(emmeans::emmeans(m, ~ temp * contact), adjust = "none"))
  expect_equal(sort(abs(r$pairs$estimate)), sort(abs(ep$estimate)), tolerance = 1e-6)
  p <- predictmeans(m, "temp:contact", mean_type = "prob")
  eq <- as.data.frame(emmeans::emmeans(m, ~ rating | temp * contact, mode = "prob"))
  eq <- eq[order(eq$temp, eq$contact, eq$rating), ]
  expect_equal(p$probs$estimate, eq$prob, tolerance = 1e-6)
  expect_equal(p$probs$std.error, eq$SE, tolerance = 1e-6)
  expect_true("ordinal_typical_subject" %in% advice(p)$id)
  expect_match(advice(p)$text[advice(p)$id == "ordinal_typical_subject"], "random-effect SD 1.1")
  expect_false("ordinal_latent_average" %in% advice(p)$id)              # judge is a random term, not averaged over
  expect_output(print(summary(p)), "within a subject \\(random effects at 0\\)")
})

test_that("clmm: a treatment between clusters gets the df of the clusters", {
  skip_on_cran()
  set.seed(3)
  d <- data.frame(g = factor(rep(1:6, each = 10)), trt = factor(rep(c("A", "B"), each = 30)))
  d$y <- cut(rnorm(6)[d$g] + rlogis(60), c(-Inf, -1, 0.5, Inf), ordered_result = TRUE)
  m <- suppressWarnings(ordinal::clmm(y ~ trt + (1 | g), data = d))
  r <- predictmeans(m, "trt", pairwise = TRUE)
  expect_equal(r$pairs$df, 4)                                               # 6 clusters - 2
  expect_output(print(summary(r)), "t tests")
})

test_that("clmm: random slopes, crossed effects, symmetric thresholds, sum contrasts; stops", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  data(wine, package = "ordinal", envir = environment())
  fits <- list(ordinal::clmm(rating ~ temp + contact + (1 + contact | judge), data = wine),
               ordinal::clmm(rating ~ temp + contact + (1 | judge) + (1 | bottle), data = wine),
               ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine, threshold = "symmetric"),
               ordinal::clmm(rating ~ temp * contact + (1 | judge), data = wine,
                             contrasts = list(temp = "contr.sum", contact = "contr.sum")))
  for (m in fits) {
    r <- predictmeans(m, "temp:contact")
    e <- as.data.frame(emmeans::emmeans(m, ~ temp * contact))
    e <- e[order(e$temp, e$contact), ]
    expect_equal(r$means$estimate, e$emmean, tolerance = 1e-6)
    expect_equal(r$means$std.error, e$SE, tolerance = 1e-6)
  }
  m <- fits[[3]]
  expect_error(predictmeans(m, "temp", bias_adjust = TRUE), "'bias_adjust' cannot be used")
  expect_error(se_ranef(m), "not covered")
  wine$x <- seq_len(nrow(wine)) / 10
  mx <- ordinal::clmm(rating ~ temp + x + (1 | judge), data = wine)
  expect_error(covariatemeans(mx, "temp", covariate = "x", trans = plogis, plot = FALSE), "'trans' cannot be used")
  mn <- suppressWarnings(ordinal::clmm(rating ~ temp + (1 | judge), nominal = ~ contact, data = wine))
  expect_error(predictmeans(mn, "temp"), "does not fit nominal effects")
  expect_warning(r <- predictmeans(m, "temp", df_method = "residual"), "using z")
  expect_true(all(is.infinite(r$means$df)))
})

test_that("vglm (parallel): the same model as clm; exact SEs on vglm's covariance; probabilities; the sign", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  ctl <- VGAM::vglm.control(epsilon = 1e-12)
  cm <- ordinal::clm(rating ~ temp * contact, data = wine)
  rc <- predictmeans(cm, "temp:contact", pairwise = TRUE)
  pc <- predictmeans(cm, "temp:contact", mean_type = "prob")$probs
  for (fam in list(VGAM::cumulative(parallel = TRUE), VGAM::cumulative(parallel = TRUE, reverse = TRUE),
                   VGAM::propodds())) {
    m <- VGAM::vglm(rating ~ temp * contact, fam, data = wine, control = ctl)
    r <- predictmeans(m, "temp:contact", pairwise = TRUE)
    expect_equal(r$means$estimate, rc$means$estimate, tolerance = 1e-6)  # the same model and MLE
    expect_equal(r$pairs$estimate, rc$pairs$estimate, tolerance = 1e-6)
    # SEs: exactly K V K' on vglm's own covariance (expected information; clm's observed: within 2%)
    lat <- predictmeans:::pm_ord_latent(m)
    K <- Kmatrix(m, "temp:contact")$K[, names(lat$coef)]
    expect_equal(r$means$std.error, unname(sqrt(diag(K %*% lat$vcov %*% t(K)))), tolerance = 1e-10)
    expect_equal(r$means$std.error, rc$means$std.error, tolerance = 0.02)
    expect_equal(predictmeans(m, "temp:contact", mean_type = "prob")$probs$estimate, pc$estimate, tolerance = 1e-6)
    expect_output(print(summary(r)), "log odds ratios of a higher class")
  }
  # other links: cloglog with reverse = TRUE is clm's loglog
  ml <- VGAM::vglm(rating ~ temp * contact, VGAM::cumulative(link = "clogloglink", parallel = TRUE, reverse = TRUE),
                   data = wine, control = ctl)
  rl <- predictmeans(ml, "temp:contact")
  expect_equal(rl$means$estimate, predictmeans(ordinal::clm(rating ~ temp * contact, data = wine, link = "loglog"),
                                                "temp:contact")$means$estimate, tolerance = 1e-6)
  expect_equal(rl$settings$response$link, "loglog")
  mp <- VGAM::vglm(rating ~ temp * contact, VGAM::cumulative(link = "probitlink", parallel = TRUE), data = wine,
                   control = ctl)
  expect_equal(predictmeans(mp, "temp:contact")$means$estimate,
               predictmeans(ordinal::clm(rating ~ temp * contact, data = wine, link = "probit"),
                            "temp:contact")$means$estimate, tolerance = 1e-6)
  # equidistant cut points (constraints on the intercepts) = clm(threshold = "equidistant")
  me <- VGAM::vglm(rating ~ temp + contact, VGAM::cumulative(parallel = TRUE), data = wine, control = ctl,
                   constraints = list("(Intercept)" = cbind(1, 0:3), temp = matrix(1, 4, 1), contact = matrix(1, 4, 1)))
  ce <- ordinal::clm(rating ~ temp + contact, data = wine, threshold = "equidistant")
  expect_equal(predictmeans(me, "temp")$means$estimate, predictmeans(ce, "temp")$means$estimate, tolerance = 1e-6)
  expect_equal(predictmeans(me, "temp", mean_type = "prob")$probs$estimate,
               predictmeans(ce, "temp", mean_type = "prob")$probs$estimate, tolerance = 1e-6)
  # acat, cratio, sratio, both directions: cell probabilities = predict(); the sign
  for (fam in list(VGAM::acat(parallel = TRUE), VGAM::acat(parallel = TRUE, reverse = TRUE),
                   VGAM::cratio(parallel = TRUE), VGAM::cratio(parallel = TRUE, reverse = TRUE),
                   VGAM::sratio(parallel = TRUE), VGAM::sratio(parallel = TRUE, reverse = TRUE))) {
    m <- VGAM::vglm(rating ~ temp * contact, fam, data = wine)
    p <- predictmeans(m, "temp:contact", mean_type = "prob")$probs
    P <- VGAM::predict(m, newdata = unique(p[c("temp", "contact")]), type = "response")
    expect_equal(p$estimate, as.vector(t(P)), tolerance = 1e-8)
    r <- predictmeans(m, "temp")
    mc <- predictmeans(m, "temp", mean_type = "mean_class")$probs
    expect_true(diff(r$means$estimate) > 0 && diff(mc$estimate) > 0)
  }
  m <- VGAM::vglm(rating ~ temp + contact, VGAM::acat(parallel = TRUE), data = wine)
  expect_output(print(summary(predictmeans(m, "temp"))), "adjacent-category log odds ratios")
})

test_that("vglm: probability SEs for a covariate with large values (analytic delta method); stops", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  wine$x <- 5000 + seq_len(nrow(wine))                       # e.g. a year
  m <- VGAM::vglm(rating ~ temp + x, VGAM::cumulative(parallel = TRUE), data = wine)
  p <- predictmeans(m, "temp", covariate = list(x = 5040), mean_type = "prob")$probs
  # analytic: P(Y <= j) = plogis(alpha_j + x'b) (vglm's sign), at cold, x = 5040
  cf <- VGAM::coef(m); V <- VGAM::vcov(m)
  a <- cf[1:4]; eta <- a + cf["x"] * 5040
  f <- dlogis(eta)
  Gc <- cbind(diag(f), 0, f * 5040)                          # d P(Y <= j) / d (alpha, tempwarm, x)
  G <- rbind(Gc, 0) - rbind(0, Gc)
  expect_equal(p$std.error[1:5], unname(sqrt(diag(G %*% V %*% t(G)))), tolerance = 1e-6)
  # multinomial and non-parallel fits: predictmeans() only (test-categorical.R)
  expect_error(Kmatrix(suppressWarnings(VGAM::vglm(rating ~ temp, VGAM::cumulative(parallel = FALSE),
                                                   data = wine)), "temp"), "vector of means")
  ms <- VGAM::vglm(rating ~ temp + VGAM::sm.scale(x), VGAM::cumulative(parallel = TRUE), data = wine)
  expect_error(predictmeans(ms, "temp"), "smart terms")
  expect_error(residplot(m), "do not yet handle models of a class response")
  mpo <- VGAM::vglm(as.numeric(rating) ~ temp, VGAM::poissonff, data = wine)
  expect_error(predictmeans(mpo, "temp"), "not poissonff")
  expect_error(residplot(mpo), "not poissonff")
})

test_that("clmm population-average class probabilities (bias_adjust = TRUE): quadrature, probit, SEs, both routes", {
  skip_on_cran()
  data(wine, package = "ordinal", envir = environment())
  m <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine)
  pr <- predictmeans(m, "temp", mean_type = "prob", bias_adjust = TRUE)$probs
  expect_true(attr(pr, "population_average"))
  s <- m$ST$judge[1, 1]
  # the latent mean of each temp (contact averaged on the latent scale), integrated over the judges
  xb <- c(0, m$beta[["tempwarm"]]) + m$beta[["contactyes"]] / 2
  ref <- unlist(lapply(xb, function(v) diff(c(0, vapply(m$alpha, function(t) stats::integrate(function(u)
    stats::plogis(t - v - u) * stats::dnorm(u, 0, s), -Inf, Inf, rel.tol = 1e-12)$value, 0), 1))))
  expect_equal(pr$estimate, unname(ref), tolerance = 1e-10)
  # delta-method SEs over (alpha, beta, log SD): vcov() of the default fit is for the log SD
  k <- 1:59                                           # Gauss-Hermite, 60 nodes (Golub-Welsch)
  J <- matrix(0, 60, 60)
  J[cbind(k, k + 1)] <- J[cbind(k + 1, k)] <- sqrt(k)
  e <- eigen(J, symmetric = TRUE)
  gh <- list(x = e$values, w = e$vectors[1, ]^2)
  f <- function(th) unlist(lapply(c(0, th[5]) + th[6] / 2, function(v)
    diff(c(0, vapply(th[1:4], function(t) sum(gh$w * stats::plogis(t - v - exp(th[7]) * gh$x)), 0), 1))))
  D <- numDeriv::jacobian(f, c(m$alpha, m$beta, log(s)))
  expect_equal(pr$std.error, sqrt(diag(D %*% stats::vcov(m) %*% t(D))), tolerance = 1e-8)
  # useMatrix = TRUE: vcov() for the SD itself; the same probabilities and SEs (to the fits' agreement)
  m2 <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine,
                      control = ordinal::clmm.control(useMatrix = TRUE))
  p2 <- predictmeans(m2, "temp", mean_type = "prob", bias_adjust = TRUE)$probs
  expect_equal(p2$std.error, pr$std.error, tolerance = 1e-4)
  # probit: exact closed form Phi((theta - x'b) / sqrt(1 + s^2))
  mp <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine, link = "probit")
  pp <- predictmeans(mp, "temp", mean_type = "prob", bias_adjust = TRUE)$probs
  sp <- mp$ST$judge[1, 1]
  xp <- c(0, mp$beta[["tempwarm"]]) + mp$beta[["contactyes"]] / 2
  expect_equal(pp$estimate, unname(unlist(lapply(xp, function(v) diff(c(0, stats::pnorm((mp$alpha - v) / sqrt(1 + sp^2)), 1))))),
               tolerance = 1e-12)
  # less extreme than for a typical subject; the advice says which it is
  pt <- predictmeans(m, "temp", mean_type = "prob")
  expect_true(max(pt$probs$estimate) > max(pr$estimate))
  pa <- predictmeans(m, "temp", mean_type = "mean_class", bias_adjust = TRUE)
  expect_match(advice(pa)$text[advice(pa)$id == "ordinal_typical_subject"], "averages over subjects")
  expect_output(print(summary(pa)), "averaged over the random effects")
  # the latent means: bias_adjust does not apply; random slopes stop
  expect_error(predictmeans(m, "temp", bias_adjust = TRUE), "population-average class probabilities")
  wine$x <- as.numeric(wine$temp)
  ms <- suppressWarnings(ordinal::clmm(rating ~ temp + contact + (1 + x | judge), data = wine))
  expect_error(predictmeans(ms, "temp", mean_type = "prob", bias_adjust = TRUE), "random slopes")
})

test_that("population-average integral: any SD and link; scalar random slopes stop; sigma refused", {
  for (lk in c("logit", "cloglog")) {
    Fd <- pm_ord_link_fd(lk)
    z <- c(-6, 0, 3)
    r <- pm_ord_pop_int(z, 5, Fd)
    ref <- vapply(z, function(z0) stats::integrate(function(u) Fd[[1]](z0 - 5 * u) * stats::dnorm(u), -Inf, Inf,
                                                   rel.tol = 1e-13, subdivisions = 2000L)$value, 0)
    expect_equal(r$F, ref, tolerance = 1e-10)
  }
  # tiny SDs (integrate() stopped with "roundoff error"): the series F(z) + s^2 f'(z) / 2, d/ds = s f'(z)
  Fl <- pm_ord_link_fd("logit")
  for (s in c(1e-7, 1e-4, 1e-3)) {
    r <- pm_ord_pop_int(c(-0.2, 0.7), s, Fl)
    fp <- stats::dlogis(c(-0.2, 0.7)) * (1 - 2 * stats::plogis(c(-0.2, 0.7)))
    expect_equal(r$F, stats::plogis(c(-0.2, 0.7)) + s^2 * fp / 2, tolerance = 1e-12)
    expect_equal(r$ds, s * fp, tolerance = 1e-6)
  }
  # the probit closed form at small and large SDs (a small SD used to lose the mass: 0/1)
  for (s in c(0.002, 0.05, 0.1, 1, 5)) {
    z <- c(-5, -1, 0, 2, 5, 40)
    expect_equal(pm_ord_pop_int(z, s, pm_ord_link_fd("probit"))$F, stats::pnorm(z / sqrt(1 + s^2)), tolerance = 1e-10,
                 info = paste("s =", s))
  }
  skip_on_cran()
  data(wine, package = "ordinal", envir = environment())
  set.seed(2)
  wine$x <- 3 * (wine$contact == "yes") + stats::rnorm(nrow(wine), 0, 0.3)
  ms <- suppressWarnings(ordinal::clmm(rating ~ temp + contact + (0 + x | judge), data = wine))
  expect_error(predictmeans(ms, "temp", mean_type = "prob", bias_adjust = TRUE), "random slopes")
  expect_false(grepl("bias_adjust", advice(predictmeans(ms, "temp", mean_type = "prob"))$text[1]))
  m <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine)
  expect_error(predictmeans(m, "temp", mean_type = "prob", bias_adjust = TRUE, sigma = 1), "'sigma' is not used")
  p <- predictmeans(m, "temp", mean_type = "prob", bias_adjust = TRUE)
  expect_equal(p$settings$bias_adjust, m$ST$judge[1, 1])
})
