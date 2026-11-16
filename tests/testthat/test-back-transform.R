## Back-transformation (R/back-transform.R): the response transformation in
## the formula, Box-Cox and Yeo-Johnson, the mean of the back-transformed
## response, and the Box-Cox check of pm_check().

bt_oats <- function() transform(as.data.frame(nlme::Oats), nitro = factor(nitro))

test_that("the transformation of the response in the formula is found, with its inverse", {
  y <- c(0.5, 2, 7, 30)
  chk <- function(expr, fwd) {
    inv <- pm_parse_trans(expr)
    expect_false(is.null(inv), label = deparse(expr))
    expect_equal(inv(fwd(y)), y, tolerance = 1e-10, label = deparse(expr))
  }
  chk(quote(log(y)), log)
  expect_identical(pm_parse_trans(quote(log(y))), exp)               # exp itself: ratios in the pairs
  chk(quote(log(y + 1)), function(y) log(y + 1))
  chk(quote(log(1 + y)), log1p)
  chk(quote(log(y - 0.1)), function(y) log(y - 0.1))
  chk(quote(log(y, 2)), log2)
  chk(quote(log10(y)), log10)
  chk(quote(log2(y)), log2)
  chk(quote(log1p(y)), log1p)
  chk(quote(sqrt(y)), sqrt)
  chk(quote(sqrt(y + 0.5)), function(y) sqrt(y + 0.5))
  chk(quote(y^0.25), function(y) y^0.25)
  chk(quote(I(y^2)), function(y) y^2)
  chk(quote(1/y), function(y) 1/y)
  chk(quote(100/(y + 1)), function(y) 100/(y + 1))
  chk(quote(pm_boxcox(y, 0.3)), function(y) pm_boxcox(y, 0.3))
  chk(quote(pm_boxcox(y, 0)), log)
  chk(quote(predictmeans::pm_yeojohnson(y, -0.5)), function(y) pm_yeojohnson(y, -0.5))
  p <- c(0.1, 0.4, 0.9)
  expect_equal(pm_parse_trans(quote(qlogis(p)))(qlogis(p)), p)
  expect_equal(pm_parse_trans(quote(asin(sqrt(p))))(asin(sqrt(p))), p)
  # not transformations of one response
  for (e in list(quote(y), quote(cbind(a, b)), quote(log(y/x)), quote(y^x), quote(scale(y)), quote(log(y) + 1)))
    expect_null(pm_parse_trans(e), label = deparse(e))
  # Yeo-Johnson takes zero and negative values
  v <- c(-3, -0.5, 0, 0.5, 4)
  for (l in c(-1, 0, 0.5, 1, 2, 2.5)) expect_equal(pm_yeojohnson_inv(l)(pm_yeojohnson(v, l)), v, tolerance = 1e-10)
  expect_equal(pm_yeojohnson(v, 1), v)
  expect_error(pm_boxcox(c(1, 0), 0.5), "positive")
})

test_that("predictmeans() back-transforms a transformed response by default; trans = I keeps the model scale", {
  O <- bt_oats()
  m <- lm(log(yield) ~ Block + nitro, data = O)
  expect_message(r <- predictmeans(m, "nitro", pairwise = TRUE, prtnum = FALSE), "The response is log\\(yield\\)")
  expect_equal(r$means$response, exp(r$means$estimate))
  expect_true("Bk_Mean" %in% names(r$mean_table$Table))
  expect_true("ratio" %in% names(r$pairs))                              # differences of logs as ratios
  expect_identical(r$settings$back_transform$estimate, "median of the original scale (bias_adjust = TRUE gives the mean)")
  expect_match(r$settings$back_transform$by, "inverse of log\\(yield\\)")
  # the same as giving trans
  r2 <- predictmeans(m, "nitro", pairwise = TRUE, trans = exp, prtnum = FALSE)
  expect_equal(r$means, r2$means)
  expect_silent(r0 <- predictmeans(m, "nitro", trans = I, prtnum = FALSE))
  expect_null(r0$means$response)
  expect_null(r0$settings$back_transform)
  # an offset in the formula; other transformations
  r1 <- suppressMessages(predictmeans(lm(log(yield + 10) ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_equal(r1$means$response, exp(r1$means$estimate) - 10)
  rs <- suppressMessages(predictmeans(lm(sqrt(yield) ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_equal(rs$means$response, rs$means$estimate^2)
  rb <- suppressMessages(predictmeans(lm(pm_boxcox(yield, 0.3) ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_equal(rb$means$response, (0.3*rb$means$estimate + 1)^(1/0.3))
  # an untransformed response, and a link: as before
  expect_silent(ru <- predictmeans(lm(yield ~ Block + nitro, data = O), "nitro", prtnum = FALSE))
  expect_null(ru$means$response)
  # mixed model, and covariatemeans()
  rl <- suppressMessages(predictmeans(lme4::lmer(log(yield) ~ nitro + (1 | Block), data = O), "nitro", prtnum = FALSE))
  expect_equal(rl$means$response, exp(rl$means$estimate))
  O$nit <- as.numeric(as.character(O$nitro))
  expect_message(cm <- covariatemeans(lm(log(yield) ~ nit * Variety, data = O), "Variety", covariate = "nit", plot = FALSE),
                 "back-transformed")
  # the back-transformed fitted values (1.x: exp(fitted + SE/2)): between the limits, at their geometric mean
  mfit <- lm(log(yield) ~ nit * Variety, data = O)
  nd <- data.frame(nit = cm$pltdf$xvar, Variety = cm$pltdf$Variety)
  expect_equal(cm$pltdf$Mean, unname(exp(predict(mfit, nd))), tolerance = 1e-8)
  expect_equal(cm$pltdf$Mean, sqrt(cm$pltdf$LL*cm$pltdf$UL), tolerance = 1e-8)
})

test_that("bias_adjust gives the mean of the back-transformed response: exact for log and sqrt", {
  O <- bt_oats()
  ml <- lm(log(yield) ~ Block + nitro, data = O)
  s2 <- sigma(ml)^2
  r <- suppressMessages(predictmeans(ml, "nitro", bias_adjust = TRUE, prtnum = FALSE))
  expect_equal(r$means$response, exp(r$means$estimate + s2/2), tolerance = 1e-9)
  expect_equal(r$means$response.std.error, exp(r$means$estimate + s2/2)*r$means$std.error, tolerance = 1e-9)
  expect_identical(r$settings$back_transform$estimate, "mean (bias adjusted)")
  ms <- lm(sqrt(yield) ~ Block + nitro, data = O)
  rs <- suppressMessages(predictmeans(ms, "nitro", bias_adjust = TRUE, prtnum = FALSE))
  expect_equal(rs$means$response, rs$means$estimate^2 + sigma(ms)^2, tolerance = 1e-9)
  # Gauss-Hermite: E f(Z) for a standard normal Z
  gh <- pm_gauss_hermite()
  expect_equal(sum(gh$w), 1)
  expect_equal(sum(gh$w*gh$z^4), 3, tolerance = 1e-10)
  expect_equal(sum(gh$w*exp(1.5*gh$z)), exp(1.5^2/2), tolerance = 1e-10)
  # a logit with a large random effect: the integral, where the second-order formula is off
  eta <- c(-1, 0.3); s <- 1.2
  bf <- pm_back_funs(stats::plogis, stats::plogis, NULL, FALSE)
  exact <- pm_bias_funs(bf, s^2)$h(eta)
  truth <- sapply(eta, function(e) integrate(function(z) plogis(e + s*z)*dnorm(z), -Inf, Inf)$value)
  expect_equal(exact, truth, tolerance = 1e-6)
  expect_gt(max(abs(pm_bias_funs(bf, s^2, "taylor")$h(eta) - truth)), 0.005)
  # a reciprocal near zero has no finite mean: the second-order formula, with a warning
  bfr <- pm_back_funs(function(x) 1/x, function(x) 1/x, NULL, FALSE)
  expect_warning(pm_bias_funs(bfr, 4)$h(0.5), "second-order formula")
})

test_that("summary() says what the back-transformed means are", {
  O <- bt_oats()
  out <- function(m, ...) capture.output(print(summary(suppressMessages(predictmeans(m, "nitro", prtnum = FALSE, ...)))))
  expect_true(any(grepl("'response' is the median of the original scale", out(lm(log(yield) ~ Block + nitro, data = O)))))
  O$cnt <- round(O$yield/10)
  expect_true(any(grepl("the inverse link: 'response' is the mean\\.", out(glm(cnt ~ Block + nitro, poisson, data = O)))))
  gm <- suppressMessages(lme4::glmer(cnt ~ nitro + (1 | Block), poisson, data = O))
  expect_true(any(grepl("'response' is the typical cluster", out(gm))))
  expect_true(any(grepl("'response' is the mean \\(bias adjusted\\)", out(gm, bias_adjust = TRUE))))
})

test_that("pm_check() gives the Box-Cox power the residuals suggest", {
  set.seed(2)
  d <- data.frame(A = factor(rep(1:4, each = 15)))
  d$y <- exp(1 + 0.3*as.numeric(d$A) + rnorm(60, sd = 0.5))
  m <- lm(y ~ A, d)
  bc <- pm_boxcox_profile(m)
  ref <- MASS::boxcox(y ~ A, data = d, lambda = seq(-3, 3, 0.01), plotit = FALSE)
  expect_equal(bc$lambda, ref$x[which.max(ref$y)])
  expect_equal(c(bc$lower, bc$upper), range(ref$x[ref$y > max(ref$y) - qchisq(0.95, 1)/2]), tolerance = 0.011)
  expect_identical(bc$name, "log(y)")
  pc <- pm_check(m)
  row <- pc[pc$check == "Transformation of the response (Box-Cox)", ]
  expect_true(row$look)
  expect_match(row$note, "log\\(y\\) is within the interval")
  # no transformation needed, a transformed response, a glm: no flag or no row
  O <- bt_oats()
  p0 <- pm_check(lm(yield ~ Block + nitro, data = O))
  expect_false(p0$look[p0$check == "Transformation of the response (Box-Cox)"])
  expect_null(pm_boxcox_profile(lm(log(yield) ~ Block + nitro, data = O)))
  expect_null(pm_boxcox_profile(glm(round(yield) ~ nitro, poisson, data = O)))
})

test_that("covariatemeans(contr = ): the difference between two lines with its own band", {
  O <- as.data.frame(nlme::Oats)
  fm <- lm(yield ~ Block + nitro * Variety, data = O)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  r <- covariatemeans(fm, "Variety", covariate = "nitro", contr = c("Golden Rain", "Victory"), covariateV = c(0, 0.3, 0.6))
  expect_s3_class(r$plt, "ggplot")
  b <- coef(fm)
  # Golden Rain is the first level: the difference is -(Victory + nitro x Variety Victory)
  expect_equal(r$pltdf$diff, -(b[["VarietyVictory"]] + c(0, 0.3, 0.6)*b[["nitro:VarietyVictory"]]), tolerance = 1e-10)
  L <- matrix(0, 3, length(b), dimnames = list(NULL, names(b)))
  L[, "VarietyVictory"] <- -1; L[, "nitro:VarietyVictory"] <- -c(0, 0.3, 0.6)
  expect_equal(r$pltdf$ses, sqrt(diag(L %*% vcov(fm) %*% t(L))), tolerance = 1e-10)
  expect_equal(r$pltdf$UL - r$pltdf$diff, qt(0.975, df.residual(fm))*r$pltdf$ses, tolerance = 1e-10)
  # by number, and a log response gives ratios
  r2 <- covariatemeans(fm, "Variety", covariate = "nitro", contr = c(1, 3), covariateV = c(0, 0.3, 0.6))
  expect_equal(r2$pltdf, r$pltdf)
  rl <- covariatemeans(lm(log(yield) ~ Block + nitro * Variety, data = O), "Variety", covariate = "nitro", contr = c(1, 3),
                       covariateV = c(0, 0.6))
  expect_equal(rl$pltdf$ratio, exp(rl$pltdf$diff))
  expect_error(covariatemeans(fm, "Variety", covariate = "nitro", contr = c(1, 1)), "two different levels")
  expect_error(covariatemeans(fm, "Variety", covariate = "nitro", contr = "Golden"), "two levels")
})

test_that("covariatemeans(contr = ) for GLMs and GLMMs: both scales equal emmeans; no aov with Error()", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  set.seed(1)
  d <- expand.grid(A = factor(c("a1", "a2", "a3")), Blk = factor(1:8), r = 1:4)
  d$x <- runif(nrow(d), 0, 2); u <- rnorm(8, sd = 0.5)
  d$cnt <- rpois(nrow(d), exp(0.2 + 0.3*as.numeric(d$A) + 0.4*d$x + u[d$Blk]))
  d$bin <- rbinom(nrow(d), 1, plogis(-1 + 0.5*as.numeric(d$A) + 0.6*d$x + u[d$Blk]))
  d$y <- rnorm(nrow(d)) + u[d$Blk]
  xs <- c(0.2, 1, 1.8)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  mods <- list(glm(cnt ~ A * x, poisson, d), glm(bin ~ A * x, binomial, d),
               suppressMessages(lme4::glmer(cnt ~ A * x + (1 | Blk), d, family = poisson)),
               suppressMessages(lme4::glmer(bin ~ A * x + (1 | Blk), d, family = binomial)),
               glmmTMB::glmmTMB(cnt ~ A * x + (1 | Blk), d, family = glmmTMB::nbinom2))
  for (m in mods) {
    r <- suppressMessages(covariatemeans(m, "A", covariate = "x", contr = c(1, 3), covariateV = xs))
    p <- r$pltdf
    rownames(p) <- NULL
    eg <- emmeans::emmeans(m, ~ A | x, at = list(x = xs), data = d)
    el <- as.data.frame(pairs(eg)); el <- el[el$contrast == "a1 - a3", ]
    er <- as.data.frame(pairs(emmeans::regrid(eg, transform = "response"))); er <- er[er$contrast == "a1 - a3", ]
    lab <- class(m)[1]
    expect_equal(p$diff, el$estimate, tolerance = 1e-8, label = lab)
    expect_equal(p$ses, el$SE, tolerance = 1e-6, label = lab)
    expect_equal(p$response.diff, er$estimate, tolerance = 1e-8, label = lab)
    expect_equal(p$response.ses, er$SE, tolerance = 1e-6, label = lab)
    expect_equal(p$response.UL - p$response.diff, qt(0.975, p$Df)*p$response.ses, tolerance = 1e-10)
    # a log link: ratios; a logit link: odds ratios; never the back-transformed limits of the difference
    if (stats::family(m)$link == "log") expect_equal(p$ratio.UL, exp(p$UL)) else expect_equal(p$odds.ratio, exp(p$diff))
    # the plot: the model scale as a ratio or an odds ratio (default); the response scale when asked
    expect_equal(ggplot2::layer_data(r$plt, 3)$y, log10(exp(p$diff)))
    rr <- suppressMessages(covariatemeans(m, "A", covariate = "x", contr = c(1, 3), covariateV = xs, contr_scale = "response"))
    expect_equal(ggplot2::layer_data(rr$plt, 3)$y, p$response.diff)
    expect_equal(rr$pltdf[names(p)], p, ignore_attr = TRUE)
  }
  # trans = I: the model scale only
  pI <- suppressMessages(covariatemeans(mods[[3]], "A", covariate = "x", contr = c(1, 3), covariateV = xs, trans = I))
  expect_null(pI$pltdf$response.diff)
  expect_equal(ggplot2::layer_data(pI$plt, 3)$y, pI$pltdf$diff)
})

test_that("an aov model with Error() and a covariate: covariatemeans() by strata equals the mixed model", {
  skip_on_cran()
  # repeated measures over time, one baseline covariate per subject
  set.seed(4)
  my <- expand.grid(time = factor(1:4), subject = factor(1:24))
  my$group <- factor(rep(c("ctl", "trt"), each = 48))
  my$baseline <- rnorm(24, 10, 2)[my$subject]
  my$y <- 5 + 0.6*my$baseline + (my$group == "trt")*as.numeric(my$time)*0.4 + rnorm(24)[my$subject] + rnorm(96, sd = 0.7)
  ma <- aov(y ~ baseline + group * time + Error(subject / time), data = my)
  ml <- lme4::lmer(y ~ baseline + group * time + (1 | subject), data = my)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  xs <- c(8, 10, 12)
  f <- function(m, ...) suppressWarnings(suppressMessages(covariatemeans(m, "group", covariate = "baseline", covariateV = xs, ...)))$pltdf
  a <- f(ma); l <- f(ml)
  expect_equal(a[c("Mean", "ses", "LL", "UL")], l[c("Mean", "ses", "LL", "UL")], tolerance = 1e-4)
  ac <- f(ma, contr = c(1, 2)); lc <- f(ml, contr = c(1, 2))
  expect_equal(ac[c("diff", "ses", "Df", "LL", "UL")], lc[c("diff", "ses", "Df", "LL", "UL")], tolerance = 1e-4)
  expect_equal(unique(round(ac$Df)), 21)                       # the subject stratum: 24 - 2 groups - 1 covariate
  # predictmeans() at a chosen baseline
  pa <- suppressMessages(predictmeans(ma, "group", covariate = list(baseline = 12), prtnum = FALSE))$means
  pl <- suppressMessages(predictmeans(ml, "group", covariate = list(baseline = 12), prtnum = FALSE))$means
  expect_equal(pa[c("estimate", "std.error")], pl[c("estimate", "std.error")], tolerance = 1e-4)
})

test_that("covariatemeans(contr = ): simultaneous band, joint test and the range where the lines differ", {
  set.seed(1)
  d <- data.frame(A = gl(3, 40, labels = c("a1", "a2", "a3")), x = runif(120, 0, 2))
  d$cnt <- rpois(120, exp(1 + 0.6*d$x*(d$A == "a3") + 0.2*d$x))
  g <- glm(cnt ~ A*x, poisson, d)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  xs <- seq(0, 2, 0.02)
  expect_message(rp <- covariatemeans(g, "A", covariate = "x", contr = c(3, 1), covariateV = xs), "a3 differs from a1")
  rs <- suppressMessages(covariatemeans(g, "A", covariate = "x", contr = c(3, 1), covariateV = xs, band = "simultaneous"))
  # pointwise: z; simultaneous: Scheffe with rank 2 (Working-Hotelling), wider
  expect_equal((rp$pltdf$UL - rp$pltdf$diff)/rp$pltdf$ses, rep(qnorm(0.975), length(xs)), tolerance = 1e-8)
  expect_equal((rs$pltdf$UL - rs$pltdf$diff)/rs$pltdf$ses, rep(sqrt(qchisq(0.95, 2)), length(xs)), tolerance = 1e-8)
  # the significant values are those whose band leaves out zero; fewer with the simultaneous band
  expect_equal(rs$pltdf$significant, rs$pltdf$LL > 0 | rs$pltdf$UL < 0)
  expect_true(all(rp$pltdf$significant[rs$pltdf$significant]))
  expect_lt(sum(rs$pltdf$significant), sum(rp$pltdf$significant))
  expect_equal(rs$significant$from, min(xs[rs$pltdf$significant]))
  expect_equal(rs$significant$direction, ">")
  # the ratio limits leave out 1 exactly where the difference's leave out 0
  expect_equal(rs$pltdf$ratio.LL > 1 | rs$pltdf$ratio.UL < 1, rs$pltdf$significant)
  # indicator columns: each set of limits includes no difference (0, or 1 for a ratio) or not
  expect_equal(rs$pltdf$includes0, rs$pltdf$LL <= 0 & rs$pltdf$UL >= 0)
  expect_equal(rs$pltdf$ratio.includes1, rs$pltdf$ratio.LL <= 1 & rs$pltdf$ratio.UL >= 1)
  expect_equal(rs$pltdf$response.includes0, rs$pltdf$response.LL <= 0 & rs$pltdf$response.UL >= 0)
  expect_equal(rs$pltdf$significant, !rs$pltdf$includes0)
  # the joint test: the Wald test that level and slope differences are both zero
  b <- coef(g); L <- rbind(c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 0, 0, 1))
  w <- drop(t(L %*% b) %*% solve(L %*% vcov(g) %*% t(L)) %*% (L %*% b))
  expect_equal(rp$test[["F"]], w/2, tolerance = 1e-8)
  expect_equal(rp$test[["numDF"]], 2)
  expect_equal(rp$test[["p.value"]], pchisq(w, 2, lower.tail = FALSE), tolerance = 1e-6)
  # the Scheffe band and the joint test agree: no significant value when the test is not significant
  O <- as.data.frame(nlme::Oats)
  fm <- lm(yield ~ Block + nitro*Variety, data = O)
  expect_message(ro <- covariatemeans(fm, "Variety", covariate = "nitro", contr = c(1, 3), band = "simultaneous"),
                 "includes no difference")
  expect_gt(ro$test[["p.value"]], 0.05)
  expect_equal(nrow(ro$significant), 0)
  expect_equal(anova(lm(yield ~ Block + nitro + I(Variety == "Marvellous") + I(Variety == "Marvellous"):nitro, O), fm)$F[2],
               ro$test[["F"]], tolerance = 1e-8)
})

test_that("covariatemeans(bias_adjust = ): the mean of the response, as predictmeans()", {
  O <- as.data.frame(nlme::Oats)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  cm <- function(...) suppressMessages(covariatemeans(...))
  # a log response: the factor exp(sigma^2/2), on the line and on its limits
  ml <- lm(log(yield) ~ Block + nitro*Variety, O)
  a <- cm(ml, "Variety", covariate = "nitro", covariateV = c(0, 0.3))$pltdf
  b <- cm(ml, "Variety", covariate = "nitro", covariateV = c(0, 0.3), bias_adjust = TRUE)$pltdf
  f <- exp(sigma(ml)^2/2)
  expect_equal(b$Mean, a$Mean*f, tolerance = 1e-8)
  expect_equal(b$LL, a$LL*f, tolerance = 1e-8)
  pm <- suppressMessages(predictmeans(ml, "Variety", covariate = c(nitro = 0.3), bias_adjust = TRUE, plot = FALSE,
                                      prtnum = FALSE))
  expect_equal(b$Mean[b$xvar == 0.3], pm$means$response, tolerance = 1e-8)
  # a square root in a mixed model: + the total variance; "taylor" is exact for a square
  mm <- lme4::lmer(sqrt(yield) ~ nitro*Variety + (1 | Block), O)
  a <- cm(mm, "Variety", covariate = "nitro", covariateV = 0.3)$pltdf
  b <- cm(mm, "Variety", covariate = "nitro", covariateV = 0.3, bias_adjust = TRUE)$pltdf
  v <- sum(as.data.frame(lme4::VarCorr(mm))$vcov)
  expect_equal(b$Mean - a$Mean, rep(v, 3), tolerance = 1e-6)
  expect_equal(cm(mm, "Variety", covariate = "nitro", covariateV = 0.3, bias_adjust = "taylor")$pltdf$Mean, b$Mean,
               tolerance = 1e-5)
  # 'sigma' given
  expect_equal(cm(ml, "Variety", covariate = "nitro", covariateV = 0.3, bias_adjust = TRUE, sigma = 0.5)$pltdf$Mean,
               a0 <- cm(ml, "Variety", covariate = "nitro", covariateV = 0.3)$pltdf$Mean*exp(0.125), tolerance = 1e-8)
  # a Poisson GLMM: the population mean, exp(random-intercept variance/2) above the typical cluster;
  # the response-scale difference of contr is that of the adjusted lines
  set.seed(2)
  d <- data.frame(A = gl(2, 60), x = runif(120, 0, 2), g = gl(12, 10))
  d$y <- rpois(120, exp(1 + 0.4*d$x*(d$A == 2) + rnorm(12, 0, 0.5)[d$g]))
  gm <- lme4::glmer(y ~ A*x + (1 | g), d, family = poisson)
  a <- cm(gm, "A", covariate = "x", covariateV = c(0.5, 1))$pltdf
  b <- cm(gm, "A", covariate = "x", covariateV = c(0.5, 1), bias_adjust = TRUE)$pltdf
  expect_equal(b$Mean/a$Mean, rep(exp(as.data.frame(lme4::VarCorr(gm))$vcov/2), 4), tolerance = 1e-6)
  r <- cm(gm, "A", covariate = "x", covariateV = c(0.5, 1), contr = c(2, 1), bias_adjust = TRUE)$pltdf
  expect_equal(r$response.diff, b$Mean[3:4] - b$Mean[1:2], tolerance = 1e-8)
  r0 <- cm(gm, "A", covariate = "x", covariateV = c(0.5, 1), contr = c(2, 1))$pltdf
  expect_equal(r$ratio, r0$ratio)                       # the ratio and the test do not change
  expect_equal(r$significant, r0$significant)
  # nothing to adjust: a warning; a wrong value: an error
  expect_warning(cm(lm(yield ~ nitro*Variety, O), "Variety", covariate = "nitro", bias_adjust = TRUE), "ignored")
  expect_error(cm(ml, "Variety", covariate = "nitro", bias_adjust = "yes"), "TRUE, FALSE")
  # aov with Error(): the stratum variance components
  fa <- aov(log(yield) ~ nitro*Variety + Error(Block/Variety), O)
  ba <- cm(fa, "Variety", covariate = "nitro", covariateV = 0.3, bias_adjust = TRUE)$pltdf$Mean
  bl <- cm(lme4::lmer(log(yield) ~ nitro*Variety + (1 | Block/Variety), O), "Variety", covariate = "nitro",
           covariateV = 0.3, bias_adjust = TRUE)$pltdf$Mean
  expect_equal(ba, bl, tolerance = 1e-3)
})
