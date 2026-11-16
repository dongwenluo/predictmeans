test_that("factors not in modelterm (incl. 'by' factors) are averaged in both parts", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date*district,
                                 smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
                                 data = W, prt = FALSE))
  overall <- semipred(sm, covariate = "construction.date", prt = FALSE)$pred_df
  by_district <- semipred(sm, "district", "construction.date", prt = FALSE)$pred_df
  avg <- aggregate(by_district$areaPerMzloty, list(cd = by_district$construction.date), mean)
  expect_equal(overall$construction.date, avg$cd)
  expect_equal(overall$areaPerMzloty, avg$x, tolerance = 1e-10)
})

test_that("a matrix covariateV gives one prediction per point for a 2-D smooth", {
  set.seed(3)
  s2 <- data.frame(u = runif(200), v = runif(200))
  s2$w <- sin(3*s2$u) + cos(3*s2$v) + rnorm(200, sd = 0.2)
  sm <- semireg(w ~ u + v, smoothZ = list(g = smZ(cbind(u, v), k = 6, type = "Ztps")),
                data = s2, prt = FALSE)
  pts <- cbind(u = c(0.2, 0.5, 0.8), v = c(0.3, 0.5, 0.7))
  p <- semipred(sm, covariate = c("u", "v"), covariateV = pts, prt = FALSE)$pred_df
  expect_equal(nrow(p), 3)
  expect_equal(as.matrix(p[c("u", "v")]), pts, ignore_attr = TRUE, tolerance = 1e-12)
  expect_true(all(p$LL < p$w & p$w < p$UL))
})

test_that("a numeric variable that is not in the prediction (here the grouping id) doesn't stop semipred()", {
  set.seed(7)
  d <- data.frame(u = runif(150, 2, 8), id = rep(1:15, each = 10))
  d$y <- sin(d$u) + rnorm(15, sd = 0.5)[d$id] + rnorm(150, sd = 0.3)
  centred <- suppressMessages(semireg(y ~ u + (1 | id), smoothZ = list(g = smZ(u)), data = d, prt = FALSE))
  plain <- suppressMessages(semireg(y ~ u + (1 | id), smoothZ = list(g = smZ(u)), data = d,
                                    ncenter = FALSE, prt = FALSE))
  expect_equal(semipred(centred, covariate = "u", prt = FALSE)$pred_df,
               semipred(plain, covariate = "u", prt = FALSE)$pred_df, tolerance = 1e-6)
})

test_that("anova() of a semireg fit tests each smooth by refitting without it (copper)", {
  skip_on_cran()
  data(copper)
  full <- suppressWarnings(semireg(log(grade) ~ xcoord + ycoord + zcoord,
    smoothZ = list(grp_z = smZ(zcoord), grp_geo = smZ(cbind(xcoord, ycoord), type = "Ztps")),
    data = copper, prt = FALSE, predict_info = FALSE))
  expect_s4_class(full$semer, "semerMod")
  a <- anova(full)
  expect_equal(a$smooth$term, c("grp_z", "grp_geo"))
  # dropping grp_geo must equal fitting with the grp_z smooth only (REML, same fixed effects)
  only_z <- suppressWarnings(semireg(log(grade) ~ xcoord + ycoord + zcoord,
    smoothZ = list(grp_z = smZ(zcoord)), data = copper, prt = FALSE, predict_info = FALSE))
  expect_equal(a$smooth$logLik.without[a$smooth$term == "grp_geo"], as.numeric(logLik(only_z$semer)),
               tolerance = 1e-6)
  r <- a$smooth[a$smooth$term == "grp_geo", ]
  expect_equal(r$LRT, 2*(as.numeric(logLik(full$semer)) - as.numeric(logLik(only_z$semer))), tolerance = 1e-6)
  expect_equal(r$p.value, 0.5*pchisq(r$LRT, 1, lower.tail = FALSE))
  expect_true(all(c("xcoord", "ycoord", "zcoord") %in% rownames(a$fixed)))
})

test_that("predictmeans() on a semireg fit: means with the smooth, as semipred(); pairs as semipred(contr =)", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date*district,
                                 smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
                                 data = W, prt = FALSE))
  r <- predictmeans(sm, "district", covariate = c(construction.date = 1980), adj = "tukey")
  expect_s3_class(r, "pm_semimeans")
  sp <- semipred(sm, "district", "construction.date", covariateV = 1980, prt = FALSE)$pred_df
  expect_equal(r$means$estimate, sp$areaPerMzloty, tolerance = 1e-10)
  expect_equal(r$means$std.error, sp$SE, tolerance = 1e-10)
  sc <- semipred(sm, "district", "construction.date", contr = c(2, 1), covariateV = 1980, prt = FALSE)$pred_df
  expect_equal(r$pairs$estimate[1], -sc$areaPerMzloty, tolerance = 1e-10)
  expect_equal(r$pairs$std.error[1], sc$SE_contr, tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(vcov(r)))), r$means$std.error)
  # default: covariates at their means
  r0 <- predictmeans(sm, "district")
  sp0 <- semipred(sm, "district", "construction.date", covariateV = mean(W$construction.date), prt = FALSE)$pred_df
  expect_equal(r0$means$estimate, sp0$areaPerMzloty, tolerance = 1e-8)
})

test_that("semipred(): two-sided limits; simultaneous bands are wider, repeat, and keep the user's RNG", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date,
                                 smoothZ = list(grp = smZ(construction.date, k = 25)), data = WarsawApts, prt = FALSE))
  p <- semipred(sm, covariate = "construction.date", prt = FALSE)$pred_df
  expect_equal((p$UL - p$areaPerMzloty)/p$SE, rep(qnorm(0.975), nrow(p)), tolerance = 1e-10)
  set.seed(9); u1 <- runif(1); set.seed(9)
  s <- semipred(sm, covariate = "construction.date", prt = FALSE, simultaneous = TRUE)$pred_df
  expect_identical(runif(1), u1)
  crit <- unique(round((s$UL - s$areaPerMzloty)/s$SE, 8))
  expect_length(crit, 1)
  expect_gt(crit, qnorm(0.975))
  expect_lt(crit, qnorm(1 - 0.05/(2*nrow(s))))      # below the Bonferroni bound
  s2 <- semipred(sm, covariate = "construction.date", prt = FALSE, simultaneous = TRUE)$pred_df
  expect_identical(s$UL, s2$UL)
})

test_that("semipred(contr = ): the difference of two curves as covariatemeans(contr = ); GLMMs on both scales", {
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ construction.date*district,
                                 smoothZ = list(grp = smZ(construction.date, k = 15, by = district, group = TRUE)),
                                 data = W, prt = FALSE))
  # a linear model: the 1.x columns and the indicators; level names in contr
  r <- semipred(sm, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"), prt = FALSE)$pred_df
  lv <- levels(W$district)
  r2 <- semipred(sm, "district", "construction.date", contr = match(c("Srodmiescie", "Mokotow"), lv), prt = FALSE)$pred_df
  expect_equal(r, r2)
  expect_equal(r$UL - r$areaPerMzloty, qnorm(0.975)*r$SE_contr, tolerance = 1e-10)
  expect_equal(r$includes0, r$LL <= 0 & r$UL >= 0)
  expect_equal(r$significant, !r$includes0)
  m <- semipred(sm, "district", "construction.date", prt = FALSE)$pred_df
  expect_equal(r$areaPerMzloty, m$areaPerMzloty[m$district == "Srodmiescie"] - m$areaPerMzloty[m$district == "Mokotow"],
               tolerance = 1e-8)
  rs <- semipred(sm, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"), prt = FALSE,
                 simultaneous = TRUE)$pred_df
  expect_true(all(rs$UL - rs$LL > r$UL - r$LL))
  expect_lte(sum(rs$significant), sum(r$significant))
  sg <- attr(rs, "significant")
  expect_equal(nrow(sg), sum(rle(rs$significant)$values))
  expect_message(semipred(sm, "district", "construction.date", contr = c(2, 1), prt = TRUE), "band, level 0.05")
  # a Poisson model: link-scale difference, ratio, and the response-scale difference by the delta method
  set.seed(3)
  n <- 300
  d <- data.frame(x = runif(n, 0, 3), g = gl(2, n/2, labels = c("a", "b")))
  d$y <- rpois(n, exp(0.5 + sin(d$x) + 0.5*(d$g == "b")*d$x/3))
  sp <- suppressWarnings(semireg(y ~ x*g, smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d,
                                 family = poisson, prt = FALSE))
  rg <- semipred(sp, "g", "x", contr = c("b", "a"), prt = FALSE)
  p <- rg$pred_df
  m <- semipred(sp, "g", "x", prt = FALSE)$pred_df
  ma <- m$y[m$g == "a"]; mb <- m$y[m$g == "b"]
  expect_equal(p$response.diff, mb - ma, tolerance = 1e-8)
  expect_equal(p$ratio, mb/ma, tolerance = 1e-8)
  expect_equal(p$diff, log(mb/ma), tolerance = 1e-8)
  expect_equal(p$ratio.UL, exp(p$diff + qnorm(0.975)*p$ses), tolerance = 1e-8)
  expect_equal(p$ratio.includes1, p$includes0)
  expect_equal(p$SE_contr, p$ses)
  # the delta-method SE of a difference of two means: between |SE_b - SE_a| and SE_b + SE_a
  sa <- m$SE_mu[m$g == "a"]; sb <- m$SE_mu[m$g == "b"]
  expect_true(all(p$response.ses <= sa + sb + 1e-8 & p$response.ses >= abs(sb - sa) - 1e-8))
  expect_match(rg$plt$labels$y, "Ratio")
  rr <- semipred(sp, "g", "x", contr = c("b", "a"), contr_scale = "response", prt = FALSE)
  expect_match(rr$plt$labels$y, "response scale")
  expect_equal(rr$pred_df, p, ignore_attr = TRUE)
})

test_that("the covariance behind semipred()'s bands is that of the fitted model, for every family", {
  # sigma^2 CovMat is the covariance of (fixed, random - their true values)
  # given the variance parameters: its fixed-effect block must equal the
  # model's own vcov() (lme4's, from the penalised least squares fit)
  set.seed(3)
  n <- 400
  d <- data.frame(x = runif(n, 0, 3), g = gl(2, n/2, labels = c("a", "b")), id = gl(20, 1, n))
  eta <- 0.3*sin(2*d$x) + 0.5*(d$g == "b")*d$x/3 + rnorm(20, 0, 0.3)[d$id]
  d$yn <- eta + rnorm(n, 0, 0.5)
  d$yp <- rpois(n, exp(0.5 + eta))
  d$yb <- rbinom(n, 1, plogis(-1.5 + 2*eta))
  q <- function(e) suppressWarnings(suppressMessages(e))
  fits <- list(
    gaussian = q(semireg(yn ~ x*g + (1 | id), smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d,
                         prt = FALSE)),
    poisson = q(semireg(yp ~ x*g + (1 | id), smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d,
                        family = poisson, prt = FALSE)),
    # a 0/1 response: before 2.0 the weights were the constant 3/pi^2 and the SEs 15-25% too small
    binary = q(semireg(yb ~ x*g + (1 | id), smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d,
                       family = binomial, prt = FALSE)))
  for (nm in names(fits)) {
    sm <- fits[[nm]]
    p <- length(lme4::fixef(sm$semer))
    V <- as.matrix(sigma(sm$semer)^2 * sm$CovMat[1:p, 1:p])
    Vm <- if (lme4::isGLMM(sm$semer)) q(as.matrix(vcov(sm$semer, use.hessian = FALSE))) else as.matrix(vcov(sm$semer))
    expect_equal(unname(V), unname(Vm), tolerance = 1e-4, label = nm)
    # the SE of the difference of two curves: (D_a - D_b) V (D_a - D_b)' with that covariance,
    # not a combination of the two curves' own SEs
    xv <- c(0.5, 1.5, 2.5)                       # semireg() centres and scales the covariate
    if (!is.null(attr(sm$data, "x_mean"))) xv <- xv - attr(sm$data, "x_mean")
    if (!is.null(attr(sm$data, "x_sd"))) xv <- xv/attr(sm$data, "x_sd")
    des <- predictmeans:::pm_semi_design(sm, "g", "x", xv, NULL, sm$data)
    D <- as.matrix(des$D_mx)
    lev <- des$pred_df$g
    L <- D[lev == "b", , drop = FALSE] - D[lev == "a", , drop = FALSE]
    r <- q(semipred(sm, "g", "x", contr = c("b", "a"), covariateV = c(0.5, 1.5, 2.5), prt = FALSE))$pred_df
    se <- sqrt(diag(L %*% as.matrix(sigma(sm$semer)^2 * des$sm_Cov) %*% t(L)))
    expect_equal(r$SE_contr, unname(se), tolerance = 1e-8, label = nm)
    m <- q(semipred(sm, "g", "x", covariateV = c(0.5, 1.5, 2.5), prt = FALSE))$pred_df
    sa <- sqrt(diag(D[lev == "a", ] %*% as.matrix(sigma(sm$semer)^2 * des$sm_Cov) %*% t(D[lev == "a", ])))
    sb <- sqrt(diag(D[lev == "b", ] %*% as.matrix(sigma(sm$semer)^2 * des$sm_Cov) %*% t(D[lev == "b", ])))
    expect_true(all(se <= sa + sb & se >= abs(sa - sb)), label = nm)
  }
})

test_that("semipred(): generalized models back-transform the limits, as predictmeans(); band = as covariatemeans()", {
  set.seed(3)
  n <- 300
  d <- data.frame(x = runif(n, 0, 3), g = gl(2, n/2, labels = c("a", "b")))
  d$y <- rpois(n, exp(0.5 + sin(d$x) + 0.5*(d$g == "b")*d$x/3))
  sp <- suppressWarnings(semireg(y ~ x*g, smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d,
                                 family = poisson, prt = FALSE))
  m <- semipred(sp, "g", "x", covariateV = 1.5, prt = FALSE)$pred_df
  pm <- suppressMessages(predictmeans(sp, "g", covariate = c(x = 1.5)))$means
  expect_equal(m$y, pm$response, tolerance = 1e-8)
  expect_equal(m$SE_mu, pm$response.std.error, tolerance = 1e-8)
  expect_equal(m$LL, pm$response.conf.low, tolerance = 1e-6)
  expect_equal(m$UL, pm$response.conf.high, tolerance = 1e-6)
  # limits are exp of the link-scale limits: positive, and not symmetric about the mean
  expect_true(all(m$LL > 0) && all(m$UL - m$y > m$y - m$LL))
  # band = "simultaneous" is simultaneous = TRUE
  expect_equal(semipred(sp, "g", "x", prt = FALSE, band = "simultaneous")$pred_df,
               semipred(sp, "g", "x", prt = FALSE, simultaneous = TRUE)$pred_df)
})

test_that("semipred(contr = ) with two covariates: the difference of two surfaces from the contrast", {
  set.seed(3)
  s2 <- data.frame(u = runif(200), v = runif(200), g = gl(2, 100, labels = c("a", "b")))
  s2$w <- sin(3*s2$u) + cos(3*s2$v) + 0.4*(s2$g == "b") + rnorm(200, sd = 0.2)
  sm <- semireg(w ~ u + v + g, smoothZ = list(s = smZ(cbind(u, v), k = 6, type = "Ztps")), data = s2, prt = FALSE)
  pts <- cbind(u = c(0.2, 0.5, 0.8), v = c(0.3, 0.5, 0.7))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  r <- semipred(sm, "g", c("u", "v"), contr = c("b", "a"), covariateV = pts, prt = FALSE)$pred_df
  # the surfaces differ by the g effect only: the same difference and SE at every point
  b <- lme4::fixef(sm$semer)[["gb"]]
  se <- sqrt(vcov(sm$semer)["gb", "gb"])
  expect_equal(r$w, rep(b, 3), tolerance = 1e-6)
  expect_equal(r$SE_contr, rep(se, 3), tolerance = 1e-4)
  expect_equal(r$UL - r$w, qnorm(0.975)*r$SE_contr, tolerance = 1e-8)
  expect_equal(r$includes0, r$LL <= 0 & r$UL >= 0)
  expect_equal(r$significant, !r$includes0)
})
