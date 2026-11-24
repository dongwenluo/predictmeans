test_that("factors not in modelterm (incl. 'by' factors) are averaged in both parts", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ s(construction.date, by = district, bs = "os", k = 15),
                                 data = W))
  overall <- semipred(sm, covariate = "construction.date")$pred_df
  by_district <- semipred(sm, "district", "construction.date")$pred_df
  avg <- aggregate(by_district$areaPerMzloty, list(cd = by_district$construction.date), mean)
  expect_equal(overall$construction.date, avg$cd)
  expect_equal(overall$areaPerMzloty, avg$x, tolerance = 1e-10)
})

test_that("a matrix covariateV gives one prediction per point for a 2-D smooth", {
  set.seed(3)
  s2 <- data.frame(u = runif(200), v = runif(200))
  s2$w <- sin(3*s2$u) + cos(3*s2$v) + rnorm(200, sd = 0.2)
  sm <- semireg(w ~ s(u, v, bs = "ztps", k = 6), data = s2)
  pts <- cbind(u = c(0.2, 0.5, 0.8), v = c(0.3, 0.5, 0.7))
  p <- semipred(sm, covariate = c("u", "v"), at = pts)$pred_df
  expect_equal(nrow(p), 3)
  expect_equal(as.matrix(p[c("u", "v")]), pts, ignore_attr = TRUE, tolerance = 1e-12)
  expect_true(all(p$LL < p$w & p$w < p$UL))
})

test_that("a numeric variable that is not in the prediction (here the grouping id) doesn't stop semipred()", {
  set.seed(7)
  d <- data.frame(u = runif(150, 2, 8), id = rep(1:15, each = 10))
  d$y <- sin(d$u) + rnorm(15, sd = 0.5)[d$id] + rnorm(150, sd = 0.3)
  fit <- suppressMessages(semireg(y ~ s(u, bs = "os") + (1 | id), data = d))
  p <- semipred(fit, covariate = "u")$pred_df
  expect_true(all(is.finite(p$y)) && all(p$LL < p$y & p$y < p$UL))
})

test_that("anova() of a semireg fit tests each smooth by refitting without it (copper)", {
  skip_on_cran()
  data(copper)
  full <- suppressWarnings(semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6),
    data = copper))
  expect_s4_class(full$semer, "semerMod")
  a <- anova(full)
  expect_equal(a$smooth$term, c("s_zcoord", "s_xcoord_ycoord"))
  # dropping grp_geo must equal fitting with the grp_z smooth only (REML, same fixed effects)
  only_z <- suppressWarnings(semireg(log(grade) ~ s(zcoord, bs = "os") + xcoord + ycoord, data = copper))
  expect_equal(a$smooth$logLik.without[a$smooth$term == "s_xcoord_ycoord"], as.numeric(logLik(only_z$semer)),
               tolerance = 1e-6)
  r <- a$smooth[a$smooth$term == "s_xcoord_ycoord", ]
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
  sm <- suppressWarnings(semireg(areaPerMzloty ~ s(construction.date, by = district, bs = "os", k = 15),
                                 data = W))
  r <- predictmeans(sm, "district", covariate = c(construction.date = 1980), adj = "tukey")
  expect_s3_class(r, "pm_semimeans")
  sp <- semipred(sm, "district", "construction.date", at = 1980)$pred_df
  expect_equal(r$means$estimate, sp$areaPerMzloty, tolerance = 1e-10)
  expect_equal(r$means$std.error, sp$SE, tolerance = 1e-10)
  sc <- semipred(sm, "district", "construction.date", contr = c(2, 1), at = 1980)$pred_df
  expect_equal(r$pairs$estimate[1], -sc$areaPerMzloty, tolerance = 1e-10)
  expect_equal(r$pairs$std.error[1], sc$SE_contr, tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(vcov(r)))), r$means$std.error)
  # default: covariates at their means
  r0 <- predictmeans(sm, "district")
  sp0 <- semipred(sm, "district", "construction.date", at = mean(W$construction.date))$pred_df
  expect_equal(r0$means$estimate, sp0$areaPerMzloty, tolerance = 1e-8)
})

test_that("semipred(): two-sided t limits with Satterthwaite df; simultaneous bands are wider, repeat, and keep the user's RNG", {
  skip_on_cran()
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  sm <- suppressWarnings(semireg(areaPerMzloty ~ s(construction.date, bs = "os", k = 25), data = WarsawApts))
  p <- semipred(sm, covariate = "construction.date")$pred_df
  expect_true(all(is.finite(p$df) & p$df > 1))
  expect_equal((p$UL - p$areaPerMzloty)/p$SE, qt(0.975, p$df), tolerance = 1e-10)
  set.seed(9); u1 <- runif(1); set.seed(9)
  s <- semipred(sm, covariate = "construction.date", band = "simultaneous")$pred_df
  expect_identical(runif(1), u1)
  expect_equal(s$df, p$df)
  crit <- (s$UL - s$areaPerMzloty)/s$SE
  # one probability for all points, each point at the t quantile of its df
  expect_lt(diff(range(pt(crit, s$df))), 1e-10)
  expect_true(all(crit > qt(0.975, s$df)))
  expect_true(all(crit < qt(1 - 0.05/(2*nrow(s)), s$df)))      # below the Bonferroni bound
  s2 <- semipred(sm, covariate = "construction.date", band = "simultaneous")$pred_df
  expect_identical(s$UL, s2$UL)
})

test_that("semipred(contr = ): the difference of two curves as covariatemeans(contr = ); GLMMs on both scales", {
  skip_if_not_installed("HRW")
  data("WarsawApts", package = "HRW", envir = environment())
  W <- WarsawApts
  W$district <- factor(W$district)
  sm <- suppressWarnings(semireg(areaPerMzloty ~ s(construction.date, by = district, bs = "os", k = 15),
                                 data = W))
  # a linear model: the 1.x columns and the indicators; level names in contr
  r <- semipred(sm, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"))$pred_df
  lv <- levels(W$district)
  r2 <- semipred(sm, "district", "construction.date", contr = match(c("Srodmiescie", "Mokotow"), lv))$pred_df
  expect_equal(r, r2)
  expect_equal(r$UL - r$areaPerMzloty, qt(0.975, r$df)*r$SE_contr, tolerance = 1e-10)
  expect_equal(r$includes0, r$LL <= 0 & r$UL >= 0)
  expect_equal(r$significant, !r$includes0)
  m <- semipred(sm, "district", "construction.date")$pred_df
  expect_equal(r$areaPerMzloty, m$areaPerMzloty[m$district == "Srodmiescie"] - m$areaPerMzloty[m$district == "Mokotow"],
               tolerance = 1e-8)
  rs <- semipred(sm, "district", "construction.date", contr = c("Srodmiescie", "Mokotow"),
                 band = "simultaneous")$pred_df
  expect_true(all(rs$UL - rs$LL > r$UL - r$LL))
  expect_lte(sum(rs$significant), sum(r$significant))
  sg <- attr(rs, "significant")
  expect_equal(nrow(sg), sum(rle(rs$significant)$values))
  expect_output(print(semipred(sm, "district", "construction.date", contr = c(2, 1))), "pointwise limits")
  # a Poisson model: link-scale difference, ratio, and the response-scale difference by the delta method
  set.seed(3)
  n <- 300
  d <- data.frame(x = runif(n, 0, 3), g = gl(2, n/2, labels = c("a", "b")))
  d$y <- rpois(n, exp(0.5 + sin(d$x) + 0.5*(d$g == "b")*d$x/3))
  sp <- suppressWarnings(semireg(y ~ s(x, by = g, bs = "os", k = 8), data = d,
                                 family = poisson))
  rg <- semipred(sp, "g", "x", contr = c("b", "a"))
  p <- rg$pred_df
  m <- semipred(sp, "g", "x")$pred_df
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
  rr <- semipred(sp, "g", "x", contr = c("b", "a"), contr_scale = "response")
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
    gaussian = q(semireg(yn ~ s(x, by = g, bs = "os", k = 8) + (1 | id), data = d)),
    poisson = q(semireg(yp ~ s(x, by = g, bs = "os", k = 8) + (1 | id), data = d,
                        family = poisson)),
    # a 0/1 response: before 2.0 the weights were the constant 3/pi^2 and the SEs 15-25% too small
    binary = q(semireg(yb ~ s(x, by = g, bs = "os", k = 8) + (1 | id), data = d,
                       family = binomial)))
  for (nm in names(fits)) {
    sm <- fits[[nm]]
    p <- length(lme4::fixef(sm$semer))
    V <- as.matrix(sigma(sm$semer)^2 * sm$CovMat[1:p, 1:p])
    Vm <- if (lme4::isGLMM(sm$semer)) q(as.matrix(vcov(sm$semer, use.hessian = FALSE))) else as.matrix(vcov(sm$semer))
    expect_equal(unname(V), unname(Vm), tolerance = 1e-4, label = nm)
    # the SE of the difference of two curves: (D_a - D_b) V (D_a - D_b)' with that covariance,
    # not a combination of the two curves' own SEs
    xv <- c(0.5, 1.5, 2.5)
    des <- predictmeans:::pm_semi_design(sm, "g", "x", xv, NULL, sm$data)
    D <- as.matrix(des$D_mx)
    lev <- des$pred_df$g
    L <- D[lev == "b", , drop = FALSE] - D[lev == "a", , drop = FALSE]
    r <- q(semipred(sm, "g", "x", contr = c("b", "a"), at = c(0.5, 1.5, 2.5)))$pred_df
    se <- sqrt(diag(L %*% as.matrix(sigma(sm$semer)^2 * des$sm_Cov) %*% t(L)))
    expect_equal(r$SE_contr, unname(se), tolerance = 1e-8, label = nm)
    m <- q(semipred(sm, "g", "x", at = c(0.5, 1.5, 2.5)))$pred_df
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
  sp <- suppressWarnings(semireg(y ~ s(x, by = g, bs = "os", k = 8), data = d,
                                 family = poisson))
  m <- semipred(sp, "g", "x", at = 1.5)$pred_df
  pm <- suppressMessages(predictmeans(sp, "g", covariate = c(x = 1.5)))$means
  expect_equal(m$y, pm$response, tolerance = 1e-8)
  expect_equal(m$SE_mu, pm$response.std.error, tolerance = 1e-8)
  expect_equal(m$LL, pm$response.conf.low, tolerance = 1e-6)
  expect_equal(m$UL, pm$response.conf.high, tolerance = 1e-6)
  # limits are exp of the link-scale limits: positive, and not symmetric about the mean
  expect_true(all(m$LL > 0) && all(m$UL - m$y > m$y - m$LL))
  # the 1.x names stop with their 2.0 spelling
  expect_error(semipred(sp, "g", "x", simultaneous = TRUE), "use band = \"simultaneous\"", fixed = TRUE)
  expect_error(semipred(sp, "g", "x", covariateV = 1), "use at = 1", fixed = TRUE)
  expect_error(semipred(sp, "g", "x", prt = TRUE), "'prt' was removed", fixed = TRUE)
  expect_error(semipred(sp, "g", "x", trellis = FALSE), "use plot_args = list(trellis", fixed = TRUE)
  expect_error(semipred(semireg = sp, "g", "x"), "use model = sp", fixed = TRUE)
  expect_error(semipred(sp, "g", "x", level = 0.1), "use conf_level = 0.9", fixed = TRUE)
})

test_that("semipred(contr = ) with two covariates: the difference of two surfaces from the contrast", {
  set.seed(3)
  s2 <- data.frame(u = runif(200), v = runif(200), g = gl(2, 100, labels = c("a", "b")))
  s2$w <- sin(3*s2$u) + cos(3*s2$v) + 0.4*(s2$g == "b") + rnorm(200, sd = 0.2)
  sm <- semireg(w ~ s(u, v, bs = "ztps", k = 6) + g, data = s2)
  pts <- cbind(u = c(0.2, 0.5, 0.8), v = c(0.3, 0.5, 0.7))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  r <- semipred(sm, "g", c("u", "v"), contr = c("b", "a"), at = pts)$pred_df
  # the surfaces differ by the g effect only: the same difference and SE at every point
  b <- lme4::fixef(sm$semer)[["gb"]]
  se <- sqrt(vcov(sm$semer)["gb", "gb"])
  expect_equal(r$w, rep(b, 3), tolerance = 1e-6)
  expect_equal(r$SE_contr, rep(se, 3), tolerance = 1e-4)
  expect_equal(r$UL - r$w, qt(0.975, r$df)*r$SE_contr, tolerance = 1e-8)
  # a difference of fixed effects only, the smoothing variance held at its
  # estimate (the default): the residual is the only variance left, df n - p;
  # with nothing held, lmer's Satterthwaite df of that coefficient
  expect_equal(r$df, rep(200 - 4, 3), tolerance = 1e-3)
  L <- matrix(as.numeric(names(lme4::fixef(sm$semer)) == "gb"), 1)
  expect_equal(pm_sr_satt_ddf(sm, L, 1:4, hold = character(0)), unname(pm_satt_ddf(L, pm_satt_parts(sm$semer))),
               tolerance = 1e-4)
  expect_equal(r$includes0, r$LL <= 0 & r$UL >= 0)
  expect_equal(r$significant, !r$includes0)
})

test_that("semipred(deriv = ): slopes and curvatures equal the exact derivatives of the O'Sullivan basis", {
  set.seed(3)
  d <- data.frame(x = runif(150, 0, 10))
  d$y <- sin(d$x/1.5) + rnorm(150, 0, 0.3)
  kn <- stats::quantile(unique(d$x), seq(0, 1, length = 12)[-c(1, 12)])
  rg <- c(1.01*min(d$x) - 0.01*max(d$x), 1.01*max(d$x) - 0.01*min(d$x))
  g <- suppressWarnings(semireg(y ~ s(x, bs = "os", xt = list(knots = kn, range = rg)), data = d))
  xs <- seq(0.5, 9.5, by = 0.5)
  # (the design of a formula smooth is the basis times a fixed matrix Tm: the
  # reparameterisation of the formula interface)
  Z0 <- ZOSull(d$x, intKnots = kn, range.x = rg)
  Tm <- solve(crossprod(Z0), crossprod(Z0, as.matrix(lme4::getME(g$semer, "Z"))))
  expect_equal(Z0 %*% Tm, as.matrix(lme4::getME(g$semer, "Z")), ignore_attr = TRUE, tolerance = 1e-8)
  for (dv in 1:2) {
    p <- semipred(g, covariate = "x", at = xs, deriv = dv)$pred_df
    Zd <- ZOSull(xs, intKnots = kn, range.x = rg, drv = dv) %*% Tm
    L <- cbind(0, if (dv == 1) 1 else 0, Zd)
    expect_equal(p$y, as.numeric(L %*% c(lme4::fixef(g$semer), g$u_lst[[1]])), tolerance = 1e-6)
    idx <- seq_len(ncol(L))
    expect_equal(p$SE, sigma(g$semer)*sqrt(rowSums((L %*% as.matrix(g$CovMat[idx, idx])) * L)), tolerance = 1e-6)
    expect_equal(p$UL - p$y, qt(0.975, p$df)*p$SE, tolerance = 1e-8)
    expect_equal(p$significant, !(p$LL <= 0 & p$UL >= 0))
  }
  # the default knots of bs = "os" (k = 10) are those given above: the same slopes
  f <- suppressWarnings(semireg(y ~ s(x, bs = "os", k = 10), data = d))
  pf <- semipred(f, covariate = "x", at = xs, deriv = 1)$pred_df
  pg <- semipred(g, covariate = "x", at = xs, deriv = 1)$pred_df
  expect_equal(pf$y, pg$y, tolerance = 1e-5)
  expect_equal(pf$SE, pg$SE, tolerance = 1e-5)
  # the ranges where the slope differs from 0 (true slope cos(x/1.5)/1.5: > 0 below 2.36, < 0 to 7.07)
  sig <- attr(semipred(f, covariate = "x", deriv = 1)$pred_df, "significant")
  expect_equal(sig$direction, c("> 0", "< 0", "> 0"))
  expect_true(sig$to[1] < 2.36 && sig$from[2] > 2.36 && sig$to[2] < 7.07 && sig$from[3] > 7.07)
  expect_output(print(semipred(f, covariate = "x", deriv = 1)), "Slope differs from 0 \\(pointwise limits\\) for x")
  expect_error(semipred(f, covariate = "x", deriv = 3), "0, 1 or 2")
  expect_error(semipred(f, covariate = "x", deriv = 1, trans = exp), "model scale")
})

test_that("semipred(deriv = , contr = ): the difference of two slopes; generalized fits on the link scale", {
  set.seed(11)
  sp <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
  sp$y <- 2 + sin(sp$t/2) + 0.08*(sp$A == "a2")*sp$t + rnorm(6, 0, 0.5)[sp$Block] + rnorm(nrow(sp), 0, 0.3)
  f <- suppressWarnings(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block), data = sp))
  a <- semipred(f, "A", "t", deriv = 1)$pred_df
  r <- semipred(f, "A", "t", deriv = 1, contr = c("a2", "a1"))$pred_df
  expect_equal(r$y, a$y[a$A == "a2"] - a$y[a$A == "a1"], tolerance = 1e-10)
  expect_equal(r$UL - r$y, qt(0.975, r$df)*r$SE_contr, tolerance = 1e-8)
  expect_output(print(semipred(f, "A", "t", deriv = 1, contr = c("a2", "a1"))), "Difference a2 - a1 of the slopes")
  # a GLMM: slopes of the linear predictor, no back-transformation
  sp$cnt <- rpois(nrow(sp), exp(1 + sin(sp$t/2)/2))
  g <- suppressWarnings(semireg(cnt ~ s(t, by = A, bs = "os", k = 6) + (1 | Block), family = poisson, data = sp))
  gp <- semipred(g, "A", "t", deriv = 1)$pred_df
  expect_null(gp$SE_mu)
  eta <- function(v) {
    des <- pm_semi_design(g, "A", "t", v, NULL, g$data)
    as.numeric(as.matrix(des$D_mx) %*% des$beta_u)
  }
  tt <- unique(gp$t)
  expect_equal(gp$cnt, (eta(tt + 1e-4) - eta(tt - 1e-4))/2e-4, tolerance = 1e-5)
})

test_that("plot_args sets the colours: band under its line, reference line, bars where a difference is clear", {
  data(Dialyzer, package = "nlme", envir = environment())
  fm <- suppressMessages(semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
                                 data = Dialyzer))
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
  p1 <- semipred(fm, "QB", "pressure", plot_args = list(palette = c("black", "grey50")))$plt
  g <- geoms(p1)
  expect_lt(match("GeomRibbon", g), match("GeomLine", g))         # the band first, the line on top
  b <- ggplot2::ggplot_build(p1)
  expect_setequal(unique(b$data[[match("GeomLine", g)]]$colour), c("black", "grey50"))
  pc <- semipred(fm, "QB", "pressure", contr = c("300", "200"),
                 plot_args = list(diff_line = "navy", ref = "orange"))$plt
  gc <- geoms(pc)
  expect_identical(pc$layers[[match("GeomLine", gc)]]$aes_params$colour, "navy")
  expect_identical(pc$layers[[match("GeomHline", gc)]]$aes_params$colour, "orange")
  expect_true("GeomRug" %in% gc)                               # 300 and 200 differ at high pressure
  expect_error(semipred(fm, "QB", "pressure", plot_args = list(colour = "red")), "can hold")
})

test_that("generalized fits: SE on the link scale next to SE_mu, the limits h(eta -+ crit SE)", {
  set.seed(5)
  d <- data.frame(x = runif(200, 0, 10))
  d$y <- rbinom(200, 1, plogis(-1 + sin(d$x / 2)))
  f <- suppressMessages(semireg(y ~ s(x, bs = "os", k = 8), family = binomial, data = d))
  p <- semipred(f, "x", "x")$pred_df
  expect_true(all(c("SE", "SE_mu") %in% names(p)))
  expect_equal(p$LL, plogis(qlogis(p$y) - qnorm(0.975) * p$SE), tolerance = 1e-10)
  expect_equal(p$SE_mu, p$y * (1 - p$y) * p$SE, tolerance = 1e-10)
})

test_that("CovMat of a Gamma semireg() fit is relative to sigma^2 (no double scaling)", {
  set.seed(2291)
  d <- data.frame(x = runif(200), g = gl(20, 10))
  d$y <- rgamma(200, shape = 2, rate = 2/exp(1 + sin(2 * pi * d$x) + rnorm(20, sd = 0.3)[d$g]))
  fit <- semireg(y ~ s(x, bs = "os") + (1 | g), data = d, family = Gamma(link = "log"))
  s2 <- sigma(fit$semer)^2
  # the fixed-effect block of CovMat times sigma^2 is lme4's covariance
  # conditional on the variance parameters, (RX'RX)^-1: lme4's working weights
  # already hold 1/sigma^2 (vcov(use.hessian = FALSE) multiplies by sigma^2
  # again, and lme4 warns that it disagrees with the Hessian)
  expect_equal(as.matrix(fit$CovMat)[1:2, 1:2] * s2, chol2inv(lme4::getME(fit$semer, "RX")),
               ignore_attr = TRUE, tolerance = 1e-6)
  # and close to vcov() from the Hessian (which also allows for the
  # uncertainty of the variance parameters)
  expect_equal(diag(as.matrix(fit$CovMat))[1:2] * s2, diag(as.matrix(vcov(fit$semer))),
               ignore_attr = TRUE, tolerance = 0.1)
  # the same model with the glmmTMB engine: covariances agree to the
  # difference of the two fits (they disagreed by a factor sigma^2 = 0.48)
  ft <- semireg_tmb(y ~ s(x, bs = "os") + (1 | g), data = d, family = Gamma(link = "log"))
  expect_equal(diag(as.matrix(fit$CovMat))[1:2] * s2, diag(as.matrix(ft$CovMat))[1:2] * sigma(ft$semer)^2,
               tolerance = 0.15)
})

test_that("ZOSull() gives a matrix for one value, and O'Sullivan predictions outside the basis range stop", {
  Z1 <- ZOSull(0.5, intKnots = c(0.25, 0.5, 0.75), range.x = c(0, 1))
  expect_equal(dim(Z1), c(1L, 5L))
  Z <- ZOSull(c(0.2, 0.5), intKnots = c(0.25, 0.5, 0.75), range.x = c(0, 1))
  expect_equal(Z1[1, ], Z[2, ])
  set.seed(5530)
  d <- data.frame(x = runif(80))
  d$y <- sin(2 * pi * d$x) + rnorm(80, sd = 0.3)
  f <- semireg(y ~ s(x, bs = "os"), data = d)
  expect_equal(nrow(semipred(f, covariate = "x", at = list(x = 0.5))$pred_df), 1)
  # outside the basis the B-splines are 0: before 2.0.0 the curve silently fell back to its linear part
  expect_snapshot(semipred(f, covariate = "x", at = list(x = 1.2)), error = TRUE)
})

test_that("the unconditional correction is computed once per fit (fit$cache) and gives the same curves", {
  set.seed(2291)
  d <- data.frame(x = runif(200), g = gl(20, 10))
  d$y <- rpois(200, exp(1 + sin(2 * pi * d$x) + rnorm(20, sd = 0.3)[d$g]))
  for (eng in c("lme4", "glmmTMB")) {
    fit <- if (eng == "lme4") semireg(y ~ s(x, bs = "os") + (1 | g), data = d, family = poisson) else
      semireg_tmb(y ~ s(x, bs = "os") + (1 | g), data = d, family = poisson)
    expect_true(is.environment(fit$cache))
    n <- 0
    real <- predictmeans:::pm_sr_uncond
    a <- with_mocked_bindings(semipred(fit, covariate = "x", unconditional = TRUE),
                              pm_sr_uncond = function(s) { n <<- n + 1; real(s) },
                              .package = "predictmeans")
    b <- with_mocked_bindings(semipred(fit, covariate = "x", unconditional = TRUE),
                              pm_sr_uncond = function(s) { n <<- n + 1; real(s) },
                              .package = "predictmeans")
    expect_identical(n, 1, label = eng)
    expect_identical(a$pred_df, b$pred_df)
    # a fit without a cache (made before 2.0.0) gives the same result
    old <- fit
    old$cache <- NULL
    expect_equal(semipred(old, covariate = "x", unconditional = TRUE)$pred_df, a$pred_df)
  }
})

test_that("semipred(contr = ): the max-|t| test of equal curves; one back-transformation with covariatemeans()", {
  set.seed(4421)
  d <- data.frame(x = rep(seq(0, 1, length = 30), 2), g = factor(rep(c("a", "b"), each = 30)))
  d$y <- sin(3 * d$x) + 0.6 * (d$g == "b") * d$x + rnorm(60, sd = 0.15)
  fl <- semireg(y ~ s(x, by = g, bs = "os", k = 6), data = d)
  r <- semipred(fl, "g", "x", contr = c("b", "a"), band = "simultaneous")
  tst <- r$settings$test
  expect_named(tst, c("statistic", "method", "points", "p.value"))
  expect_equal(tst$statistic, max(abs(r$pairs$statistic)))
  # the p-value is the smallest level at which the simultaneous band leaves out 0
  expect_equal(tst$p.value, min(r$pairs$p.adj))
  expect_lt(tst$p.value, 0.01)
  d$y0 <- sin(3 * d$x) + rnorm(60, sd = 0.15)
  f0 <- semireg(y0 ~ s(x, by = g, bs = "os", k = 6), data = d)
  expect_gt(semipred(f0, "g", "x", contr = c("b", "a"))$settings$test$p.value, 0.05)
  expect_match(paste(capture.output(print(r)), collapse = " "), "Curves equal: max |t|", fixed = TRUE)

  # a response transformed in the formula is back-transformed as covariatemeans() does
  d$z <- exp(1 + sin(3 * d$x) + rnorm(60, sd = 0.15))
  flog <- semireg(log(z) ~ s(x, by = g, bs = "os", k = 6), data = d)
  s1 <- semipred(flog, "g", "x")
  expect_equal(s1$means$response, exp(s1$means$estimate))
  expect_equal(s1$means$response, semipred(flog, "g", "x", trans = exp)$means$response)
  expect_null(semipred(flog, "g", "x", trans = I)$means$response)
  # bias_adjust: the residual variance only (the smooths are the curve, not variation)
  s3 <- semipred(flog, "g", "x", bias_adjust = TRUE)
  expect_equal(s3$means$response, s1$means$response * exp(sigma(flog$semer)^2 / 2))
  expect_match(s3$settings$back_transform$estimate, "bias adjusted")
})
