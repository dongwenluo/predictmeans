## covariatemeans() results in the layout of predictmeans() (R/curves.R)

test_that("covariatemeans(): classed result, means in the predictmeans layout, no plot unless asked", {
  O <- as.data.frame(nlme::Oats)
  fm <- lm(yield ~ Block + nitro * Variety, data = O)
  xs <- c(0, 0.3, 0.6)
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  r <- covariatemeans(fm, "Variety", covariate = "nitro", at = xs)
  expect_null(grDevices::dev.list())                     # nothing drawn
  expect_s3_class(r, c("covariatemeans", "pm_curves"))
  expect_named(r, c("means", "pairs", "settings", "plt", "pred_df"))
  m <- r$means
  expect_named(m, c("Variety", "nitro", "estimate", "std.error", "df", "conf.low", "conf.high"))
  # the same values as the 1.x table of the plot
  expect_equal(m$estimate, r$pred_df$Mean)
  expect_equal(m$std.error, r$pred_df$ses)
  expect_equal(m$conf.low, r$pred_df$LL)
  expect_equal(m$conf.high, r$pred_df$UL)
  expect_equal(m$nitro, r$pred_df$xvar)
  # and as predictmeans() at the same covariate values
  pm <- predictmeans(fm, "Variety", covariate = list(nitro = 0.3), plot = FALSE)$means
  m3 <- m[m$nitro == 0.3, ]
  expect_equal(m3$estimate, pm$estimate[match(m3$Variety, pm$Variety)], tolerance = 1e-10)
  expect_equal(m3$std.error, pm$std.error[match(m3$Variety, pm$Variety)], tolerance = 1e-10)
  expect_equal(m3$df, pm$df[match(m3$Variety, pm$Variety)])
  expect_identical(as.data.frame(r), m)
  if (requireNamespace("generics", quietly = TRUE)) expect_identical(generics::tidy(r), m)
  expect_error(as.data.frame(r, what = "pairs"), "no 'pairs' table")
  expect_s3_class(ggplot2::autoplot(r), "ggplot")
  expect_equal(r$settings$df_method, "residual")
  expect_null(r$settings$back_transform)
  expect_output(print(r), "Curves of 'Variety' along nitro")
  # plot = TRUE draws, as 1.x did
  grDevices::pdf(NULL)
  n0 <- length(grDevices::dev.list())
  on.exit(grDevices::dev.off())
  r2 <- covariatemeans(fm, "Variety", covariate = "nitro", at = xs, plot = TRUE)
  expect_equal(r2$means, m)
})

test_that("covariatemeans(): back-transformed curves have response columns, as predictmeans()", {
  skip_if_not_installed("lme4")
  set.seed(3)
  d <- expand.grid(x = seq(0, 10, length = 11), A = factor(c("a1", "a2", "a3")), site = factor(1:8))
  d$y <- rpois(nrow(d), exp(0.5 + 0.1*d$x + 0.2*(d$A == "a2") + rnorm(8, 0, 0.3)[d$site]))
  g <- glm(y ~ A*x, poisson, data = d)
  r <- suppressMessages(covariatemeans(g, "A", covariate = "x", at = c(2, 5)))
  m <- r$means
  expect_true(all(c("response", "response.std.error", "response.conf.low", "response.conf.high") %in% names(m)))
  expect_equal(m$response, exp(m$estimate))
  expect_equal(m$response.std.error, exp(m$estimate)*m$std.error)
  expect_equal(m$response.conf.low, exp(m$conf.low))
  expect_equal(m$response, r$pred_df$Mean)
  expect_equal(m$response.conf.high, r$pred_df$UL)
  expect_equal(r$settings$back_transform$by, "the inverse link")
  pm <- suppressMessages(predictmeans(g, "A", covariate = list(x = 5), plot = FALSE))$means
  m5 <- m[m$x == 5, ]
  expect_equal(m5$response.std.error, pm$response.std.error[match(m5$A, pm$A)], tolerance = 1e-8)
  expect_equal(m5$response.conf.low, pm$response.conf.low[match(m5$A, pm$A)], tolerance = 1e-8)
})

test_that("covariatemeans(contr = ): pairs in the predictmeans layout; Scheffe p.adj agrees with its limits", {
  set.seed(1)
  d <- data.frame(A = gl(3, 40, labels = c("a1", "a2", "a3")), x = runif(120, 0, 2))
  d$cnt <- rpois(120, exp(1 + 0.6*d$x*(d$A == "a3") + 0.2*d$x))
  g <- glm(cnt ~ A*x, poisson, d)
  xs <- seq(0, 2, 0.1)
  rp <- covariatemeans(g, "A", covariate = "x", contr = c(3, 1), at = xs)
  rs <- covariatemeans(g, "A", covariate = "x", contr = c(3, 1), at = xs, band = "simultaneous")
  expect_null(rp$means)
  p <- rp$pairs
  expect_named(p, c("level1", "level2", "x", "estimate", "std.error", "df", "statistic", "p.value", "conf.low",
                    "conf.high", "ratio", "ratio.conf.low", "ratio.conf.high", "response", "response.std.error",
                    "response.conf.low", "response.conf.high"))
  expect_equal(unique(p$level1), "a3")
  expect_equal(unique(p$level2), "a1")
  expect_equal(p$estimate, rp$pred_df$diff)
  expect_equal(p$conf.low, rp$pred_df$LL)
  expect_equal(p$ratio.conf.high, rp$pred_df$ratio.UL)
  expect_equal(p$response.conf.low, rp$pred_df$response.LL)
  expect_equal(p$p.value, 2*pnorm(-abs(p$estimate/p$std.error)))
  expect_identical(as.data.frame(rp), p)                 # the default table: pairs
  # simultaneous: conf.* stay pointwise; adj.conf.* are the Scheffe band of pred_df, and p.adj < 0.05 exactly
  # where they leave out 0
  s <- rs$pairs
  expect_equal(s$conf.low, p$conf.low)
  expect_equal(s$adj.conf.low, rs$pred_df$LL)
  expect_equal(s$adj.conf.high, rs$pred_df$UL)
  expect_equal(s$p.adj, pchisq(s$statistic^2, 2, lower.tail = FALSE))
  expect_equal(s$p.adj < 0.05, rs$pred_df$significant)
  expect_equal(rs$settings$band, "simultaneous")
  expect_equal(rs$settings$contrast, c("a3", "a1"))
  expect_equal(unlist(rs$settings$test), setNames(rs$test, c("df1", "df2", "F.ratio", "p.value")))
  expect_output(print(rs), "Difference a3 - a1 of 'A' along x")
})

test_that("semipred(): the layout of predictmeans(); the simultaneous band in adj.conf.*, with p.adj that agrees", {
  data(Dialyzer, package = "nlme")
  fm <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject), data = Dialyzer)
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  out <- semipred(fm, "QB", "pressure", at = c(0.5, 1, 1.5, 2))
  expect_null(grDevices::dev.list())                     # nothing drawn
  expect_s3_class(out, c("semipred", "pm_curves"))
  m <- out$means
  expect_named(m, c("QB", "pressure", "estimate", "std.error", "df", "conf.low", "conf.high"))
  expect_equal(m$estimate, out$pred_df$rate)
  expect_equal(m$std.error, out$pred_df$SE)
  expect_equal(m$df, out$pred_df$df)
  expect_equal(m$conf.low, out$pred_df$LL)
  expect_identical(as.data.frame(out), m)
  expect_match(out$settings$df_method, "Satterthwaite")
  expect_false(out$settings$unconditional)
  expect_s3_class(ggplot2::autoplot(out), "ggplot")
  expect_output(print(out), "Curves of 'QB' along pressure")
  # simultaneous: conf.* stay pointwise, adj.conf.* are the band of pred_df
  os <- semipred(fm, "QB", "pressure", band = "simultaneous")
  op <- semipred(fm, "QB", "pressure")
  expect_equal(os$means$conf.low, op$means$conf.low)
  expect_equal(os$means$adj.conf.low, os$pred_df$LL)
  expect_true(all(os$means$adj.conf.high > os$means$conf.high))
  # the difference of two curves: p.adj < 0.05 exactly where the band leaves out 0
  d <- semipred(fm, "QB", "pressure", contr = c("300", "200"), band = "simultaneous")
  p <- d$pairs
  expect_named(p, c("level1", "level2", "pressure", "estimate", "std.error", "df", "statistic", "p.value", "p.adj",
                    "conf.low", "conf.high", "adj.conf.low", "adj.conf.high"))
  expect_equal(unique(p$level1), "300")
  expect_equal(p$estimate, d$pred_df$rate)
  expect_equal(p$std.error, d$pred_df$SE_contr)
  expect_equal(p$adj.conf.high, d$pred_df$UL)
  expect_equal(p$p.adj < 0.05, d$pred_df$significant)
  expect_true(all(p$p.adj >= p$p.value))
  expect_equal(d$significant, attr(d$pred_df, "significant"))
  expect_identical(as.data.frame(d), p)
  expect_equal(d$settings$contrast, c("300", "200"))
  # slopes, with the test of slope = 0
  s <- semipred(fm, "QB", "pressure", deriv = 1, band = "simultaneous")
  expect_null(s$means)
  expect_named(s$slopes, c("QB", "pressure", "estimate", "std.error", "df", "statistic", "p.value", "p.adj",
                           "conf.low", "conf.high", "adj.conf.low", "adj.conf.high"))
  expect_equal(s$slopes$adj.conf.low, s$pred_df$LL)
  expect_equal(s$slopes$p.adj < 0.05, s$pred_df$significant)
  expect_equal(s$settings$deriv, 1)
  expect_identical(as.data.frame(s), s$slopes)
  expect_output(print(s), "Slopes of 'QB'")
})

test_that("semipred(): generalized fits have response columns; surfaces keep the boundary", {
  set.seed(4)
  db <- expand.grid(x = seq(0, 1, length = 25), g = factor(c("a", "b")), cl = factor(1:12))
  db$y <- rbinom(nrow(db), 1, plogis(-1 + 2*sin(pi*db$x) + 0.8*(db$g == "b")*db$x + rnorm(12, 0, 0.5)[db$cl]))
  fb <- semireg(y ~ s(x, by = g, bs = "os", k = 6) + (1 | cl), family = binomial, data = db)
  b <- semipred(fb, "g", "x", at = c(0.2, 0.5, 0.8))
  m <- b$means
  expect_equal(m$response, plogis(m$estimate))
  expect_equal(m$response, b$pred_df$y)
  expect_equal(m$response.conf.low, b$pred_df$LL)
  expect_equal(m$response.std.error, b$pred_df$SE_mu)
  expect_equal(b$settings$back_transform$estimate, "typical cluster (random effects at zero)")
  expect_true(b$settings$unconditional)
  bc <- semipred(fb, "g", "x", at = c(0.2, 0.5, 0.8), contr = c("b", "a"))
  expect_equal(bc$pairs$odds.ratio, bc$pred_df$odds.ratio)
  expect_equal(bc$pairs$odds.ratio.conf.low, bc$pred_df$odds.ratio.LL)
  expect_equal(bc$pairs$response, bc$pred_df$response.diff)
  expect_equal(bc$pairs$response.conf.high, bc$pred_df$response.UL)
  # a surface inside semiboundary(): the rows of pred_df
  data(copper)
  fcu <- semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6), data = copper)
  s2 <- semipred(fcu, covariate = c("xcoord", "ycoord"), trans = exp, boundary = TRUE)
  expect_equal(nrow(s2$means), nrow(s2$pred_df))
  expect_equal(s2$means$response, s2$pred_df$Trans_y)
  expect_equal(s2$means$response.std.error, s2$means$response*s2$means$std.error, tolerance = 1e-6)
  expect_true(s2$settings$boundary)
  expect_output(print(s2), "Surface over xcoord, ycoord")
})
