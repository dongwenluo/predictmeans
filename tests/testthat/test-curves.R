## covariatemeans() results in the layout of predictmeans() (R/curves.R)

test_that("covariatemeans(): classed result, means in the predictmeans layout, no plot unless asked", {
  O <- as.data.frame(nlme::Oats)
  fm <- lm(yield ~ Block + nitro * Variety, data = O)
  xs <- c(0, 0.3, 0.6)
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  r <- covariatemeans(fm, "Variety", covariate = "nitro", covariateV = xs)
  expect_null(grDevices::dev.list())                     # nothing drawn
  expect_s3_class(r, c("covariatemeans", "pm_curves"))
  expect_named(r, c("means", "pairs", "settings", "plt", "pltdf"))
  m <- r$means
  expect_named(m, c("Variety", "nitro", "estimate", "std.error", "df", "conf.low", "conf.high"))
  # the same values as the 1.x table of the plot
  expect_equal(m$estimate, r$pltdf$Mean)
  expect_equal(m$std.error, r$pltdf$ses)
  expect_equal(m$conf.low, r$pltdf$LL)
  expect_equal(m$conf.high, r$pltdf$UL)
  expect_equal(m$nitro, r$pltdf$xvar)
  # and as predictmeans() at the same covariate values
  pm <- predictmeans(fm, "Variety", covariate = list(nitro = 0.3), plot = FALSE)$means
  m3 <- m[m$nitro == 0.3, ]
  expect_equal(m3$estimate, pm$estimate[match(m3$Variety, pm$Variety)], tolerance = 1e-10)
  expect_equal(m3$std.error, pm$std.error[match(m3$Variety, pm$Variety)], tolerance = 1e-10)
  expect_equal(m3$df, pm$df[match(m3$Variety, pm$Variety)])
  expect_identical(as.data.frame(r), m)
  if (requireNamespace("generics", quietly = TRUE)) expect_identical(generics::tidy(r), m)
  expect_error(as.data.frame(r, what = "pairs"), "no pairs")
  expect_s3_class(ggplot2::autoplot(r), "ggplot")
  expect_equal(r$settings$df_method, "residual")
  expect_null(r$settings$back_transform)
  expect_output(print(r), "Curves of 'Variety' along nitro")
  # plot = TRUE draws, as 1.x did
  grDevices::pdf(NULL)
  n0 <- length(grDevices::dev.list())
  on.exit(grDevices::dev.off())
  r2 <- covariatemeans(fm, "Variety", covariate = "nitro", covariateV = xs, plot = TRUE)
  expect_equal(r2$means, m)
})

test_that("covariatemeans(): back-transformed curves have response columns, as predictmeans()", {
  skip_if_not_installed("lme4")
  set.seed(3)
  d <- expand.grid(x = seq(0, 10, length = 11), A = factor(c("a1", "a2", "a3")), site = factor(1:8))
  d$y <- rpois(nrow(d), exp(0.5 + 0.1*d$x + 0.2*(d$A == "a2") + rnorm(8, 0, 0.3)[d$site]))
  g <- glm(y ~ A*x, poisson, data = d)
  r <- suppressMessages(covariatemeans(g, "A", covariate = "x", covariateV = c(2, 5)))
  m <- r$means
  expect_true(all(c("response", "response.std.error", "response.conf.low", "response.conf.high") %in% names(m)))
  expect_equal(m$response, exp(m$estimate))
  expect_equal(m$response.std.error, exp(m$estimate)*m$std.error)
  expect_equal(m$response.conf.low, exp(m$conf.low))
  expect_equal(m$response, r$pltdf$Mean)
  expect_equal(m$response.conf.high, r$pltdf$UL)
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
  rp <- covariatemeans(g, "A", covariate = "x", contr = c(3, 1), covariateV = xs)
  rs <- covariatemeans(g, "A", covariate = "x", contr = c(3, 1), covariateV = xs, band = "simultaneous")
  expect_null(rp$means)
  p <- rp$pairs
  expect_named(p, c("level1", "level2", "x", "estimate", "std.error", "df", "statistic", "p.value", "conf.low",
                    "conf.high", "ratio", "ratio.conf.low", "ratio.conf.high", "response", "response.std.error",
                    "response.conf.low", "response.conf.high"))
  expect_equal(unique(p$level1), "a3")
  expect_equal(unique(p$level2), "a1")
  expect_equal(p$estimate, rp$pltdf$diff)
  expect_equal(p$conf.low, rp$pltdf$LL)
  expect_equal(p$ratio.conf.high, rp$pltdf$ratio.UL)
  expect_equal(p$response.conf.low, rp$pltdf$response.LL)
  expect_equal(p$p.value, 2*pnorm(-abs(p$estimate/p$std.error)))
  expect_identical(as.data.frame(rp), p)                 # the default table: pairs
  # simultaneous: conf.* stay pointwise; adj.conf.* are the Scheffe band of pltdf, and p.adj < 0.05 exactly
  # where they leave out 0
  s <- rs$pairs
  expect_equal(s$conf.low, p$conf.low)
  expect_equal(s$adj.conf.low, rs$pltdf$LL)
  expect_equal(s$adj.conf.high, rs$pltdf$UL)
  expect_equal(s$p.adj, pchisq(s$statistic^2, 2, lower.tail = FALSE))
  expect_equal(s$p.adj < 0.05, rs$pltdf$significant)
  expect_equal(rs$settings$band, "simultaneous")
  expect_equal(rs$settings$contrast, c("a3", "a1"))
  expect_equal(rs$settings$test, rs$test)
  expect_output(print(rs), "Difference a3 - a1 of 'A' along x")
})
