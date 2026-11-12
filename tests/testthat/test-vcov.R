## vcov = : robust, clustered or given covariance of the coefficients

oats_het <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  set.seed(1)
  d$y2 <- d$yield + c(0, 0, 0, 40)[as.integer(d$nitro)] * rnorm(nrow(d))
  d
}

test_that("HC0-HC3 for lm equal car::hccm()", {
  skip_if_not_installed("car")
  d <- oats_het()
  fm <- lm(y2 ~ Block + nitro + Variety, data = d)
  for (ty in c("HC0", "HC1", "HC2", "HC3")) {
    expect_equal(unname(pm_vcov_sandwich(fm, ty)), unname(car::hccm(fm, type = tolower(ty))), tolerance = 1e-10)
  }
})

test_that("HC3 for glm by the sandwich package's definitions (bread, meat, dispersion)", {
  set.seed(2)
  d <- data.frame(A = factor(rep(1:3, each = 10)), x = rnorm(30))
  d$y <- rgamma(30, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x + (d$A == "3")))
  g <- glm(y ~ A + x, family = Gamma(link = "log"), data = d)
  X <- model.matrix(g)
  w <- weights(g, "working")
  wres <- residuals(g, "working") * w
  disp <- sum(wres^2) / sum(w)
  ef <- wres * X / disp
  n <- nrow(X)
  bread <- summary(g)$cov.unscaled * n * disp
  h <- hatvalues(g)
  meat <- crossprod(ef / (1 - h)) / n
  expect_equal(unname(pm_vcov_sandwich(g, "HC3")), unname(bread %*% meat %*% bread / n), tolerance = 1e-8)
})

test_that("vcov = ~cluster: CR1 by hand, with missing responses", {
  d <- oats_het()
  d$y2[c(3, 40)] <- NA
  fm <- lm(y2 ~ nitro + Variety, data = d)
  V <- pm_vcov_resolve(fm, ~ Block, names(coef(fm)))
  dd <- d[!is.na(d$y2), ]
  X <- model.matrix(fm)
  B <- solve(crossprod(X))
  e <- residuals(fm)
  M <- Reduce(`+`, lapply(split(seq_len(nrow(X)), dd$Block), function(i) tcrossprod(crossprod(X[i, ], e[i]))))
  G <- nlevels(dd$Block); n <- nrow(X); p <- ncol(X)
  expect_equal(unname(V[, ]), unname(G / (G - 1) * (n - 1) / (n - p) * B %*% M %*% B), tolerance = 1e-10)
  expect_identical(attr(V, "label"), "CR1, clusters: Block")
})

test_that("predictmeans(), contrastmeans() and pm_slopes() with vcov = \"HC3\" equal emmeans with that covariance", {
  skip_if_not_installed("car")
  skip_if_not_installed("emmeans")
  d <- oats_het()
  fm <- lm(y2 ~ Block + nitro + Variety, data = d)
  Vh <- car::hccm(fm, type = "hc3")
  pm <- predictmeans(fm, "nitro", vcov = "HC3", plot = FALSE, prtnum = FALSE)
  em <- summary(emmeans::emmeans(fm, ~ nitro, vcov. = Vh))
  expect_equal(pm$means$std.error, em$SE, tolerance = 1e-8)
  ep <- summary(pairs(emmeans::emmeans(fm, ~ nitro, vcov. = Vh), reverse = TRUE, adjust = "none"))
  expect_equal(pm$pairs$std.error, ep$SE, tolerance = 1e-8)
  expect_equal(pm$pairs$p.value, ep$p.value, tolerance = 1e-8)
  expect_identical(pm$settings$vcov, "HC3")
  expect_false(pm$settings$vcov_adjust)
  # contrasts and their joint Wald F
  cn <- rbind(c(-1, 1, 0, 0), c(0, -1, 1, 0))
  ct <- contrastmeans(fm, "nitro", ctrmatrix = cn, vcov = "HC3", joint = TRUE)
  K <- Kmatrix(fm, "nitro")$K
  L <- cn %*% K
  b <- L %*% coef(fm)
  expect_equal(unname(ct[[1]][, "Std. Error"]), round(sqrt(diag(L %*% Vh %*% t(L))), 4))
  Fh <- drop(t(b) %*% solve(L %*% Vh %*% t(L), b)) / 2
  expect_equal(ct[[3]]$F.ratio, round(Fh, 4))
  expect_match(attr(ct[[1]], "Note"), "HC3")
  # slopes
  d$N <- as.numeric(as.character(d$nitro))
  fs <- lm(y2 ~ N * Variety + Block, data = d)
  sl <- pm_slopes(fs, "N", "Variety", vcov = "HC3")
  et <- summary(emmeans::emtrends(fs, ~ Variety, var = "N", vcov. = car::hccm(fs, type = "hc3")))
  expect_equal(sl$slopes$std.error, et$SE, tolerance = 1e-6)
})

test_that("a given matrix or function; lmer keeps its df; errors and permlist", {
  skip_on_cran()
  d <- oats_het()
  m <- lme4::lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = d)
  p0 <- predictmeans(m, "nitro", vcov_adjust = FALSE, plot = FALSE, prtnum = FALSE)
  pk <- predictmeans(m, "nitro", plot = FALSE, prtnum = FALSE)
  p1 <- predictmeans(m, "nitro", vcov = as.matrix(vcov(m)) * 1.2, plot = FALSE, prtnum = FALSE)
  expect_equal(p1$means$std.error, p0$means$std.error * sqrt(1.2))
  expect_equal(p1$means$df, pk$means$df)
  expect_identical(p1$settings$vcov, "user matrix")
  p2 <- predictmeans(m, "nitro", vcov = function(x) as.matrix(vcov(x)) * 2, plot = FALSE, prtnum = FALSE)
  expect_equal(p2$means$std.error, p0$means$std.error * sqrt(2))
  expect_error(predictmeans(m, "nitro", vcov = "HC3", plot = FALSE, prtnum = FALSE), "lm, aov and glm")
  expect_error(predictmeans(m, "nitro", vcov = TRUE, plot = FALSE, prtnum = FALSE), "vcov_adjust")
  expect_error(predictmeans(m, "nitro", vcov = diag(2), plot = FALSE, prtnum = FALSE), "matrix")
  fm <- lm(yield ~ nitro, data = d)
  expect_error(predictmeans(fm, "nitro", vcov = "HC5", plot = FALSE, prtnum = FALSE), "not available")
  r <- permmodels(fm, nperm = 19, seed = 1, prt = FALSE)
  expect_warning(pp <- predictmeans(fm, "nitro", vcov = "HC3", permlist = r, plot = FALSE, prtnum = FALSE),
                 "not used with 'permlist'")
  expect_identical(pp$settings$vcov, "model")
})
