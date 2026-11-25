# residvariogram(): estimator against Genstat and brute force, layout checks,
# envelope from simulations, advice, pm_check() hook.

vario_field <- function(phi = c(0, 0), seed = 1) {
  set.seed(seed)
  d <- expand.grid(row = 1:12, col = 1:10)
  d$block <- factor((d$row - 1) %/% 3 + 1)
  d$variety <- factor(unlist(lapply(1:4, function(b) sample(rep(1:15, 2)))))
  Rr <- phi[1]^abs(outer(1:12, 1:12, "-"))
  Rc <- phi[2]^abs(outer(1:10, 1:10, "-"))
  E <- t(chol(kronecker(Rc, Rr) + diag(0.3, 120))) %*% rnorm(120)
  d$y <- 10 + as.numeric(d$variety) / 5 + as.vector(E)
  d
}

## brute force: every pair, absolute row and column lags
vario_brute <- function(r, c, e, robust = FALSE) {
  ij <- which(upper.tri(diag(length(e))), arr.ind = TRUE)
  dr <- abs(r[ij[, 1]] - r[ij[, 2]]); dc <- abs(c[ij[, 1]] - c[ij[, 2]])
  dz <- e[ij[, 1]] - e[ij[, 2]]
  key <- paste(dr, dc)
  N <- tapply(dz, key, length)
  g <- if (robust) tapply(sqrt(abs(dz)), key, mean)^4 / (2 * (0.457 + 0.494 / N)) else tapply(dz^2, key, mean) / 2
  list(g = g, N = N)
}

test_that("the distance variogram reproduces Genstat FVARIOGRAM, Example 8.3.1 (Brooms Barn)", {
  b <- read.csv(test_path("_data", "variogram", "broomsbarn.csv"), strip.white = TRUE, na.strings = "*")
  b <- b[!is.na(b$K), ]
  ref <- read.table(test_path("_data", "variogram", "genstat-example-8.3.1.txt"), header = TRUE)
  v <- residvariogram(log10(b$K), coords = c("East", "North"), data = b, directions = c(0, 45, 90, 135),
                      segment = 45, step = 1, max_lag = 13, min_pairs = 1, nsim = 0, plot = FALSE)
  expect_equal(nrow(v), 48L)
  v <- v[order(v$direction, v$lag), ]
  ref <- ref[order(ref$dir, ref$lag), ]
  expect_equal(v$n_pairs, ref$n)
  expect_lt(max(abs(v$gamma - ref$gamma)), 1e-5)                     # Genstat prints 5 decimals
  expect_lt(max(abs(v$lag - ref$lag)), 5e-4)                         # 3 decimals
})

test_that("the grid variogram equals the sum over all pairs, with missing plots and both diagonals", {
  set.seed(3)
  d <- expand.grid(row = 1:7, col = 1:6)[-c(4, 11, 30), ]
  e <- rnorm(nrow(d))
  bf <- vario_brute(d$row, d$col, e)
  v <- residvariogram(e, row = "row", col = "col", data = d, max_lag = c(6, 5), min_pairs = 1, nsim = 0,
                      plot = FALSE)
  k <- paste(v$lag_row, v$lag_col)
  expect_setequal(k, names(bf$g))
  expect_equal(v$gamma, as.vector(bf$g[k]), tolerance = 1e-12)
  expect_equal(v$n_pairs, as.vector(bf$N[k]))
  vr <- residvariogram(e, row = "row", col = "col", data = d, max_lag = c(6, 5), min_pairs = 1, nsim = 0,
                       robust = TRUE, plot = FALSE)
  bfr <- vario_brute(d$row, d$col, e, robust = TRUE)
  expect_equal(vr$gamma, as.vector(bfr$g[paste(vr$lag_row, vr$lag_col)]), tolerance = 1e-12)
})

test_that("distances: haversine km, and the class rule [k, k + 1) x step", {
  # one degree of latitude: 111.195 km on the mean Earth radius
  d <- data.frame(lon = c(0, 0, 0), lat = c(0, 1, 2), z = c(1, 2, 4))
  v <- residvariogram(d$z, coords = c("lon", "lat"), data = d, lonlat = TRUE, step = 100, max_lag = 300,
                      min_pairs = 1, nsim = 0, plot = FALSE)
  expect_equal(v$lag, c(6371.0088 * pi / 180, 2 * 6371.0088 * pi / 180), tolerance = 1e-9)
  expect_equal(v$gamma, c(mean(c(1, 4)) / 2, 9 / 2))
})

test_that("positions that repeat (rows numbered within blocks) stop; with 'by' each trial has its own", {
  d <- vario_field()
  fm <- lm(y ~ block + variety, data = d)
  d$row_in_block <- (d$row - 1) %% 3 + 1
  expect_error(residvariogram(fm, row = "row_in_block", col = "col", data = d, nsim = 0, plot = FALSE),
               "not within blocks")
  expect_error(residvariogram(fm, row = "row", col = "variety", nsim = 0, plot = FALSE), "share row")
  expect_error(residvariogram(fm, row = "row", col = "col", coords = c("row", "col"), plot = FALSE), "either")
  d2 <- rbind(transform(d, trial = "A"), transform(d, trial = "B", y = y + rnorm(120)))
  fm2 <- lm(y ~ trial * (block + variety), data = d2)
  v <- residvariogram(fm2, row = "row", col = "col", by = "trial", nsim = 0, plot = FALSE)
  expect_equal(unique(v$trial), c("A", "B"))
  expect_error(residvariogram(fm2, row = "row", col = "col", nsim = 0, plot = FALSE), "with 'by'")
})

test_that("lm: the simulated residuals are those of refits to data simulated from the model", {
  d <- vario_field()
  fm <- lm(y ~ block + variety, data = d)
  S <- predictmeans:::pm_vario_sims(fm, rstandard(fm), factor(rep(1, 120)), rep(TRUE, 120), 2, seed = 5, ncore = 1)
  z <- predictmeans:::pm_with_seed(5, matrix(rnorm(240), 120))
  for (b in 1:2) {
    r <- rstandard(lm(y ~ block + variety, data = transform(d, y = fitted(fm) + z[, b])))
    expect_equal(S[b, ] / sd(S[b, ]), unname(r / sd(r)), tolerance = 1e-10)   # the scale is removed per trial
  }
})

test_that("the envelope holds the null; spatial correlation gives the advice (lm)", {
  v0 <- residvariogram(lm(y ~ block + variety, data = vario_field(seed = 2)), row = "row", col = "col",
                       nsim = 39, plot = FALSE)
  expect_s3_class(v0, "pm_variogram")
  expect_equal(c(max(v0$lag_row), max(v0$lag_col)), c(6, 5))            # half the field by default
  expect_gt(attr(v0, "test")$p.value, 0.05)
  expect_equal(nrow(advice(v0)), 0L)
  v1 <- residvariogram(lm(y ~ block + variety, data = vario_field(c(0.8, 0.8), seed = 2)), row = "row",
                       col = "col", nsim = 39, plot = FALSE)
  expect_lt(attr(v1, "test")$p.value, 0.05)
  expect_true("variogram_spatial" %in% advice(v1)$id)
  expect_output(print(v1), "Global envelope test")
  expect_s3_class(autoplot(v1), "ggplot")
  expect_s3_class(autoplot(v1, type = "surface"), "ggplot")
})

test_that("mixed models: simulate and refit; the result does not depend on ncore", {
  skip_on_cran()
  d <- vario_field(c(0.6, 0.7))
  fm <- suppressMessages(lmer(y ~ variety + (1 | block), data = d))
  v1 <- residvariogram(fm, row = "row", col = "col", nsim = 19, plot = FALSE)
  v2 <- residvariogram(fm, row = "row", col = "col", nsim = 19, ncore = 2, plot = FALSE)
  expect_identical(v1, v2)
  expect_equal(attr(v1, "nsim"), 19L)
  fg <- nlme::gls(y ~ variety, correlation = nlme::corExp(form = ~ row + col), data = d)
  vg <- residvariogram(fg, row = "row", col = "col", nsim = 19, plot = FALSE)
  expect_equal(attr(vg, "nsim"), 19L)
})

test_that("pm_check() adds the spatial row only when given a layout", {
  skip_on_cran()
  d <- vario_field(c(0.8, 0.8), seed = 2)
  fm <- lm(y ~ block + variety, data = d)
  expect_false(any(grepl("Spatial", pm_check(fm)$check)))
  pc <- pm_check(fm, row = "row", col = "col")
  k <- pc$check == "Spatial pattern (residual variogram)"
  expect_true(any(k))
  expect_true(pc$look[k])
})
