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

test_that("neighbours less alike (competition) get their own advice, not 'longer lags only'", {
  d <- vario_field(seed = 4)
  set.seed(4)
  z <- matrix(rnorm(130), 13, 10)
  d$y <- 10 + as.numeric(d$variety) / 5 + as.vector(z[-1, ] - 0.8 * z[-13, ])     # MA(1) down the columns, negative
  v <- residvariogram(lm(y ~ block + variety, data = d), row = "row", col = "col", nsim = 39, plot = FALSE)
  expect_lt(attr(v, "test")$p.value, 0.05)
  expect_true("variogram_neighbours" %in% advice(v)$id)
  expect_false("variogram_trend" %in% advice(v)$id)
})

test_that("several trials: p adjusted (Holm), advice and pm_check use it", {
  skip_on_cran()
  d <- rbind(transform(vario_field(c(0.8, 0.8), seed = 2), trial = "A"), transform(vario_field(seed = 3), trial = "B"))
  fm <- lm(y ~ trial * (block + variety), data = d)
  v <- residvariogram(fm, row = "row", col = "col", by = "trial", plot = FALSE)        # nsim 99, as pm_check()
  tst <- attr(v, "test")
  expect_equal(tst$p.adj, p.adjust(tst$p.value, "holm"))
  expect_output(print(v), "adjusted")
  pc <- pm_check(fm, row = "row", col = "col", by = "trial")
  k <- pc$check == "Spatial pattern (residual variogram)"
  expect_equal(pc$p.value[k], min(tst$p.adj))
  expect_equal(pc$look[k], any(grepl("^variogram_", advice(v)$id) & advice(v)$level == "caution"))
})

test_that("clear stops and notes: too few pairs, small nsim, zero weights, no envelope; na.exclude works", {
  d <- vario_field()
  d4 <- d[d$row <= 4 & d$col <= 4, ]
  expect_error(residvariogram(lm(y ~ 1, data = d4), row = "row", col = "col", nsim = 0, plot = FALSE),
               "No lag has at least 20 pairs")
  fm <- lm(y ~ block + variety, data = d)
  expect_error(residvariogram(fm, row = "row", col = "col", nsim = 9, plot = FALSE), "at least 19")
  expect_error(residvariogram(fm, row = "row", col = "col", max_lag = 0, plot = FALSE), "at least 1")
  dw <- transform(d, w = c(0, rep(1, 119)))
  expect_error(residvariogram(lm(y ~ block + variety, data = dw, weights = w), row = "row", col = "col",
                              nsim = 0, plot = FALSE), "weight 0")
  dn <- d
  dn$y[c(5, 40)] <- NA
  fe <- lm(y ~ block + variety, data = dn, na.action = na.exclude)
  fo <- lm(y ~ block + variety, data = dn)
  expect_equal(residvariogram(fe, row = "row", col = "col", nsim = 19, plot = FALSE),
               residvariogram(fo, row = "row", col = "col", nsim = 19, plot = FALSE))
  dq <- transform(d, n = rpois(120, 5))
  expect_message(v <- residvariogram(glm(n ~ variety, quasipoisson, data = dq), row = "row", col = "col",
                                     nsim = 19, plot = FALSE), "quasi family")
  expect_equal(attr(v, "nsim"), 0)
  pc <- suppressMessages(pm_check(glm(n ~ variety, quasipoisson, data = dq), row = "row", col = "col"))
  expect_equal(pc$status[pc$check == "Spatial pattern (residual variogram)"], "i")
})

test_that("an incomplete grid gets a note; directions with no pairs are dropped; antimeridian", {
  d <- vario_field()[-seq(3, 120, by = 7), ]                     # 18 of 120 plots missing
  v <- residvariogram(lm(y ~ block + variety, data = d), row = "row", col = "col", nsim = 0, plot = FALSE)
  expect_true("variogram_incomplete" %in% advice(v)$id)
  p <- data.frame(x = 1:30, y = 0, z = rnorm(30))
  v <- residvariogram(p$z, coords = c("x", "y"), data = p, directions = c(0, 90), segment = 10, step = 1,
                      min_pairs = 1, nsim = 0, plot = FALSE)
  expect_equal(unique(v$direction), 0)
  a <- data.frame(lon = c(-179.9, 179.9, 0), lat = c(0, 0.2, 60), z = c(0, 1, 2))
  v <- residvariogram(a$z, coords = c("lon", "lat"), data = a, lonlat = TRUE, directions = c(0, 135),
                      segment = 20, step = 50, min_pairs = 1, nsim = 0, plot = FALSE)
  expect_equal(v$direction, 135)
})

test_that("the grid variogram reproduces Genstat FVARIOGRAM on the Slate Hall residuals", {
  # dev/genstat/variogram-slatehall.gen, run 2026-11-25: Genstat keeps the two diagonals apart (signed row
  # lag); the grid lag (dr, dc) pools them, so the reference is their pair-weighted mean
  p <- function(f) test_path("_data", "variogram", f)
  sh <- read.csv(p("slatehall.csv"))
  sh <- sh[order(sh$fieldrow, sh$fieldcolumn), ]
  r <- read.csv(p("slatehall-genstat-residuals.csv"))
  g <- read.csv(p("slatehall-genstat-grid.csv"))
  fm <- lm(yield ~ factor(replicates) + factor(variety), data = sh)
  expect_lt(max(abs(residuals(fm) - r$res1)), 1e-8)                     # Genstat prints 8 decimals
  g$lag_row <- abs(g$i - 6)
  g$lag_col <- g$j - 1
  ref <- do.call(rbind, lapply(split(g, list(g$model, g$lag_row, g$lag_col), drop = TRUE), function(z)
    data.frame(model = z$model[1], lag_row = z$lag_row[1], lag_col = z$lag_col[1],
               gamma = if (nrow(z) == 1 || z$lag_col[1] == 0) z$gamma[1] else sum(z$gamma * z$n) / sum(z$n),
               n = if (nrow(z) == 1 || z$lag_col[1] == 0) z$n[1] else sum(z$n))))
  for (m in c("res1", "res2")) {                                          # regression; AR1 x AR1 REML
    v <- residvariogram(r[[m]], row = "fieldrow", col = "fieldcolumn", data = r, max_lag = c(5, 7), min_pairs = 1,
                        nsim = 0, plot = FALSE)
    x <- merge(v, ref[ref$model == m, ], by = c("lag_row", "lag_col"))
    expect_equal(nrow(x), 47L)
    expect_equal(x$n_pairs, x$n)
    expect_lt(max(abs(x$gamma.x - x$gamma.y)), 1e-7)
  }
  # from the fit: residuals scaled by their SD (equal leverages here), so gamma / var(residuals)
  v <- residvariogram(fm, row = "fieldrow", col = "fieldcolumn", data = sh, max_lag = c(5, 7), min_pairs = 1,
                      nsim = 0, plot = FALSE)
  x <- merge(v, ref[ref$model == "res1", ], by = c("lag_row", "lag_col"))
  expect_equal(x$gamma.x, x$gamma.y / var(r$res1), tolerance = 1e-6)
})
