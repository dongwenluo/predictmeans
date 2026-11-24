## Means over empty cells keep the 1.x rule (the average of the observed
## cells), but two means that average over different cells are compared on
## the cells both have (R/empty-cells.R; dev/review-defaults.md: comparing
## the averages rejected 8-78% of true nulls in dev/sim-empty-cells.R, the
## common cells 4-6%).

ec_factorial <- function(drop) {
  d <- expand.grid(rep = 1:3, B = factor(paste0("B", 1:4)), A = factor(paste0("A", 1:3)))
  d <- d[!paste(d$A, d$B, sep = ":") %in% drop, ]
  set.seed(7)
  d$y <- 10 + 2 * as.integer(d$B) + rnorm(nrow(d))
  d
}

test_that("lm: pairs of means over empty cells use the cells both means have", {
  d <- ec_factorial("A1:B4")
  fm <- lm(y ~ A * B, data = d)
  cm <- tapply(d$y, list(d$A, d$B), mean)
  expect_warning(p <- predictmeans(fm, "A", pairwise = TRUE, plot = FALSE),
                 "No data for A:B A1:B4. The mean of A = A1 averages only its observed cells")
  expect_equal(p$means$estimate, unname(rowMeans(cm, na.rm = TRUE)), tolerance = 1e-8)     # means: 1.x
  pr <- p$pairs
  expect_equal(pr$estimate[pr$level1 == "A1" & pr$level2 == "A2"], mean(cm["A1", 1:3] - cm["A2", 1:3]), tolerance = 1e-8)
  expect_equal(pr$estimate[pr$level1 == "A2" & pr$level2 == "A3"], mean(cm["A2", ] - cm["A3", ]), tolerance = 1e-8)
  expect_identical(pr$over, c("B B1, B2, B3", "B B1, B2, B3", NA))
  # SE: sigma^2 / 3 replicates, two cells per difference, averaged over 3 levels
  s2 <- summary(fm)$sigma^2
  expect_equal(pr$std.error[1], sqrt(s2 * 2 / 3 / 3), tolerance = 1e-8)
  # the same in contrastmeans() and, for the whole term, pm_joint()
  cmt <- suppressWarnings(suppressMessages(contrastmeans(fm, "A", ctrmatrix = rbind(c(1, -1, 0)))))[[1]]
  expect_equal(unname(cmt[1, "Estimate"]), round(pr$estimate[1], 4))
  jt <- suppressWarnings(suppressMessages(pm_joint(fm, "A")))
  expect_identical(jt$over, "B B1, B2, B3")
  Fv <- anova(lm(y ~ A * B, data = d[d$B != "B4", ]))
  expect_equal(jt$F.ratio, {
    sub <- d[d$B != "B4", ]
    m <- lm(y ~ A * B, data = sub)
    K <- Kmatrix(m, "A")$K
    L <- rbind(K[1, ] - K[3, ], K[2, ] - K[3, ])
    b <- coef(m)
    drop(t(L %*% b) %*% solve(L %*% vcov(m) %*% t(L) * s2 / summary(m)$sigma^2, L %*% b)) / 2
  }, tolerance = 1e-8)
})

test_that("no cell in common: the pair is NA", {
  d <- ec_factorial(c("A1:B3", "A1:B4", "A3:B1", "A3:B2"))
  fm <- lm(y ~ A * B, data = d)
  p <- suppressWarnings(predictmeans(fm, "A", pairwise = TRUE, plot = FALSE))
  r <- p$pairs[p$pairs$level1 == "A1" & p$pairs$level2 == "A3", ]
  expect_true(is.na(r$estimate) && is.na(r$std.error))
  expect_identical(r$over, "none")
  expect_identical(p$pairs$over[p$pairs$level1 == "A1" & p$pairs$level2 == "A2"], "B B1, B2")
  expect_true(all(p$means$estimable))
})

test_that("no warning and no change where no mean averages over an empty cell", {
  d <- ec_factorial("A1:B4")
  fm <- lm(y ~ A * B, data = d)
  expect_no_warning(p <- predictmeans(fm, "A:B", pairwise = TRUE, plot = FALSE))
  expect_equal(nrow(p$means), 11)
  expect_null(p$pairs$over)
  # an additive model: the means average over all cells, estimable, no warning
  fa <- lm(y ~ A + B, data = d)
  expect_no_warning(pa <- predictmeans(fa, "A", plot = FALSE))
  expect_null(pa$pairs$over)
})

test_that("lmer split plot with a lost cell: the comparisons on the varieties both rates have", {
  O <- as.data.frame(nlme::Oats)
  O$N <- factor(O$nitro)
  O$Variety <- factor(O$Variety, ordered = FALSE)
  Om <- subset(O, !(Variety == "Victory" & N == "0.6"))
  fm <- suppressMessages(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = Om))
  expect_warning(p <- predictmeans(fm, "N", pairwise = TRUE, plot = FALSE), "N = 0.6")
  K <- Kmatrix(fm, "N:Variety")$K
  gm <- c("Golden Rain", "Marvellous")
  L <- colMeans(K[paste0("0.6:", gm), ] - K[paste0("0:", gm), ])
  r <- p$pairs[p$pairs$level1 == "0" & p$pairs$level2 == "0.6", ]
  expect_equal(r$estimate, -sum(L * lme4::fixef(fm)), tolerance = 1e-8)
  expect_equal(r$estimate, -42.5, tolerance = 1e-3)
  expect_identical(r$over, "Variety Golden Rain, Marvellous")
  # pairs of the complete levels are unchanged
  expect_equal(p$pairs$estimate[p$pairs$level1 == "0" & p$pairs$level2 == "0.2"], -19.5, tolerance = 1e-8)
  expect_true(is.na(p$pairs$over[p$pairs$level1 == "0" & p$pairs$level2 == "0.2"]))
})

test_that("a combination without data is not an empty cell when its factors do not interact (Sample_date/Farm)", {
  # farms within dates (date/farm: a farm's effect differs by date), treatments
  # crossed with dates; no data for treatment B on farm F2 at date d2. Farm
  # and treatment do not interact in the model, so its cells are those of
  # date:farm and of date:trt, which all have data: the mean of d2:B over both
  # farms can be estimated, as the mean of d2:A.
  set.seed(11)
  d <- expand.grid(rep = 1:3, trt = factor(c("A", "B")), farm = factor(c("F1", "F2")), date = factor(c("d1", "d2")))
  d <- d[!(d$date == "d2" & d$farm == "F2" & d$trt == "B"), ]
  d$y <- 5 + as.integer(d$farm) * as.integer(d$date) + 0.5 * (d$trt == "B") + rnorm(nrow(d), 0, 0.3)
  fm <- lm(y ~ date/farm + date * trt, data = d)
  expect_no_warning(p <- predictmeans(fm, "date:trt", pairwise = TRUE, plot = FALSE))
  pr <- p$pairs
  expect_null(pr$over)
  pred <- function(dt, f, t) unname(predict(fm, data.frame(date = dt, farm = f, trt = t)))
  over_farms <- function(dt, t) mean(c(pred(dt, "F1", t), pred(dt, "F2", t)))
  expect_equal(p$means$estimate, c(over_farms("d1", "A"), over_farms("d1", "B"), over_farms("d2", "A"),
                                   over_farms("d2", "B")), tolerance = 1e-8)
  # treatments at d2: the same difference on either farm
  r <- pr[pr$level1 == "d2:A" & pr$level2 == "d2:B", ]
  expect_equal(r$estimate, pred("d2", "F1", "A") - pred("d2", "F1", "B"), tolerance = 1e-8)
  r <- pr[pr$level1 == "d1:B" & pr$level2 == "d2:B", ]
  expect_equal(r$estimate, over_farms("d1", "B") - over_farms("d2", "B"), tolerance = 1e-8)
  # with a farm x treatment interaction in the model it is an empty cell
  fm3 <- lm(y ~ date/farm * trt, data = d)
  expect_warning(p3 <- predictmeans(fm3, "date:trt", pairwise = TRUE, plot = FALSE), "d2:F2:B")
  expect_equal(p3$means$estimate[4], unname(predict(fm3, data.frame(date = "d2", farm = "F1", trt = "B"))),
               tolerance = 1e-8)
  r <- p3$pairs[p3$pairs$level1 == "d2:A" & p3$pairs$level2 == "d2:B", ]
  expect_identical(r$over, "farm F1")
})

test_that("the empty-cell warning names the means by their factor levels, once", {
  set.seed(2)
  d <- expand.grid(B = factor(c("b1", "b2")), C = factor(c("c1", "c2", "c3")), rep = 1:3)
  d$x <- rnorm(nrow(d)); d$y <- as.numeric(d$B) + d$x + rnorm(nrow(d))
  e <- d[!(d$B == "b2" & d$C == "c1"), ]
  w <- character()
  # (B and C must interact for b2:c1 to be an empty cell)
  withCallingHandlers(invisible(pm_slopes(lm(y ~ x*B + B*C, data = e), "x", "B")),
                      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 1)
  expect_match(w, "The mean of B = b2 averages", fixed = TRUE)
  expect_warning(covariatemeans(lm(y ~ x + B*C, data = e), "B", covariate = "x", plot = FALSE),
                 "The mean of B = b2 averages", fixed = TRUE)
})

test_that("factors confounded within another factor (each pasture on one farm): the combinations that exist", {
  # dates x farms, each pasture on one farm at a date; Farm and Pasture do not
  # interact in the model and are not nested over the whole data (pasture P1
  # is on farm F1 early and on F3 late), so the mean over the Farm x Pasture
  # cross cannot be estimated. The means average over the combinations that
  # exist, with a warning naming the two factors, and equal those of the model
  # without the aliased factor and of the model with the nesting written.
  set.seed(5)
  d <- rbind(expand.grid(date = c("d1", "d2"), past = c("P1", "P2", "P3", "P4"), trt = c("A", "B"), rep = 1:3),
             expand.grid(date = "d3", past = c("P1", "P4"), trt = c("A", "B"), rep = 1:3))
  d$farm <- ifelse(d$date == "d3", ifelse(d$past == "P1", "F3", "F4"), ifelse(d$past == "P4", "F2", "F1"))
  for (v in c("date", "past", "trt", "farm")) d[[v]] <- factor(d[[v]])
  d$y <- as.integer(d$date) + as.integer(d$past)/2 + 0.5*(d$trt == "B") + rnorm(nrow(d), 0, 0.3)
  fm <- lm(y ~ date/farm + date/(trt + past)^2, data = d)
  expect_true(anyNA(coef(fm)))
  expect_warning(p <- predictmeans(fm, "date:trt", pairwise = TRUE, plot = FALSE),
                 "'farm' and 'past' are confounded in the data.*farm/past")
  expect_false(anyNA(p$means$estimate))
  f0 <- lm(y ~ date/(trt + past)^2, data = d)                    # without the aliased factor
  p0 <- suppressWarnings(predictmeans(f0, "date:trt", pairwise = TRUE, plot = FALSE))
  expect_equal(p$means$estimate, p0$means$estimate, tolerance = 1e-8)
  expect_equal(p$means$std.error, p0$means$std.error, tolerance = 1e-8)
  ok <- !is.na(p0$pairs$estimate)
  expect_equal(p$pairs$estimate[ok], p0$pairs$estimate[ok], tolerance = 1e-8)
  # the cell means averaged with equal weights over the pastures at each date
  cm <- aggregate(y ~ date + trt + past, d, mean)
  expect_equal(p$means$estimate, aggregate(y ~ date + trt, cm, mean)$y[order(rep(1:3, 2), rep(1:2, each = 3))],
               tolerance = 1e-8)
  # weights = "cells": no warning, the combinations weighted by their numbers
  expect_no_warning(pw <- suppressMessages(predictmeans(fm, "date:trt", weights = "cells", plot = FALSE)))
  expect_false(anyNA(pw$means$estimate))
})

test_that("means over empty cells are labelled, and nothing read off them contradicts their pair tests", {
  skip_on_cran()
  q <- function(e) suppressWarnings(suppressMessages(e))
  Oats <- as.data.frame(nlme::Oats); Oats$N <- factor(Oats$nitro)
  Om <- subset(Oats, !(Variety == "Victory" & N == "0.6"))
  fm <- q(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = Om))
  pm <- q(predictmeans(fm, "Variety", pairwise = TRUE, plot = FALSE))
  m <- pm$means
  # the cells each mean averages, as 'over' of the pairs
  expect_equal(m$over, c(NA, NA, "N 0, 0.2, 0.4"))
  expect_equal(which(names(m) == "over"), which(names(m) == "estimate") - 1)
  # no comparison interval for it: its pairs are on common cells
  expect_true(all(is.na(m[3, c("comparison.low", "comparison.high")])))
  expect_false(anyNA(m$comparison.low[1:2]))
  cmp <- pm$settings$comparison
  expect_equal(cmp$left_out, "Victory")
  expect_equal(cmp$pairs, 1)
  # the multiplier from the one pair left: overlap = its test
  p <- pm$pairs[1, ]
  expect_equal(cmp$mult, qt(0.975, p$df)*p$std.error/sum(m$std.error[1:2]), tolerance = 1e-10)
  # letter_ci: letters from the p-values instead (Marvellous vs Victory: p = 0.07)
  expect_message(pl <- suppressWarnings(predictmeans(fm, "Variety", letter_ci = TRUE, plot = FALSE)),
                 "letters are from the p-values")
  expect_equal(pl$means$group, pm$means$group)
  expect_equal(pm$means$group, c("A", "A", "A"))
  expect_match(attr(pm$LSD, "Note"), "compared on common cells")
  # the plot: an open point, no comparison bar, a caption
  g <- ggplot2::autoplot(pm, interval = "comparison")
  b <- ggplot2::ggplot_build(g)
  expect_equal(pm_layer(g, "GeomPoint")$shape, c(16, 16, 1))
  expect_true(all(is.na(pm_layer(g, "GeomErrorbar")[3, c("ymin", "ymax")])))
  expect_match(g$labels$caption, "Open points")
  expect_equal(nrow(pm_layer(ggplot2::autoplot(pm, interval = "ci"), "GeomErrorbar")), 3)   # confidence limits: all
  # complete data: no 'over' column, no caption, intervals for all
  pf <- q(predictmeans(q(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = Oats)), "Variety",
                       pairwise = TRUE, plot = FALSE))
  expect_null(pf$means$over)
  expect_length(pf$settings$comparison$left_out, 0)
  expect_null(ggplot2::autoplot(pf, interval = "comparison")$labels$caption)
})

test_that("confounded factors: means labelled; intervals kept where the pairs are the differences of the means", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  set.seed(9)
  dc <- expand.grid(C = factor(c("c1", "c2")), vol = factor(c("v1", "v2", "v3", "v4")), rep = 1:3)
  dc$drop <- factor(ifelse(dc$vol == "v3", "fine", "coarse"))
  dc$y <- as.numeric(dc$vol) + as.numeric(dc$C) + rnorm(nrow(dc))
  pv <- q(predictmeans(lm(y ~ C*vol + drop, data = dc), "vol", pairwise = TRUE, plot = FALSE))
  expect_equal(pv$means$over, c("drop coarse", "drop coarse", "drop fine", "drop coarse"))
  # each pair is the difference of the shown means: the intervals stay
  m <- setNames(pv$means$estimate, pv$means$vol)
  expect_equal(abs(pv$pairs$estimate), unname(abs(m[pv$pairs$level1] - m[pv$pairs$level2])), tolerance = 1e-10)
  expect_length(pv$settings$comparison$left_out, 0)
  expect_false(anyNA(pv$means$comparison.low))
  expect_null(attr(pv$LSD, "Note"))
})
