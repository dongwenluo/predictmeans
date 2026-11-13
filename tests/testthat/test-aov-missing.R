## aov() with Error() and missing values: estimated as Genstat's ANOVA does
## (the residual sum of squares of the lowest stratum made smallest), the
## completed design analysed by strata, exact SEDs; a message recommends REML.

am_oats <- function() {
  o <- as.data.frame(nlme::Oats)
  o$Variety <- factor(gsub(" ", "", o$Variety))
  o$Nitro <- factor(sprintf("%.1f", o$nitro))
  o$Block <- factor(as.character(o$Block), levels = c("I", "II", "III", "IV", "V", "VI"))
  o
}

test_that("one sub-plot missing: Genstat's missing value, stratum mean squares and means (Genstat run 2, pattern C)", {
  G <- read.csv(test_path("_data", "genstat", "genstat-reml.csv.gz"), colClasses = "character")
  g <- G[G$analysis == "run2:YC:anova", ]
  o <- am_oats()
  d <- droplevels(o[!(o$Block == "II" & o$Variety == "Victory" & o$Nitro == "0.6"), ])
  a <- aov(yield ~ Variety * Nitro + Error(Block/Variety), data = d)
  expect_message(sa <- pm_aov_strata(pm_aov_carrier(a)), "REML")
  expect_equal(attr(sa$strata, "missing")$estimate, as.numeric(g$value[g$table == "missing_value"]), tolerance = 0.05 / 113)
  expect_equal(sa$strata$df[-1], c(5, 10, 44))
  res <- strsplit(g$value[g$table == "aov_table" & g$term == "Residual"], "\\|")
  expect_equal(sa$strata$ms[3], as.numeric(res[[1]][3]), tolerance = 0.05 / 568)        # 568.2
  expect_equal(sa$strata$ms[4], as.numeric(utils::tail(res[[2]], 1)), tolerance = 0.05 / 178)   # 178.7 (after "(1)")
  for (tm in c("Variety", "Nitro", "Variety:Nitro")) {
    p <- suppressMessages(predictmeans(a, tm, plot = FALSE, prtnum = FALSE))
    gm <- g[g$table == "aov_mean" & g$term == sub(":", ".", tm), ]
    lab <- do.call(paste, c(lapply(p$means[strsplit(tm, ":")[[1]]], as.character), sep = ":"))
    expect_equal(round(p$means$estimate[match(gm$key, lab)], 1), as.numeric(gm$value), tolerance = 1e-8, label = tm)
  }
  # pairs that do not involve the missing plot's cell: Genstat's SEDs and df
  p <- suppressMessages(predictmeans(a, "Variety:Nitro", plot = FALSE, prtnum = FALSE))$pairs
  free <- !grepl("Victory:0.6", p$level1) & !grepl("Victory:0.6", p$level2)
  v1 <- sub(":.*", "", p$level1); v2 <- sub(":.*", "", p$level2)
  expect_equal(round(unique(round(p$std.error[free & v1 == v2], 6)), 2), 7.72)
  expect_equal(round(unique(round(p$std.error[free & v1 != v2], 6)), 2), 9.59)
  expect_equal(round(unique(round(p$df[free & v1 != v2], 4)), 2), 31.42)
  # the data with the row kept and its response NA give the same
  d2 <- o
  d2$yield[o$Block == "II" & o$Variety == "Victory" & o$Nitro == "0.6"] <- NA
  a2 <- aov(yield ~ Variety * Nitro + Error(Block/Variety), data = d2)
  p2 <- suppressMessages(predictmeans(a2, "Variety:Nitro", plot = FALSE, prtnum = FALSE))
  expect_equal(p2$means$estimate, suppressMessages(predictmeans(a, "Variety:Nitro", plot = FALSE,
                                                                prtnum = FALSE))$means$estimate, tolerance = 1e-10)
})

test_that("randomized blocks with one value missing: Yates' formula", {
  set.seed(4)
  d <- expand.grid(trt = factor(1:5), blk = factor(1:4))
  d$y <- 10 + as.integer(d$trt) + rnorm(4)[d$blk] + rnorm(20)
  y <- d$y; y[7] <- NA; d$y <- y
  sa <- suppressMessages(pm_aov_strata(pm_aov_carrier(aov(y ~ trt + Error(blk), data = d))))
  tt <- 5; bb <- 4
  T <- sum(y[d$trt == d$trt[7]], na.rm = TRUE); B <- sum(y[d$blk == d$blk[7]], na.rm = TRUE); G <- sum(y, na.rm = TRUE)
  expect_equal(attr(sa$strata, "missing")$estimate, (tt * T + bb * B - G) / ((tt - 1) * (bb - 1)), tolerance = 1e-8)
  expect_equal(sa$strata$df[sa$strata$stratum == "Within"], 11)                  # 12 - 1 missing value
})

test_that("a whole plot missing: the next stratum up determines its total", {
  o <- am_oats()
  o$yield[(o$Block == "I" & o$Variety == "GoldenRain")] <- NA
  a <- suppressWarnings(aov(yield ~ Variety * Nitro + Error(Block/Variety), data = o))
  sa <- suppressMessages(pm_aov_strata(pm_aov_carrier(a)))
  expect_equal(nrow(attr(sa$strata, "missing")), 4)
  expect_equal(sa$strata$df[-1], c(5, 9, 42))
})

test_that("complete orthogonal designs: no message; designs that are not orthogonal: a message", {
  o <- am_oats()
  expect_no_message(pm_aov_strata(pm_aov_carrier(aov(yield ~ Variety * Nitro + Error(Block/Variety), data = o))))
  set.seed(2)
  d <- data.frame(blk = factor(rep(1:4, c(6, 6, 5, 7))), trt = factor(c(rep(1:3, 2), rep(1:3, 2), 1:3, 1, 2,
                                                                       rep(1:3, 2), 1)))
  d$y <- rnorm(nrow(d))
  expect_message(pm_aov_strata(pm_aov_carrier(aov(y ~ trt + Error(blk), data = d))), "not orthogonal")
})
