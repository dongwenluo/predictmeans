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
  expect_warning(p <- predictmeans(fm, "A", pairwise = TRUE, plot = FALSE, prtnum = FALSE),
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
  p <- suppressWarnings(predictmeans(fm, "A", pairwise = TRUE, plot = FALSE, prtnum = FALSE))
  r <- p$pairs[p$pairs$level1 == "A1" & p$pairs$level2 == "A3", ]
  expect_true(is.na(r$estimate) && is.na(r$std.error))
  expect_identical(r$over, "none")
  expect_identical(p$pairs$over[p$pairs$level1 == "A1" & p$pairs$level2 == "A2"], "B B1, B2")
  expect_true(all(p$means$estimable))
})

test_that("no warning and no change where no mean averages over an empty cell", {
  d <- ec_factorial("A1:B4")
  fm <- lm(y ~ A * B, data = d)
  expect_no_warning(p <- predictmeans(fm, "A:B", pairwise = TRUE, plot = FALSE, prtnum = FALSE))
  expect_equal(nrow(p$means), 11)
  expect_null(p$pairs$over)
  # an additive model: the means average over all cells, estimable, no warning
  fa <- lm(y ~ A + B, data = d)
  expect_no_warning(pa <- predictmeans(fa, "A", plot = FALSE, prtnum = FALSE))
  expect_null(pa$pairs$over)
})

test_that("lmer split plot with a lost cell: the comparisons on the varieties both rates have", {
  O <- as.data.frame(nlme::Oats)
  O$N <- factor(O$nitro)
  O$Variety <- factor(O$Variety, ordered = FALSE)
  Om <- subset(O, !(Variety == "Victory" & N == "0.6"))
  fm <- suppressMessages(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = Om))
  expect_warning(p <- predictmeans(fm, "N", pairwise = TRUE, plot = FALSE, prtnum = FALSE), "N = 0.6")
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
