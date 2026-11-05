# Direct tests of the pure functions behind predictmeans()'s stages.

test_that("pair index and contrast matrix follow column order", {
  pp <- predictmeans:::pm_pair_index(4)
  expect_equal(pp$i, c(2, 3, 3, 4, 4, 4))
  expect_equal(pp$j, c(1, 1, 2, 1, 2, 3))
  CM <- predictmeans:::pm_pair_matrix(3)
  expect_equal(CM, rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1)))
  expect_length(predictmeans:::pm_pair_index(1)$i, 0)
})

test_that("label matrix keeps levels with ':' intact", {
  lab <- data.frame(A = factor(c("a:1", "b")), B = factor(c("x", "y:z")), x = 1:2)
  expect_equal(predictmeans:::pm_label_matrix(lab, 2), matrix(c("a:1", "b", "x", "y:z"), 2))
})

test_that("mean table: t limits, or mean +/- 2 SE at 0.05 with permutation", {
  mt <- data.frame(pm = c(10, 12), ses = c(1, 2), Df = c(5, 5), A = factor(c("a", "b")))
  r <- predictmeans:::pm_mean_table(mt, "A", perm = FALSE, slevel = 0.1, ndecimal = 4)
  expect_equal(r$table[["LL(90%)"]], round(c(10, 12) - qt(0.95, 5) * c(1, 2), 4))
  rp <- predictmeans:::pm_mean_table(mt, "A", perm = TRUE, slevel = 0.1, ndecimal = 4)
  expect_equal(rp$slevel, 0.05)
  expect_equal(rp$table[["LL(95%)"]], c(8, 8))
  expect_true(all(is.na(rp$table$Df)))
})

test_that("compare: SEDs and nested factors drop atvar", {
  d <- data.frame(A = factor(rep(c("a1", "a2"), each = 6)),
                  B = factor(rep(c("b1", "b2", "b3", "b4"), each = 3)))
  set.seed(1); d$y <- rnorm(12)
  m <- lm(y ~ A/B, d)                      # B nested in A
  k <- suppressWarnings(Kmatrix(m, "A:B"))
  mp <- predictmeans:::mymodelparm(m)
  ok <- table(d$A, d$B)[cbind(k$fctnames$A, k$fctnames$B)] > 0   # cells with data
  K <- k$K[ok, mp$estimable, drop = FALSE]
  r <- predictmeans:::pm_compare_sed(K, k$K, k$fctnames[ok, ], k$fctnames, c("A", "B"),
                                     mp$coef, mp$vcov, "none", 0.05, atvar = "A")
  expect_equal(unname(r$SED.out[1:2]), rep(sqrt(2 * sigma(m)^2 / 3), 2))
  expect_null(r$atvar)                     # no pairs within a level of B
})

test_that("LSDs: t quantile, or 2 x SED with permutation", {
  sed <- c(Max.SED = 2, Min.SED = 1, Aveg.SED = 1.5)
  l <- predictmeans:::pm_compare_lsd(sed, NULL, "A", 0.05, 10, 0.05, 4, perm = FALSE)
  expect_equal(as.vector(l), round(qt(0.975, 10) * c(2, 1, 1.5), 5))
  lp <- predictmeans:::pm_compare_lsd(sed, NULL, "A", 0.05, NULL, 0.05, 4, perm = TRUE)
  expect_equal(as.vector(lp), c(4, 2, 3))
})

test_that("means that cannot be estimated are NA; estimable differences kept", {
  set.seed(1)
  d <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), r = 1:4)
  d$x <- as.numeric(d$A == "a3") * 2 + as.numeric(d$B == "b2")   # aliased with A + B
  d$y <- rnorm(nrow(d))
  m <- lm(y ~ A + B + x, d)
  expect_message(r <- predictmeans(m, "A", pairwise = TRUE, plot = FALSE, prtnum = FALSE),
                 "cannot be estimated")
  expect_true(all(is.na(r$means$estimate)) && all(!r$means$estimable))
  expect_true(all(is.na(r$`Predicted Means`)))
  # x = 2*(A == a3) + (B == b2): a1 vs a2 is estimable; a3 vs the others,
  # at a fixed x, depends on the dropped x coefficient and is not
  expect_equal(r$pairs$estimate[1], unname(-coef(m)["Aa2"]))
  expect_true(all(is.na(r$pairs$estimate[2:3])))
  pdf(NULL); on.exit(dev.off())
  expect_message(predictmeans(m, "A", prtnum = FALSE, prtplt = FALSE), "no mean plots")
})
