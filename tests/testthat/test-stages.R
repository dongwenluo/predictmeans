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

test_that("mean table: t limits; with permutation the model's t limits at the user's level (ci_df)", {
  mt <- data.frame(pm = c(10, 12), ses = c(1, 2), Df = c(5, 5), A = factor(c("a", "b")))
  r <- predictmeans:::pm_mean_table(mt, "A", perm = FALSE, slevel = 0.1, ndecimal = 4)
  expect_equal(r$table[["LL(90%)"]], round(c(10, 12) - qt(0.95, 5) * c(1, 2), 4))
  # before 2.0.0: mean +/- 2 SE labelled 95% whatever the level
  rp <- predictmeans:::pm_mean_table(mt, "A", perm = TRUE, slevel = 0.1, ndecimal = 4, ci_df = c(8, 8))
  expect_equal(rp$ci_level, 0.1)
  expect_equal(rp$table[["LL(90%)"]], round(c(10, 12) - qt(0.95, 8) * c(1, 2), 4))
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
                                     mp$coef, mp$vcov, atvar = "A")
  expect_equal(unname(r$SED.out[1:2]), rep(sqrt(2 * sigma(m)^2 / 3), 2))
  expect_null(r$atvar)                     # no pairs within a level of B
})

test_that("LSDs: SED x the critical value, or 2 x SED with permutation", {
  sed <- c(Max.SED = 2, Min.SED = 1, Aveg.SED = 1.5)
  l <- predictmeans:::pm_compare_lsd(sed, NULL, "A", qt(0.975, 10), 10, 0.05, 4, perm = FALSE)
  expect_equal(as.vector(l), round(qt(0.975, 10) * c(2, 1, 1.5), 5))
  lp <- predictmeans:::pm_compare_lsd(sed, NULL, "A", NA, NULL, 0.05, 4, perm = TRUE)
  expect_equal(as.vector(lp), c(4, 2, 3))
})

test_that("the printed LSDs go with the adjusted p-values: the critical value of 'adj' and its family", {
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fm <- lm(yield ~ Block + nitro*Variety, data = O)
  df <- df.residual(fm)
  q <- function(...) predictmeans(fm, ..., plot = FALSE)
  sed <- q("nitro")[["Standard Error of Differences"]][["Aveg.SED"]]
  lsd <- function(a) q("nitro", adj = a)$LSD[["Aveg.LSD"]]
  expect_equal(lsd("none"), round(qt(0.975, df)*sed, 5))
  expect_equal(lsd("tukey"), round(qtukey(0.95, 4, df)/sqrt(2)*sed, 5))            # the HSD: 4 means
  expect_equal(lsd("bonferroni"), round(qt(1 - 0.05/(2*6), df)*sed, 5))            # 6 pairs
  expect_equal(lsd("sidak"), round(qt(1 - (1 - 0.95^(1/6))/2, df)*sed, 5))
  expect_equal(lsd("scheffe"), round(sqrt(3*qf(0.95, 3, df))*sed, 5))
  # step-wise methods have no critical value: not adjusted, and a note says so
  expect_equal(lsd("holm"), lsd("none"))
  expect_match(attr(q("nitro", adj = "holm")$LSD, "Note"), "not adjusted")
  expect_match(attr(q("nitro", adj = "tukey")$LSD, "Note"), "adjusted by the 'tukey' method")
  expect_null(attr(q("nitro")$LSD, "Note"))
  # a difference is beyond the LSD exactly when its adjusted p-value is below the level (equal SEDs)
  for (a in c("tukey", "bonferroni", "sidak", "scheffe")) {
    r <- q("nitro", adj = a)
    expect_identical(abs(r$pairs$estimate) > r$LSD[["Aveg.LSD"]], r$pairs$p.adj < 0.05, label = a)
    expect_equal(r$LSD[["Aveg.LSD"]], round((r$pairs$adj.conf.high[1] - r$pairs$adj.conf.low[1])/2, 5), label = a)
  }
  # with 'atvar' the family is the pairs within a level, as for the p-values (was: all 66 pairs)
  ra <- q("nitro:Variety", atvar = "Variety", adj = "bonferroni")
  sed_a <- ra[["Standard Error of Differences"]][["Aveg.SED"]]
  expect_equal(ra$LSD[["Aveg.LSD"]], round(qt(1 - 0.05/(2*6), df)*sed_a, 5))
  rt <- q("nitro:Variety", atvar = "Variety", adj = "tukey")
  expect_equal(rt$LSD[["Aveg.LSD"]], round(qtukey(0.95, 4, df)/sqrt(2)*sed_a, 5))
  # each pair's LSD goes with its adjusted p-value (the 1.x LSD matrix is not printed)
  rb <- q("nitro", adj = "bonferroni")
  expect_equal(rb$pairs$adj.conf.high[1] - rb$pairs$estimate[1], qt(1 - 0.05/12, df)*rb$pairs$std.error[1], tolerance = 1e-8)
  expect_null(rb[["Pairwise LSDs"]])
})

test_that("means that cannot be estimated are NA; estimable differences kept", {
  set.seed(1)
  d <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), r = 1:4)
  d$x <- as.numeric(d$A == "a3") * 2 + as.numeric(d$B == "b2")   # aliased with A + B
  d$y <- rnorm(nrow(d))
  m <- lm(y ~ A + B + x, d)
  expect_message(r <- predictmeans(m, "A", pairwise = TRUE, plot = FALSE),
                 "cannot be estimated")
  expect_true(all(is.na(r$means$estimate)) && all(!r$means$estimable))
  expect_true(all(is.na(r$`Predicted Means`)))
  # x = 2*(A == a3) + (B == b2): a1 vs a2 is estimable; a3 vs the others,
  # at a fixed x, depends on the dropped x coefficient and is not
  expect_equal(r$pairs$estimate[1], unname(-coef(m)["Aa2"]))
  expect_true(all(is.na(r$pairs$estimate[2:3])))
})

test_that("predictmeans() draws nothing; the option predictmeans.plot is ignored", {
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fm <- lm(yield ~ Block + nitro*Variety, data = O)
  f <- tempfile(fileext = ".pdf")
  pages <- function(expr) {
    grDevices::pdf(f); on.exit(grDevices::dev.off())
    force(expr)
    grDevices::dev.off(); on.exit()
    # PDF page objects, read as bytes (the file is not text in a UTF-8 locale)
    length(grepRaw("/Type /Page[^s]", readBin(f, "raw", file.size(f)), all = TRUE))
  }
  expect_equal(pages(predictmeans(fm, "nitro")), 0)
  op <- options(predictmeans.plot = TRUE); on.exit(options(op), add = TRUE)
  expect_equal(pages(predictmeans(fm, "nitro")), 0)
  expect_gt(pages(print(plot(predictmeans(fm, "nitro")))), 0)
})
