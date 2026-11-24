oats_pw <- function() transform(as.data.frame(nlme::Oats), nitro = factor(nitro))

test_that("pairwise only sets what is printed; adj, atvar and permlist ask for the comparisons", {
  fm <- lm(yield ~ nitro * Variety, data = oats_pw())
  # NULL (default): shown when comparisons are asked for, as 1.x
  d <- predictmeans(fm, "nitro")
  expect_false("Pairwise p-value" %in% names(d))
  expect_false(d$settings$adjust != "none")
  t <- predictmeans(fm, "nitro", adj = "tukey")
  expect_true("Pairwise p-value" %in% names(t))
  # FALSE hides the matrices and letters, the pairs are the same
  f <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", pairwise = FALSE)
  a <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  expect_false(any(c("Golden Rain", "p_valueMatrix") %in% names(f)))
  expect_null(f$means$group)
  expect_equal(f$pairs, a$pairs)
  expect_equal(f$LSD, a$LSD)
  expect_identical(names(f$mean_table)[1:2], c("Variety", "nitro"))
  # TRUE without adj: the unadjusted p-values and letters
  p <- predictmeans(fm, "nitro", pairwise = TRUE)
  expect_true(all(c("Pairwise p-value", "p_valueMatrix") %in% names(p)))
  expect_false("Pairwise LSDs" %in% names(p))                    # each pair's LSD is pairs$lsd
  # many means: adj compares all the pairs even with pairwise = FALSE
  op <- options(predictmeans.max_means = 5)
  on.exit(options(op))
  expect_equal(nrow(predictmeans(fm, "nitro:Variety", adj = "holm", pairwise = FALSE)$pairs), 66)
  expect_null(suppressMessages(predictmeans(fm, "nitro:Variety"))$pairs)
  expect_error(predictmeans(fm, "nitro", pairwise = "yes"), "TRUE, FALSE or NULL")
})

test_that("one LSD when all are equal; the range, a note and advice when they differ", {
  d <- oats_pw()
  eq <- predictmeans(lm(yield ~ nitro * Variety, data = d), "nitro", adj = "tukey")
  expect_true(attr(eq$LSD, "equal"))
  expect_false(grepl("LSDs differ between pairs", attr(eq$LSD, "Note")))
  expect_false("lsd_varies" %in% advice(eq)$id)
  set.seed(4471)
  fu <- lm(yield ~ nitro * Variety, data = d[-sample(nrow(d), 10), ])
  un <- predictmeans(fu, "nitro:Variety", adj = "tukey")
  expect_false(attr(un$LSD, "equal"))
  expect_match(attr(un$LSD, "Note"), "The LSDs differ between pairs")
  expect_match(advice(un)$text[advice(un)$id == "lsd_varies"], "^2 of 66 pairs are judged differently")
  out <- capture.output(print(summary(un)))
  expect_true(any(grepl("they differ", out)))
  # no comparisons asked for: no advice, the summary is still there
  expect_false("lsd_varies" %in% advice(predictmeans(fu, "nitro:Variety"))$id)
})
