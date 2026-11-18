## Several hundred means (lines of a breeding trial): SEDs and LSDs are
## summarised from the covariance of the means, without a row per pair.

test_that("the SED and LSD summaries without a pair table equal those from the pairs", {
  set.seed(11)
  d <- expand.grid(A = factor(1:4), B = factor(1:3), Block = factor(1:8))
  d$y <- as.numeric(d$A) + rnorm(8)[d$Block] + rnorm(nrow(d))
  d$cnt <- rpois(nrow(d), 4)
  d <- d[-c(3, 17, 40), ]                        # unequal SEDs
  mods <- list(lm = lm(y ~ A*B + Block, d),
               lmer = lme4::lmer(y ~ A*B + (1|Block), d),
               lme = nlme::lme(y ~ A*B, random = ~1|Block, data = d),
               glmer = suppressMessages(lme4::glmer(cnt ~ A*B + (1|Block), d, family = poisson)))
  for (nm in names(mods)) for (term in c("A", "A:B")) {
    full <- suppressMessages(predictmeans(mods[[nm]], term))
    op <- options(predictmeans.max_means = 2)
    expect_message(short <- predictmeans(mods[[nm]], term), "pairs are not tabulated")
    options(op)
    expect_null(short$pairs)
    expect_false(is.null(full$pairs))
    for (part in c("Predicted Means", "Standard Error of Means", "Standard Error of Differences", "LSD", "mean_table"))
      expect_equal(short[[part]], full[[part]], tolerance = 1e-8, info = paste(nm, term, part))
    # (comparison intervals need the pairs)
    expect_equal(short$means, full$means[setdiff(names(full$means), c("comparison.low", "comparison.high"))],
                 tolerance = 1e-8, info = paste(nm, term))
  }
})

test_that("nested factors: no pairs within the levels of the outer factor, as from the pairs", {
  d <- data.frame(G = factor(rep(c("a", "b"), each = 12)), T = factor(rep(1:6, each = 4)))
  set.seed(4); d$y <- rnorm(24) + as.numeric(d$T)
  m <- lm(y ~ G/T, d)
  full <- suppressMessages(suppressWarnings(predictmeans(m, "G:T")))
  op <- options(predictmeans.max_means = 2); on.exit(options(op))
  short <- suppressMessages(suppressWarnings(predictmeans(m, "G:T")))
  expect_equal(short[["Standard Error of Differences"]], full[["Standard Error of Differences"]], tolerance = 1e-8)
  expect_equal(short$LSD, full$LSD, tolerance = 1e-8)
})

test_that("300 lines: summaries in a second; pairs when asked, without letter groups", {
  skip_on_cran()
  set.seed(5)
  d <- expand.grid(Line = factor(sprintf("L%03d", 1:300)), Rep = factor(1:2))
  d$y <- rnorm(300, sd = 0.7)[d$Line] + rnorm(2)[d$Rep] + rnorm(nrow(d))
  m <- lm(y ~ Line + Rep, d)
  tm <- system.time(expect_message(pm <- predictmeans(m, "Line"), "44850 pairs"))[["elapsed"]]
  expect_lt(tm, 5)
  expect_null(pm$pairs)
  expect_equal(nrow(pm$means), 300)
  # balanced: one SED, sqrt(2 s2/2)
  expect_equal(unname(pm[["Standard Error of Differences"]]), rep(summary(m)$sigma, 3), tolerance = 1e-8)
  expect_equal(unname(pm$LSD[1]), round(qt(0.975, df.residual(m))*summary(m)$sigma, 5), tolerance = 1e-8)
  out <- capture.output(print(summary(pm)))
  expect_true(any(grepl("^SED: min", out)) && any(grepl("^LSD \\(5%, 299 df\\)", out)))
  # the half-LSD bars of the plot use the summary
  expect_s3_class(ggplot2::autoplot(pm, interval = "lsd"), "ggplot")
  expect_error(ggplot2::autoplot(pm, type = "pvalue"), "pairwise = TRUE")
  d2 <- droplevels(d[d$Line %in% levels(d$Line)[1:120], ])
  expect_message(pp <- predictmeans(lm(y ~ Line + Rep, d2), "Line", pairwise = TRUE),
                 "No letter groups for 120 means")
  expect_equal(nrow(pp$pairs), 120*119/2)
  expect_true(all(is.na(pp$mean_table$Table$LetterGrp)))
})

test_that("PMplot() and autoplot(type = 'pvalue') draw the same square level plot", {
  set.seed(2013)
  pmatrix <- matrix(0, 8, 8)
  pmatrix[lower.tri(pmatrix)] <- runif(28)
  g <- PMplot(pmatrix)
  expect_s3_class(g, "ggplot")
  b <- ggplot2::ggplot_build(g)
  expect_equal(nrow(b$data[[1]]), 28)                          # one tile per pair
  expect_equal(range(b$data[[1]]$x), c(1, 7))                  # the empty last column and first row are left out
  expect_equal(sort(unique(b$data[[1]]$y)), -(8:2))            # (reversed y axis: first level at the top)
  expect_equal(g$coordinates$ratio, 1)                 # square tiles
  expect_equal(nrow(ggplot2::ggplot_build(PMplot(pmatrix, values = TRUE))$data[[2]]), 28)
  expect_equal(nlevels(PMplot(pmatrix, conf.level = 0.9)$data$band), 2)
  expect_equal(nlevels(g$data$band), 4)
  expect_message(expect_null(PMplot(matrix(NA_real_, 3, 3))), "no p-values")
  expect_error(PMplot(matrix(0, 2, 3)), "square")
  # a list: blocks down the diagonal
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fm <- lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = O)
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "BH")
  gl <- PMplot(pm$p_valueMatrix)
  ga <- ggplot2::autoplot(pm, type = "pvalue")
  expect_equal(nrow(gl$data), 18)
  expect_equal(gl$data[c("row", "col", "band")], ga$data[c("row", "col", "band")])
  expect_equal(sort(ga$data$p), sort(pm$pairs$p.adj[!is.na(pm$pairs$by)]))
  expect_length(ga$layers, 3)                                  # tiles and the lines between the blocks
  # predictmeans() draws it when asked
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  predictmeans(fm, "nitro:Variety", pairwise = TRUE, plot = "pval")
  grDevices::dev.off()
  expect_gt(length(grepRaw("/Type /Page[^s]", readBin(f, "raw", file.size(f)), all = TRUE)), 0)
})
