oats_fac <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}

test_that("autoplot() draws the means with one LSD bar per atvar panel", {
  skip_on_cran()
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
  p <- ggplot2::autoplot(pm)
  expect_s3_class(p, "ggplot")
  expect_equal(sort(pm_layer(p, "GeomPoint")$y), sort(pm$means$estimate))
  # panels are the levels of atvar
  expect_s3_class(p$facet, "FacetWrap")
  # one LSD bar per panel: the average Tukey LSD of the pairs within it
  bl <- pm_layer(p, "GeomErrorbar")
  expect_equal(nrow(bl), 3)
  w <- pm$pairs[!is.na(pm$pairs$by), ]
  lsd <- tapply((w$adj.conf.high - w$adj.conf.low)/2, w$by, mean)
  expect_equal(bl$ymax - bl$ymin, as.vector(lsd[levels(oats_fac()$Variety)]), tolerance = 1e-8)
  expect_match(gsub("\n", " ", p$labels$title), "average LSD \\(5%, tukey\\)")
  expect_null(p$labels$caption)                                # equal LSDs: no range note
  # SED bar; confidence limits on every mean; none
  sd <- pm_layer(ggplot2::autoplot(pm, interval = "sed"), "GeomErrorbar")
  expect_equal(sd$ymax - sd$ymin, rep(mean(w$std.error), 3), tolerance = 1e-8)
  ci <- pm_layer(ggplot2::autoplot(pm, interval = "ci"), "GeomErrorbar")
  expect_equal(sort(ci$ymax), sort(pm$means$conf.high))
  expect_false(any(vapply(ggplot2::autoplot(pm, interval = "none")$layers,
                          function(l) inherits(l$geom, "GeomErrorbar"), logical(1))))
})

test_that("without atvar or pairwise: one average LSD (or SED) for all the means; unequal LSDs are noted", {
  d <- oats_fac()
  fm <- lm(yield ~ nitro*Variety, data = d)
  pm <- predictmeans(fm, "nitro:Variety")
  b <- pm_layer(ggplot2::autoplot(pm), "GeomErrorbar")
  expect_equal(nrow(b), 1)
  expect_equal(b$ymax - b$ymin, pm$LSD[["Aveg.LSD"]], tolerance = 1e-4)
  set.seed(4471)
  fu <- lm(yield ~ nitro*Variety, data = d[-sample(nrow(d), 10), ])
  g <- ggplot2::autoplot(predictmeans(fu, "nitro:Variety"))
  expect_match(g$labels$caption, "LSDs range from")
  # many means, pairs not tabulated: the summary's average LSD
  op <- options(predictmeans.max_means = 5)
  on.exit(options(op))
  big <- suppressMessages(predictmeans(fm, "nitro:Variety"))
  expect_null(big$pairs)
  expect_equal(diff(unlist(pm_layer(ggplot2::autoplot(big), "GeomErrorbar")[c("ymin", "ymax")])),
               big$LSD[["Aveg.LSD"]], tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("x, group and panel place the factors; type = 'bar'; lines = FALSE", {
  fm <- lm(yield ~ nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro:Variety")
  g <- ggplot2::autoplot(pm, x = "Variety", group = "nitro")
  expect_identical(g$labels$x, "Variety")
  expect_match(deparse(g$mapping$colour), "nitro")
  expect_error(ggplot2::autoplot(pm, x = "Varety"), "must name factors")
  gb <- ggplot2::autoplot(pm, type = "bar")
  expect_equal(sort(pm_layer(gb, "GeomCol")$y), sort(pm$means$estimate))
  expect_false(any(vapply(ggplot2::autoplot(pm, lines = FALSE)$layers, function(l) inherits(l$geom, "GeomLine"),
                          logical(1))))
  expect_identical(pm_plot_layout(c("A", "B", "C"), NULL), list(x = "A", group = "B", panel = "C"))
  expect_identical(pm_plot_layout(c("A", "B", "C"), c("C", "B")), list(x = "A", group = NULL, panel = c("C", "B")))
  expect_s3_class(ggplot2::autoplot(predictmeans(lm(yield ~ nitro*Variety*Block, oats_fac()),
                                                 "nitro:Variety:Block", atvar = c("Variety", "Block"))
                                    )$facet, "FacetGrid")
})

test_that("type = 'ci': horizontal limits, on the response scale when back-transformed, with the data", {
  d <- oats_fac()
  fm <- lm(log(yield) ~ nitro*Variety, data = d)
  pm <- suppressMessages(predictmeans(fm, "nitro:Variety"))
  g <- ggplot2::autoplot(pm, type = "ci", data = TRUE)
  eb <- pm_layer(g, "GeomErrorbar")
  expect_equal(sort(eb$xmax), sort(pm$means$response.conf.high))
  expect_equal(sort(pm_layer(g, "GeomPoint", 1)$x), sort(d$yield))          # the data, original scale
  expect_identical(g$labels$x, "Median of yield")
  gm <- ggplot2::autoplot(pm, data = TRUE)                                  # model scale: log(yield)
  expect_equal(sort(pm_layer(gm, "GeomPoint", 1)$y), sort(log(d$yield)))
  expect_error(ggplot2::autoplot(pm, type = "response"), "type = \"ci\"")
  fl <- lm(yield ~ nitro, data = d)
  expect_error(ggplot2::autoplot(predictmeans(fl, "nitro"), scale = "response"), "not back-transformed")
  # a GLM: the data only on the response scale
  data("Drug", package = "predictmeans", envir = environment())
  fb <- glm(cbind(R, N - R) ~ Drug + log10Dose, data = Drug, family = binomial)
  pb <- predictmeans(fb, "Drug")
  expect_message(ggplot2::autoplot(pb, data = TRUE), "not shown on the link scale")
  gc <- ggplot2::autoplot(pb, type = "ci", data = TRUE)
  expect_equal(sort(pm_layer(gc, "GeomPoint", 1)$x), sort(Drug$R / Drug$N))
})

test_that("scales: fixed puts every LSD bar at one place; free_y at each panel's lowest point", {
  d <- oats_fac()
  pm <- predictmeans(lm(yield ~ nitro*Variety, data = d), "nitro:Variety", atvar = "Variety")
  fx <- pm_layer(ggplot2::autoplot(pm), "GeomErrorbar")
  expect_equal(length(unique(fx$ymin)), 1)
  expect_equal(unique(fx$ymin), min(pm$means$estimate))
  gf <- ggplot2::autoplot(pm, scales = "free_y")
  fr <- pm_layer(gf, "GeomErrorbar")
  low <- tapply(pm$means$estimate, pm$means$Variety, min)
  expect_equal(sort(fr$ymin), sort(as.vector(low)))
  expect_identical(gf$facet$params$free$y, TRUE)
  expect_match(gf$labels$caption, "within its own panel")
  expect_error(ggplot2::autoplot(pm, scales = "free_x"))
})

test_that("results keep the observations for data = TRUE, not the model", {
  d <- oats_fac()
  fm <- lm(log(yield) ~ nitro * Variety, data = d)
  pm <- suppressMessages(predictmeans(fm, "nitro:Variety"))
  expect_null(attr(pm, "model"))
  pd <- attr(pm, "plot_data")
  expect_identical(names(pd), c("nitro", "Variety", ".y", ".y_orig"))
  expect_equal(pd$.y_orig, d$yield)
  # small: far below the model when saved
  sz <- function(o) { f <- tempfile(); saveRDS(o, f); on.exit(unlink(f)); file.size(f) }
  expect_lt(sz(pm), sz(fm) / 3)
  # untransformed response: one column
  expect_null(attr(predictmeans(lm(yield ~ nitro, data = d), "nitro"), "plot_data")$.y_orig)
  # without it: data = TRUE asks for the data; a data frame still works
  attr(pm, "plot_data") <- NULL
  expect_error(ggplot2::autoplot(pm, data = TRUE), "autoplot\\(x, data = df\\)")
  expect_s3_class(ggplot2::autoplot(pm, data = pm_data(fm)), "ggplot")
})

test_that("plot() draws autoplot()", {
  skip_on_cran()
  pm <- predictmeans(lm(yield ~ nitro, data = oats_fac()), "nitro")
  pdf(NULL)
  on.exit(dev.off())
  expect_s3_class(plot(pm), "ggplot")
  pv <- predictmeans(lm(yield ~ nitro * Variety, data = oats_fac()), "nitro:Variety")
  expect_identical(plot(pv, xvar = "Variety", group = "nitro")$labels$x, "Variety")
})

test_that("autoplot(type = 'pvalue') shows one tile per comparison", {
  fm <- lm(yield ~ Block + nitro*Variety, data = oats_fac())
  pm <- predictmeans(fm, "nitro", adj = "BH", plot = FALSE)
  tiles <- ggplot2::layer_data(ggplot2::autoplot(pm, type = "pvalue"), 1)
  expect_equal(nrow(tiles), nrow(pm$pairs))
  pa <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "BH", plot = FALSE)
  pp <- ggplot2::autoplot(pa, type = "pvalue")
  expect_equal(nrow(ggplot2::layer_data(pp, 1)), 3 * 6)
  # one square panel per level of 'atvar', titled by the level
  expect_s3_class(pp$facet, "FacetWrap")
  expect_setequal(levels(pp$data$panel), levels(oats_fac()$Variety))
  expect_identical(pp$theme$aspect.ratio, 1)
  # every band is in the legend, in the order of the p-values
  fs <- ggplot2::ggplot_build(pp)$plot$scales$get_scales("fill")
  expect_identical(fs$get_limits(), pm_p_bands(0.05)$labels)
})

test_that("PMplot(): legend in p-value order with every band; named atvar panels; legend inside", {
  set.seed(2013)
  P <- matrix(0, 8, 8); P[lower.tri(P)] <- runif(28); P[2, 1] <- 0.07
  g <- PMplot(P)
  ld <- ggplot2::layer_data(g, 1)
  b <- pm_p_bands(0.05)
  band <- b$labels[findInterval(P[lower.tri(P)], b$breaks, left.open = TRUE) + 1L]
  # the colour of each tile is that of its band in the legend
  d <- g$data[order(g$data$col, g$data$row), ]
  expect_identical(unname(stats::setNames(b$colours, b$labels)[as.character(d$band)]), b$colours[match(band, b$labels)])
  expect_setequal(ld$fill, unique(b$colours[match(band, b$labels)]))
  expect_identical(g$theme$legend.justification, c(1, 1))
  # small triangles: legend outside, so it cannot cover tiles or leave the plot
  expect_identical(PMplot(P[1:4, 1:4])$theme$legend.position, "right")
  Oats$nitro <- factor(Oats$nitro)
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
  pa <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "BH", plot = FALSE)
  expect_named(pa$p_valueMatrix, levels(Oats$Variety))
  ga <- PMplot(pa$p_valueMatrix)
  expect_identical(levels(ga$data$panel), levels(Oats$Variety))
  expect_equal(nrow(ga$data), 3 * 6)
})

test_that("summary() describes the result", {
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats_fac())
  pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
  out <- capture.output(print(summary(pm)))
  expect_true(any(grepl("Predicted means for 'nitro:Variety'.*Kenward-Roger", out)))
  expect_true(any(grepl("18 pairwise comparisons within each level of Variety", out)))
  expect_true(any(grepl("^SED: min", out)))
  expect_true(any(grepl("^LSD \\(5%, tukey\\): [0-9.]+ \\(equal for all pairs\\)", out)))
  expect_true(any(grepl("^Comparison intervals: mean \\+/- ", out)))
  expect_true(any(grepl("adjusted by 'tukey'", out)))
  expect_true(any(grepl("Means sharing a letter", out)))
  perm <- suppressWarnings(permmodels(fm, method = "raw", nperm = 9, ncore = 1, seed = 1))
  pp <- predictmeans(fm, "nitro", permlist = perm, plot = FALSE, ncore = 1)
  expect_true(any(grepl("permutation \\(9 permutations\\)", capture.output(print(summary(pp))))))
})
