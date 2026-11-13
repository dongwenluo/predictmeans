## Denominator df by grouping level ("inner-outer", R/df-innerouter.R): equal
## to nlme's fixDF for nested random effects; the default for GLMMs.

io_oats <- function() {
  o <- as.data.frame(nlme::Oats)
  o$N <- factor(o$nitro)
  o$Variety <- factor(o$Variety, ordered = FALSE)
  o$Block <- factor(o$Block, ordered = FALSE)
  o
}

test_that("pm_inner_outer_df() equals nlme's fixDF (coefficients and terms) for nested random effects", {
  o <- io_oats()
  set.seed(3)
  d <- expand.grid(rep = 1:3, plot = 1:4, block = 1:4, site = 1:5)
  d$site <- factor(d$site); d$block <- factor(paste(d$site, d$block)); d$plot <- factor(paste(d$block, d$plot))
  d$S <- factor(as.integer(d$site) %% 2); d$xb <- rnorm(nlevels(d$block))[d$block]
  d$P <- factor(rep(rep(1:4, each = 3), 20)); d$xr <- rnorm(nrow(d)); d$y <- rnorm(nrow(d))
  cases <- list(
    list(y ~ Variety * N, "(1 | Block/Variety)", ~ 1 | Block/Variety, o),
    list(y ~ Variety * nitro, "(1 | Block/Variety)", ~ 1 | Block/Variety, o),
    list(y ~ Variety * N, "(1 | Block/Variety)", ~ 1 | Block/Variety, o[-c(1:4, 30, 55), ]),
    list(y ~ S + xb + P + xr + S:P, "(1 | site/block/plot)", ~ 1 | site/block/plot, d),
    list(y ~ S * xb + P * xr, "(1 | site/block/plot)", ~ 1 | site/block/plot, d[-c(5, 77, 130), ]))
  for (k in seq_along(cases)) {
    cs <- cases[[k]]
    dd <- cs[[4]]
    if (!"y" %in% names(dd)) dd$y <- dd$yield
    m <- suppressMessages(lme4::lmer(stats::as.formula(paste(deparse(cs[[1]]), "+", cs[[2]])), data = dd))
    e <- nlme::lme(cs[[1]], random = cs[[3]], data = dd)
    io <- pm_inner_outer_df(m)
    expect_equal(as.vector(io), as.vector(e$fixDF$X[names(io)]), info = paste("case", k))
    expect_equal(unname(attr(io, "terms")), unname(e$fixDF$terms[names(attr(io, "terms"))]), info = paste("case", k))
  }
})

io_clusters <- function() {
  set.seed(4)
  d <- expand.grid(unit = 1:3, W = factor(1:2), clus = factor(1:8))
  d$Trt <- factor(ifelse(as.integer(d$clus) <= 4, "a", "b"))
  d$r <- rbinom(nrow(d), 10, plogis(-0.5 + 0.5 * rnorm(8)[d$clus]))
  d
}

test_that("GLMMs: inner-outer by default; a treatment on whole clusters gets the cluster df", {
  d <- io_clusters()
  g <- lme4::glmer(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  p <- predictmeans(g, "Trt:W", plot = FALSE, prtnum = FALSE)
  expect_identical(p$settings$df_method, "inner-outer")
  # 8 clusters, 2 treatments: 6 df between clusters; within clusters 48 - 8 - 2 = 38
  between <- p$pairs$level1 == "a:1" & p$pairs$level2 == "b:1"
  within <- p$pairs$level1 == "a:1" & p$pairs$level2 == "a:2"
  expect_equal(p$pairs$df[between], 6)
  expect_equal(p$pairs$df[within], 38)
  expect_length(unique(p$means$df), 1)                    # the means of a term share one df
  expect_equal(pm_joint(g, "Trt")$df2, 6)
  expect_equal(unname(contrastmeans(g, "Trt", ctrmatrix = rbind(c(1, -1)))[[1]][, "df"]), 6)
  # residual df (predictmeans 1.x) on request
  r <- predictmeans(g, "Trt", df_method = "residual", plot = FALSE, prtnum = FALSE)
  expect_equal(r$means$df, rep(unname(summary(g)$AICtab["df.resid"]), 2))
  expect_identical(r$settings$df_method, "residual")
  # Kenward-Roger is not available for GLMMs
  expect_warning(k <- predictmeans(g, "Trt", df_method = "Kenward-Roger", plot = FALSE, prtnum = FALSE), "GLMMs")
  expect_identical(k$settings$df_method, "inner-outer")
  # containment and between-within work for glmer too
  expect_equal(unique(predictmeans(g, "Trt", df_method = "between-within", plot = FALSE, prtnum = FALSE)$pairs$df), 6)
})

test_that("glmmTMB gets the same inner-outer df as glmer", {
  skip_if_not_installed("glmmTMB")
  d <- io_clusters()
  g <- lme4::glmer(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  t <- glmmTMB::glmmTMB(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  expect_equal(unname(pm_inner_outer_df(t)), unname(pm_inner_outer_df(g)))
  expect_equal(predictmeans(t, "Trt:W", plot = FALSE, prtnum = FALSE)$pairs$df,
               predictmeans(g, "Trt:W", plot = FALSE, prtnum = FALSE)$pairs$df)
})

test_that("crossed random effects: a treatment on columns gets columns minus treatments", {
  set.seed(5)
  d <- expand.grid(W = factor(1:2), Row = factor(1:4), Col = factor(1:12))
  d$T <- factor((as.integer(d$Col) - 1) %% 3 + 1)
  d$y <- rpois(nrow(d), exp(1.5 + 0.3 * rnorm(4)[d$Row] + 0.3 * rnorm(12)[d$Col]))
  g <- lme4::glmer(y ~ T * W + (1 | Row) + (1 | Col), family = poisson, data = d)
  io <- pm_inner_outer_df(g)
  expect_equal(unname(io[c("T2", "T3")]), c(9, 9))                      # 12 columns - 3 treatments
  # W varies within every row and column: N - rank([X, rows, columns]) = 96 - (6 + 4 + 12 - 1 - 2) = 77
  expect_equal(unname(io["W2"]), 96 - qr(cbind(lme4::getME(g, "X"), model.matrix(~ Row - 1, d), model.matrix(~ Col - 1, d)))$rank)
})
