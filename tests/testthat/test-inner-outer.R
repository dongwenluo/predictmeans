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
  skip_on_cran()
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

io_with_options <- function(opts, code) {
  op <- options(opts)
  on.exit(options(op))
  code
}

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
  p <- predictmeans(g, "Trt:W", plot = FALSE)
  expect_identical(p$settings$df_method, "inner-outer")
  # 8 clusters, 2 treatments: 6 df between clusters; within clusters the trials count,
  # 480 - 8 - 2 = 470 (rows, 48 - 8 - 2 = 38, held 0-2.6% at 5%: dev/sim-df-io-trials.R)
  between <- p$pairs$level1 == "a:1" & p$pairs$level2 == "b:1"
  within <- p$pairs$level1 == "a:1" & p$pairs$level2 == "a:2"
  expect_equal(p$pairs$df[between], 6)
  expect_equal(p$pairs$df[within], 470)
  io_with_options(list(predictmeans.inner_outer_count = "rows"), {   # SAS GLIMMIX events/trials
    expect_equal(unname(pm_inner_outer_df(g)[c("(Intercept)", "Trtb", "W2")]), c(38, 6, 38))
  })
  expect_length(unique(p$means$df), 1)                    # the means of a term share one df
  expect_equal(pm_joint(g, "Trt")$df2, 6)
  expect_equal(unname(contrastmeans(g, "Trt", ctrmatrix = rbind(c(1, -1)))[[1]][, "df"]), 6)
  # residual df (predictmeans 1.x) on request
  r <- predictmeans(g, "Trt", df_method = "residual", plot = FALSE)
  expect_equal(r$means$df, rep(unname(summary(g)$AICtab["df.resid"]), 2))
  expect_identical(r$settings$df_method, "residual")
  # Kenward-Roger is not available for GLMMs
  expect_warning(k <- predictmeans(g, "Trt", df_method = "Kenward-Roger", plot = FALSE), "GLMMs")
  expect_identical(k$settings$df_method, "inner-outer")
  # containment and between-within work for glmer too
  expect_equal(unique(predictmeans(g, "Trt", df_method = "between-within", plot = FALSE)$pairs$df), 6)
})

test_that("glmmTMB gets the same inner-outer df as glmer", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  d <- io_clusters()
  g <- lme4::glmer(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  t <- glmmTMB::glmmTMB(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  expect_equal(unname(pm_inner_outer_df(t)), unname(pm_inner_outer_df(g)))
  expect_equal(predictmeans(t, "Trt:W", plot = FALSE)$pairs$df,
               predictmeans(g, "Trt:W", plot = FALSE)$pairs$df)
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

test_that("inner-outer df: aggregated binomial counts get the df of one row per trial", {
  skip_on_cran()
  d <- io_clusters()
  g <- lme4::glmer(cbind(r, 10 - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  gp <- lme4::glmer(r / 10 ~ Trt * W + (1 | clus), weights = rep(10, nrow(d)), family = binomial, data = d)
  b <- d[rep(seq_len(nrow(d)), each = 10), ]
  b$y <- unlist(lapply(d$r, function(r) rep(1:0, c(r, 10 - r))))
  gb <- lme4::glmer(y ~ Trt * W + (1 | clus), family = binomial, data = b)
  expect_equal(pm_inner_outer_df(g), pm_inner_outer_df(gb), ignore_attr = c("count", "trials_df"))
  expect_equal(pm_inner_outer_df(gp), pm_inner_outer_df(gb), ignore_attr = c("count", "trials_df"))
  # an observation-level random effect (overdispersion) makes the rows a stratum: effects within clusters get 38
  d$obs <- factor(seq_len(nrow(d)))
  go <- suppressMessages(lme4::glmer(cbind(r, 10 - r) ~ Trt * W + (1 | clus) + (1 | obs), family = binomial, data = d))
  expect_equal(attr(pm_inner_outer_df(go), "terms"), c(Trt = 6, W = 38, "Trt:W" = 38))
  expect_error(io_with_options(list(predictmeans.inner_outer_count = "trials"), pm_inner_outer_df(g)),
               "\"cases\" or \"rows\"")
})

test_that("inner-outer cases: classes, zero trials, weights that are not counts, settings, advice", {
  skip_on_cran()
  d <- io_clusters()
  d$n <- 10
  b <- d[rep(seq_len(nrow(d)), d$n), ]
  b$y <- unlist(lapply(seq_len(nrow(d)), function(i) rep(1:0, c(d$r[i], d$n[i] - d$r[i]))))
  gb <- lme4::glmer(y ~ Trt * W + (1 | clus), family = binomial, data = b)
  ref <- attr(pm_inner_outer_df(gb), "terms")
  # glmmTMB: cbind and proportion + weights
  skip_if_not_installed("glmmTMB")
  t1 <- glmmTMB::glmmTMB(cbind(r, n - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  t2 <- glmmTMB::glmmTMB(r / n ~ Trt * W + (1 | clus), weights = n, family = binomial, data = d)
  expect_equal(attr(pm_inner_outer_df(t1), "terms"), ref)
  expect_equal(attr(pm_inner_outer_df(t2), "terms"), ref)
  # responses that are not counts still count rows: lmer and Poisson glmer with prior weights
  set.seed(1)
  d$yn <- rnorm(nrow(d)); d$cnt <- rpois(nrow(d), 3); d$pw <- runif(nrow(d), 1, 3)
  lw <- lme4::lmer(yn ~ Trt * W + (1 | clus), weights = pw, data = d)
  pw <- suppressWarnings(lme4::glmer(cnt ~ Trt * W + (1 | clus), weights = round(pw), family = poisson, data = d))
  expect_equal(attr(pm_inner_outer_df(lw), "count"), "rows")
  expect_equal(attr(pm_inner_outer_df(pw), "count"), "rows")
  expect_equal(attr(pm_inner_outer_df(lw), "terms")[["W"]], 38)
  # a cluster with no trials: left out, as in the fit to one row per trial
  d0 <- d; d0$n[d0$clus == 1] <- 0; d0$r[d0$clus == 1] <- 0
  g0 <- suppressWarnings(lme4::glmer(cbind(r, n - r) ~ Trt * W + (1 | clus), family = binomial, data = d0))
  gb0 <- lme4::glmer(y ~ Trt * W + (1 | clus), family = binomial, data = b[b$clus != 1, ])
  expect_equal(attr(pm_inner_outer_df(g0), "terms"), attr(pm_inner_outer_df(gb0), "terms"))
  # binomial weights that are not whole numbers are not trials (sampling or precision weights): rows
  gh <- suppressWarnings(lme4::glmer(r / n ~ Trt * W + (1 | clus), weights = n + 0.5, family = binomial, data = d))
  expect_equal(attr(pm_inner_outer_df(gh), "terms")[["W"]], 38)
  expect_identical(attr(pm_inner_outer_df(gh), "count"), "rows")
  # the result says what was counted
  g <- lme4::glmer(cbind(r, n - r) ~ Trt * W + (1 | clus), family = binomial, data = d)
  expect_identical(as.vector(predictmeans(g, "W")$settings$df_cases), "trials")
  expect_identical(as.vector(io_with_options(list(predictmeans.inner_outer_count = "rows"), predictmeans(g, "W"))$settings$df_cases),
                   "rows")
  expect_true(is.na(predictmeans(lw, "W")$settings$df_cases))
  expect_false("overdispersion_trials" %in% advice(predictmeans(g, "W"))$id)
  # overdispersed rows: the advice says to add an observation-level random effect
  set.seed(7)
  dd <- expand.grid(trt = factor(1:4), blk = factor(1:6)); dd$m <- 25
  dd$y <- rbinom(nrow(dd), 25, plogis(rnorm(6)[dd$blk] + rnorm(nrow(dd), 0, 0.8)))
  go <- lme4::glmer(cbind(y, m - y) ~ trt + (1 | blk), family = binomial, data = dd)
  expect_true("overdispersion_trials" %in% advice(predictmeans(go, "trt"))$id)
  cm <- contrastmeans(go, "trt", ctrmatrix = rbind(c(1, -1, 0, 0)))
  expect_true("overdispersion_trials" %in% advice(cm)$id)
  expect_identical(as.vector(cm$settings$df_cases), "trials")
  # the same data with an observation-level random effect: effects within blocks get the rows' df; label rows
  dd$obs <- factor(seq_len(nrow(dd)))
  gr <- lme4::glmer(cbind(y, m - y) ~ trt + (1 | blk) + (1 | obs), family = binomial, data = dd)
  pr <- predictmeans(gr, "trt")
  expect_identical(as.vector(pr$settings$df_cases), "rows")
  expect_equal(unique(pr$means$df), 15)
  expect_false("overdispersion_trials" %in% advice(pr)$id)
  # glmmTMB with a proportion and weights: the same dispersion check as glmer
  tp <- glmmTMB::glmmTMB(I(y / m) ~ trt + (1 | blk), weights = m, family = binomial, data = dd)
  expect_identical(pm_check_dispersion(tp)$kind, "tested")
  # a treatment on whole clusters: its tests use the clusters' df, so no advice about trials
  set.seed(8)
  dc <- expand.grid(rep = 1:2, clus = factor(1:12)); dc$m <- 25
  dc$trt <- factor(as.integer(dc$clus) %% 2)
  dc$y <- rbinom(nrow(dc), 25, plogis(rnorm(12)[dc$clus] + rnorm(nrow(dc), 0, 0.8)))
  gc <- lme4::glmer(cbind(y, m - y) ~ trt + (1 | clus), family = binomial, data = dc)
  pc <- predictmeans(gc, "trt", pairwise = TRUE)
  expect_equal(pc$pairs$df, 10)
  expect_false("overdispersion_trials" %in% advice(pc)$id)
  # df given by the user: not inner-outer, so nothing counted
  expect_true(is.na(predictmeans(go, "trt", df_method = 5)$settings$df_cases))
  # binary data with sampling weights: no dispersion test
  set.seed(3)
  db <- data.frame(y = rbinom(120, 1, 0.4), g = factor(rep(1:12, 10)), w = runif(120, 1, 3))
  gw <- suppressWarnings(lme4::glmer(y ~ 1 + (1 | g), weights = w, family = binomial, data = db))
  expect_identical(pm_check_dispersion(gw)$kind, "binary")
})

test_that("inner-outer cases: clmm with case weights = one row per case", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  dat <- ordinal::wine
  agg <- stats::aggregate(list(w = rep(1, nrow(dat))), dat[c("rating", "temp", "contact", "judge")], sum)
  m1 <- ordinal::clmm(rating ~ temp * contact + (1 | judge), data = dat)
  m2 <- ordinal::clmm(rating ~ temp * contact + (1 | judge), weights = w, data = agg)
  expect_equal(pm_inner_outer_df(m2), pm_inner_outer_df(m1), ignore_attr = c("count", "trials_df"))
  expect_identical(attr(pm_inner_outer_df(m2), "count"), "case weights")
})
