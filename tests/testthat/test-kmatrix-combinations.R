## Kmatrix(): covariate values by name, any combination of fixed-effect
## variables, averaging by groups of interacting variables.

km_data <- function() {
  set.seed(42)
  d <- expand.grid(A = factor(c("a1", "a2", "a3")), B = factor(c("b1", "b2")), C = factor(1:4), Blk = factor(1:3))
  d$x <- rnorm(nrow(d)); d$z <- runif(nrow(d), 1, 3); d$y <- rnorm(nrow(d)) + d$x
  d
}

test_that("covariate values go by name; a wrong number of unnamed values is an error", {
  d <- km_data()
  m <- lm(y ~ A + B + x + z, d)
  f <- function(...) Kmatrix(m, "A", ...)$fctnames
  expect_equal(f(covariate = c(z = 2, x = 0.5))[1, c("x", "z")], data.frame(x = 0.5, z = 2), ignore_attr = TRUE)
  expect_equal(f(covariate = c(0.5, 2))[1, c("x", "z")], data.frame(x = 0.5, z = 2), ignore_attr = TRUE)
  expect_equal(f(covariate = c(x = 0.5))$z[1], mean(d$z))            # the others at their means
  expect_error(f(covariate = 0.5), "1 value\\(s\\) and the model has 2 covariate")
  expect_error(f(covariate = c(1, 2, 3)), "3 value")
  expect_error(f(covariate = c(w = 1)), "Unknown covariate")
  expect_error(f(covariate = "w"), "Unknown covariate")
  # a list: several values give rows, one value fixes the covariate
  k <- Kmatrix(m, "A", covariate = list(x = c(-1, 1), z = 2.5))
  expect_equal(nrow(k$K), 6)
  expect_equal(unique(k$fctnames$z), 2.5)
  expect_equal(k$K[, "x"], rep(c(-1, 1), 3), ignore_attr = TRUE)
  expect_equal(k$K, Kmatrix(m, "A", "x", at = c(-1, 1))$K + outer(rep(1, 6), (2.5 - mean(d$z)) * (colnames(k$K) == "z")),
               ignore_attr = TRUE)
  # covariateV as a named list: by name, not by position
  k1 <- Kmatrix(m, "A", c("x", "z"), at = list(z = c(1.5, 2.5), x = c(-1, 1)))
  k2 <- Kmatrix(m, "A", c("x", "z"), at = list(c(-1, 1), c(1.5, 2.5)))
  expect_equal(k1$K, k2$K)
  # joint points
  kj <- Kmatrix(m, "A", covariate = data.frame(x = c(-1, 0), z = c(1, 3)))
  expect_equal(nrow(kj$K), 6)
  expect_equal(unname(kj$K[1:2, c("x", "z")]), cbind(c(-1, 0), c(1, 3)))
})

test_that("any combination of fixed-effect variables gives rows, model term or not", {
  d <- km_data()
  m <- lm(y ~ A + B + C + x, d)
  k <- Kmatrix(m, "A:B")
  expect_equal(nrow(k$K), 6)
  nd <- transform(expand.grid(B = levels(d$B), A = levels(d$A))[2:1], x = mean(d$x))
  # the additive prediction averaged over C
  p <- sapply(seq_len(6), function(i) mean(predict(m, data.frame(A = nd$A[i], B = nd$B[i], x = nd$x[i], C = levels(d$C)))))
  expect_equal(as.vector(k$K %*% coef(m)), p, tolerance = 1e-10)
  # a covariate in the term: its mean, or the values given
  expect_equal(nrow(Kmatrix(m, "A:x")$K), 3)
  expect_equal(Kmatrix(m, "A:x")$K, Kmatrix(m, "A")$K, ignore_attr = TRUE)
  expect_equal(nrow(Kmatrix(m, "A:x", list(x = c(0, 1, 2)))$K), 9)
  expect_equal(nrow(Kmatrix(m, "x")$K), 1)
  expect_error(Kmatrix(m, "Q"), "'Q' is not in the fixed effects of the model. Their variables are: A, B, C, x")
  expect_error(Kmatrix(m, "a"), "Did you mean 'A'")
  # an aov model with Error() directly
  ao <- aov(y ~ A * B + Error(Blk), d)
  expect_equal(nrow(Kmatrix(ao, "A:B")$K), 6)
})

test_that("each group of interacting variables is averaged on its own grid: the same K as the full grid", {
  d <- km_data()
  mods <- list(lm(y ~ A + B + C + Blk + x + z, d), lm(y ~ A * B + C * Blk + x, d), lm(y ~ A * x + B + poly(z, 2), d),
               lme4::lmer(y ~ A * B + C + x + (1 | Blk), d),
               mgcv::gam(y ~ A * B + C + s(x) + s(Blk, bs = "re"), data = d, method = "REML"))
  for (m in mods) for (w in c("equal", "proportional", "outer", "cells", "flat")) {
    a <- Kmatrix(m, "A", weights = w)
    b <- Kmatrix(m, "A", weights = w, refgrid = TRUE)
    expect_null(a$refgrid)
    expect_s3_class(b$refgrid, "pm_refgrid")
    expect_equal(a$K, b$K, tolerance = 1e-10, info = paste(class(m)[1], w))
    expect_equal(a$fctnames, b$fctnames, tolerance = 1e-10)
  }
  expect_equal(Kmatrix(mods[[1]], "A", cov_reduce = "average")$K, Kmatrix(mods[[1]], "A", cov_reduce = "average", refgrid = TRUE)$K,
               tolerance = 1e-10)
  # many levels of factors that do not interact: no grid of all combinations
  set.seed(1)
  g <- expand.grid(Line = factor(sprintf("L%03d", 1:200)), Rep = factor(1:3)); g$Block <- factor(paste(g$Rep, sample(rep(1:20, 10))))
  g$Site <- factor(sample(1:4, nrow(g), TRUE)); g$y <- rnorm(nrow(g))
  mg <- lm(y ~ Line + Block + Site, g)                        # 200 x 60 x 4 = 48,000 combinations
  # (in these data the lines of a block are the same in every replicate, so Line and Block are
  # confounded: the means are over the combinations that exist, with a warning)
  tm <- system.time(k <- suppressWarnings(Kmatrix(mg, "Line")))[["elapsed"]]
  expect_lt(tm, 3)
  expect_equal(unname(k$K[2, "LineL002"]), 1)
  expect_equal(unname(k$K[1, grep("^Site", colnames(k$K))]), rep(0.25, 3))
})

test_that("predictmeans(): a combination that is not a term, its df, and combinations without data", {
  d <- km_data()
  d <- d[!(d$A == "a3" & d$B == "b2"), ]                      # no data for a3:b2
  m <- lm(y ~ A + B + C, d)
  expect_message(p <- predictmeans(m, "A:B"), "is not a term of the model")
  # the model is additive, so the mean of a3:b2 can be estimated and is shown
  expect_equal(nrow(p$means), 6)
  expect_equal(p$means$n[6], 0)
  expect_equal(p$means$estimate[6], mean(predict(m, data.frame(A = "a3", B = "b2", C = levels(d$C)))), tolerance = 1e-10)
  expect_false(anyNA(p[["Predicted Means"]]))
  # with the interaction in the model it is an empty cell and is left out
  mi <- lm(y ~ A * B + C, d)
  pi <- predictmeans(mi, "A:B")
  expect_equal(nrow(pi$means), 5)
  expect_equal(sum(is.na(pi[["Predicted Means"]])), 1)
  # df by term (lme): a combination gets the smallest df of its factors' terms
  oats <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  fl <- nlme::lme(yield ~ nitro + Variety, random = ~ 1 | Block/Variety, data = oats)
  q <- function(term) suppressMessages(predictmeans(fl, term, df_method = "containment", pairwise = TRUE))
  dv <- q("Variety")$pairs$df[1]; dn <- q("nitro")$pairs$df[1]
  expect_lt(dv, dn)
  expect_equal(unique(q("nitro:Variety")$pairs$df), dv)
  # covariates: set by name; a covariate in the term is not tabulated by predictmeans()
  mx <- lm(y ~ A + x + z, d)
  expect_equal(predictmeans(mx, "A", covariate = list(z = 2, x = 1))$means$estimate,
               unname(predict(mx, data.frame(A = levels(d$A), x = 1, z = 2))), tolerance = 1e-10)
  expect_error(predictmeans(mx, "A:x"), "is a covariate")
  expect_error(predictmeans(mx, "A", covariate = list(x = c(0, 1))), "one mean for each combination")
})

test_that("nesting the formula does not say (confounded factors): the combinations that exist, with a warning", {
  set.seed(1)
  dn <- data.frame(G = factor(rep(c("ctl", "trt"), c(4, 12))), T = factor(rep(c("c", "t1", "t2", "t3"), each = 4)))
  dn$y <- rnorm(16) + as.numeric(dn$T)
  cell <- tapply(dn$y, dn$T, mean)
  # declared (G/T): the means of the combinations that exist
  pg <- suppressMessages(predictmeans(lm(y ~ G/T, dn), "G"))
  expect_equal(pg$means$estimate, c(cell[["c"]], mean(cell[c("t1", "t2", "t3")])), tolerance = 1e-10)
  # not declared (G + T): G and T cannot be told apart, so the mean over all their combinations does
  # not exist; the means are those of the combinations that exist (as G/T), and a warning says so
  # (decision of 2026-11-16; NA before)
  expect_warning(pa <- predictmeans(lm(y ~ G + T, dn), "T"),
                 "'G' and 'T' are confounded in the data.*G/T in the formula")
  expect_equal(pa$means$estimate, as.vector(cell), tolerance = 1e-10)
  pa2 <- suppressWarnings(predictmeans(lm(y ~ G + T, dn), "G"))
  expect_equal(pa2$means$estimate, pg$means$estimate, tolerance = 1e-10)
  expect_equal(pa2$means$std.error, pg$means$std.error, tolerance = 1e-10)
  pc <- suppressMessages(predictmeans(lm(y ~ G + T, dn), "T", weights = "cells"))
  expect_equal(pc$means$estimate, as.vector(cell), tolerance = 1e-10)
  # declared nesting also with weights other than "equal"
  kp <- Kmatrix(lm(y ~ G/T, dn), "G", weights = "proportional")
  expect_true(all(attr(kp$K, "estimable")))
})

test_that("cells are those of each interaction term; one observation per cell is no exception", {
  set.seed(3)
  d <- expand.grid(A = factor(1:3), B = factor(1:3)); d$y <- rnorm(9)
  d1 <- d[-9, ]                                               # unreplicated, cell 3:3 empty
  expect_warning(p1 <- predictmeans(lm(y ~ A * B, d1), "A"), "No data for A:B 3:3")
  expect_equal(p1$means$estimate[3], mean(d1$y[d1$A == "3"]), tolerance = 1e-10)
  # two separate interactions: a combination of all four factors without data is not an empty cell
  d4 <- expand.grid(A = factor(1:2), B = factor(1:2), C = factor(1:2), D = factor(1:2), r = 1:2); d4$y <- rnorm(32)
  d4 <- d4[!(d4$A == "2" & d4$B == "2" & d4$C == "2" & d4$D == "2"), ]
  expect_no_warning(k <- Kmatrix(lm(y ~ A * B + C * D, d4), "A"))
  expect_null(k$empty_cells)
  expect_true(all(attr(k$K, "estimable")))
  # a factor that only interacts with a covariate makes no cells
  d5 <- expand.grid(A = factor(1:2), B = factor(1:3), r = 1:3); d5$x <- rnorm(18); d5$y <- rnorm(18)
  d5 <- d5[!(d5$A == "2" & d5$B == "3"), ]
  expect_no_warning(k5 <- Kmatrix(lm(y ~ A * x + B, d5), "A"))
  expect_equal(unname(k5$K[2, c("B2", "B3")]), c(1/3, 1/3))
})

test_that("a new class needs only coef(), vcov() and terms(): contrasts come from pm_contrasts()", {
  d <- km_data()
  fit <- structure(list(fit = lm(y ~ A * B + x, d, contrasts = list(A = "contr.sum")), call = quote(fit_my(y ~ A * B + x, d))),
                   class = "mymodel_k")
  registerS3method("coef", "mymodel_k", function(object, ...) coef(object$fit))
  registerS3method("vcov", "mymodel_k", function(object, ...) vcov(object$fit))
  registerS3method("terms", "mymodel_k", function(x, ...) terms(x$fit))
  registerS3method("pm_contrasts", "mymodel_k", function(model, ...) model$fit$contrasts,
                   envir = asNamespace("predictmeans"))
  registerS3method("pm_data", "mymodel_k", function(model, ...) d, envir = asNamespace("predictmeans"))
  expect_equal(Kmatrix(fit, "A")$K, Kmatrix(fit$fit, "A")$K)
  expect_identical(pm_contrasts(fit$fit), fit$fit$contrasts)
})
