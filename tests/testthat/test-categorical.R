## Class probabilities of multinomial and non-parallel fits (R/categorical.R,
## dev/review-vglm.md B3b): nnet::multinom, VGAM::vglm(multinomial()) and
## non-parallel ordinal vglm fits. References: emmeans on multinom (mode =
## "prob"; regrid(transform = "logit") for the logit pairs), the same model by
## vglm, the binomial glm for two classes, predict() and a numerical delta
## method for vglm.

housing_fits <- function() {
  data(housing, package = "MASS", envir = environment())
  m <- nnet::multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = housing, trace = FALSE, reltol = 1e-12)
  list(d = housing, m = m)
}
## our pair estimate is level1 - level2 (as emmeans' "A - B"): a key for both
pair_key <- function(a, b, by) paste(sub(":.*", "", a), sub(":.*", "", b), by)

test_that("multinom class probabilities, logit and probability pairs equal emmeans", {
  skip_on_cran()
  skip_if_not_installed("nnet")
  skip_if_not_installed("emmeans")
  h <- housing_fits()
  r <- predictmeans(h$m, "Infl")
  e <- emmeans::emmeans(h$m, ~ Infl | Sat, mode = "prob")
  se <- as.data.frame(summary(e))
  k <- match(paste(r$means$Infl, r$means$Sat), paste(se$Infl, se$Sat))
  expect_false(anyNA(k))
  expect_equal(r$means$response, se$prob[k], tolerance = 1e-8)
  expect_equal(r$means$response.std.error, se$SE[k], tolerance = 1e-8)
  expect_equal(r$means$estimate, qlogis(se$prob[k]), tolerance = 1e-8)
  # limits on the logit scale: inside (0, 1)
  expect_true(all(r$means$response.conf.low > 0 & r$means$response.conf.high < 1))
  # pairs within each class only: log odds ratios of the class against the others
  expect_true(all(sub(".*:", "", r$pairs$level1) == sub(".*:", "", r$pairs$level2)))
  expect_equal(nrow(r$pairs), 9L)
  pl <- as.data.frame(summary(pairs(emmeans::regrid(e, transform = "logit"))))
  kp <- match(pair_key(r$pairs$level1, r$pairs$level2, r$pairs$by),
              paste(sub(" - .*", "", pl$contrast), sub(".* - ", "", pl$contrast), pl$Sat))
  expect_false(anyNA(kp))
  expect_equal(r$pairs$estimate, pl$estimate[kp], tolerance = 1e-8)
  expect_equal(r$pairs$std.error, pl$SE[kp], tolerance = 1e-8)
  expect_equal(r$pairs$odds.ratio, exp(pl$estimate[kp]), tolerance = 1e-8)
  expect_equal(r$settings$df_method, "asymptotic (z)")
  expect_equal(r$settings$model_class, "multinom")
  expect_true("class_prob_average" %in% advice(r)$id)
  # mean_type = "prob_diff": differences of probabilities, as emmeans' pairs
  r2 <- predictmeans(h$m, "Infl", mean_type = "prob_diff")
  expect_equal(r2$means$estimate, se$prob[k], tolerance = 1e-8)
  pp <- as.data.frame(summary(pairs(e)))
  kq <- match(pair_key(r2$pairs$level1, r2$pairs$level2, r2$pairs$by),
              paste(sub(" - .*", "", pp$contrast), sub(".* - ", "", pp$contrast), pp$Sat))
  expect_equal(r2$pairs$estimate, pp$estimate[kq], tolerance = 1e-8)
  expect_equal(r2$pairs$std.error, pp$SE[kq], tolerance = 1e-8)
  # letters per class, Tukey within each class (3 pairs per family)
  r3 <- predictmeans(h$m, "Infl", adj = "tukey")
  expect_false(anyNA(r3$means$group))
  expect_equal(as.vector(table(r3$pairs$by)[c("Low", "Medium", "High")]), rep(3L, 3))
})

test_that("multinom linear predictors (mean_type = 'link') are K B, with the reference class in the advice", {
  skip_if_not_installed("nnet")
  h <- housing_fits()
  r <- predictmeans(h$m, "Infl", mean_type = "link")
  # equal-weight average of the design over Type and Cont at each level of Infl
  g <- expand.grid(Infl = levels(h$d$Infl), Type = levels(h$d$Type), Cont = levels(h$d$Cont))
  X <- model.matrix(~ Infl + Type + Cont, g)
  K <- rowsum(X, g$Infl)[levels(h$d$Infl), ] / 8
  B <- coef(h$m)
  ref <- c(K %*% t(B))
  V <- vcov(h$m)
  KB <- kronecker(diag(2), K)                         # predictor-major, as vcov()'s order
  expect_equal(r$means$estimate, ref, tolerance = 1e-10)
  expect_equal(r$means$std.error, sqrt(diag(KB %*% V %*% t(KB))), tolerance = 1e-10)
  expect_equal(levels(r$means$predictor), c("log(P[Medium]/P[Low])", "log(P[High]/P[Low])"))
  expect_match(advice(r)$text[advice(r)$id == "class_reference"], "reference class 'Low'")
})

test_that("vglm multinomial equals multinom; two classes equal the binomial glm", {
  skip_on_cran()
  skip_if_not_installed("nnet")
  skip_if_not_installed("VGAM")
  h <- housing_fits()
  housing <- h$d
  v <- suppressWarnings(VGAM::vglm(Sat ~ Infl + Type + Cont, VGAM::multinomial, weights = Freq, data = housing))
  rm <- suppressMessages(predictmeans(h$m, "Infl:Cont"))
  rv <- suppressMessages(predictmeans(v, "Infl:Cont"))
  expect_equal(rv$means$response, rm$means$response, tolerance = 1e-5)
  expect_equal(rv$means$response.std.error, rm$means$response.std.error, tolerance = 1e-4)
  expect_equal(rv$pairs$estimate, rm$pairs$estimate, tolerance = 1e-5)
  rl <- predictmeans(v, "Infl", mean_type = "link")
  expect_equal(levels(rl$means$predictor), c("log(P[Low]/P[High])", "log(P[Medium]/P[High])"))   # VGAM: the last class
  expect_match(advice(rl)$text[advice(rl)$id == "class_reference"], "reference class 'High'")
  # two classes: the class "TRUE" is the binomial glm's mean (cells of the model: no averaging)
  housing$hi <- factor(housing$Sat == "High")
  mb <- nnet::multinom(hi ~ Infl + Type + Cont, weights = Freq, data = housing, trace = FALSE, reltol = 1e-14)
  gb <- glm(hi ~ Infl + Type + Cont, weights = Freq, data = housing, family = binomial)
  rb <- suppressMessages(predictmeans(mb, "Infl:Type:Cont"))
  rg <- suppressMessages(predictmeans(gb, "Infl:Type:Cont"))
  t2 <- rb$means[rb$means$hi == "TRUE", ]
  kk <- match(paste(t2$Infl, t2$Type, t2$Cont), paste(rg$means$Infl, rg$means$Type, rg$means$Cont))
  expect_equal(t2$estimate, rg$means$estimate[kk], tolerance = 1e-5)
  expect_equal(t2$std.error, rg$means$std.error[kk], tolerance = 1e-5)
})

test_that("non-parallel cumulative vglm: probabilities, mean class and predictors against predict() and a numerical delta method", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  data(housing, package = "MASS", envir = environment())
  v <- VGAM::vglm(Sat ~ Infl + Cont, VGAM::cumulative(parallel = FALSE ~ Infl), weights = Freq, data = housing)
  expect_true(predictmeans:::pm_cat_is(v))
  nd <- expand.grid(Infl = levels(housing$Infl), Cont = levels(housing$Cont))
  key <- paste(nd$Infl, nd$Cont)
  probs <- function(th) { v2 <- v; v2@coefficients <- setNames(th, names(VGAM::coef(v))); VGAM::predict(v2, nd, type = "response") }
  th <- VGAM::coef(v)
  V <- VGAM::vcov(v)
  P <- probs(th)
  G <- numDeriv::jacobian(function(t) c(probs(t)), th)          # rows: cell within class
  r <- predictmeans(v, "Infl:Cont", mean_type = "prob_diff")
  k <- match(paste(r$means$Infl, r$means$Cont), key) + (as.integer(r$means$Sat) - 1) * nrow(nd)
  expect_equal(r$means$estimate, c(P)[k], tolerance = 1e-8)
  expect_equal(r$means$std.error, sqrt(diag(G %*% V %*% t(G)))[k], tolerance = 1e-6)
  mc <- predictmeans(v, "Infl:Cont", mean_type = "mean_class")
  s <- kronecker(t(1:3), diag(nrow(nd)))                         # sum_j j P_j for each cell
  km <- match(paste(mc$means$Infl, mc$means$Cont), key)
  expect_equal(mc$means$estimate, c(P %*% 1:3)[km], tolerance = 1e-8)
  expect_equal(mc$means$std.error, sqrt(diag(s %*% G %*% V %*% t(G) %*% t(s)))[km], tolerance = 1e-6)
  lk <- predictmeans(v, "Infl:Cont", mean_type = "link")
  pr <- VGAM::predict(v, nd, se.fit = TRUE)
  kl <- match(paste(lk$means$Infl, lk$means$Cont), key) + (as.integer(lk$means$predictor) - 1) * nrow(nd)
  expect_equal(lk$means$estimate, c(pr$fitted.values)[kl], tolerance = 1e-8)
  expect_equal(lk$means$std.error, c(pr$se.fit)[kl], tolerance = 1e-8)
})

test_that("separation and probabilities near 0 or 1 are flagged", {
  skip_on_cran()
  skip_if_not_installed("nnet")
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  m <- nnet::multinom(rating ~ temp + contact, data = wine, trace = FALSE)   # no warm 1s, no cold 5s
  r <- predictmeans(m, "temp")
  a <- advice(r)
  expect_true("class_separation" %in% a$id)
  expect_match(a$text[a$id == "class_separation"], "warm:1.*cold:5|cold:5.*warm:1")
  expect_true("class_prob_near_bound" %in% advice(predictmeans(m, "temp", mean_type = "prob_diff"))$id)
})

test_that("multinomial and non-parallel fits stop clearly where they are not handled", {
  skip_on_cran()
  skip_if_not_installed("nnet")
  data(housing, package = "MASS", envir = environment())
  m <- nnet::multinom(Sat ~ Infl + Cont, weights = Freq, data = housing, trace = FALSE)
  expect_error(Kmatrix(m, "Infl"), "vector of means")
  expect_error(contrastmeans(m, "Infl", ctrmatrix = rbind(c(1, -1, 0))), "vector of means")
  expect_error(covariatemeans(m, "Infl", covariate = "Cont"), "vector of means")
  expect_error(residplot(m), "class response")
  expect_message(expect_null(pm_plot_data(predictmeans(m, "Infl"), TRUE, "Infl", "response")), "observations are classes")
  md <- nnet::multinom(Sat ~ Infl + Cont, weights = Freq, data = housing, trace = FALSE, decay = 0.1)
  expect_error(predictmeans(md, "Infl"), "penalised")
  expect_error(predictmeans(m, "Infl", mean_type = "conditional"), "does not apply")
  expect_error(predictmeans(m, "Infl", mean_type = "mean_class"), "unordered")
  expect_error(predictmeans(m, "Infl", trans = exp), "'trans' cannot be used")
  expect_warning(predictmeans(m, "Infl", df_method = "residual"), "does not apply")
  expect_error(predictmeans(lm(Freq ~ Infl, housing), "Infl", mean_type = "link"), "multinomial and non-parallel")
  # a user's atvar: pairs within (class, level), never across classes
  r <- suppressMessages(predictmeans(m, "Infl:Cont", atvar = "Cont"))
  expect_true(all(sub(".*:", "", r$pairs$level1) == sub(".*:", "", r$pairs$level2)))
  expect_setequal(unique(na.omit(r$pairs$by)), c(outer(levels(housing$Sat), levels(housing$Cont), paste, sep = " : ")))
})

test_that("crossing cumulative curves of a non-parallel fit stop, not negative probabilities", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  set.seed(1)
  d <- data.frame(x = runif(300, -1, 1), f = factor(sample(c("a", "b"), 300, TRUE)))
  u <- runif(300)
  d$y <- ordered(ifelse(u < plogis(-1 + 3 * d$x), 1, ifelse(u < pmax(plogis(-1 + 3 * d$x), plogis(1 - 3 * d$x)), 2, 3)))
  v <- suppressWarnings(VGAM::vglm(y ~ x + f, VGAM::cumulative(parallel = FALSE), data = d))
  expect_error(predictmeans(v, "f", covariate = list(x = 1), mean_type = "prob_diff"),
               "class '2' has probability -0.3.* at .*x = 1")
  expect_error(predictmeans(v, "f", covariate = list(x = 1)), "cross")
  # only some means cross (an additive 2 x 2 fit extrapolated to the cell without data): those are NA
  pr <- list(c(0.269, 0.354, 0.377), c(0.5, 0.122, 0.378), c(0.5, 0.122, 0.378))
  cells <- data.frame(A = factor(c(0, 1, 0)), B = factor(c(0, 0, 1)))
  d2 <- do.call(rbind, lapply(1:3, function(i) data.frame(cells[rep(i, 3), ], y = factor(1:3), n = round(1000 * pr[[i]]))))
  d2$y <- ordered(d2$y)
  v2 <- VGAM::vglm(y ~ A + B, VGAM::cumulative(parallel = FALSE), weights = n, data = d2)
  r2 <- predictmeans(v2, "A", pairwise = TRUE)
  expect_equal(r2$means$estimable, rep(c(TRUE, FALSE), 3))
  expect_true(all(is.na(r2$means$response[!r2$means$estimable])))
  expect_true(all(is.na(r2$pairs$estimate)))
  a <- advice(r2)
  expect_true("class_crossing" %in% a$id)
  expect_false("nonestimable_means" %in% a$id)
  expect_s3_class(predictmeans(v, "f", covariate = list(x = 0)), "predictmeans")
})

test_that("review of B3b: aliasing stops, settings, plots, advice and other families", {
  skip_on_cran()
  skip_if_not_installed("nnet")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("emmeans")
  # an aliased design (T nested in G): nnet keeps arbitrary coefficients -> a stop, not arbitrary means
  set.seed(2)
  d <- data.frame(G = factor(rep(1:3, each = 40)), B = factor(rep(1:2, 60)))
  d$T <- factor(as.integer(d$G) * 2 - rep(0:1, 60))
  d$y <- factor(sample(c("a", "b", "c"), 120, TRUE))
  ma <- nnet::multinom(y ~ G + T + B, data = d, trace = FALSE)
  expect_error(predictmeans(ma, "G"), "rank deficient \\(aliased")
  data(housing, package = "MASS", envir = environment())
  m <- nnet::multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = housing, trace = FALSE)
  # settings: the probabilities are averaged on the response scale; the fit's family
  r <- predictmeans(m, "Infl")
  expect_equal(r$settings$average, "response")
  expect_equal(r$settings$response$family, "multinomial")
  expect_equal(predictmeans(m, "Infl", average = "response")$means, r$means)
  # SEDs and LSDs only within a class
  expect_equal(colnames(attr(r[["Standard Error of Differences"]], "For the Same Level of Factor")), "Sat")
  # no "not a term: the differences are the same" message for probabilities
  expect_silent(suppressWarnings(predictmeans(m, "Infl:Cont")))
  expect_message(predictmeans(m, "Infl:Cont", mean_type = "link"), "not a term")
  # p-value plot
  expect_s3_class(ggplot2::autoplot(predictmeans(m, "Infl", pairwise = TRUE), type = "pvalue"), "ggplot")
  # advice: factors averaged and the weights; a covariate held at a value is not averaged over
  a <- advice(predictmeans(m, "Infl", weights = "proportional"))
  expect_match(a$text[a$id == "class_prob_average"], "over Type, Cont .*weights = \"proportional\"")
  housing$x <- as.numeric(housing$Type)
  mx <- nnet::multinom(Sat ~ Infl + x, weights = Freq, data = housing, trace = FALSE)
  expect_false("class_prob_average" %in% advice(predictmeans(mx, "Infl", covariate = list(x = 2)))$id)
  # a covariate and proportional weights: = emmeans
  rx <- predictmeans(mx, "Infl", covariate = list(x = 2), mean_type = "prob_diff")
  ex <- as.data.frame(summary(emmeans::emmeans(mx, ~ Infl | Sat, at = list(x = 2), mode = "prob")))
  k <- match(paste(rx$means$Infl, rx$means$Sat), paste(ex$Infl, ex$Sat))
  expect_equal(rx$means$estimate, ex$prob[k], tolerance = 1e-8)
  expect_equal(rx$means$std.error, ex$SE[k], tolerance = 1e-8)
  # proportional weights on individual rows (with frequency weights emmeans weights by them, predictmeans by rows)
  hi <- housing[rep(seq_len(nrow(housing)), housing$Freq), ]
  hi <- hi[-(1:200), ]                                         # unequal margins
  mi <- nnet::multinom(Sat ~ Infl + Type + Cont, data = hi, trace = FALSE)
  rp <- predictmeans(mi, "Infl", mean_type = "prob_diff", weights = "proportional")
  ep <- as.data.frame(summary(emmeans::emmeans(mi, ~ Infl | Sat, mode = "prob", weights = "proportional")))
  k <- match(paste(rp$means$Infl, rp$means$Sat), paste(ep$Infl, ep$Sat))
  expect_equal(rp$means$estimate, ep$prob[k], tolerance = 1e-8)
  # a user atvar with Tukey: families of the levels within each (class, level)
  rt <- predictmeans(m, "Infl:Cont", atvar = "Cont", adj = "tukey")
  expect_equal(as.vector(table(rt$pairs$by[!is.na(rt$pairs$by)])), rep(3L, 6))
  # vglm multinomial with refLevel = 2, partly parallel; cratio and sratio non-parallel: = predict()
  nd <- expand.grid(Infl = levels(housing$Infl), Cont = levels(housing$Cont))
  key <- paste(nd$Infl, nd$Cont)
  fams <- list(suppressWarnings(VGAM::vglm(Sat ~ Infl + Cont, VGAM::multinomial(refLevel = 2, parallel = TRUE ~ Cont),
                                            weights = Freq, data = housing)),
               VGAM::vglm(Sat ~ Infl + Cont, VGAM::cratio(parallel = FALSE ~ Infl), weights = Freq, data = housing),
               VGAM::vglm(Sat ~ Infl + Cont, VGAM::sratio(parallel = FALSE ~ Cont), weights = Freq, data = housing))
  for (v in fams) {
    P <- VGAM::predict(v, nd, type = "response")
    rv <- predictmeans(v, "Infl:Cont", mean_type = "prob_diff")
    kv <- match(paste(rv$means$Infl, rv$means$Cont), key) + (as.integer(rv$means$Sat) - 1) * nrow(nd)
    expect_equal(rv$means$estimate, c(P)[kv], tolerance = 1e-8)
    th <- VGAM::coef(v)
    G <- numDeriv::jacobian(function(t) { v2 <- v; v2@coefficients <- setNames(t, names(th))
                                          c(VGAM::predict(v2, nd, type = "response")) }, th)
    expect_equal(rv$means$std.error, sqrt(diag(G %*% VGAM::vcov(v) %*% t(G)))[kv], tolerance = 1e-6)
  }
})

test_that("proportional weights sum frequency weights: aggregated data = one row per case = emmeans", {
  skip_on_cran()
  skip_if_not_installed("nnet")
  skip_if_not_installed("emmeans")
  data(housing, package = "MASS", envir = environment())
  h <- housing[!(housing$Type == "Tower" & housing$Cont == "Low"), ]          # unequal margins
  hi <- h[rep(seq_len(nrow(h)), h$Freq), ]
  m <- nnet::multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = h, trace = FALSE, reltol = 1e-12)
  mi <- nnet::multinom(Sat ~ Infl + Type + Cont, data = hi, trace = FALSE, reltol = 1e-12)
  r <- predictmeans(m, "Infl", mean_type = "prob_diff", weights = "proportional")
  ri <- predictmeans(mi, "Infl", mean_type = "prob_diff", weights = "proportional")
  expect_equal(r$means$estimate, ri$means$estimate, tolerance = 1e-6)
  expect_equal(r$settings$weights_by, "frequency weights")
  expect_equal(ri$settings$weights_by, "rows")
  e <- as.data.frame(emmeans::emmeans(m, ~ Infl | Sat, mode = "prob", weights = "proportional"))
  k <- match(paste(r$means$Infl, r$means$Sat), paste(e$Infl, e$Sat))
  expect_equal(r$means$estimate, e$prob[k], tolerance = 1e-8)
  h$hi <- as.numeric(h$Sat == "High")
  hi$hi <- as.numeric(hi$Sat == "High")
  g <- glm(hi ~ Infl + Type + Cont, weights = Freq, data = h, family = binomial)
  gi <- glm(hi ~ Infl + Type + Cont, data = hi, family = binomial)
  expect_equal(predictmeans(g, "Infl", weights = "proportional")$means$estimate,
               predictmeans(gi, "Infl", weights = "proportional")$means$estimate, tolerance = 1e-6)
  # an lm's prior weights are precisions: rows are counted
  expect_null(predictmeans:::pm_freq_weights(lm(Freq ~ Infl, weights = rep(1:2, length.out = nrow(h)), data = h),
                                             nrow(h)))
})

test_that("case counts follow the response, not the class (review of the weights rule)", {
  skip_on_cran()
  skip_if_not_installed("emmeans")
  skip_if_not_installed("glmmTMB")
  set.seed(5)
  d <- expand.grid(A = factor(1:3), B = factor(1:4), rep = 1:2)
  d$n <- sample(2:40, nrow(d), TRUE)
  d$s <- rbinom(nrow(d), d$n, 0.4)
  d$id <- factor(seq_len(nrow(d)) %% 6)
  pr <- function(m) predictmeans(m, "A", weights = "proportional")$means$estimate
  gb <- glm(cbind(s, n - s) ~ A + B, data = d, family = binomial)
  gq <- glm(cbind(s, n - s) ~ A + B, data = d, family = quasibinomial)
  expect_equal(pr(gq), pr(gb), tolerance = 1e-10)
  e <- as.data.frame(emmeans::emmeans(gb, ~ A, weights = "proportional"))
  expect_equal(pr(gb), e$emmean, tolerance = 1e-8)
  mr <- lme4::glmer(cbind(s, n - s) ~ A + B + (1 | id), data = d, family = binomial)
  mt <- glmmTMB::glmmTMB(cbind(s, n - s) ~ A + B + (1 | id), data = d, family = binomial)
  kr <- Kmatrix(mr, "A", weights = "proportional")$K
  kt <- Kmatrix(mt, "A", weights = "proportional")$K
  expect_equal(c(kr), c(kt[, colnames(kr)]), tolerance = 1e-12)
  expect_equal(c(Kmatrix(gb, "A", weights = "proportional")$K), c(kr), tolerance = 1e-12)
  # glmer with na.exclude: unpadded weights; polr without its model frame
  d2 <- d
  d2$s[3] <- NA
  mx <- lme4::glmer(cbind(s, n - s) ~ A + B + (1 | id), data = d2, family = binomial, na.action = na.exclude)
  expect_silent(kx <- Kmatrix(mx, "A", weights = "proportional")$K)
  g2 <- glm(cbind(s, n - s) ~ A + B, data = d2, family = binomial)
  expect_equal(c(kx), c(Kmatrix(g2, "A", weights = "proportional")$K), tolerance = 1e-12)
  data(housing, package = "MASS", envir = environment())
  h <- housing[!(housing$Type == "Tower" & housing$Cont == "Low"), ]
  p0 <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = h, Hess = TRUE, model = FALSE)
  p1 <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = h, Hess = TRUE)
  expect_equal(predictmeans(p0, "Infl", weights = "proportional")$means$estimate,
               predictmeans(p1, "Infl", weights = "proportional")$means$estimate, tolerance = 1e-10)
  expect_equal(predictmeans(p0, "Infl", weights = "proportional")$settings$weights_by, "frequency weights")
})

test_that("class means: a response named y, summary with NA means, semireg weights", {
  skip_on_cran()
  skip_if_not_installed("VGAM")
  pr <- list(c(0.269, 0.354, 0.377), c(0.5, 0.122, 0.378), c(0.5, 0.122, 0.378))
  cells <- data.frame(A = factor(c(0, 1, 0)), B = factor(c(0, 0, 1)))
  d2 <- do.call(rbind, lapply(1:3, function(i) data.frame(cells[rep(i, 3), ], y = factor(1:3), n = round(1000 * pr[[i]]))))
  d2$y <- ordered(d2$y)
  v2 <- VGAM::vglm(y ~ A + B, VGAM::cumulative(parallel = FALSE), weights = n, data = d2)
  r <- predictmeans(v2, "A", pairwise = TRUE)
  expect_true("y_class" %in% names(r$means))                       # not the plot's column y
  p <- ggplot2::autoplot(r)
  expect_equal(length(unique(ggplot2::ggplot_build(p)$layout$layout$PANEL)), 3L)
  out <- capture.output(print(summary(r)))
  expect_true(any(grepl("3 pairwise comparisons within each level of y_class \\(3 not estimable\\)", out)))
  expect_false(any(grepl("SED: min NA", out)))
})

## ---- clm with nominal or scale effects (B4 G1) ----

clm_check <- function(fit, term) {
  m <- predictmeans(fit, term)$means
  vars <- strsplit(term, ":")[[1]]
  nd <- m[m$rating == levels(m$rating)[1], vars, drop = FALSE]
  pr <- predict(fit, newdata = nd, type = "prob", se.fit = TRUE)
  ij <- cbind(match(do.call(paste, m[vars]), do.call(paste, nd)), as.integer(m$rating))
  expect_equal(m$response, unname(pr$fit[ij]), tolerance = 1e-8)
  expect_equal(m$std.error * m$response * (1 - m$response), unname(pr$se.fit[ij]), tolerance = 1e-8)  # logit SE
}

test_that("clm nominal and scale effects: class probabilities and SEs = predict.clm()", {
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  clm_check(ordinal::clm(rating ~ temp, nominal = ~ contact, data = wine), "temp:contact")
  skip_on_cran()
  clm_check(ordinal::clm(rating ~ temp + contact, scale = ~ temp, data = wine), "temp:contact")
  clm_check(ordinal::clm(rating ~ temp, nominal = ~ contact, scale = ~ temp, data = wine, link = "probit"),
            "temp:contact")
  clm_check(ordinal::clm(rating ~ temp, nominal = ~ contact, data = wine, threshold = "equidistant"), "temp:contact")
})

test_that("clm nominal = vglm(parallel = FALSE ~ z); averages over scale variables = delta method", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  skip_if_not_installed("VGAM")
  data(wine, package = "ordinal", envir = environment())
  a <- ordinal::clm(rating ~ temp, nominal = ~ contact, data = wine)
  v <- VGAM::vglm(rating ~ temp + contact, VGAM::cumulative(parallel = FALSE ~ contact), data = wine)
  # same model (logLik), so the same means; the SEs differ by up to 2% (vglm: expected information)
  expect_equal(predictmeans(a, "temp")$means$estimate, predictmeans(v, "temp")$means$estimate, tolerance = 1e-5)
  la <- predictmeans(a, "temp:contact", mean_type = "link")$means
  expect_equal(la$estimate, predictmeans(v, "temp:contact", mean_type = "link")$means$estimate, tolerance = 1e-5)
  expect_identical(levels(la$predictor)[1], "logit(P(Y <= 1))")
  # scale over the judges: probabilities averaged over judge and contact, numDeriv delta method on coef()
  s <- ordinal::clm(rating ~ temp + contact, scale = ~ judge, data = wine)
  pm <- predictmeans(s, "temp", mean_type = "prob_diff")$means
  f <- function(th) {
    na <- length(s$alpha); nb <- length(s$beta)
    g <- expand.grid(temp = levels(wine$temp), contact = levels(wine$contact), judge = levels(wine$judge))
    X <- model.matrix(~ temp + contact, g)[, -1]; S <- model.matrix(~ judge, g)[, -1]
    cp <- stats::plogis(outer(-drop(X %*% th[na + seq_len(nb)]), th[seq_len(na)], "+") /
                        exp(drop(S %*% th[-seq_len(na + nb)])))
    P <- cbind(cp, 1) - cbind(0, cp)
    as.vector(rowsum(P, g$temp) / (nrow(g) / 2))
  }
  D <- numDeriv::jacobian(f, stats::coef(s))
  expect_equal(pm$estimate, f(stats::coef(s)), tolerance = 1e-10)
  expect_equal(pm$std.error, sqrt(diag(D %*% stats::vcov(s) %*% t(D))), tolerance = 1e-6)
  mc <- predictmeans(s, "temp", mean_type = "mean_class")$means
  expect_equal(mc$estimate, as.vector(matrix(f(stats::coef(s)), 2) %*% (1:5)), tolerance = 1e-10)
  # refusals: the link scale with scale effects, aliased terms, other functions
  expect_error(predictmeans(s, "temp", mean_type = "link"), "not linear in the coefficients")
  expect_error(predictmeans(ordinal::clm(rating ~ temp + contact, nominal = ~ contact, data = wine), "temp"),
               "aliased")
  expect_error(contrastmeans(a, "temp", ctrmatrix = rbind(c(1, -1))), "clm nominal or scale effects")
  # a factor of the nominal formula only is a term (no "not a term" note)
  expect_silent(predictmeans(a, "contact", mean_type = "link"))
  # model = FALSE with a missing value in a scale variable: case weights match the rows of the fit
  w2 <- wine
  w2$s <- w2$judge
  w2$s[3] <- NA
  w2$f <- rep(1:2, length.out = nrow(w2))
  s1 <- ordinal::clm(rating ~ temp + contact, scale = ~ s, weights = f, data = w2)
  s0 <- ordinal::clm(rating ~ temp + contact, scale = ~ s, weights = f, data = w2, model = FALSE)
  expect_silent(r0 <- predictmeans(s0, "temp", weights = "proportional"))
  expect_equal(r0$means$estimate, predictmeans(s1, "temp", weights = "proportional")$means$estimate, tolerance = 1e-10)
})

test_that("clm nominal effects: crossing curves give NA means or stop, as for vglm", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  set.seed(1)
  d <- data.frame(x = runif(300, -1, 1), f = factor(sample(c("a", "b"), 300, TRUE)))
  u <- runif(300)
  d$y <- ordered(ifelse(u < plogis(-1 + 3 * d$x), 1, ifelse(u < pmax(plogis(-1 + 3 * d$x), plogis(1 - 3 * d$x)), 2, 3)))
  m <- ordinal::clm(y ~ f, nominal = ~ x, data = d)
  expect_error(predictmeans(m, "f", covariate = list(x = 1)), "cross")
  expect_s3_class(predictmeans(m, "f", covariate = list(x = 0)), "predictmeans")
  pr <- list(c(0.269, 0.354, 0.377), c(0.5, 0.122, 0.378), c(0.5, 0.122, 0.378))
  cells <- data.frame(A = factor(c(0, 1, 0)), B = factor(c(0, 0, 1)))
  d2 <- do.call(rbind, lapply(1:3, function(i) data.frame(cells[rep(i, 3), ], y = factor(1:3), n = round(1000 * pr[[i]]))))
  d2$y <- ordered(d2$y)
  expect_warning(m2 <- ordinal::clm(y ~ 1, nominal = ~ A + B, weights = n, data = d2), "thresholds are increasing")
  r2 <- predictmeans(m2, "A", pairwise = TRUE)
  expect_equal(r2$means$estimable, rep(c(TRUE, FALSE), 3))
  expect_true("class_crossing" %in% advice(r2)$id)
})
