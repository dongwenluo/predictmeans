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
