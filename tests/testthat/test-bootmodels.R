test_that("bootmodels(): bootstrap p-values, repeatable, usable as permlist", {
  skip_on_cran()
  data(Oats, package = "nlme", envir = environment())
  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
  set.seed(4); u1 <- runif(1); set.seed(4)
  b1 <- bootmodels(fm, nsim = 19, seed = 1)
  expect_identical(runif(1), u1)
  b2 <- bootmodels(fm, nsim = 19, seed = 1)
  expect_identical(b1$ANOVA, b2$ANOVA)
  expect_identical(attr(b1, "resampling"), "parametric bootstrap")
  expect_length(b1$permlist, 19)
  # nitro has a very large effect: no simulated F reaches the observed one
  expect_equal(b1$ANOVA["nitro", "Boot_p_value"], round(1/20, 4))
  expect_true(all(b1$ANOVA$Boot_p_value >= 1/20 - 1e-8 & b1$ANOVA$Boot_p_value <= 1))
  r <- predictmeans(fm, "nitro", permlist = b1, plot = FALSE, ncore = 1)
  expect_true(all(r$pairs$p.method == "bootstrap"))
  b3 <- bootmodels(fm, nsim = 9, type = "III", seed = 1)
  expect_true(is.na(b3$ANOVA["(Intercept)", "Boot_p_value"]))
})

test_that("bootmodels() works for glmer and refuses other models", {
  skip_on_cran()
  g <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  b <- bootmodels(g, nsim = 9, seed = 2)
  expect_equal(nrow(b$COEFFICENTS), 4)
})

test_that("bootmodels() for lm, aov, glm, gls, lme and glmmTMB; lm p-values approach the exact F tests", {
  skip_on_cran()
  O <- transform(as.data.frame(nlme::Oats), nitro = factor(nitro))
  # lm: F is pivotal, so the bootstrap p-values approach those of the F tests
  fl <- lm(yield ~ Block + nitro*Variety, O)
  b <- bootmodels(fl, nsim = 999, seed = 1)
  pF <- anova(fl)[["Pr(>F)"]]
  ok <- !is.na(pF) & pF > 0.01
  expect_equal(b$ANOVA$Boot_p_value[ok], pF[ok], tolerance = 0.05)
  expect_true(is.na(b$ANOVA["Residuals", "Boot_p_value"]))
  expect_equal(bootmodels(aov(yield ~ Block + nitro*Variety, O), nsim = 19, seed = 1)$ANOVA,
               bootmodels(fl, nsim = 19, seed = 1)$ANOVA)
  # glm: Poisson (deviance), binomial with counts, negative binomial; usable as permlist
  gp <- glm(breaks ~ wool*tension, poisson, warpbreaks)
  bg <- bootmodels(gp, nsim = 19, seed = 1)
  expect_equal(rownames(bg$ANOVA), c("wool", "tension", "wool:tension"))
  r <- predictmeans(gp, "tension", permlist = bg, plot = FALSE)
  expect_true(all(r$pairs$p.method == "bootstrap"))
  expect_length(bootmodels(glm(cbind(incidence, size - incidence) ~ period, binomial, lme4::cbpp), nsim = 9, seed = 1)$permlist, 9)
  expect_length(bootmodels(MASS::glm.nb(breaks ~ wool*tension, warpbreaks), nsim = 9, seed = 1)$permlist, 9)
  # lme and gls: refits to responses simulated from the null model's covariance
  bl <- bootmodels(nlme::lme(yield ~ nitro*Variety, random = ~1 | Block/Variety, O), nsim = 19, seed = 1)
  expect_equal(bl$ANOVA["nitro", "Boot_p_value"], round(1/20, 4))
  bs <- bootmodels(nlme::gls(distance ~ age*Sex, nlme::Orthodont, correlation = nlme::corCompSymm(form = ~1 | Subject)),
                   nsim = 9, seed = 1)
  expect_true(is.na(bs$ANOVA["(Intercept)", "Boot_p_value"]))
  # glmmTMB: no Type I tests, Type II used
  skip_if_not_installed("glmmTMB")
  tm <- glmmTMB::glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = lme4::cbpp)
  expect_message(bt <- bootmodels(tm, nsim = 9, seed = 1), "Type II")
  expect_equal(rownames(bt$ANOVA), "period")
  # refused: quasi families and aov with Error()
  expect_error(bootmodels(glm(breaks ~ wool, quasipoisson, warpbreaks), nsim = 5), "quasi family")
  expect_error(bootmodels(aov(yield ~ nitro + Error(Block), O), nsim = 5), "Error\\(\\)")
})

test_that("bootmodels(): lme and gls responses are simulated with the null model's covariance of the data", {
  g <- nlme::gls(distance ~ age*Sex, nlme::Orthodont, correlation = nlme::corCompSymm(form = ~1 | Subject))
  bp <- predictmeans:::pm_boot_parts(g)
  Y <- do.call(cbind, predictmeans:::pm_with_seed(1, bp$simulate(4000)))
  rho <- coef(bp$null_fit$modelStruct$corStruct, unconstrained = FALSE)[[1]]
  s2 <- bp$null_fit$sigma^2
  expect_equal(mean(apply(Y, 1, var)), s2, tolerance = 0.03)
  # two measurements of one subject, and of two subjects
  expect_equal(cov(Y[1, ], Y[2, ]), rho*s2, tolerance = 0.06)
  expect_equal(cov(Y[1, ], Y[5, ]), 0, tolerance = 0.06*s2)
  expect_equal(rowMeans(Y), as.numeric(fitted(bp$null_fit)), tolerance = 0.02)
})

test_that("clmm: parametric bootstrap and whole-cluster permutation of a between-cluster treatment", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  set.seed(5)
  d <- expand.grid(rep = 1:3, W = factor(1:2), clus = factor(1:8))
  d$trt <- factor(ifelse(as.integer(d$clus) <= 4, "a", "b"))
  eta <- rnorm(8)[d$clus] + 0.8 * (d$trt == "b") + 0.5 * (d$W == 2)
  d$y <- cut(eta + rlogis(nrow(d)), c(-Inf, -1, 0.5, 2, Inf), labels = 1:4, ordered_result = TRUE)
  m <- ordinal::clmm(y ~ trt * W + (1 | clus), data = d)
  b <- suppressMessages(bootmodels(m, nsim = 19, seed = 1))
  expect_equal(rownames(b$ANOVA), c("trt", "W", "trt:W"))
  expect_true(all(is.finite(b$ANOVA$Boot_p_value)))
  expect_identical(suppressMessages(bootmodels(m, nsim = 19, seed = 1))$ANOVA, b$ANOVA)    # reproducible
  # the bootstrap p-values reach the latent pairs
  p <- predictmeans(m, "trt", permlist = b, pairwise = TRUE)
  expect_identical(p$pairs$p.method, "bootstrap")
  # treatments of whole clusters relabelled: only trt is tested; 70 allocations, all used (exact)
  pr <- suppressMessages(permmodels(m, nperm = 99))
  pp <- suppressWarnings(as.numeric(pr$ANOVA$Perm_p_value))
  expect_true(is.finite(pp[1]))
  expect_true(all(is.na(pp[2:3])))
  expect_match(pr$notes$anova, "70 allocations, all used \\(exact p-value; 35 refits")
  # reference: the exact p-value by hand, refitting all 70 allocations of the 8 clusters (4 + 4)
  st <- function(dd) suppressWarnings(pm_anova(ordinal::clmm(y ~ trt * W + (1 | clus), data = dd), type = "II",
                                               test_statistic = "Chisq")["trt", "Chisq"])
  obs <- st(d)
  al <- utils::combn(8, 4)
  sims <- apply(al, 2, function(a) {
    dd <- d
    dd$trt <- factor(ifelse(as.integer(dd$clus) %in% a, "a", "b"))
    st(dd)
  })
  expect_equal(pp[1], round(mean(sims >= obs * (1 - 1e-4)), 4))
  # the same treatment coded 0/1: the same exact p-value (two values: one refit per relabelling too)
  d01 <- d
  d01$trt <- as.numeric(d$trt == "b")
  m01 <- ordinal::clmm(y ~ trt * W + (1 | clus), data = d01)
  expect_equal(as.numeric(suppressMessages(permmodels(m01, nperm = 99))$ANOVA["trt", "Perm_p_value"]), pp[1])
  # a cluster-level baseline covariate stays with its cluster (not a stratum); random allocations reproducible
  d$x <- stats::rnorm(8)[d$clus]
  mx <- ordinal::clmm(y ~ trt + x + (1 | clus), data = d)
  px <- suppressMessages(permmodels(mx, nperm = 19, seed = 3))
  expect_match(px$notes$anova, "trt: relabelled among the levels of clus; 70 allocations, 19 random ones")
  expect_identical(suppressMessages(permmodels(mx, nperm = 19, seed = 3))$ANOVA, px$ANOVA)
  # a farm-level treatment under (1 | farm/plot): farms are the units
  set.seed(4)
  fd <- expand.grid(sub = 1:3, plot = factor(1:2), farm = factor(1:6))
  fd$trt <- factor(ifelse(as.integer(fd$farm) <= 3, "a", "b"))
  fd$y <- cut(rnorm(6)[fd$farm] + 0.7 * (fd$trt == "b") + rlogis(nrow(fd)), c(-Inf, -0.5, 0.8, Inf), labels = 1:3,
              ordered_result = TRUE)
  mf <- suppressWarnings(ordinal::clmm(y ~ trt + (1 | farm/plot), data = fd))
  pf <- suppressMessages(permmodels(mf, nperm = 99))
  expect_match(pf$notes$anova, "trt: relabelled among the levels of farm; 20 allocations")
  expect_true(as.numeric(pf$ANOVA["trt", "Perm_p_value"]) >= 0.1)
  # an RCBD with subsamples: trt is randomised within blocks (2^3 = 8 allocations, smallest p 2/8)
  set.seed(2)
  r <- expand.grid(sub = 1:4, trt = factor(c("a", "b")), block = factor(1:3))
  r$plot <- factor(paste(r$block, r$trt))
  r$y <- cut(rnorm(3)[r$block] + rnorm(6)[r$plot] + 0.5 * (r$trt == "b") + rlogis(nrow(r)), c(-Inf, -0.5, 0.8, Inf),
             labels = 1:3, ordered_result = TRUE)
  mr <- suppressWarnings(ordinal::clmm(y ~ block + trt + (1 | plot), data = r))
  prr <- suppressMessages(permmodels(mr, nperm = 99))
  expect_match(prr$notes$anova, "trt: relabelled among the levels of plot within block; 8 allocations")
  expect_true(as.numeric(prr$ANOVA["trt", "Perm_p_value"]) >= 0.25)
  expect_error(suppressMessages(permmodels(m, perms = permindex(d, nsim = 19, seed = 1))), "leave 'perms' out")
  # a permlist of permmodels() does not reach the pairs (only cluster-level terms are tested)
  expect_error(predictmeans(m, "W", permlist = pr), "must come from bootmodels")
})

test_that("clmm with its random-effect SD at 0: covariance of the other parameters, the SD known", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  data(wine, package = "ordinal", envir = environment())
  m <- ordinal::clmm(rating ~ temp * contact + (1 | judge), data = wine)
  bp <- pm_boot_parts(m)
  ys <- pm_with_seed(1, bp$simulate(99))
  r <- NULL
  for (y in ys) {                                        # a refit whose vcov() fails (SD at 0, log-SD route)
    f <- bp$refit(y)
    if (inherits(tryCatch(stats::vcov(f), error = function(e) e), "error")) {
      r <- f
      break
    }
  }
  skip_if(is.null(r), "no refit with a failing vcov()")
  V <- pm_ord_vcov(r)
  expect_true(isTRUE(attr(V, "sd_fixed")))
  # reference: with the SD at 0 the model is clm()'s, with the same covariance
  dd <- wine
  dd$.pm_y <- y
  ce <- ordinal::clm(.pm_y ~ temp * contact, data = dd)
  expect_equal(unname(V[names(stats::coef(ce)), names(stats::coef(ce))]), unname(stats::vcov(ce)), tolerance = 1e-4)
  p <- predictmeans(r, "temp")
  expect_true("singular_fit" %in% advice(p)$id)             # the SD at 0, treated as known, is said
})
