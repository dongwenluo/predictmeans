## varcomp(): every model class, one rule for SEs (delta method, observed
## information) and intervals (Satterthwaite's chi-square for variances).
## References to SAS and Genstat are in test-sas-reference.R and
## test-genstat-reference.R.

vc_oats <- function() transform(as.data.frame(nlme::Oats), nitro = factor(nitro))

test_that("lmer, lme, aov with Error() and glmmTMB give the same components of a balanced split plot", {
  skip_on_cran()
  O <- vc_oats()
  a <- varcomp(lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = O))
  b <- varcomp(nlme::lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = O))
  d <- suppressMessages(varcomp(aov(yield ~ nitro*Variety + Error(Block/Variety), data = O)))
  e <- varcomp(glmmTMB::glmmTMB(yield ~ nitro*Variety + (1|Block/Variety), data = O, REML = TRUE))
  expect_identical(rownames(a), rownames(e))
  expect_identical(rownames(b), c("Block.(Intercept)", "Variety.(Intercept)", "residual"))
  expect_identical(rownames(d), c("Block", "Block:Variety", "residual"))
  ref <- a[c(2, 1, 3), ]
  for (x in list(b, d)) expect_equal(unname(as.matrix(x)), unname(as.matrix(ref)), tolerance = 1e-4)
  expect_equal(unname(as.matrix(e)), unname(as.matrix(a)), tolerance = 1e-4)
  expect_identical(attr(a, "ci_method"), "Satterthwaite")
})

test_that("lmer SEs: the delta method on (theta, sigma) equals the information of R/kr.R", {
  skip_on_cran()
  O <- vc_oats()
  m <- lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = O[-c(3, 17, 40), ])
  W <- attr(pm_vcov_kr(m, "observed"), "W")
  expect_equal(varcomp(m)$SE, sqrt(diag(W)), tolerance = 1e-4)
  W2 <- attr(pm_vcov_kr(m, "expected"), "W")
  expect_equal(varcomp(m, information = "expected")$SE, round(sqrt(diag(W2)), 4))
})

test_that("Satterthwaite's interval: nu = 2 (est/SE)^2; the exact chi-square interval for lm", {
  O <- vc_oats()
  m <- lm(yield ~ Block + nitro*Variety, data = O)
  v <- varcomp(m)
  s2 <- summary(m)$sigma^2
  df <- m$df.residual
  expect_equal(v$vcov, round(s2, 4))
  expect_equal(v$SE, round(s2 * sqrt(2/df), 4))
  expect_equal(unname(unlist(v[1, 3:4])), round(df * s2 / qchisq(c(0.975, 0.025), df), 4))
  vw <- varcomp(m, ci_method = "Wald")
  expect_equal(unname(unlist(vw[1, 3:4])), round(s2 * exp(c(-1, 1) * qnorm(0.975) * sqrt(2/df)), 4))
  expect_equal(colnames(varcomp(m, ci = FALSE)), c("vcov", "SE"))
  expect_equal(varcomp(glm(yield ~ Block + nitro*Variety, data = O))$vcov, v$vcov)
})

test_that("glm: Pearson dispersion, fixed dispersion, negative binomial theta", {
  O <- vc_oats()
  g <- glm(yield ~ Block + nitro, family = Gamma(link = "log"), data = O)
  expect_equal(varcomp(g)["dispersion", "vcov"], round(summary(g)$dispersion, 4))
  p <- varcomp(glm(round(yield/10) ~ nitro, family = poisson, data = O))
  expect_equal(p$vcov, 1)
  expect_true(is.na(p$SE))
  expect_match(attr(p, "note"), "fixed at 1")
  set.seed(1)
  d <- data.frame(x = rep(1:4, 25))
  d$y <- MASS::rnegbin(100, mu = exp(1 + d$x/3), theta = 2)
  nb <- MASS::glm.nb(y ~ x, data = d)
  vn <- varcomp(nb)
  expect_equal(c(vn$vcov, vn$SE), round(c(nb$theta, nb$SE.theta), 4))
})

test_that("aov with Error(): ANOVA estimates of the stratum variances, also with missing values", {
  skip_on_cran()
  O <- vc_oats()
  O$yield[c(5, 30)] <- NA
  v <- suppressMessages(varcomp(aov(yield ~ nitro*Variety + Error(Block/Variety), data = O)))
  expect_identical(rownames(v), c("Block", "Block:Variety", "residual"))
  # the components solve E(MS) = C sigma2 for the stratum mean squares of predictmeans' analysis
  st <- suppressMessages(pm_aov_strata(pm_aov_carrier(aov(yield ~ nitro*Variety + Error(Block/Variety), data = O))))$strata
  expect_equal(v["residual", "vcov"], round(st$ms[st$stratum == "Within"], 4))
  expect_true(all(v$vcov > 0))
  # an incomplete-block design (not orthogonal) gives a component too
  set.seed(4)
  ib <- data.frame(blk = factor(rep(1:6, each = 3)), trt = factor(c(1,2,3, 1,2,4, 1,3,4, 2,3,4, 1,2,3, 2,3,4)))
  ib$y <- rnorm(18) + as.numeric(ib$trt) + rnorm(6)[ib$blk]
  vi <- suppressMessages(varcomp(aov(y ~ trt + Error(blk), data = ib)))
  expect_identical(rownames(vi), c("blk", "residual"))
})

test_that("complex structures: random slopes, || terms, pdDiag, corSymm, AR(1), varIdent, varPower", {
  skip_on_cran()
  data(Orthodont, package = "nlme")
  a <- varcomp(lme4::lmer(distance ~ age + (age|Subject), data = Orthodont))
  b <- varcomp(nlme::lme(distance ~ age, random = ~age|Subject, data = Orthodont))
  expect_identical(rownames(a), c("Subject.(Intercept)", "Subject.age.(Intercept)", "Subject.age", "residual"))
  expect_identical(rownames(b), rownames(a))
  expect_equal(a$SE, b$SE, tolerance = 2e-3)
  # the covariance: the interval of the correlation times the two SDs
  expect_true(a[2, 3] < a[2, 1] && a[2, 1] < a[2, 4])
  dd <- varcomp(lme4::lmer(distance ~ age + (age || Subject), data = Orthodont))
  expect_identical(rownames(dd), c("Subject.(Intercept)", "Subject.age", "residual"))
  pd <- varcomp(nlme::lme(distance ~ age, random = list(Subject = nlme::pdDiag(~age)), data = Orthodont))
  expect_identical(rownames(pd), rownames(dd))
  expect_equal(pd$vcov, dd$vcov, tolerance = 1e-3)
  un <- varcomp(nlme::gls(distance ~ factor(age), data = Orthodont,
                          correlation = nlme::corSymm(form = ~1|Subject), weights = nlme::varIdent(form = ~1|age)))
  expect_identical(rownames(un)[1:6], paste0("cor.", c("1.2", "1.3", "1.4", "2.3", "2.4", "3.4")))
  expect_true(all(un[1:6, 3] < un[1:6, 1] & un[1:6, 1] < un[1:6, 4] & un[1:6, 4] < 1))
  ar <- varcomp(nlme::lme(distance ~ age, random = ~1|Subject, correlation = nlme::corAR1(),
                          weights = nlme::varPower(), data = Orthodont))
  expect_identical(rownames(ar), c("Subject.(Intercept)", "cor.Phi", "var.power", "residual"))
  # varPower's power: a Wald interval on its own scale
  expect_equal(ar["var.power", 4] - ar["var.power", "vcov"], ar["var.power", "vcov"] - ar["var.power", 3],
               tolerance = 1e-3)
})

test_that("glmer: several grouping factors, nAGQ = 0, the same as glmmTMB", {
  skip_on_cran()
  data(cbpp, package = "lme4")
  g <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp)
  v <- varcomp(g)
  t <- varcomp(glmmTMB::glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial,
                                data = cbpp))
  expect_equal(v$vcov, t$vcov, tolerance = 1e-3)
  expect_equal(v$SE, t$SE, tolerance = 0.02)
  cbpp$obs <- factor(seq_len(nrow(cbpp)))
  g2 <- suppressMessages(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd) + (1 | obs),
                                     family = binomial, data = cbpp))
  expect_identical(rownames(varcomp(g2)), c("obs.(Intercept)", "herd.(Intercept)"))
  g0 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp,
                    nAGQ = 0)
  v0 <- varcomp(g0)
  expect_true(is.finite(v0$SE) && v0$SE > 0)
})

test_that("glmmTMB structured terms: variance and correlation parameters", {
  skip_on_cran()
  data(Orthodont, package = "nlme")
  m <- glmmTMB::glmmTMB(distance ~ age + ar1(factor(age) + 0 | Subject), data = Orthodont)
  v <- varcomp(m)
  expect_identical(rownames(v), c("Subject.var", "cor.Subject", "residual"))
  R <- stats::cov2cor(glmmTMB::VarCorr(m)$cond$Subject)
  expect_equal(v["cor.Subject", "vcov"], round(R[2, 1], 4))
})

test_that("gam and semireg fits; information and ci_method choices", {
  skip_on_cran()
  set.seed(1)
  d <- data.frame(x = runif(200), g = factor(rep(1:20, 10)))
  d$y <- sin(6*d$x) + rnorm(20)[d$g] + rnorm(200, sd = 0.3)
  v <- varcomp(mgcv::gam(y ~ s(x) + s(g, bs = "re"), data = d, method = "REML"))
  expect_identical(rownames(v), c("s(x)", "s(g)", "residual"))
  expect_true(all(v$SE > 0))
  # gam.vcomp()'s shapes: matrix or vector (mgcv < 1.9-3), a list holding them (later)
  m <- cbind(std.dev = c(a = 1, scale = 2), lower = c(0.5, 1.5), upper = c(2, 3))
  for (x in list(m, list(all = NULL, vc = m)))
    expect_identical(pm_gam_vcomp_matrix(x), m)
  for (x in list(c(a = 1, scale = 2), list(vc = c(a = 1, scale = 2), all = c(a = 1, b = 1, scale = 2))))
    expect_identical(pm_gam_vcomp_matrix(x), cbind(std.dev = c(a = 1, scale = 2), lower = NA_real_, upper = NA_real_))
  data(Dialyzer)
  fm <- semireg(rate ~ pressure*QB, smoothZ = list(qb = smZ(pressure, k = 4, by = QB, group = TRUE)),
                data = Dialyzer, prt = FALSE)
  expect_identical(rownames(varcomp(fm)), c("qb.200.(Intercept)", "qb.300.(Intercept)", "residual"))
  m <- lme4::lmer(yield ~ nitro + (1|Block/Variety), data = vc_oats(), REML = FALSE)
  expect_message(varcomp(m, information = "average"), "REML fits")
  expect_message(varcomp(lm(yield ~ nitro, data = vc_oats()), ci_method = "profile"), "lmer and glmer")
})

test_that("edge cases: a variance at zero (Dyestuff2), crossed effects (Penicillin), one-way (Dyestuff)", {
  skip_on_cran()
  data(Dyestuff, Dyestuff2, Penicillin, package = "lme4")
  b1 <- suppressMessages(varcomp(lme4::lmer(Yield ~ 1 + (1|Batch), data = Dyestuff2)))
  b2 <- varcomp(nlme::lme(Yield ~ 1, random = ~1|Batch, data = Dyestuff2))
  for (b in list(b1, b2)) {
    expect_equal(b["Batch.(Intercept)", "vcov"], 0)
    expect_true(is.na(b["Batch.(Intercept)", "SE"]) && is.na(b["Batch.(Intercept)", 4]))
    expect_equal(b["Batch.(Intercept)", 3], 0)
  }
  expect_equal(b1["residual", ], b2["residual", ], tolerance = 1e-4)
  # the ANOVA estimate is negative there
  ba <- suppressMessages(varcomp(aov(Yield ~ 1 + Error(Batch), data = Dyestuff2)))
  expect_lt(ba["Batch", "vcov"], 0)
  # balanced designs: REML = ANOVA, with the same SEs and intervals
  p1 <- varcomp(lme4::lmer(diameter ~ 1 + (1|plate) + (1|sample), data = Penicillin))
  p2 <- suppressMessages(varcomp(aov(diameter ~ 1 + Error(plate + sample), data = Penicillin)))
  expect_equal(unname(as.matrix(p1)), unname(as.matrix(p2)), tolerance = 1e-3)
  d1 <- varcomp(lme4::lmer(Yield ~ 1 + (1|Batch), data = Dyestuff))
  d2 <- suppressMessages(varcomp(aov(Yield ~ 1 + Error(Batch), data = Dyestuff)))
  expect_equal(unname(as.matrix(d1)), unname(as.matrix(d2)), tolerance = 1e-4)
})
