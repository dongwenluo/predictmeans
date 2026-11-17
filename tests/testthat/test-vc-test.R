## The engine of the tests of variance components (R/vc-test.R, step 1):
## each varcomp() row's hypothesis, its refit and the likelihood ratio
## statistic, checked against the reduced models fitted by hand (and
## lmerTest::ranova() where it applies).

q <- function(e) suppressWarnings(suppressMessages(e))
ll <- function(m) as.numeric(logLik(m))
sleep <- function() {
  data(sleepstudy, package = "lme4", envir = environment())
  sleepstudy$g <- factor(rep(1:6, length.out = nrow(sleepstudy)))
  sleepstudy$tf <- factor(sleepstudy$Days)
  sleepstudy
}

test_that("lmer: drop a term, reduce (x | g) to (1 | g), covariance 0; rows as varcomp()", {
  d <- sleep()
  m <- q(lme4::lmer(Reaction ~ Days + (Days | Subject) + (1 | g), d))
  tt <- q(pm_vc_test(m))
  expect_setequal(tt$row, rownames(varcomp(m, ci = FALSE)))
  lr <- stats::setNames(tt$LRT, tt$row)
  expect_equal(lr[["Subject.(Intercept)"]], 2*(ll(m) - ll(q(lme4::lmer(Reaction ~ Days + (1 | g), d)))), tolerance = 1e-6)
  expect_equal(lr[["Subject.Days"]], 2*(ll(m) - ll(q(lme4::lmer(Reaction ~ Days + (1 | Subject) + (1 | g), d)))),
               tolerance = 1e-6)
  expect_equal(lr[["Subject.Days.(Intercept)"]],
               2*(ll(m) - ll(q(lme4::lmer(Reaction ~ Days + (Days || Subject) + (1 | g), d)))), tolerance = 1e-4)
  expect_equal(lr[["g.(Intercept)"]], 2*(ll(m) - ll(q(lme4::lmer(Reaction ~ Days + (Days | Subject), d)))),
               tolerance = 1e-6)
  kind <- stats::setNames(tt$kind, tt$row)
  expect_equal(unname(kind[c("Subject.(Intercept)", "Subject.Days", "Subject.Days.(Intercept)", "g.(Intercept)",
                             "residual")]), c("drop_term", "reduce", "zero_cov", "drop", "none"))
  skip_if_not_installed("lmerTest")
  rt <- lmerTest::ranova(lmerTest::as_lmerModLmerTest(m))
  expect_equal(lr[["Subject.Days"]], rt$LRT[2], tolerance = 1e-6)
})

test_that("glmer: dropping the only random term refits a glm; nested terms (1 | Block/Variety)", {
  data(cbpp, package = "lme4", envir = environment())
  g1 <- q(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp))
  tt <- q(pm_vc_test(g1))
  expect_equal(tt$row, "herd.(Intercept)")       # no residual row for a binomial model
  expect_equal(tt$LRT, 2*(ll(g1) - ll(glm(cbind(incidence, size - incidence) ~ period, binomial, cbpp))),
               tolerance = 1e-6)
  data(Oats, package = "nlme", envir = environment())
  Oats <- as.data.frame(Oats)
  o1 <- q(lme4::lmer(yield ~ Variety*factor(nitro) + (1 | Block/Variety), Oats))
  to <- q(pm_vc_test(o1))
  lr <- stats::setNames(to$LRT, to$row)
  expect_equal(lr[["Block.(Intercept)"]],
               2*(ll(o1) - ll(q(lme4::lmer(yield ~ Variety*factor(nitro) + (1 | Variety:Block), Oats)))), tolerance = 1e-6)
  expect_equal(lr[["Variety:Block.(Intercept)"]],
               2*(ll(o1) - ll(q(lme4::lmer(yield ~ Variety*factor(nitro) + (1 | Block), Oats)))), tolerance = 1e-6)
})

test_that("glmmTMB: unstructured terms, ar1() correlation; residual ar1 (dispformula = ~0) variances not tested", {
  d <- sleep()
  obj <- function(m) m$fit$objective
  t1 <- q(glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject) + (1 | g), d))
  tt <- q(pm_vc_test(t1))
  lr <- stats::setNames(tt$LRT, tt$row)
  expect_equal(lr[["Subject.Days"]], 2*(obj(q(glmmTMB::glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | g), d))) - obj(t1)),
               tolerance = 1e-5)
  expect_equal(lr[["Subject.Days.(Intercept)"]],
               2*(obj(q(glmmTMB::glmmTMB(Reaction ~ Days + diag(Days | Subject) + (1 | g), d))) - obj(t1)), tolerance = 1e-4)
  t2 <- q(glmmTMB::glmmTMB(Reaction ~ Days + ar1(tf + 0 | Subject), dispformula = ~0, data = d))
  ta <- q(pm_vc_test(t2))
  k <- stats::setNames(ta$kind, ta$row)
  expect_equal(unname(k[c("Subject.var", "cor.Subject")]), c("none", "zero_cor"))
  expect_equal(ta$LRT[ta$row == "cor.Subject"],
               2*(obj(q(glmmTMB::glmmTMB(Reaction ~ Days + homdiag(tf + 0 | Subject), dispformula = ~0, data = d))) - obj(t2)),
               tolerance = 1e-5)
})

test_that("semireg: smooth blocks as anova(), other random terms by refit; semireg_tmb with ar1()", {
  d <- sleep()
  s1 <- q(semireg(Reaction ~ s(Days, bs = "os", k = 4) + (Days | Subject), data = d, prt = FALSE))
  tt <- q(pm_vc_test(s1))
  lr <- stats::setNames(tt$LRT, tt$row)
  sm <- anova(s1)$smooth
  expect_equal(lr[["s_Days.(Intercept)"]], sm$LRT[sm$term == "s_Days"], tolerance = 1e-6)
  red <- q(pm_sr_refit(s1, random = "(1 | Subject)"))
  expect_equal(lr[["Subject.Days"]], 2*(ll(s1$semer) - ll(red)), tolerance = 1e-6)
  expect_equal(tt$kind[tt$row == "Subject.Days"], "reduce")
  s2 <- q(semireg_tmb(Reaction ~ s(Days, bs = "os", k = 4) + ar1(tf + 0 | Subject), dispformula = ~0, data = d,
                      prt = FALSE))
  t2 <- q(pm_vc_test(s2))
  expect_true(all(c("cor.Subject", "s_Days.(Intercept)") %in% t2$row))
  expect_true(is.finite(t2$LRT[t2$row == "cor.Subject"]))
  expect_error(pm_vc_test(lm(Reaction ~ Days, d)), "lmer, glmer, glmmTMB and semireg")
})

test_that("likelihood ratio statistics equal SAS PROC MIXED / GLIMMIX and Genstat REML (real runs, 2026-11-17)", {
  # dev/sas/sas-results-2026-11-17-vctest.txt (pm-reference-vctest.sas),
  # dev/genstat/genstat-results-2026-11-17-run4.txt (pm-reference-4.gen):
  # differences of -2 (restricted) log-likelihoods of full and reduced fits
  data(Orthodont, package = "nlme", envir = environment())
  O <- as.data.frame(Orthodont)
  O$occ <- factor(O$age)
  O$agec <- O$age - 11                       # Genstat centres covariates
  lr <- function(m) { t <- q(pm_vc_test(m)); stats::setNames(t$LRT, t$row) }
  a <- lr(q(lme4::lmer(distance ~ age + (age | Subject), O)))
  expect_equal(a[["Subject.(Intercept)"]], 509.16949644 - 442.63668588, tolerance = 1e-6)
  expect_equal(a[["Subject.age"]], 447.00251560 - 442.63668588, tolerance = 1e-6)
  expect_equal(a[["Subject.age.(Intercept)"]], 443.31458017 - 442.63668588, tolerance = 1e-5)
  data(Oats, package = "nlme", envir = environment())
  Oa <- as.data.frame(Oats)
  Oa$N <- factor(Oa$nitro)
  b <- lr(q(lme4::lmer(yield ~ Variety*N + (1 | Block/Variety), Oa)))
  expect_equal(b[["Variety:Block.(Intercept)"]], 536.68996670 - 529.02850701, tolerance = 1e-6)
  expect_equal(b[["Block.(Intercept)"]], 534.00671108 - 529.02850701, tolerance = 1e-6)
  expect_equal(b[["Variety:Block.(Intercept)"]], 426.42 - 418.76, tolerance = 2e-3)   # Genstat, 2 decimals
  data(cbpp, package = "lme4", envir = environment())
  g <- lr(q(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp)))
  expect_equal(g[["herd.(Intercept)"]], 2*99.029199489 - 184.05256373, tolerance = 1e-4)   # Laplace, GLIMMIX
  # random intercept + AR(1) residuals; the AR(1) correlation 0
  ar <- lr(q(glmmTMB::glmmTMB(distance ~ age + (1 | Subject) + ar1(occ + 0 | Subject), dispformula = ~0, data = O,
                              REML = TRUE)))
  expect_equal(ar[["cor.Subject.1"]], 447.00251560 - 446.92550572, tolerance = 1e-4)
  # Genstat: intercept + uncorrelated slope of centred age; the slope variance 0
  gs <- lr(q(lme4::lmer(distance ~ age + (1 | Subject) + (0 + agec | Subject), O)))
  expect_equal(gs[["Subject.agec"]], 252.19 - 250.17, tolerance = 5e-3)
})

test_that("p-values: mixtures equal SAS GLIMMIX COVTEST; the exact restricted LRT has the known mass at 0", {
  # dev/sas/sas-results-2026-11-17-vctest.txt (4 decimals)
  data(Orthodont, package = "nlme", envir = environment())
  O <- as.data.frame(Orthodont)
  a <- q(pm_vc_test(lme4::lmer(distance ~ age + (age | Subject), O)))
  pa <- stats::setNames(a$p.value, a$row)
  expect_equal(pa[["Subject.age"]], 0.0747, tolerance = 6e-5/0.0747)          # 0.5 chi2(1) + 0.5 chi2(2)
  expect_equal(a$test[a$row == "Subject.age"], "0.5 chi2(1) + 0.5 chi2(2)")
  expect_equal(a$test[a$row == "Subject.(Intercept)"], "chi2(3), conservative")
  data(Oats, package = "nlme", envir = environment())
  Oa <- as.data.frame(Oats)
  Oa$N <- factor(Oa$nitro)
  m <- q(lme4::lmer(yield ~ Variety*N + (1 | Block/Variety), Oa))
  b <- q(pm_vc_test(m, method = "mixture"))
  pb <- stats::setNames(b$p.value, b$row)
  expect_equal(pb[["Block.(Intercept)"]], 0.0128, tolerance = 6e-5/0.0128)
  expect_equal(pb[["Variety:Block.(Intercept)"]], 0.0028, tolerance = 6e-5/0.0028)
  # with two random effects the default is the mixture (the exact distribution is exact only for one)
  expect_equal(q(pm_vc_test(m))$p.value, b$p.value)
  # one random effect: the exact restricted LRT, reproducible (fixed seed)
  m1 <- q(lme4::lmer(yield ~ Variety*N + (1 | Block), Oa))
  e1 <- q(pm_vc_test(m1))
  e2 <- q(pm_vc_test(m1))
  expect_identical(e1$p.value, e2$p.value)
  expect_match(e1$test[1], "exact restricted LRT")
  expect_true(e1$p.value[1] > 0 && e1$p.value[1] < 1)
  # balanced one-way design: P(RLRT = 0) = P(F(K - 1, K(J - 1)) < 1) under the null
  K <- 10; J <- 5
  Z <- stats::model.matrix(~ g - 1, data.frame(g = factor(rep(1:K, each = J))))
  null <- pm_rlrt_null(matrix(1, K*J, 1), Z, nsim = 20000)
  expect_equal(mean(null < 1e-8), stats::pf(1, K - 1, K*(J - 1)), tolerance = 4*sqrt(0.25/20000)/0.54)
  # GLMMs: the mixture
  data(cbpp, package = "lme4", envir = environment())
  g <- q(pm_vc_test(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = cbpp)))
  expect_equal(g$test, "0.5 chi2(0) + 0.5 chi2(1)")
  expect_lt(g$p.value, 1e-4)
})

test_that("varcomp(test = TRUE): LRT and p.value by row, hypotheses printed, one warning with guidance", {
  data(Orthodont, package = "nlme", envir = environment())
  O <- as.data.frame(Orthodont)
  m <- q(lme4::lmer(distance ~ age + (age | Subject), O))
  v <- q(varcomp(m, test = TRUE))
  tt <- q(pm_vc_test(m))
  expect_s3_class(v, "pm_varcomp")
  expect_equal(v$p.value, signif(tt$p.value[match(rownames(v), tt$row)], 4))
  expect_true(is.na(v["residual", "p.value"]))
  expect_equal(attr(v, "tests")$row, rownames(v)[rownames(v) %in% tt$row])
  out <- capture.output(print(v))
  expect_true(any(grepl("(age | Subject) -> (1 | Subject); 0.5 chi2(1) + 0.5 chi2(2)", out, fixed = TRUE)))
  # without test = TRUE the table is as before
  expect_false("p.value" %in% names(q(varcomp(m))))
  # few levels: a warning that names the row and how to check
  set.seed(1)
  d <- expand.grid(g = factor(1:4), r = 1:6)
  d$y <- rbinom(nrow(d), 1, 0.4)
  g <- q(lme4::glmer(y ~ 1 + (1 | g), family = binomial, data = d))
  expect_warning(varcomp(g, test = TRUE), "g.\\(Intercept\\) \\(4 levels.*permlmer")
  # classes without random effects: the table without tests, and a message
  expect_message(r <- varcomp(lm(distance ~ age, O), test = TRUE), "lmer, glmer, glmmTMB and semireg")
  expect_false("p.value" %in% names(r))
})

test_that("anova() of a semireg fit uses the same tests as varcomp(test = TRUE); one smooth: exact", {
  set.seed(2)
  d <- data.frame(x = runif(80))
  d$y <- sin(3*d$x) + rnorm(80, 0, 0.3)
  d$g <- gl(8, 10)
  fit <- q(semireg(y ~ s(x, bs = "os", k = 8) + (1 | g), data = d, prt = FALSE))
  a <- anova(fit)$smooth
  v <- q(varcomp(fit, test = TRUE))
  expect_equal(signif(a$p.value, 4), v["s_x.(Intercept)", "p.value"])
  expect_equal(a$test, "0.5 chi2(0) + 0.5 chi2(1)")       # with (1 | g): the mixture
  one <- q(semireg(y ~ s(x, bs = "os", k = 8), data = d, prt = FALSE))
  expect_match(anova(one)$smooth$test, "exact restricted LRT")
})
