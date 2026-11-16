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
