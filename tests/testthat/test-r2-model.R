test_that("R2_model(): Piepho's values as numbers (R2_glmm(), removed in 2.0, printed them in percent)", {
  fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  g <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  for (m in list(fm, g)) {
    v <- R2_model(m)
    expect_s3_class(v, "pm_r2")
    expect_equal(100*unclass(v), pm_r2_piepho(m), ignore_attr = TRUE)
  }
  expect_error(R2_glmm(fm), "use R2_model(model)", fixed = TRUE, class = "defunctError")
  expect_error(ci_mcp(1, 2), "letter_ci = TRUE", fixed = TRUE, class = "defunctError")
})

test_that("R2_model(method = 'nakagawa'): published sleepstudy values; random slopes as Johnson (2014)", {
  fm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  v <- R2_model(fm, "nakagawa")
  # performance::r2_nakagawa(): marginal 0.279, conditional 0.799
  expect_equal(round(unclass(v), 3), c(Marginal = 0.279, Conditional = 0.799), ignore_attr = TRUE)
  X <- getME(fm, "X")
  S <- VarCorr(fm)$Subject
  vr <- mean(rowSums((X %*% S) * X))
  expect_equal(attr(v, "variances")[["random"]], vr, tolerance = 1e-10)
  expect_equal(attr(v, "variances")[["observation"]], sigma(fm)^2)
  # same model fitted by glmmTMB
  ft <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = TRUE)
  expect_equal(unclass(R2_model(ft, "nakagawa")), unclass(v), tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("R2_model(method = 'nakagawa'): latent variance for binomial, error for other families", {
  g <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), data = cbpp, family = binomial)
  v <- R2_model(g, "nakagawa")
  expect_equal(attr(v, "variances")[["observation"]], pi^2/3)
  vf <- var(drop(getME(g, "X") %*% fixef(g)))
  vr <- VarCorr(g)$herd[1, 1]
  expect_equal(unclass(v)[["Marginal"]], vf/(vf + vr + pi^2/3), tolerance = 1e-10)
  gp <- glmer(incidence ~ period + (1 | herd), data = cbpp, family = poisson)
  expect_error(R2_model(gp, "nakagawa"), "gaussian and binomial")
})

test_that("R2_model(): every model class; equal to the known value of each", {
  u <- function(x) unclass(x)[seq_along(x)]
  Oats <- transform(nlme::Oats, nitro = factor(nitro), Block = factor(Block, ordered = FALSE))
  l <- lm(yield ~ nitro*Variety, Oats)
  adj <- summary(l)$adj.r.squared
  # lm, aov, gaussian glm, gls and gam without smooths: the adjusted R2
  expect_equal(u(R2_model(l)), c(Total = adj, Fixed = adj, Random = 0), tolerance = 1e-10)
  expect_equal(u(R2_model(aov(yield ~ nitro*Variety, Oats))), u(R2_model(l)), tolerance = 1e-10)
  expect_equal(u(R2_model(glm(yield ~ nitro*Variety, data = Oats))), u(R2_model(l)), tolerance = 1e-10)
  expect_equal(u(R2_model(nlme::gls(yield ~ nitro*Variety, Oats))), u(R2_model(l)), tolerance = 1e-6)
  expect_equal(u(R2_model(mgcv::gam(yield ~ nitro*Variety, data = Oats))), u(R2_model(l)), tolerance = 1e-8)
  expect_output(print(R2_model(l)), "Adjusted R2 \\(Piepho\\)")
  # lme and aov with Error(): as lmer
  fm <- lme4::lmer(yield ~ nitro*Variety + (1 | Block/Variety), Oats)
  r <- u(R2_model(fm))
  fl <- u(R2_model(nlme::lme(yield ~ nitro*Variety, random = ~1 | Block/Variety, Oats)))
  expect_equal(unname(fl[1:3]), unname(r[1:3]), tolerance = 1e-3)
  expect_equal(unname(fl[c("Block", "Variety")]), unname(r[c("Block", "Variety:Block")]), tolerance = 1e-3)
  fa <- u(R2_model(aov(yield ~ nitro*Variety + Error(Block/Variety), Oats)))
  expect_equal(unname(fa[1:3]), unname(r[1:3]), tolerance = 1e-3)
  expect_equal(unname(fa[c("Block", "Block:Variety")]), unname(r[c("Block", "Variety:Block")]), tolerance = 1e-3)
  # random slopes
  s1 <- u(R2_model(lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)))
  s2 <- u(R2_model(nlme::lme(Reaction ~ Days, random = ~Days | Subject, lme4::sleepstudy)))
  expect_equal(unname(s2), unname(s1), tolerance = 1e-3)
  # gls with compound symmetry: the fixed share of the random-intercept model
  data(Orthodont, package = "nlme")
  g1 <- u(R2_model(nlme::gls(distance ~ age*Sex, Orthodont, correlation = nlme::corCompSymm(form = ~1 | Subject))))
  g2 <- u(R2_model(nlme::lme(distance ~ age*Sex, random = ~1 | Subject, Orthodont)))
  expect_equal(g1[["Fixed"]], g2[["Fixed"]], tolerance = 1e-4)
  expect_equal(g1[["Random"]], 0)
  g4 <- R2_model(nlme::gls(distance ~ age*Sex, Orthodont, correlation = nlme::corAR1(form = ~1 | Subject),
                           weights = nlme::varIdent(form = ~1 | Sex)))
  expect_true(g4[["Fixed"]] > 0 && g4[["Fixed"]] < 1)
  # gam: a random-effect smooth as the lmer random effect
  m1 <- u(R2_model(mgcv::gam(yield ~ nitro*Variety + s(Block, bs = "re"), data = Oats, method = "REML")))
  m0 <- u(R2_model(lme4::lmer(yield ~ nitro*Variety + (1 | Block), Oats)))
  expect_equal(unname(m1), unname(m0), tolerance = 1e-3)
  m2 <- mgcv::gam(accel ~ s(times), data = MASS::mcycle, method = "REML")
  expect_equal(R2_model(m2)[["Fixed"]], summary(m2)$r.sq, tolerance = 0.02)
})

test_that("R2_model(): glm as the glmmTMB fit without random effects; binary on the latent scale", {
  u <- function(x) unclass(x)[seq_along(x)]
  p1 <- glm(breaks ~ wool*tension, poisson, warpbreaks)
  p2 <- glmmTMB::glmmTMB(breaks ~ wool*tension, family = poisson, data = warpbreaks)
  expect_equal(unname(u(R2_model(p1))), unname(u(R2_model(p2))[1:3]), tolerance = 1e-4)
  # overdispersion takes the R2 down
  expect_lt(R2_model(glm(breaks ~ wool*tension, quasipoisson, warpbreaks))[["Fixed"]], R2_model(p1)[["Fixed"]])
  expect_true(is.finite(R2_model(MASS::glm.nb(breaks ~ wool*tension, warpbreaks))[["Fixed"]]))
  b1 <- glm(am ~ wt, binomial, mtcars)
  X <- model.matrix(b1); b <- coef(b1); n <- nrow(X)
  XtPX <- crossprod(X) - tcrossprod(colSums(X))/n
  th <- (drop(b %*% XtPX %*% b) - sum(XtPX * vcov(b1)))/(n - 1)
  expect_equal(R2_model(b1)[["Fixed"]], th/(th + pi^2/3), tolerance = 1e-10)
  expect_error(R2_model(b1, "nakagawa"), "lmer, glmer, glmmTMB and lme")
  expect_error(R2_model(lm(cbind(mpg, wt) ~ hp, mtcars)), "several responses")
})

test_that("R2_model(method = 'nakagawa') for lme equals lmer; variance functions give the mean residual variance", {
  u <- function(x) unclass(x)[seq_along(x)]
  f1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  f2 <- nlme::lme(Reaction ~ Days, random = ~Days | Subject, lme4::sleepstudy)
  expect_equal(u(R2_model(f2, "nakagawa")), u(R2_model(f1, "nakagawa")), tolerance = 1e-3)
  expect_equal(unname(u(R2_model(f2, "nakagawa"))), c(0.279, 0.799), tolerance = 0.01)
  O <- transform(nlme::Oats, nitro = factor(nitro), Block = factor(Block, ordered = FALSE))
  g1 <- lme4::lmer(yield ~ nitro*Variety + (1 | Block/Variety), O)
  g2 <- nlme::lme(yield ~ nitro*Variety, random = ~1 | Block/Variety, O)
  expect_equal(u(R2_model(g2, "nakagawa")), u(R2_model(g1, "nakagawa")), tolerance = 1e-3)
  # a variance function: the observation-level variance is the mean of sigma^2 w_i^-2
  data(Orthodont, package = "nlme")
  h <- nlme::lme(distance ~ age*Sex, random = ~1 | Subject, Orthodont, weights = nlme::varIdent(form = ~1 | Sex))
  v <- attr(R2_model(h, "nakagawa"), "variances")
  w <- nlme::varWeights(h$modelStruct$varStruct)
  expect_equal(v[["observation"]], mean(h$sigma^2/w^2), tolerance = 1e-8)
})

test_that("Piepho R2 reproduces Table 2 of Piepho (2023) (SAS GLIMMIX values, %)", {
  # published to 2 decimals; the paper's egg random share (41.80) is its total minus
  # its fixed share, so it is checked to 0.01 instead of 0.005
  r2 <- function(m, ...) 100 * unclass(R2_model(m, ...))[c("Fixed", "Random", "Total")]
  e <- new.env()
  data(Body_length, Poisson_ind, package = "predictmeans", envir = e)
  load(test_path("fixtures", "Binary_color.RData"), envir = e)
  body <- lme4::lmer(BodyL ~ Treatment + Habitat + Sex + (1 | Population/Container), data = e$Body_length)
  expect_lt(max(abs(r2(body) - c(Fixed = 40.09, Random = 33.30, Total = 73.39))), 0.005)
  egg <- lme4::glmer(Egg ~ Treatment + Habitat + (1 | Population/Container/individual), family = poisson,
                     data = e$Poisson_ind)
  expect_lt(max(abs(r2(egg, over_disp = TRUE) - c(Fixed = 9.13, Random = 41.80, Total = 50.93))), 0.01)
  colour <- lme4::glmer(Colour ~ Treatment + Habitat + (1 | Population/Container), family = binomial,
                        data = e$Binary_color)
  # (GLMMs to 0.01: the glmer estimates move with lme4's version, by 0.006 for lme4 1.1-35)
  expect_lt(max(abs(r2(colour) - c(Fixed = 7.46, Random = 21.99, Total = 29.46))), 0.01)
})
