oats_se <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  d
}

test_that("conditional SEs are lme4's condVar (1.x numbers)", {
  d <- oats_se()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = d)
  pv <- lme4::ranef(fm, condVar = TRUE)
  l <- as.list(se_ranef(fm, type = "conditional"))
  for (g in names(pv)) expect_equal(unname(l[[g]][, 1]), sqrt(attr(pv[[g]], "postVar")[1, 1, ]), tolerance = 1e-8)
  # correlated slopes
  f2 <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  c2 <- as.list(se_ranef(f2, "Subject", type = "conditional"))
  pv2 <- attr(lme4::ranef(f2, condVar = TRUE)$Subject, "postVar")
  expect_equal(unname(c2), t(sapply(1:18, function(j) sqrt(diag(pv2[, , j])))), tolerance = 1e-8)
  expect_equal(colnames(c2), c("(Intercept)", "Days"))
})

test_that("prediction SEs are the inverse of the mixed model equations", {
  d <- oats_se()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = d)
  X <- as.matrix(lme4::getME(fm, "X"))
  Z <- as.matrix(lme4::getME(fm, "Z"))
  vc <- as.data.frame(lme4::VarCorr(fm))
  G <- diag(rep(vc$vcov[1:2], c(18, 6)))
  s2 <- vc$vcov[3]
  W <- cbind(X, Z)
  M <- crossprod(W) / s2 + rbind(matrix(0, ncol(X), ncol(W)), cbind(matrix(0, ncol(Z), ncol(X)), solve(G)))
  C <- solve(M)
  se_dense <- sqrt(diag(C)[-seq_len(ncol(X))])
  r <- se_ranef(fm)
  expect_equal(r$std.error, unname(se_dense), tolerance = 1e-8)
  rc <- se_ranef(fm, type = "conditional")
  expect_true(all(r$std.error >= rc$std.error - 1e-10))
  expect_s3_class(r, "pm_ranef")
  expect_named(as.data.frame(r), c("group", "level", "term", "estimate", "std.error", "conf.low", "conf.high"))
  expect_equal(generics::tidy(r), as.data.frame(r))
  expect_equal(attr(r, "settings")$type, "prediction")
  expect_output(print(r), "prediction error SE")
})

test_that("lme and lmer agree; lme with corAR1 runs", {
  d <- oats_se()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = d)
  m <- nlme::lme(yield ~ nitro * Variety, random = ~ 1 | Block / Variety, data = d)
  a <- se_ranef(fm, predict = TRUE)
  b <- se_ranef(m, predict = TRUE)
  key <- function(x, g) x[x$group == g, ]
  expect_equal(key(b, "Block")$std.error, key(a, "Block")$std.error, tolerance = 1e-4)
  expect_equal(sort(key(b, "Variety")$prediction), sort(key(a, "Variety:Block")$prediction), tolerance = 1e-4)
  expect_equal(key(b, "Block")$estimate, unname(nlme::ranef(m)$Block[, 1]), tolerance = 1e-6)
  expect_equal(se_ranef(m, "Block", type = "conditional")$std.error,
               se_ranef(fm, "Block", type = "conditional")$std.error, tolerance = 1e-4)
  d2 <- lme4::sleepstudy
  m2 <- nlme::lme(Reaction ~ Days, random = ~ 1 | Subject, correlation = nlme::corAR1(), data = d2)
  r2 <- se_ranef(m2, type = "prediction")
  expect_equal(r2$estimate, unname(nlme::ranef(m2)[, 1]), tolerance = 1e-6)
  expect_true(all(r2$std.error >= se_ranef(m2, type = "conditional")$std.error - 1e-10))
})

test_that("predict adds the fixed part: overall mean for Block, the Variety mean for Block:Variety", {
  d <- oats_se()
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = d)
  pm <- predictmeans(fm, "Variety", vcov_adjust = FALSE)$means
  r <- se_ranef(fm, predict = TRUE)
  bl <- r[r$group == "Block", ]
  expect_equal(bl$prediction - bl$estimate, rep(mean(pm$estimate), 6), tolerance = 1e-8)
  bv <- r[r$group == "Variety:Block", ]
  v <- sub(":.*", "", bv$level)
  expect_equal(bv$prediction - bv$estimate, pm$estimate[match(v, pm$Variety)], tolerance = 1e-8)
  expect_true(all(bv$prediction.std.error > 0))
  expect_error(se_ranef(fm, "Plot"), "\"Variety:Block\", \"Block\"")
  expect_equal(dim(as.list(se_ranef(fm, "Block"))), c(6L, 1L))
})

test_that("crossed terms, zero variance, glmer, glmmTMB, semireg, and models without random effects", {
  d <- oats_se()
  fc <- lme4::lmer(yield ~ nitro + (1 | Block) + (1 | Variety), data = d)
  rc <- se_ranef(fc)
  expect_setequal(unique(rc$group), c("Block", "Variety"))
  expect_true("few_re_levels" %in% advice(rc)$id)
  set.seed(1)
  d$y0 <- d$yield + rnorm(72, 0, 0.01) - ave(d$yield, d$Block)
  f0 <- suppressMessages(lme4::lmer(y0 ~ Variety + (1 | Block), data = d))
  r0 <- se_ranef(f0)
  expect_true("variance_at_zero" %in% advice(r0)$id)
  expect_equal(r0$std.error, rep(0, 6), tolerance = 1e-6)
  # glmer: link and response
  gm <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = lme4::cbpp)
  g <- se_ranef(gm, predict = TRUE)
  pv <- attr(lme4::ranef(gm, condVar = TRUE)$herd, "postVar")[1, 1, ]
  expect_equal(se_ranef(gm, type = "conditional")$std.error, sqrt(pv), tolerance = 1e-6)
  expect_true(all(g$std.error >= sqrt(pv) - 1e-10))
  expect_equal(g$response, plogis(g$prediction))
  expect_true(all(c("response.std.error", "response.conf.low", "response.conf.high") %in% names(g)))
  expect_match(attr(g, "settings")$scale, "logit")
  # glmmTMB: us and ar1
  t1 <- glmmTMB::glmmTMB(yield ~ nitro * Variety + (1 | Block / Variety), data = d, REML = TRUE)
  a1 <- se_ranef(t1)
  a0 <- se_ranef(fc <- lme4::lmer(yield ~ nitro * Variety + (1 | Block / Variety), data = d))
  expect_equal(sort(a1$std.error), sort(a0$std.error), tolerance = 1e-3)
  s <- lme4::sleepstudy
  s$t <- factor(s$Days)
  t2 <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), data = s)
  t3 <- glmmTMB::glmmTMB(Reaction ~ Days + ar1(t + 0 | Subject), data = s)
  r2 <- se_ranef(t2)
  r3 <- se_ranef(t3)
  expect_equal(nrow(r2), 36)
  expect_equal(nrow(r3), 180)
  expect_true(all(is.finite(r3$std.error) & r3$std.error > 0))
  cv <- attr(glmmTMB::ranef(t2)$cond$Subject, "condVar")
  expect_equal(as.list(se_ranef(t2, type = "conditional"))$Subject[, 2], sqrt(cv[2, 2, ]), ignore_attr = TRUE, tolerance = 1e-6)
  # semireg
  d$n0 <- d$nitro0 <- as.numeric(as.character(d$nitro))
  sm <- semireg(yield ~ Variety + s(n0, k = 4) + (1 | Block), data = d)
  rs <- se_ranef(sm, "Block", predict = TRUE)
  expect_equal(rs$estimate, unname(lme4::ranef(sm$semer)$Block[, 1]), tolerance = 1e-6)
  expect_true(all(is.finite(rs$prediction.std.error)))
  # no random effects
  expect_error(se_ranef(nlme::gls(yield ~ Variety, data = d)), "random effects")
  expect_error(se_ranef(lm(yield ~ Variety, data = d)), "random effects")
  expect_error(se_ranef(aov(yield ~ Variety + Error(Block), data = d)), "strata")
})

test_that("prediction SEs equal SAS PROC MIXED Std Err Pred (unbalanced oats, dev/sas/se-ranef.sas)", {
  sas <- readRDS(test_path("_data", "se-ranef-sas.rds"))
  O <- data.frame(nlme::Oats)
  O$Block <- as.character(O$Block)
  O$Variety <- sub("Golden Rain", "GoldenRain", as.character(O$Variety))
  O$nitro <- format(O$nitro, nsmall = 1)
  drop <- with(O, (Block == "I" & Variety == "Marvellous") | (Block == "III" & Variety == "Victory" & nitro == "0.4") |
                  (Block == "V" & Variety == "GoldenRain" & nitro == "0.0") | (Block == "VI" & Variety == "Victory" & nitro == "0.6"))
  D <- O[!drop, c("Block", "Variety", "nitro", "yield")]
  D[] <- lapply(D, function(v) if (is.character(v)) factor(v) else v)
  # at SAS's REML variance estimates (Block 284.80, Block:Variety 110.97,
  # residual 171.54), so only the SE formula is compared, not the optimizers
  fm <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = D,
                   start = list(theta = sqrt(c(110.97, 284.80) / 171.54)), control = lme4::lmerControl(optimizer = NULL))
  r <- as.data.frame(se_ranef(fm))
  m <- merge(sas, data.frame(group = r$group, level = r$level, R_se = r$std.error, R_est = r$estimate),
             by = c("group", "level"))
  expect_equal(nrow(m), 23)
  expect_lt(max(abs(m$R_se - m$se)), 2e-4)     # SAS prints 4 decimals
  expect_lt(max(abs(m$R_est - m$est)), 5e-4)
  # at lme4's own REML optimum the variances differ in the 4th digit from SAS's; SEs agree to 2e-4 relative
  fr <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = D)
  rr <- as.data.frame(se_ranef(fr))
  m2 <- merge(sas, data.frame(group = rr$group, level = rr$level, R_se = rr$std.error), by = c("group", "level"))
  expect_lt(max(abs(m2$R_se / m2$se - 1)), 2e-4)
})
