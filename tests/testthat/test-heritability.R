# heritability() (R/heritability.R). References: the forms agree exactly for
# balanced data (theory); Cullis from Henderson's mixed model equations by hand;
# Piepho from emmeans's unadjusted covariance of the BLUEs; the standard form
# from the variance components. (inti::H2cal() 0.7.5 gave 0.837 for the balanced
# standard and Piepho forms, and 0.797 for the unbalanced Piepho form with its
# Kenward-Roger covariance; dev/NEWS-2.0-development.md.)

h2_data <- function() {
  set.seed(11)
  d <- expand.grid(gen = factor(sprintf("G%02d", 1:20)), blk = factor(1:4))
  g <- rnorm(20, 0, 2); b <- rnorm(4, 0, 1)
  d$y <- 10 + g[d$gen] + b[d$blk] + rnorm(nrow(d), 0, 2)
  d
}

test_that("balanced data: the four forms are equal; the standard form from the components", {
  skip_on_cran()
  d <- h2_data()
  m <- lme4::lmer(y ~ (1 | blk) + (1 | gen), d)
  h <- heritability(m, "gen", method = "all")
  expect_equal(h$estimate, rep(h$estimate[1], 4), tolerance = 1e-6)
  vc <- as.data.frame(lme4::VarCorr(m))
  s2g <- vc$vcov[vc$grp == "gen"]
  expect_equal(h$estimate[1], s2g / (s2g + sigma(m)^2 / 4), tolerance = 1e-10)
  expect_true(h$std.error[1] > 0 && h$conf.low[1] < h$estimate[1] && h$conf.high[1] > h$estimate[1])
  expect_s3_class(h, "pm_heritability")
  expect_output(print(h), "line-mean basis")
})

test_that("unbalanced data: Cullis = Henderson's equations by hand; Piepho = emmeans; Oakey = Cullis", {
  skip_on_cran()
  d <- h2_data()
  set.seed(2)
  du <- d[-sample(nrow(d), 15), ]
  m <- lme4::lmer(y ~ (1 | blk) + (1 | gen), du)
  h <- heritability(m, "gen", method = "all")
  X <- lme4::getME(m, "X"); Z <- as.matrix(lme4::getME(m, "Z")); s2 <- sigma(m)^2
  vc <- as.data.frame(lme4::VarCorr(m))
  q <- vapply(lme4::getME(m, "Ztlist"), nrow, 1)
  G <- diag(rep(vc$vcov[match(sub("\\..*", "", names(q)), vc$grp)], q))
  C <- solve(rbind(cbind(crossprod(X), crossprod(X, Z)), cbind(crossprod(Z, X), crossprod(Z) + s2 * solve(G)))) * s2
  gi <- ncol(X) + which(rownames(lme4::getME(m, "Zt")) %in% levels(du$gen))
  Cg <- C[gi, gi]; k <- nrow(Cg)
  vd <- 2 * (k * sum(diag(Cg)) - sum(Cg)) / (k * (k - 1))
  s2g <- vc$vcov[vc$grp == "gen"]
  expect_equal(h$estimate[h$method == "cullis"], 1 - vd / (2 * s2g), tolerance = 1e-8)
  expect_equal(h$estimate[h$method == "oakey"], h$estimate[h$method == "cullis"], tolerance = 1e-8)
  skip_if_not_installed("emmeans")
  mf <- lme4::lmer(y ~ 0 + gen + (1 | blk), du)
  se2 <- summary(graphics::pairs(emmeans::emmeans(mf, "gen", lmer.df = "satterthwaite")))$SE^2
  expect_equal(h$estimate[h$method == "piepho"], s2g / (s2g + mean(se2) / 2), tolerance = 1e-6)
})

test_that("several environments: the standard form divides by environments and plots per genotype", {
  skip_on_cran()
  set.seed(5)
  d <- expand.grid(gen = factor(1:15), env = factor(1:3), rep = factor(1:2))
  d$y <- rnorm(15, 0, 2)[d$gen] + rnorm(3, 0, 3)[d$env] + rnorm(45, 0, 1)[interaction(d$gen, d$env)] +
    rnorm(nrow(d))
  m <- lme4::lmer(y ~ env + (1 | env:rep) + (1 | gen) + (1 | gen:env), d)
  h <- heritability(m, "gen", method = c("standard"))
  vc <- as.data.frame(lme4::VarCorr(m))
  v <- function(g) vc$vcov[vc$grp == g]
  expect_equal(h$estimate, v("gen") / (v("gen") + v("gen:env") / 3 + sigma(m)^2 / 6), tolerance = 1e-10)
  expect_equal(unname(attr(h, "components")$n_terms), 3)
})

test_that("clear stops: no random genotype, other classes", {
  d <- h2_data()
  m <- suppressMessages(lme4::lmer(y ~ gen + (1 | blk), d))
  expect_error(heritability(m, "gen"), "random intercept")
  expect_error(heritability(lm(y ~ gen, d), "gen"), "lmer fits")
})
