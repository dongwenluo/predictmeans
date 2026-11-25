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
  h <- heritability(m, "gen", method = "all", nsim = 2000, seed = 1)
  expect_equal(h$estimate[1:4], rep(h$estimate[1], 4), tolerance = 1e-6)          # standard to Oakey
  expect_equal(heritability(m, "gen", method = "piepho", piepho = "refit")$estimate, h$estimate[1], tolerance = 1e-6)
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
  pr <- heritability(m, "gen", method = "piepho", piepho = "refit")$estimate      # re-estimated components
  expect_equal(pr, s2g / (s2g + mean(se2) / 2), tolerance = 1e-6)
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

test_that("family (narrow sense, single plant) and repeatability: the formulas, delta-method SEs, advice", {
  skip_on_cran()
  set.seed(3)
  f <- expand.grid(tree = 1:4, fam = factor(1:30), blk = factor(1:6))
  f$y <- rnorm(30, 0, 1)[f$fam] + rnorm(6, 0, 1)[f$blk] + rnorm(180, 0, 0.7)[interaction(f$fam, f$blk)] +
    rnorm(nrow(f), 0, 2)
  m <- lme4::lmer(y ~ (1 | blk) + (1 | fam) + (1 | fam:blk), f)
  vc <- as.data.frame(lme4::VarCorr(m)); v <- function(g) vc$vcov[vc$grp == g]
  p <- v("fam") + v("fam:blk") + sigma(m)^2                 # blocks are not in the single-plant variance
  hf <- heritability(m, "fam", method = "family", relationship = 1/4)
  expect_equal(hf$estimate, 4 * v("fam") / p, tolerance = 1e-10)
  expect_true(hf$std.error > 0 && hf$conf.low < hf$estimate && hf$conf.high > hf$estimate)
  expect_equal(heritability(m, "fam", method = "repeatability")$estimate, v("fam") / p, tolerance = 1e-10)
  expect_error(heritability(m, "fam", method = "family"), "relationship")
  # delta-method SE of the family form = 4 x that of the repeatability (a linear change of the same ratio)
  expect_equal(hf$std.error, 4 * heritability(m, "fam", method = "repeatability")$std.error, tolerance = 1e-6)
  expect_true("h2_above_one" %in% advice(heritability(m, "fam", method = "family", relationship = 0.02))$id)
})

test_that("bootstrap SEs of the mixed-model forms: the same for any ncore; advice without them", {
  skip_on_cran()
  d <- h2_data()
  set.seed(2)
  du <- d[-sample(nrow(d), 15), ]
  m <- lme4::lmer(y ~ (1 | blk) + (1 | gen), du)
  a <- heritability(m, "gen", method = "all", nboot = 20, nsim = 500, seed = 1)
  b <- heritability(m, "gen", method = "all", nboot = 20, nsim = 500, seed = 1, ncore = 2)
  expect_identical(a$std.error, b$std.error)
  expect_true(all(a$std.error > 0))
  expect_true(all(a$conf.low < a$estimate & a$conf.high > a$estimate))
  expect_equal(attr(a, "nboot"), 20)
  expect_true(all(c("h2_standard_unbalanced") %in% advice(a)$id))
  expect_true("h2_no_se" %in% advice(heritability(m, "gen"))$id)
})

test_that("Piepho and Möhring (2007), Example 1 (oats alpha design, John and Williams 1995): published values", {
  skip_on_cran()
  d <- read.csv(test_path("_data", "john-alpha.csv"), stringsAsFactors = TRUE)
  m <- lme4::lmer(yield ~ rep + (1 | rep:block) + (1 | gen), d)
  h <- heritability(m, "gen", method = "all", nsim = 2000, seed = 1)
  expect_equal(round(heritability(m, "gen", method = "piepho", piepho = "refit")$estimate, 3), 0.803)  # their Eq. 19
  expect_equal(round(h$estimate[h$method == "cullis"], 3), 0.809)     # Cullis, their Eq. 20
  expect_true("h2_standard_unbalanced" %in% advice(h)$id)               # incomplete blocks
  m2 <- lme4::lmer(yield ~ rep + (1 | gen), d)                          # blocks dropped: RCB
  expect_equal(round(heritability(m2, "gen", method = "all", nsim = 2000, seed = 1)$estimate[1:4], 3), rep(0.780, 4))
})

test_that("the standard, Cullis and Oakey forms equal those of the heritable package", {
  skip_on_cran()
  skip_if_not_installed("heritable")
  d <- read.csv(test_path("_data", "john-alpha.csv"), stringsAsFactors = TRUE)
  m <- lme4::lmer(yield ~ rep + (1 | rep:block) + (1 | gen), d)
  h <- heritability(m, "gen", method = c("all"), nsim = 2000, seed = 1)
  for (k in c("Cullis", "Oakey")) {
    expect_equal(h$estimate[h$method == tolower(k)], as.numeric(heritable::H2(m, target = "gen", method = k)),
                 tolerance = 1e-6)
  }
  m2 <- lme4::lmer(yield ~ rep + (1 | gen), d)
  expect_equal(heritability(m2, "gen", method = "standard")$estimate,
               as.numeric(heritable::H2(m2, target = "gen", method = "Standard")), tolerance = 1e-6)
})

test_that("Schmidt et al. (2019), Examples 1 and 2 (oats, john.alpha): every published value of their Figure 1", {
  skip_on_cran()
  d <- read.csv(test_path("_data", "john-alpha.csv"), stringsAsFactors = TRUE)
  ex1 <- heritability(lme4::lmer(yield ~ (1 | gen), d), "gen", method = "all", nsim = 1e5, seed = 1)
  expect_equal(round(ex1$estimate[1:4], 2), rep(0.58, 4))                     # standard, Cullis, Piepho, Oakey
  expect_equal(round(ex1$estimate[ex1$method == "reliability"], 3), 0.556)
  expect_equal(round(ex1$estimate[ex1$method == "simulated"], 3), 0.578)
  ex2 <- heritability(lme4::lmer(yield ~ rep + (1 | rep:block) + (1 | gen), d), "gen", method = "all",
                      nsim = 1e5, seed = 1)
  e <- setNames(round(ex2$estimate, 3), ex2$method)
  expect_equal(unname(e[c("cullis", "piepho", "oakey", "reliability", "simulated")]),
               c(0.809, 0.808, 0.809, 0.775, 0.804))
  expect_equal(round(ex2$estimate[ex2$method == "standard"], 2), 0.84)
})

