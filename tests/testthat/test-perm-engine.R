## The shared permutation engine (R/perm-engine.R).

test_that("permutation indices: chunks give the same stream, seeds leave the user's stream alone", {
  s <- pm_perm_stream(7, seed = 3)
  chunked <- cbind(s(2), s(3))
  whole <- pm_perm_stream(7, seed = 3)(5)
  expect_identical(chunked, whole)
  expect_identical(dim(whole), c(7L, 5L))
  expect_true(all(apply(whole, 2, function(p) setequal(p, 1:7))))
  # same draws as set.seed() + sample.int(), the 1.x code
  set.seed(3)
  expect_identical(whole, vapply(1:5, function(i) sample.int(7), integer(7)))
  # with a seed the user's stream is not touched
  set.seed(10); a <- runif(1)
  set.seed(10); invisible(pm_perm_stream(7, seed = 3)(4)); b <- runif(1)
  expect_identical(a, b)
  # without one it draws from the user's stream
  set.seed(10); x <- pm_perm_stream(7)(2)
  set.seed(10); expect_identical(x, vapply(1:2, function(i) sample.int(7), integer(7)))
  expect_identical(dim(pm_perm_stream(1, seed = 1)(3)), c(1L, 3L))
})

test_that("p-values: +1 rule, exact, rounding, NAs", {
  sims <- rbind(c(1, 2, 3, 4), c(0, 0, 0, 5))
  expect_equal(pm_perm_p(c(2.5, 1), sims), c(3, 2) / 5)
  expect_equal(pm_perm_p(c(2.5, 1), sims, exact = TRUE), c(2, 1) / 4)
  # vector of simulated values for one statistic
  expect_equal(pm_perm_p(2.5, c(1, 2, 3, 4)), 3 / 5)
  # noise below the rounding digits counts as a tie
  expect_equal(pm_perm_p(1e-9, c(0, 0, 0), digits = 4), 1)
  expect_equal(pm_perm_p(1e-9, c(0, 0, 0), digits = 12), 1 / 4)
  # NA simulated values are left out; an NA observed value gives NA
  expect_equal(pm_perm_p(2.5, c(1, NA, 3, 4)), 3 / 4)
  expect_identical(pm_perm_p(c(NA, 1), sims)[1], NA_real_)
  # no simulated values: p = 1
  expect_equal(pm_perm_p(1, numeric(0)), 1)
})

test_that("refits: errors become failures and are counted", {
  res <- pm_perm_map(1:5, function(k) if (k %% 2 == 0) stop("boom") else k, ncore = 1)
  d <- pm_perm_drop(res)
  expect_identical(unlist(d$res), c(1L, 3L, 5L))
  expect_identical(d$failures, 2L)
  expect_identical(pm_perm_failnote(0), "")
  expect_match(pm_perm_failnote(2), "2 refits failed")
  skip_on_os("windows")
  skip_on_cran()
  expect_identical(pm_perm_map(1:4, function(k) k^2, ncore = 2), as.list((1:4)^2))
})

test_that("|t| of contrasts from a list of refits", {
  refits <- list(list(coef = c(1, 2, 3), vcov = diag(3)),
                 list(coef = c(-1, 0, 2), vcov = diag(c(1, 4, 9))))
  rK <- rbind(c(1, -1, 0), c(0, 1, -1))
  t1 <- abs(rK %*% refits[[1]]$coef) / sqrt(diag(rK %*% refits[[1]]$vcov %*% t(rK)))
  t2 <- abs(rK %*% refits[[2]]$coef) / sqrt(diag(rK %*% refits[[2]]$vcov %*% t(rK)))
  expect_equal(pm_perm_tstats(refits, rK), cbind(as.vector(t1), as.vector(t2)))
  expect_identical(dim(pm_perm_tstats(refits, rK[1, , drop = FALSE])), c(1L, 2L))
})

test_that("permmodels(): seeds, response variables, failures", {
  data(Drug, package = "predictmeans", envir = environment())
  m1 <- glm(cbind(R, N - R) ~ Drug + log10Dose, binomial, Drug)
  Drug$F <- Drug$N - Drug$R
  m2 <- glm(cbind(R, F) ~ Drug + log10Dose, binomial, Drug)
  r1 <- permmodels(m1, nperm = 19, seed = 1, prt = FALSE, ncore = 1)
  r2 <- permmodels(m2, nperm = 19, seed = 1, prt = FALSE, ncore = 1)
  # both columns of cbind(R, F) are permuted together (1.x permuted only F)
  expect_identical(r1$ANOVA$Perm_p_value, r2$ANOVA$Perm_p_value)
  expect_null(attr(r1, "failures"))
  # a transformed response is permuted on its variable, with 'data'
  Oats <- as.data.frame(nlme::Oats)
  m3 <- lm(log(yield) ~ nitro + Variety, Oats)
  r3 <- permmodels(m3, nperm = 19, seed = 1, prt = FALSE, ncore = 1, data = Oats)
  expect_true(all(as.numeric(r3$ANOVA[1:2, "Perm_p_value"]) > 0))
  expect_error(permmodels(m3, nperm = 9, prt = FALSE, ncore = 1), "response variable")
  expect_error(permmodels(lm(cbind(yield, nitro) ~ Variety, Oats), nperm = 9, prt = FALSE, ncore = 1),
               "one response")
  # a seed leaves the user's random number stream as it was
  set.seed(5); a <- runif(1)
  set.seed(5); invisible(permmodels(m3, nperm = 9, seed = 2, prt = FALSE, ncore = 1, data = Oats)); b <- runif(1)
  expect_identical(a, b)
})

test_that("permlmer() and bootmodels() keep the user's stream with a seed", {
  skip_on_cran()
  fm1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  fm2 <- lme4::lmer(Reaction ~ Days + (Days || Subject), lme4::sleepstudy)
  set.seed(5); a <- runif(1)
  set.seed(5); p1 <- permlmer(fm1, fm2, nperm = 9, ncore = 1, seed = 3); b <- runif(1)
  expect_identical(a, b)
  expect_identical(p1$`Perm-p`[2], permlmer(fm1, fm2, nperm = 9, ncore = 1, seed = 3)$`Perm-p`[2])
  expect_null(attr(p1, "failures"))
})

test_that("refits work for models fitted in a function with formula and family variables", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fit_lmer <- function(f, d) lme4::lmer(f, data = d)
  m <- fit_lmer(yield ~ nitro + (1 | Block), Oats)
  r <- permmodels(m, nperm = 9, prt = FALSE, ncore = 1, seed = 1)
  expect_identical(length(r$permlist), 9L)
  data(Drug, package = "predictmeans", envir = environment())
  fit_glm <- function(f, fam, d) glm(f, family = fam, data = d)
  g <- fit_glm(cbind(R, N - R) ~ Drug + log10Dose, binomial, Drug)
  rg <- permmodels(g, nperm = 9, prt = FALSE, ncore = 1, seed = 1)
  expect_identical(length(rg$permlist), 9L)
})
