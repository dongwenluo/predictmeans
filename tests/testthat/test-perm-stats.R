## permstats() and permcombine(): uses of the permuted statistics (plan step 14)

cars_data <- function() {
  d <- mtcars
  d$cyl <- factor(d$cyl)
  d
}

test_that("permstats(): parametric and permutation p-values, effect sizes, size of the parametric test", {
  d <- cars_data()
  fm <- lm(mpg ~ wt + cyl, data = d)
  set.seed(7)
  P <- replicate(199, sample(nrow(d)))
  r <- permmodels(fm, perms = P, type = "II", prt = FALSE)
  s <- permstats(r)
  expect_equal(rownames(s), c("wt", "cyl", "coef:wt", "coef:cyl6", "coef:cyl8"))
  dr <- drop1(fm, test = "F")
  expect_equal(s[c("wt", "cyl"), "p.param"], dr[c("wt", "cyl"), "Pr(>F)"])
  expect_equal(s[c("coef:wt", "coef:cyl6"), "p.param"], unname(coef(summary(fm))[c("wt", "cyl6"), "Pr(>|t|)"]))
  expect_equal(s[c("wt", "cyl"), "p.perm"], as.numeric(r$ANOVA[c("wt", "cyl"), "Perm_p_value"]))
  expect_equal(s$df2, rep(df.residual(fm), 5))
  # by hand for cyl: Freedman-Lane F from the model without cyl
  red <- lm(mpg ~ wt, data = d)
  Fs <- apply(P, 2, function(i) {
    ys <- fitted(red) + residuals(red)[i]
    anova(lm(ys ~ wt, data = d), lm(ys ~ wt + cyl, data = d))$F[2]
  })
  Fo <- s["cyl", "statistic"]
  expect_equal(s["cyl", "ses"], (Fo - mean(Fs)) / sd(Fs), tolerance = 1e-6)
  expect_equal(s["cyl", "param.size"], mean(pf(Fs, 2, df.residual(fm), lower.tail = FALSE) <= 0.05), tolerance = 1e-6)
  expect_equal(permstats(r, level = 0.5)["cyl", "param.size"],
               mean(pf(Fs, 2, df.residual(fm), lower.tail = FALSE) <= 0.5), tolerance = 1e-6)
  expect_error(permstats(permmodels(fm, nperm = 9, method = "raw", prt = FALSE)), "Freedman-Lane")
})

test_that("permstats(): df of the parametric tests of mixed models", {
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fl <- lme4::lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Oats)
  s <- permstats(permmodels(fl, nperm = 19, type = "III", seed = 1, prt = FALSE, refit = FALSE))
  at <- anova(lmerTest::as_lmerModLmerTest(fl), type = "III")
  expect_equal(s[c("nitro", "Variety"), "df2"], at[c("nitro", "Variety"), "DenDF"], tolerance = 1e-4)
  ct <- coef(summary(lmerTest::as_lmerModLmerTest(fl)))
  expect_equal(s["coef:nitro0.2", "df2"], ct["nitro0.2", "df"], tolerance = 1e-4)
  fe <- nlme::lme(yield ~ nitro + Variety, random = ~ 1 | Block/Variety, data = Oats)
  se <- permstats(permmodels(fe, nperm = 19, seed = 1, prt = FALSE, refit = FALSE))
  expect_equal(se[c("nitro", "Variety"), "df2"], anova(fe)[c("nitro", "Variety"), "denDF"])
  # refits: the p-values of the table
  rr <- permmodels(fl, nperm = 9, seed = 1, prt = FALSE)
  expect_equal(permstats(rr)[c("nitro", "Variety"), "p.perm"], rr$ANOVA[c("nitro", "Variety"), "Perm_p_value"])
})

test_that("permcombine(): Fisher and Tippett combination and min-p adjustment by hand", {
  d <- cars_data()
  set.seed(8)
  P <- replicate(99, sample(nrow(d)))
  r1 <- permmodels(lm(mpg ~ wt + cyl, data = d), perms = P, type = "II", prt = FALSE)
  r2 <- permmodels(lm(qsec ~ wt + cyl, data = d), perms = P, type = "II", prt = FALSE)
  cf <- permcombine(mpg = r1, qsec = r2)
  expect_equal(rownames(cf), rownames(permstats(r1)))
  expect_equal(cf["cyl", "mpg"], as.numeric(r1$ANOVA["cyl", "Perm_p_value"]))
  # by hand for cyl
  S <- sapply(list(r1, r2), function(r) {
    st <- r$perm$stats
    i <- match("cyl", st$names)
    round(c(st$obs[i], st$T[i, ]), 6)
  })
  Pm <- apply(S, 2, function(v) sapply(v, function(x) mean(v >= x)))
  C <- -2 * rowSums(log(Pm))
  expect_equal(cf["cyl", "combined"], (sum(C[-1] >= C[1]) + 1) / 100)
  minp <- apply(Pm[-1, ], 1, min)
  expect_equal(cf["cyl", "adj.mpg"], (sum(minp <= Pm[1, 1]) + 1) / 100)
  expect_true(all(cf$adj.mpg >= cf$mpg - 1e-12) && all(cf$adj.qsec >= cf$qsec - 1e-12))
  ct <- permcombine(list(mpg = r1, qsec = r2), method = "tippett", tests = "cyl")
  expect_equal(ct["cyl", "combined"], (sum(-apply(Pm[-1, ], 1, min) >= -min(Pm[1, ])) + 1) / 100)
  # different permutations are refused
  r3 <- permmodels(lm(qsec ~ wt + cyl, data = d), nperm = 99, seed = 1, type = "II", prt = FALSE)
  expect_error(permcombine(r1, r3), "same permutations")
  expect_error(permcombine(r1, r2, tests = "carb"), "not in every result")
})
