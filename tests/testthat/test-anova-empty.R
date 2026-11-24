## Type III tests with empty cells (aliased coefficients): hypotheses on the
## cells with data (pm_type3_map_cells(), R/anova.R); the same contrasts as
## dev/sas/pm-reference-emptycells.sas.

oats_empty <- function() {
  O <- as.data.frame(nlme::Oats); O$N <- factor(O$nitro)
  list(full = O, empty = subset(O, !(Variety == "Victory" & N == "0.6")))
}
# the contrasts of the SAS script, on the 12 N:Variety means (Victory:0.6 not used)
sas_contrasts <- function(levs) {
  cell <- function(v, n) as.numeric(levs == paste(n, v, sep = ":"))
  G <- function(n) cell("Golden Rain", n); M <- function(n) cell("Marvellous", n); V <- function(n) cell("Victory", n)
  list(Variety = rbind(G("0") + G("0.2") + G("0.4") - V("0") - V("0.2") - V("0.4"),
                       M("0") + M("0.2") + M("0.4") - V("0") - V("0.2") - V("0.4")),
       N = rbind(G("0") - G("0.6") + M("0") - M("0.6"), G("0.2") - G("0.6") + M("0.2") - M("0.6"),
                 G("0.4") - G("0.6") + M("0.4") - M("0.6")),
       NxV = rbind(G("0") - G("0.6") - M("0") + M("0.6"), G("0.2") - G("0.6") - M("0.2") + M("0.6"),
                   G("0.4") - G("0.6") - M("0.4") + M("0.6"), G("0") - G("0.4") - V("0") + V("0.4"),
                   G("0.2") - G("0.4") - V("0.2") + V("0.4")))
}

test_that("anova() of an lmer fit with an empty cell: each term on the cells with data, = pm_joint()", {
  skip_on_cran()
  q <- function(e) suppressWarnings(suppressMessages(e))
  d <- oats_empty()
  fm <- q(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = d$empty))
  a <- anova(fm)
  expect_equal(rownames(a), c("N", "Variety", "N:Variety"))
  expect_equal(a$Df, c(3, 2, 5))                              # 6 - 1 interaction df: one empty cell
  expect_match(paste(attr(a, "heading"), collapse = " "), "Empty cells")
  levs <- rownames(q(Kmatrix(fm, "N:Variety"))$K)
  cs <- sas_contrasts(levs)
  for (k in 1:3) {
    term <- c("Variety", "N", "N:Variety")[k]
    jt <- q(pm_joint(fm, "N:Variety", ctrmatrix = cs[[k]]))
    expect_equal(a[term, "F"], jt$F.ratio, tolerance = 1e-8)
    expect_equal(a[term, "Df.res"], jt$df2, tolerance = 1e-6)
    expect_equal(a[term, "Df"], jt$df1)
  }
  # main effects: the common-cell tests of pm_joint() on the term's own means
  expect_equal(a["Variety", "F"], q(pm_joint(fm, "Variety"))$F.ratio, tolerance = 1e-8)
  expect_equal(a["N", "F"], q(pm_joint(fm, "N"))$F.ratio, tolerance = 1e-8)
  # Satterthwaite too
  s <- anova(fm, ddf = "Satterthwaite")
  expect_equal(s$Df, c(3, 2, 5))
  expect_true(all(is.finite(s$F)))
})

test_that("the hypotheses on the cells equal the sum-to-zero Type III ones when no cell is empty", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  ff <- q(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = oats_empty()$full))
  a <- anova(ff)
  m <- pm_type3_map_cells(ff, names(lme4::fixef(ff)))
  PhiA <- pm_vcov_kr(ff)
  for (j in 1:3) {
    kr <- pm_kr_ftest(ff, m$hyp[[j]], PhiA)
    expect_equal(kr[["stat"]], a[j + 1, "F"], tolerance = 1e-8)
    expect_equal(kr[["ddf"]], a[j + 1, "Df.res"], tolerance = 1e-6)
  }
})

test_that("lm with an empty cell: Type III interaction = the model comparison; permmodels and glm", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  e <- oats_empty()$empty
  fl <- lm(yield ~ Block + N * Variety, data = e)
  a <- pm_anova(fl, "III")
  cmpr <- stats::anova(lm(yield ~ Block + N + Variety, data = e), fl)
  expect_equal(a["N:Variety", "Sum Sq"], cmpr[2, "Sum of Sq"], tolerance = 1e-8)
  expect_equal(a["N:Variety", "Df"], 5)
  expect_equal(a["Variety", "F value"], q(pm_joint(fl, "Variety"))$F.ratio, tolerance = 1e-8)
  pr <- q(permmodels(fl, nperm = 19, type = "III", seed = 1, ncore = 1))
  expect_true(all(c("N", "Variety", "N:Variety") %in% rownames(pr$ANOVA)))
  # glm: Wald tests on the cells; likelihood ratio tests need a complete design
  g <- glm(round(yield) ~ N * Variety, family = poisson, data = e)
  expect_equal(pm_anova(g, "III", test_statistic = "Wald")$Df, c(3, 2, 5))
  expect_error(pm_anova(g, "III"), "complete design")
})

test_that("empty cell: equal to the SAS run (dev/sas/sas-results-2026-11-18-emptycells.txt)", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  e <- oats_empty()$empty
  # PROC GLM, Block fixed (= SAS Type IV for the factorial model)
  a <- pm_anova(lm(yield ~ Block + N * Variety, data = e), "III")
  expect_equal(a[c("Variety", "N", "N:Variety"), "Sum Sq"], c(1628.11111, 12077.72917, 296.22917), tolerance = 1e-8)
  expect_equal(a["Residuals", "Sum Sq"], 11955.75758, tolerance = 1e-8)
  expect_equal(round(a[c("Variety", "N", "N:Variety"), "F value"], 2), c(3.40, 16.84, 0.25))
  # PROC MIXED, DDFM=KR and Satterthwaite (the same values at SAS's printed
  # digits); mixed-model references to 2%: lme4's optimum (Block 205.14,
  # Block:Variety 79.84, residual 179.80) and SAS's (205.18, 79.76, 179.82)
  # have the same -2 REML log-likelihood, 484.539
  fm <- q(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = e))
  for (ddf in c("Kenward-Roger", "Satterthwaite")) {
    m <- anova(fm, ddf = ddf)
    expect_equal(m[c("Variety", "N", "N:Variety"), "F"], c(1.94, 22.39, 0.33), tolerance = 0.02)
    expect_equal(m[c("Variety", "N", "N:Variety"), "Df.res"], c(11.9, 40.1, 40.1), tolerance = 0.02)
    expect_equal(m[c("Variety", "N:Variety"), "Pr(>F)"], c(0.1864, 0.8922), tolerance = 0.02)
  }
})

test_that("empty cell: the interaction test equals Genstat's REML Wald test (dev/genstat run 2, pattern B)", {
  e <- oats_empty()$empty
  fm <- suppressWarnings(suppressMessages(lme4::lmer(yield ~ N * Variety + (1 | Block/Variety), data = e)))
  m <- anova(fm)
  # Genstat: Variety.Nitro F 0.33 on 5 and 40.0 df, p 0.892 (dropping it from the full model); Genstat gives
  # only sequential tests for the main effects (marginality), so SAS is the reference for those
  expect_equal(m["N:Variety", "F"], 0.33, tolerance = 0.02)
  expect_equal(m["N:Variety", "Df"], 5)
  expect_equal(m["N:Variety", "Df.res"], 40.0, tolerance = 0.02)
  expect_equal(m["N:Variety", "Pr(>F)"], 0.892, tolerance = 0.02)
})
