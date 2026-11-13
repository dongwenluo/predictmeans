## Examples from "Experimental Design and Analysis for Tree Improvement"
## (Williams, Matheson and Harwood; data also in package eda4treeR), with the
## Genstat and SAS runs the maintainer supplied (2026-11-13; the output files
## are named in each test). Edge cases: nested treatments, a split plot with
## comparisons across strata, an incomplete site-by-seedlot table, families
## within a provenance as random effects. Genstat's printed values have 1-3
## decimals, SAS's 4-6.

tree_data <- function(f) read.csv(test_path("_data", "tree", f), stringsAsFactors = TRUE)

test_that("Exam 4.3, treatments nested in control/treated: Genstat's means and SEDs (TABC4_1)", {
  d <- tree_data("ex4.3.csv")
  d$percent <- d$count * 100 / 25
  d$treat <- factor(d$treat, levels = c("control", "nick", "bw&s", "bw1min"))
  d$contcomp <- factor(ifelse(d$treat == "control", "control", "Treated"), levels = c("control", "Treated"))
  d[c("repl", "seedlot")] <- lapply(d[c("repl", "seedlot")], factor)
  fm <- aov(percent ~ repl + contcomp + seedlot + treat/contcomp + contcomp/seedlot + treat/contcomp/seedlot, data = d)
  pm <- function(term) suppressWarnings(predictmeans(fm, term, plot = FALSE, prtnum = FALSE))
  # Genstat: averages over the treatments that exist, equal weights (SAS and emmeans give the control half the
  # weight of a seedlot mean: 38.67 for seedlot 18211)
  s <- pm("seedlot")
  expect_equal(round(s$means$estimate, 1), c(58.0, 52.3, 49.0, 40.7, 48.7, 59.7))
  expect_equal(unique(round(s$pairs$std.error, 2)), 4.03)
  expect_equal(unique(s$pairs$df), 46)
  cc <- pm("contcomp")
  expect_equal(round(cc$means$estimate, 1), c(2.0, 67.9))
  expect_equal(round(cc$pairs$std.error, 2), 2.68)
  ct <- pm("contcomp:treat")
  expect_equal(round(ct$means$estimate, 1), c(2.0, 56.9, 65.8, 80.9))
  expect_equal(unique(round(ct$pairs$std.error, 2)), 3.29)
  cs <- pm("contcomp:seedlot")
  expect_equal(round(cs$means$estimate, 1),
               c(0.0, 1.3, 6.7, 2.7, 0.0, 1.3, 77.3, 69.3, 63.1, 53.3, 64.9, 79.1))
  expect_equal(sort(unique(round(cs$pairs$std.error, 2))), c(4.65, 6.58, 8.05))   # max.rep, max-min, min.rep
})

test_that("Exam 4.4, split plot: Genstat's means, SEDs and df across strata (TABC4_3)", {
  d <- tree_data("ex4.4.csv")
  d[c("repl", "mainpl", "irrig", "fert", "seedlot")] <- lapply(d[c("repl", "mainpl", "irrig", "fert", "seedlot")], factor)
  m <- lme4::lmer(ht ~ irrig * fert * seedlot + (1 | repl) + (1 | repl:mainpl), data = d)
  pm <- function(term) predictmeans(m, term, plot = FALSE, prtnum = FALSE)
  ir <- pm("irrig")
  expect_equal(round(ir$means$estimate, 2), c(9.95, 10.06))
  expect_equal(round(ir$pairs$std.error, 3), 0.296)
  expect_equal(ir$pairs$df, 3, tolerance = 1e-3)
  sl <- pm("seedlot")
  expect_equal(round(sl$means$estimate, 2), c(8.29, 10.18, 11.40, 10.15))
  expect_equal(unique(round(sl$pairs$std.error, 3)), 0.410)
  expect_equal(unique(round(sl$pairs$df, 2)), 12)
  # comparisons across strata: Genstat 0.583 on 14.69 df (Satterthwaite); within an irrigation level 0.579 on 12
  is <- pm("irrig:seedlot")
  expect_equal(round(is$means$estimate, 2), c(8.02, 10.06, 11.65, 10.06, 8.55, 10.30, 11.16, 10.24), tolerance = 0.006)
  same <- substr(is$pairs$level1, 1, 4) == substr(is$pairs$level2, 1, 4)
  expect_equal(unique(round(is$pairs$std.error[same], 3)), 0.579)
  expect_equal(unique(round(is$pairs$std.error[!same], 3)), 0.583)
  expect_equal(unique(round(is$pairs$df[same], 2)), 12)
  expect_equal(unique(round(is$pairs$df[!same], 2)), 14.69)
})

test_that("Exam 5.2, incomplete site-by-seedlot table: SAS LS-means and SEs (TABD5_2)", {
  e <- tree_data("ex5.2.csv")
  e$site <- factor(e$site); e$seedlot <- factor(e$seedlot)
  f <- lm(ht ~ site + seedlot, data = e)
  expect_equal(nobs(f), 192)                       # 30 of the 222 site-by-seedlot cells are empty
  ref <- read.csv(test_path("_data", "tree", "ex5.2-sas-lsmeans.csv"), colClasses = c("character", "character",
                                                                                     "numeric", "numeric"))
  for (eff in c("site", "seedlot")) {
    p <- predictmeans(f, eff, plot = FALSE, prtnum = FALSE)$means
    r <- ref[ref$effect == eff, ]
    i <- match(r$level, as.character(p[[eff]]))
    expect_false(anyNA(i))
    expect_true(all(p$estimable))
    expect_equal(p$estimate[i], r$estimate, tolerance = 1e-7, info = eff)
    expect_equal(p$std.error[i], r$std.error, tolerance = 1e-6, info = eff)
  }
})

test_that("Exam 6.2, families within a provenance (random): SAS PROC MIXED (TABD6_1)", {
  d <- tree_data("ex6.2a.csv")
  d <- droplevels(d[!d$prov %in% c("Qld", "Sabah"), ])
  d$repl <- factor(d$repl); d$family <- factor(d$family)
  m <- lme4::lmer(mean ~ repl + (1 | family), data = d)
  vc <- as.data.frame(lme4::VarCorr(m))
  expect_equal(vc$vcov, c(0.2584, 1.1667), tolerance = 2e-4)       # family(prov), Residual
  # SAS's default df here are containment (120), with the unadjusted covariance
  p <- predictmeans(m, "repl", df_method = "containment", plot = FALSE, prtnum = FALSE)
  expect_equal(p$means$estimate, c(14.3883, 14.2327, 13.7788, 14.4402), tolerance = 5e-5)
  expect_equal(p$means$std.error, rep(0.1864, 4), tolerance = 5e-4)
  expect_equal(unique(p$means$df), 120)
  # the provenance mean: 14.2100, SE 0.1158 on 40 df (the families' containment df)
  mi <- lme4::lmer(mean ~ 1 + repl + (1 | family), data = d)
  K <- matrix(c(1, rep(1/4, 3)), 1)
  expect_equal(drop(K %*% lme4::fixef(mi)), 14.2100, tolerance = 5e-5)
  expect_equal(sqrt(drop(K %*% as.matrix(vcov(mi)) %*% t(K))), 0.1158, tolerance = 5e-4)
})
