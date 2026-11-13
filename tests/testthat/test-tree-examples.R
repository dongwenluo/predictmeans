## Examples from "Experimental Design and Analysis for Tree Improvement"
## (Williams, Matheson and Harwood; data also in package eda4treeR), with the
## Genstat and SAS runs the maintainer supplied (2026-11-13; the output files
## are named in each test). Edge cases: nested treatments, a split plot with
## comparisons across strata, an incomplete site-by-seedlot table, families
## within a provenance as random effects. Genstat's printed values have 1-3
## decimals, SAS's 4-6.

tree_data <- function(f) read.csv(test_path("_data", "tree", f), stringsAsFactors = TRUE)

## Agreement with values printed to 'digits' decimals: within half a unit of
## the last printed digit (a tie such as 3.3425, printed 3.343, rounds either
## way with the last bits of floating point, which differ by platform).
expect_printed <- function(x, ref, digits, info = NULL) {
  ref <- rep_len(ref, length(x))
  ok <- abs(x - ref) <= 0.5 * 10^-digits + 1e-9
  expect(length(x) > 0 && all(ok), paste0("not within half a unit of the last printed digit: ",
                                          paste(format(x[!ok], digits = 10), "vs", ref[!ok], collapse = "; ")),
         info = info)
}
## the same for a set of printed values: each value is one of them, and each of them occurs
expect_printed_set <- function(x, refs, digits) {
  near <- outer(x, refs, function(a, b) abs(a - b) <= 0.5 * 10^-digits + 1e-9)
  expect(all(rowSums(near) > 0) && all(colSums(near) > 0),
         paste0("values ", paste(format(sort(unique(x)), digits = 10), collapse = ", "), " do not match ",
                paste(refs, collapse = ", ")))
}

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
  expect_printed(s$means$estimate, c(58.0, 52.3, 49.0, 40.7, 48.7, 59.7), 1)
  expect_printed(s$pairs$std.error, 4.03, 2)
  expect_equal(unique(s$pairs$df), 46)
  cc <- pm("contcomp")
  expect_printed(cc$means$estimate, c(2.0, 67.9), 1)
  expect_printed(cc$pairs$std.error, 2.68, 2)
  ct <- pm("contcomp:treat")
  expect_printed(ct$means$estimate, c(2.0, 56.9, 65.8, 80.9), 1)
  expect_printed(ct$pairs$std.error, 3.29, 2)
  cs <- pm("contcomp:seedlot")
  expect_printed(cs$means$estimate, c(0.0, 1.3, 6.7, 2.7, 0.0, 1.3, 77.3, 69.3, 63.1, 53.3, 64.9, 79.1), 1)
  expect_printed_set(cs$pairs$std.error, c(4.65, 6.58, 8.05), 2)   # max.rep, max-min, min.rep
})

test_that("Exam 4.4, split plot: Genstat's means, SEDs and df across strata (TABC4_3)", {
  d <- tree_data("ex4.4.csv")
  d[c("repl", "mainpl", "irrig", "fert", "seedlot")] <- lapply(d[c("repl", "mainpl", "irrig", "fert", "seedlot")], factor)
  m <- lme4::lmer(ht ~ irrig * fert * seedlot + (1 | repl) + (1 | repl:mainpl), data = d)
  pm <- function(term) predictmeans(m, term, plot = FALSE, prtnum = FALSE)
  ir <- pm("irrig")
  expect_printed(ir$means$estimate, c(9.95, 10.06), 2)
  expect_printed(ir$pairs$std.error, 0.296, 3)
  expect_equal(ir$pairs$df, 3, tolerance = 1e-3)
  sl <- pm("seedlot")
  expect_printed(sl$means$estimate, c(8.29, 10.18, 11.40, 10.15), 2)
  expect_printed(sl$pairs$std.error, 0.410, 3)
  expect_printed(sl$pairs$df, 12, 2)
  # comparisons across strata: Genstat 0.583 on 14.69 df (Satterthwaite); within an irrigation level 0.579 on 12
  is <- pm("irrig:seedlot")
  expect_printed(is$means$estimate, c(8.02, 10.06, 11.65, 10.06, 8.55, 10.30, 11.16, 10.24), 2)
  same <- substr(is$pairs$level1, 1, 4) == substr(is$pairs$level2, 1, 4)
  expect_printed(is$pairs$std.error[same], 0.579, 3)
  expect_printed(is$pairs$std.error[!same], 0.583, 3)
  expect_printed(is$pairs$df[same], 12, 2)
  expect_printed(is$pairs$df[!same], 14.69, 2)
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

test_that("Exam 3.1, trees within plots (multistratum aov): Genstat's means and SEDs on the plot stratum (TABC3_1)", {
  d <- tree_data("ex3.1a.csv")
  d$repl <- factor(d$repl); d$plot <- factor(d$plot)
  ref <- list(ht = list(m = c(10.29, 7.10, 5.51, 4.94, 2.73), sed = 0.861),
              dgl = list(m = c(18.87, 13.99, 11.43, 17.59, 5.95), sed = 1.074))
  for (v in names(ref)) {
    # dgl: the plot-stratum mean square (9.23) is below the tree-stratum one (19.47), a negative plot variance
    # by ANOVA; analysed by strata the seedlots are still compared with the plot stratum (REML, which fixes that
    # variance at zero, gave SED 1.5 on 75 df)
    fm <- aov(stats::as.formula(paste(v, "~ seedlot + Error(repl/plot)")), data = d)
    p <- suppressWarnings(predictmeans(fm, "seedlot", plot = FALSE, prtnum = FALSE))
    expect_printed(p$means$estimate, ref[[v]]$m, 2, info = v)
    expect_printed(p$pairs$std.error, ref[[v]]$sed, 3, info = v)
    expect_printed(p$pairs$df, 4, 2, info = v)            # the plot stratum
    expect_identical(p$settings$df_method, "multistratum")
  }
})

test_that("Exam 4.4 as aov(... + Error(repl/mainpl)): Genstat's SEDs and df by strata, no mixed model (TABC4_3)", {
  d <- tree_data("ex4.4.csv")
  d[c("repl", "mainpl")] <- lapply(d[c("repl", "mainpl")], factor)
  fa <- aov(ht ~ irrig * fert * seedlot + Error(repl/mainpl), data = d)
  pm <- function(term) predictmeans(fa, term, plot = FALSE, prtnum = FALSE)
  ir <- pm("irrig")
  expect_printed(ir$means$estimate, c(9.95, 10.06), 2)
  expect_printed(ir$pairs$std.error, 0.296, 3)
  expect_equal(ir$pairs$df, 3, tolerance = 1e-8)
  sl <- pm("seedlot")
  expect_printed(sl$pairs$std.error, 0.410, 3)
  expect_equal(unique(round(sl$pairs$df, 8)), 12)
  is <- pm("irrig:seedlot")
  same <- substr(is$pairs$level1, 1, 4) == substr(is$pairs$level2, 1, 4)
  expect_printed(is$pairs$std.error[same], 0.579, 3)
  expect_printed(is$pairs$std.error[!same], 0.583, 3)
  expect_equal(unique(round(is$pairs$df[same], 8)), 12)
  expect_printed(is$pairs$df[!same], 14.69, 2)
  # the ANOVA F tests of the strata
  a <- summary(fa)
  wp <- a[["Error: repl:mainpl"]][[1]]
  expect_equal(pm_joint(fa, "fert")$F.ratio, wp["fert", "F value"], tolerance = 1e-8)
  expect_equal(pm_joint(fa, "fert")$df2, 3, tolerance = 1e-8)
  sp <- a[["Error: Within"]][[1]]
  expect_equal(pm_joint(fa, "seedlot")$F.ratio, sp[grep("^seedlot", rownames(sp))[1], "F value"], tolerance = 1e-8)
  expect_identical(attr(pm_joint(fa, "seedlot"), "df_method"), "multistratum")
  expect_identical(attr(contrastmeans(fa, "seedlot", ctrmatrix = rbind(c(1, -1, 0, 0)))[[1]], "Note") != "", TRUE)
})

test_that("Exam 8.1, inoculation (on whole replicates) x provenances within countries: Genstat by strata (TABC8_1)", {
  d <- tree_data("ex8.1a.csv")
  d$repl <- factor(d$repl); d$prov <- factor(d$prov)
  fa <- aov(dbh ~ inoc * (country/prov) + Error(repl), data = d)
  pm <- function(term) expect_no_warning(predictmeans(fa, term, plot = FALSE, prtnum = FALSE))   # nesting: no warning
  i <- pm("inoc")
  expect_printed(i$means$estimate, c(3.183, 3.625), 3)
  expect_printed(i$pairs$std.error, 0.1306, 4)                         # replicate stratum, 2 df
  expect_equal(i$pairs$df, 2, tolerance = 1e-8)
  cn <- pm("country")
  expect_printed(cn$means$estimate[1:5], c(2.631, 3.343, 3.686, 2.498, 2.612), 3)
  sed <- function(p, a, b) p$pairs$std.error[(p$pairs$level1 == a & p$pairs$level2 == b) | (p$pairs$level1 == b & p$pairs$level2 == a)]
  expect_printed(sed(cn, "Benin", "Guam"), 0.5503, 4)                 # 4 and 4 trees: Genstat's min.rep
  expect_printed(sed(cn, "Benin", "Malaysia"), 0.4102, 4)             # 4 and 36: max-min
  expect_equal(unique(round(cn$pairs$df, 8)), 116)
  ic <- pm("inoc:country")
  pdf <- function(a, b) ic$pairs$df[(ic$pairs$level1 == a & ic$pairs$level2 == b) | (ic$pairs$level1 == b & ic$pairs$level2 == a)]
  expect_printed(sed(ic, "1 week:Benin", "7 weeks:Guam"), 0.7826, 4)  # across inoculation: both strata
  expect_printed(pdf("1 week:Benin", "7 weeks:Guam"), 117.16, 2)
  expect_printed(sed(ic, "1 week:Benin", "7 weeks:Malaysia"), 0.5859, 4)
  expect_printed(pdf("1 week:Benin", "7 weeks:Malaysia"), 110.85, 2)
  expect_printed(sed(ic, "1 week:Malaysia", "7 weeks:Malaysia"), 0.2722, 4)
  expect_printed(pdf("1 week:Malaysia", "7 weeks:Malaysia"), 31.61, 2)
  expect_printed(sed(ic, "1 week:Benin", "1 week:Guam"), 0.7782, 4)   # same inoculation: 116 df
  expect_printed(pdf("1 week:Benin", "1 week:Guam"), 116, 2)
  cp <- pm("country:prov")
  expect_equal(nrow(cp$means), 59)
  expect_printed(cp$pairs$std.error, 0.5503, 4)
})

test_that("Exam 8.2, clone trial with standards (REML; random rows, columns and clones): SAS DDFM=KR and Genstat (TABD8_5, TABC8_5)", {
  d <- tree_data("ex8.2.csv")
  for (v in c("repl", "row", "column", "standard", "clone")) d[[v]] <- factor(d[[v]])
  d$contcompv <- d$contcompf                                   # clone effects for clones only (Genstat contcompv.clone)
  d$contcompf <- factor(d$contcompf)
  m <- suppressMessages(lme4::lmer(dbhmean ~ repl + column + contcompf + contcompf:standard + (1 | repl:row) +
                                     (1 | repl:column) + (0 + contcompv | clone), data = d))
  vc <- as.data.frame(lme4::VarCorr(m))
  vv <- setNames(vc$vcov, vc$grp)
  expect_equal(unname(vv[c("repl:row", "repl:column", "clone", "Residual")]), c(0.08023, 0.05284, 0.4950, 0.3992),
               tolerance = 2e-3)
  p <- expect_no_warning(predictmeans(m, "contcompf:standard", plot = FALSE, prtnum = FALSE))
  expect_equal(p$means$estimate, c(7.7027, 6.0842, 6.5544, 8.9710, 8.9076), tolerance = 2e-5)   # 4 standards, clones
  expect_equal(p$means$std.error, c(0.3146, 0.3125, 0.3134, 0.3137, 0.1168), tolerance = 5e-4)
  expect_equal(p$means$df, c(221, 224, 224, 223, 63.6), tolerance = 4e-3)
  cf <- predictmeans(m, "contcompf", plot = FALSE, prtnum = FALSE)
  expect_equal(cf$means$estimate, c(7.3281, 8.9076), tolerance = 2e-5)       # standards averaged with equal weights
  expect_equal(cf$means$std.error, c(0.1648, 0.1168), tolerance = 5e-4)
  expect_equal(cf$means$df, c(215, 63.6), tolerance = 4e-3)
  # Kenward-Roger F tests (SAS Type 1; replicates first, standards within the standard group last)
  expect_equal(unlist(pm_joint(m, "repl")[c("F.ratio", "df2")]), c(F.ratio = 2.05, df2 = 26.5), tolerance = 5e-3)
  st <- rbind(c(1, -1, 0, 0, 0), c(1, 0, -1, 0, 0), c(1, 0, 0, -1, 0))
  colnames(st) <- paste(p$means$contcompf, p$means$standard, sep = ":")    # the five existing means
  js <- pm_joint(m, "contcompf:standard", ctrmatrix = st)
  expect_equal(c(js$F.ratio, js$df2), c(16.86, 207), tolerance = 5e-3)
})
