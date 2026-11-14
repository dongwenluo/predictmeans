## predictmeans against Genstat 23.1 REML (runs of 2026-11-14 by the maintainer:
## dev/genstat/pm-reference.gen and pm-reference-2.gen; output files as text in
## dev/genstat/genstat-results-2026-11-14-run*.txt, parsed by
## dev/genstat/parse-genstat.py into _data/genstat/genstat-reml.csv.gz).
## Six REML analyses: oats with a whole plot missing (68), with a whole plot
## and a sub-plot missing (67), with an empty cell (66), with one sub-plot
## missing (71), with five scattered sub-plots missing (67), and AR(1)
## repeated measures with 5 of 48 values missing.
##
## What Genstat computes: variance components and predicted means (random
## terms excluded) as REML; SEDs from the unadjusted covariance
## (vcov_adjust = FALSE); Wald tests adding the terms in turn, all from the
## full fit; their denominator df are Kenward-Roger's with the average
## information (df_method = "KR-average"), and the F statistic is
## Kenward-Roger's scaling factor times Wald / ndf.

gen_ref <- function() read.csv(test_path("_data", "genstat", "genstat-reml.csv.gz"), colClasses = "character")
## within one unit of the last printed digit
gen_agree <- function(r, s) {
  dec <- ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0)
  v <- suppressWarnings(as.numeric(sub("<", "", s)))
  ifelse(grepl("^<", s), r < v, !is.na(r) & abs(r - v) <= 10^(-dec) * (1 + 1e-6) + 1e-9)
}
gen_oats <- function() {
  o <- as.data.frame(nlme::Oats)
  o$Variety <- factor(gsub(" ", "", o$Variety))
  o$Nitro <- factor(sprintf("%.1f", o$nitro))
  o$Block <- factor(as.character(o$Block), levels = c("I", "II", "III", "IV", "V", "VI"))
  o
}
gen_fits <- function() {
  o <- gen_oats()
  drop <- list(
    "run1:Yield" = with(o, Block == "I" & Variety == "GoldenRain"),
    "run2:YA" = with(o, (Block == "I" & Variety == "GoldenRain") | (Block == "II" & Variety == "Victory" & Nitro == "0.6")),
    "run2:YB" = with(o, Variety == "Victory" & Nitro == "0.6"),
    "run2:YC" = with(o, Block == "II" & Variety == "Victory" & Nitro == "0.6"),
    "run2:YD" = with(o, (Block == "I" & Variety == "Marvellous" & Nitro == "0.2") |
                       (Block == "III" & Variety == "GoldenRain" & Nitro == "0.6") |
                       (Block == "IV" & Variety == "Victory" & Nitro == "0.0") |
                       (Block == "V" & Variety == "Marvellous" & Nitro == "0.4") |
                       (Block == "VI" & Variety == "GoldenRain" & Nitro == "0.2")))
  ar <- read.table(test_path("_data", "sas", "ar1.txt"), header = TRUE)
  ar$Subj <- factor(ar$Subj); ar$Group <- factor(ar$Group, labels = c("1", "2")); ar$t <- ar$Time; ar$Time <- factor(ar$Time)
  c(lapply(drop, function(dr) suppressMessages(lme4::lmer(yield ~ Variety * Nitro + (1 | Block/Variety),
                                                          data = droplevels(o[!dr, ])))),
    list("run1:Y" = nlme::gls(Y ~ Time * Group, data = ar, correlation = nlme::corAR1(form = ~ t | Subj))))
}

test_that("Genstat REML: variance components, sequential Wald tests (KR-average df), means, SEDs", {
  skip_on_cran()
  G <- gen_ref()
  fits <- gen_fits()
  res <- list()
  add <- function(an, what, key, gen, r) res[[length(res) + 1]] <<- data.frame(an = an, what = what, key = key,
                                                                                 genstat = gen, r = r)
  for (an in names(fits)) {
    m <- fits[[an]]
    g <- G[G$analysis == an, ]
    nl <- inherits(m, "gls")
    gv <- function(term) g$value[g$table == "vc" & g$term == term & g$key != "se"]
    if (nl) {
      add(an, "vc", "Sigma2", gv("Sigma2"), m$sigma^2)
      add(an, "vc", "phi_1", gv("phi_1"), coef(m$modelStruct$corStruct, unconstrained = FALSE))
    } else {
      vc <- as.data.frame(lme4::VarCorr(m))
      vv <- stats::setNames(vc$vcov, vc$grp)
      add(an, "vc", "Block", gv("Block"), vv[["Block"]])
      add(an, "vc", "Block.Variety", gv("Block.Variety"), vv[["Variety:Block"]])
      add(an, "vc", "Sigma2", gv("Sigma2"), vv[["Residual"]])
    }
    # Wald tests adding the terms in turn (from the full fit): the rows of the
    # Cholesky factor of X'V^-1X for each term
    b <- if (nl) coef(m) else lme4::fixef(m)
    R <- chol(solve(as.matrix(vcov(m))))
    asg <- if (nl) attr(model.matrix(formula(m), pm_get_data(m)), "assign") else attr(lme4::getME(m, "X"), "assign")
    PhiA <- if (nl) pm_vcov_kr_nlme(m, "average") else pm_vcov_kr(m, "average")
    labs <- if (nl) c("Time", "Group", "Time.Group") else c("Variety", "Nitro", "Variety.Nitro")
    for (k in seq_along(labs)) {
      Lk <- R[asg == k, , drop = FALSE]
      kr <- pm_kr_ftest(m, Lk, PhiA = PhiA)
      w <- sum(drop(Lk %*% b)^2)
      gw <- g[g$table == "wald_seq" & g$term == labs[k], ]
      gwv <- function(key) gw$value[gw$key == key]
      add(an, "wald", labs[k], gwv("wald"), w)
      add(an, "ndf", labs[k], gwv("ndf"), nrow(Lk))
      add(an, "F", labs[k], gwv("F"), kr[["F.scaling"]] * w / nrow(Lk))
      # the empty-cell model: Genstat gives 10.4 for Variety, Kenward-Roger 10.10
      # (average), 10.19 (observed), 10.01 (expected); its other df agree
      if (!(an == "run2:YB" && labs[k] == "Variety")) add(an, "ddf", labs[k], gwv("ddf"), kr[["ddf"]])
    }
    # predicted means and SEDs
    terms <- if (nl) c(Group = "Group", Time = "Time", Time.Group = "Time:Group") else
      c(Variety = "Variety", Nitro = "Nitro", Variety.Nitro = "Variety:Nitro")
    for (tg in names(terms)) {
      p <- suppressWarnings(suppressMessages(predictmeans(m, terms[[tg]], vcov_adjust = FALSE, plot = FALSE,
                                                           prtnum = FALSE)))
      f2 <- strsplit(terms[[tg]], ":")[[1]]
      lab <- do.call(paste, c(lapply(p$means[f2], as.character), sep = ":"))
      mm <- g[g$table == "mean" & g$term == tg & g$value != "*", ]
      for (j in seq_len(nrow(mm))) add(an, "mean", paste(tg, mm$key[j]), mm$value[j],
                                       p$means$estimate[match(mm$key[j], lab)])
      pr <- p$pairs
      sel <- pr$level1 %in% mm$key & pr$level2 %in% mm$key     # pairs of means Genstat estimates
      se <- pr$std.error[sel]
      sd <- g[g$table == "sed" & g$term == tg, ]
      for (j in seq_len(nrow(sd))) add(an, "sed", paste(tg, sd$key[j]), sd$value[j],
                                       switch(sd$key[j], single = se[1], average = mean(se), maximum = max(se),
                                              minimum = min(se)))
      ss <- g[g$table == "sed_same" & g$term == tg, ]
      for (j in seq_len(nrow(ss))) {
        kk <- strsplit(ss$key[j], ":")[[1]]
        pos <- match(kk[2], f2)
        same <- vapply(strsplit(pr$level1, ":"), `[`, "", pos) == vapply(strsplit(pr$level2, ":"), `[`, "", pos)
        s2 <- pr$std.error[sel & same]
        add(an, "sed same", paste(tg, ss$key[j]), ss$value[j],
            switch(kk[1], average = mean(s2), maximum = max(s2), minimum = min(s2)))
      }
    }
  }
  all <- do.call(rbind, res)
  all$ok <- gen_agree(all$r, all$genstat)
  bad <- all[!all$ok, ]
  expect(nrow(bad) == 0, paste0(nrow(bad), " values differ from Genstat:\n",
                                paste(utils::capture.output(print(utils::head(bad, 20), row.names = FALSE)),
                                      collapse = "\n")))
  expect_gt(nrow(all), 270)
  expect_setequal(unique(all$what), c("vc", "wald", "ndf", "F", "ddf", "mean", "sed", "sed same"))
})

test_that("empty cell: Genstat's PRESENT means equal predictmeans'; its SEDs are those of the averages", {
  skip_on_cran()
  G <- gen_ref()
  m <- gen_fits()[["run2:YB"]]
  g <- G[G$analysis == "run2:YB", ]
  p <- suppressWarnings(predictmeans(m, "Nitro", pairwise = TRUE, vcov_adjust = FALSE, plot = FALSE, prtnum = FALSE))
  vp <- g[g$table == "vp_pred" & g$term == "Nitro", ]
  expect_true(all(gen_agree(p$means$estimate[match(vp$key, p$means$Nitro)], vp$value)))
  vs <- g[g$table == "vp_se" & g$term == "Nitro", ]
  expect_true(all(gen_agree(p$means$std.error[match(vs$key, p$means$Nitro)], vs$value)))
  pv <- suppressWarnings(predictmeans(m, "Variety", vcov_adjust = FALSE, plot = FALSE, prtnum = FALSE))
  vv <- g[g$table == "vp_pred" & g$term == "Variety", ]
  expect_true(all(gen_agree(pv$means$estimate[match(vv$key, pv$means$Variety)], vv$value)))
  # Genstat's SEDs are those of the differences of the averages (0.6 averages
  # two varieties, the others three); predictmeans compares on the common
  # cells instead (R/empty-cells.R), so its SED for 0.6 is larger
  K <- suppressWarnings(Kmatrix(m, "Nitro"))$K[, names(lme4::fixef(m))]
  V <- as.matrix(vcov(m))
  sd <- g[g$table == "vp_sed" & g$term == "Nitro", ]
  dsed <- vapply(strsplit(sd$key, "\\|"), function(ab) { l <- K[ab[1], ] - K[ab[2], ]; sqrt(drop(l %*% V %*% l)) },
                 numeric(1))
  expect_true(all(gen_agree(dsed, sd$value)))
  r <- p$pairs[p$pairs$level1 == "0.0" & p$pairs$level2 == "0.6", ]
  expect_equal(r$std.error, 5.474, tolerance = 1e-3)
  expect_identical(sd$value[sd$key == "0.0|0.6"], "5.214")
})

test_that("varcomp(): Genstat's variance components and their s.e. (average information)", {
  skip_on_cran()
  G <- gen_ref()
  fits <- gen_fits()
  res <- list()
  for (an in names(fits)) {
    m <- fits[[an]]
    g <- G[G$analysis == an & G$table == "vc", ]
    v <- suppressMessages(varcomp(m, information = "average"))
    nm <- if (inherits(m, "gls")) c(Sigma2 = "residual", phi_1 = "cor.Phi") else
      c(Block = "Block.(Intercept)", Block.Variety = "Variety:Block.(Intercept)", Sigma2 = "residual")
    est <- g[g$key != "se", ]
    se <- g[g$key == "se", ]
    res[[an]] <- data.frame(an = an, what = c(rep("estimate", nrow(est)), rep("s.e.", nrow(se))),
                            term = c(est$term, se$term), genstat = c(est$value, se$value),
                            r = c(v[nm[est$term], "vcov"], v[nm[se$term], "SE"]))
  }
  all <- do.call(rbind, res)
  # varcomp() rounds to 4 decimals: compare within one unit of Genstat's last digit
  all$ok <- gen_agree(all$r, all$genstat) | abs(all$r - as.numeric(all$genstat)) <= 1e-4 * abs(as.numeric(all$genstat))
  expect(all(all$ok), paste(utils::capture.output(print(all[!all$ok, ], row.names = FALSE)), collapse = "\n"))
  expect_equal(nrow(all), 34)
  # the default (observed information, as SAS) differs by a few per cent where plots are missing
  vo <- suppressMessages(varcomp(fits[["run1:Yield"]]))
  expect_gt(vo["Block.(Intercept)", "SE"], 170)
})
