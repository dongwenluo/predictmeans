# Regression tests for bugs found while splitting predictmeans() into stages.

oats_lmer <- function() {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  suppressMessages(lmerTest::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats))
}

test_that("permlist works with atvar", {
  m <- oats_lmer()
  perm <- suppressWarnings(permmodels(m, nperm = 9, prt = FALSE, ncore = 1, seed = 1))
  r <- predictmeans(m, "nitro:Variety", permlist = perm, atvar = "Variety", ncore = 1,
                    plot = FALSE, prtnum = FALSE)
  expect_equal(nrow(r$means), 12)
  expect_true(all(is.na(r$means$df)))
  expect_true(all(r$pairs$p.method == "permutation"))
  # p-values within each variety are adjusted and present
  within <- r$pairs[!is.na(r$pairs$by), ]
  expect_equal(nrow(within), 3 * 6)
  expect_false(anyNA(within$p.adj))
})

test_that("a user-supplied Df is used for every pair of an lmer model", {
  m <- oats_lmer()
  r <- predictmeans(m, "nitro", Df = 20, adj = "BH", plot = FALSE, prtnum = FALSE)
  pv <- r$`Pairwise p-value`
  tv <- pv[upper.tri(pv)]
  expect_equal(t(pv)[upper.tri(pv)], round(p.adjust(2 * pt(-abs(tv), 20), "BH"), 4),
               tolerance = 1e-4)
  expect_true(all(r$means$df == 20))
  r_t <- predictmeans(m, "nitro", Df = 20, adj = "tukey", plot = FALSE, prtnum = FALSE)
  pt_ <- r_t$`Pairwise p-value`
  expect_false(anyNA(t(pt_)[upper.tri(pt_)]))
  expect_equal(t(pt_)[upper.tri(pt_)],
               round(ptukey(sqrt(2) * abs(tv), 4, 20, lower.tail = FALSE), 4), tolerance = 1e-4)
})

test_that("factor levels containing ':' give the same results (B4)", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  Oats <- Oats[-c(1, 2, 30), ]                 # unbalanced: SEDs differ
  O2 <- Oats
  levels(O2$Variety) <- c("Golden:Rain", "Marv:1", "Vic")
  fit <- function(d) suppressMessages(lme4::lmer(yield ~ nitro*Variety + (1|Block/Variety), data = d))
  m1 <- fit(Oats); m2 <- fit(O2)
  strip <- function(x) {
    if (is.list(x) && !is.data.frame(x)) return(lapply(unname(x), strip))
    if (is.data.frame(x)) { x[] <- lapply(x, function(v) if (is.factor(v) || is.character(v)) NULL else v); return(unname(as.list(x))) }
    if (is.numeric(x)) { a <- attributes(x); x <- as.vector(x); if (!is.null(a[["For the Same Level of Factor"]])) attr(x, "s") <- unname(a[["For the Same Level of Factor"]]); return(x) }
    NULL
  }
  for (args in list(list(modelterm = "nitro:Variety", pairwise = TRUE),
                    list(modelterm = "nitro:Variety", atvar = "Variety", adj = "BH"),
                    list(modelterm = "nitro:Variety", atvar = "nitro", adj = "tukey"))) {
    r1 <- expect_silent(do.call(predictmeans, c(list(m1), args, plot = FALSE, prtnum = FALSE)))
    r2 <- expect_silent(do.call(predictmeans, c(list(m2), args, plot = FALSE, prtnum = FALSE)))
    keep <- !grepl("Plot|plot|settings|^Golden|^Marv|^Vic", names(r1))
    expect_equal(strip(r2[!grepl("Plot|plot|settings|^Golden|^Marv|^Vic", names(r2))]), strip(r1[keep]))
    expect_equal(r2$pairs[c("estimate", "std.error", "df", "p.value", "p.adj")],
                 r1$pairs[c("estimate", "std.error", "df", "p.value", "p.adj")])
    expect_equal(r2$means$group, r1$means$group)
  }
  pdf(NULL); on.exit(dev.off())
  expect_no_error(predictmeans(m2, "nitro:Variety", atvar = "Variety", adj = "BH", prtnum = FALSE, prtplt = FALSE))
})

test_that("atvar groups with one or no observed level get letters without error", {
  set.seed(1)
  d <- expand.grid(repl = factor(1:3), seedlot = factor(1:4),
                   treat = factor(c("control", "nick", "bwss"), levels = c("control", "nick", "bwss")))
  d$contcomp <- factor(ifelse(d$treat == "control", "control", "Treated"))
  d$y <- as.numeric(d$treat) * 10 + as.numeric(d$seedlot) + rnorm(nrow(d))
  fm <- aov(y ~ repl + contcomp + seedlot + treat/contcomp + contcomp/seedlot + treat/contcomp/seedlot, data = d)
  pm <- suppressWarnings(predictmeans(fm, "contcomp:seedlot:treat", atvar = c("treat", "contcomp"),
                                      adj = "BH", plot = FALSE, prtnum = FALSE))
  expect_equal(nrow(pm$mean_table$Table), 12)
  expect_false(anyNA(pm$mean_table$Table$LetterGrp))
  # p-values go to the groups with data: none for the empty control : Treated
  # group, a full set for nick : Treated (group order differs from row order)
  num <- function(tab) suppressWarnings(as.numeric(tab[lower.tri(tab)]))
  expect_true(all(is.na(num(pm[["control : Treated"]]))))
  expect_false(anyNA(num(pm[["nick : Treated"]])))
})

test_that("atvar groups with a single observed level get a letter", {
  set.seed(2)
  d <- expand.grid(date = factor(c("d1", "d2", "d3")), trt = factor(c("a", "b")), rep = 1:4)
  d <- d[!(d$date == "d2" & d$trt == "a"), ]          # d2 has only trt b
  d$y <- as.numeric(d$trt) + rnorm(nrow(d))
  fm <- lm(y ~ date * trt, data = d)
  pm <- suppressWarnings(predictmeans(fm, "date:trt", atvar = "date", adj = "BH",
                                      plot = FALSE, prtnum = FALSE))
  expect_equal(nrow(pm$mean_table$Table), 5)
  expect_false(anyNA(pm$mean_table$Table$LetterGrp))
})

test_that("deprecated bar = TRUE keeps the other plots (mean plot not NULL)", {
  Oats <- as.data.frame(nlme::Oats)
  Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro * Variety, data = Oats)
  pdf(NULL); on.exit(dev.off())
  expect_warning(pm <- predictmeans(fm, "nitro", adj = "BH", bar = TRUE, prtplt = FALSE, prtnum = FALSE),
                 'which.plots = c\\("mean", "bar", "pval", "back"\\)')
  expect_s3_class(pm$predictmeansPlot$meanPlot, "ggplot")
  expect_s3_class(pm$predictmeansBarPlot, "ggplot")
  pm2 <- suppressWarnings(predictmeans(fm, "nitro", adj = "BH", bar = FALSE, prtplt = FALSE, prtnum = FALSE))
  expect_s3_class(pm2$predictmeansPlot$meanPlot, "ggplot")
  expect_null(pm2$predictmeansBarPlot)
})

test_that("confounded factors: non-estimable means and differences are NA, no error", {
  # vol 400 occurs only with droplet fine (as in a spray trial)
  set.seed(3)
  d <- data.frame(droplet = factor(rep(c("medium", "fine"), c(40, 10)), levels = c("medium", "fine")),
                  vol = factor(c(rep(c("130", "200"), 20), rep("400", 10))),
                  conc = factor(rep(c("none", "half"), 25)))
  d$y <- 5 + 0.3 * (d$conc == "half") - 0.4 * (d$droplet == "fine") + rnorm(50, sd = 0.2)
  fm <- lm(y ~ droplet + conc + vol, data = d)
  pv <- suppressMessages(predictmeans(fm, "vol", adj = "BH", plot = FALSE, prtnum = FALSE))
  expect_true(all(is.na(pv$means$estimate)))
  pr <- pv$pairs
  expect_false(is.na(pr$estimate[pr$level1 == "130" & pr$level2 == "200"]))
  expect_true(all(is.na(pr$estimate[pr$level2 == "400" | pr$level1 == "400"])))
  expect_equal(sum(!is.na(pv$mean_table$Table$LetterGrp)), 2)
  pd <- suppressMessages(predictmeans(fm, "droplet", adj = "BH", plot = FALSE, prtnum = FALSE))
  expect_true(is.na(pd$pairs$estimate))
  pc <- suppressMessages(predictmeans(fm, "conc", adj = "BH", plot = FALSE, prtnum = FALSE))
  expect_false(anyNA(pc$pairs$estimate))
  # weights = "cells" averages over the observed cells only
  pw <- suppressMessages(predictmeans(fm, "vol", weights = "cells", plot = FALSE, prtnum = FALSE))
  expect_false(anyNA(pw$means$estimate))
})

test_that("letters: one helper handles NA p-values overall and within atvar groups", {
  P <- matrix(c(1, .01, NA, .01, 1, NA, NA, NA, 1), 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  L <- predictmeans:::pm_letter_groups(P, c(TRUE, TRUE, TRUE), 0.05)
  expect_true(is.na(L[["c"]]) && L[["a"]] != L[["b"]])
  # a level alone because it is the only one with data gets "A"
  expect_equal(unname(predictmeans:::pm_letter_groups(P, c(FALSE, TRUE, FALSE), 0.05)), c(NA, "A", NA))
  # a level left alone by non-estimable comparisons gets NA
  P2 <- matrix(c(1, NA, NA, 1), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_true(all(is.na(predictmeans:::pm_letter_groups(P2, c(TRUE, TRUE), 0.05))))

  # within atvar groups: vol 400 is confounded with droplet fine, so within
  # each conc group the 400 comparisons are not estimable but 130 vs 200 is
  set.seed(3)
  d <- data.frame(droplet = factor(rep(c("medium", "fine"), c(40, 10)), levels = c("medium", "fine")),
                  vol = factor(c(rep(c("130", "200"), 20), rep("400", 10))),
                  conc = factor(c(rep(c("none", "none", "half", "half"), 10), rep(c("none", "half"), 5))))
  d$y <- 5 + 0.3 * (d$conc == "half") - 0.4 * (d$droplet == "fine") + rnorm(50, sd = 0.2)
  fm <- lm(y ~ droplet + conc * vol, data = d)
  pm <- suppressMessages(predictmeans(fm, "conc:vol", atvar = "conc", adj = "BH", plot = FALSE, prtnum = FALSE))
  tab <- pm$mean_table$Table
  expect_equal(nrow(tab), 6)
  expect_true(all(is.na(tab$LetterGrp[tab$vol == "400"])))
  expect_false(anyNA(tab$LetterGrp[tab$vol != "400"]))
})
