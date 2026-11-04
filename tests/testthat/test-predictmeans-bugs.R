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
})
