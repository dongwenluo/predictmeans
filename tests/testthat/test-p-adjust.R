# "sidak", "scheffe" and "mvt" adjustments, checked against emmeans.

oats_unbal <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro, labels = c("n0", "n2", "n4", "n6"))
  set.seed(7)
  d[-sample(nrow(d), 9), ]
}
pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b))

test_that("sidak, scheffe and mvt pairs match emmeans (overall and within atvar)", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("mvtnorm")
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats_unbal())
  for (a in c("sidak", "scheffe", "mvt")) {
    tol <- if (a == "mvt") 1e-3 else 1e-8
    pm <- predictmeans(fm, "nitro", adj = a, df_method = "KR-expected", plot = FALSE)
    em <- summary(pairs(emmeans::emmeans(fm, ~nitro, lmer.df = "kenward-roger"), adjust = a))
    lv <- do.call(rbind, strsplit(as.character(em$contrast), " - "))
    i <- match(pair_key(pm$pairs$level1, pm$pairs$level2), pair_key(lv[, 1], lv[, 2]))
    expect_equal(pm$pairs$p.adj, em$p.value[i], tolerance = tol)
    # the legacy p-value matrix holds the same adjusted p-values
    P <- pm$`Pairwise p-value`
    expect_equal(sort(P[lower.tri(P)]), sort(round(pm$pairs$p.adj, 4)))

    pa <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = a, df_method = "KR-expected", plot = FALSE)
    ea <- summary(pairs(emmeans::emmeans(fm, ~nitro | Variety, lmer.df = "kenward-roger"), adjust = a))
    lv <- do.call(rbind, strsplit(as.character(ea$contrast), " - "))
    s <- pa$pairs[!is.na(pa$pairs$by), ]
    j <- match(paste(s$by, pair_key(sub(":.*", "", s$level1), sub(":.*", "", s$level2))),
               paste(ea$Variety, pair_key(lv[, 1], lv[, 2])))
    expect_equal(s$p.adj, ea$p.value[j], tolerance = tol)
  }
})

test_that("contrastmeans(): sidak and scheffe (rank of ctrmatrix) match emmeans", {
  skip_if_not_installed("emmeans")
  fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = oats_unbal())
  cm <- rbind(c(-1, 1, 0, 0), c(0, -1, 1, 0), c(0, 0, -1, 1))
  for (a in c("sidak", "scheffe")) {
    r <- contrastmeans(fm, "nitro", ctrmatrix = cm, adj = a, df_method = "KR-expected")[[1]]
    e <- summary(emmeans::contrast(emmeans::emmeans(fm, ~nitro, lmer.df = "kenward-roger"),
                                   method = list(a = cm[1, ], b = cm[2, ], c = cm[3, ]), adjust = a))
    expect_equal(unname(r[, "Pr(>|t|)"]), round(e$p.value, 4))
  }
})

test_that("mvt repeats exactly and leaves the user's random numbers alone", {
  skip_if_not_installed("mvtnorm")
  d <- oats_unbal()   # (not in the call: predictmeans() re-evaluates the data argument)
  fm <- lm(yield ~ Block + nitro*Variety, data = d)
  set.seed(5); u1 <- runif(1)
  set.seed(5)
  p1 <- predictmeans(fm, "nitro", adj = "mvt", plot = FALSE)$pairs$p.adj
  u2 <- runif(1)
  p2 <- predictmeans(fm, "nitro", adj = "mvt", plot = FALSE)$pairs$p.adj
  expect_identical(u1, u2)
  expect_identical(p1, p2)
})

test_that("with permutation p-values sidak works; scheffe and mvt stop", {
  d <- oats_unbal()
  fm <- nlme::lme(yield ~ nitro*Variety, random = ~1|Block/Variety, data = d)
  pl <- permmodels(fm, method = "raw", nperm = 30, ncore = 1, seed = 1, prt = FALSE)
  r <- predictmeans(fm, "nitro", permlist = pl, adj = "sidak", plot = FALSE, ncore = 1)
  m <- sum(!is.na(r$pairs$p.value))
  expect_equal(r$pairs$p.adj, 1 - (1 - r$pairs$p.value)^m)
  expect_error(predictmeans(fm, "nitro", permlist = pl, adj = "scheffe", plot = FALSE, ncore = 1), "not available")
})
