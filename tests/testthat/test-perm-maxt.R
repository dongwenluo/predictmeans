## Westfall-Young max-|t| adjustment of permutation p-values (adj = "tukey",
## "mvt": single step; "westfall": step down)

test_that("pm_perm_maxT() follows Westfall and Young", {
  set.seed(1)
  T <- abs(matrix(rnorm(3 * 99), 3))
  obs <- c(2.5, 0.4, 1.8)
  mx <- apply(T, 2, max)
  ss <- (vapply(obs, function(v) sum(mx >= v), 0) + 1) / 100
  expect_equal(pm_perm_maxT(obs, T), ss)
  # step down: in decreasing order of |t|, the maximum over the remaining rows
  o <- order(obs, decreasing = TRUE)
  q <- c((sum(apply(T[o, , drop = FALSE], 2, max) >= obs[o[1]]) + 1) / 100,
         (sum(apply(T[o[2:3], , drop = FALSE], 2, max) >= obs[o[2]]) + 1) / 100,
         (sum(T[o[3], ] >= obs[o[3]]) + 1) / 100)
  sd <- numeric(3); sd[o] <- cummax(q)
  expect_equal(pm_perm_maxT(obs, T, stepdown = TRUE), sd)
  expect_true(all(pm_perm_maxT(obs, T, stepdown = TRUE) <= pm_perm_maxT(obs, T) + 1e-12))
  expect_identical(pm_perm_maxT(c(NA, 1), T[1:2, ])[1], NA_real_)
})

test_that("predictmeans(): permutation Tukey and Westfall with a permlist", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fm <- lm(yield ~ Block + nitro * Variety, data = Oats)
  r <- permmodels(fm, nperm = 999, seed = 1, prt = FALSE)
  tk <- predictmeans(fm, "nitro", permlist = r, adj = "tukey", plot = FALSE, prtnum = FALSE, ncore = 1)
  wy <- predictmeans(fm, "nitro", permlist = r, adj = "westfall", plot = FALSE, prtnum = FALSE, ncore = 1)
  pt <- predictmeans(fm, "nitro", adj = "tukey", plot = FALSE, prtnum = FALSE)
  # the single-step max-|t| is the permutation version of Tukey's test
  expect_lt(max(abs(tk$pairs$p.adj - pt$pairs$p.adj)), 0.02)
  expect_true(all(wy$pairs$p.adj <= tk$pairs$p.adj + 1e-12))
  expect_true(all(wy$pairs$p.adj >= wy$pairs$p.value - 1e-12))
  # within each atvar group
  pa <- predictmeans(fm, "nitro:Variety", atvar = "Variety", permlist = r, adj = "westfall", plot = FALSE,
                     prtnum = FALSE, ncore = 1)
  w <- !is.na(pa$pairs$by)
  expect_true(all(!is.na(pa$pairs$p.adj[w])))
  expect_true(all(pa$pairs$p.adj[w] >= pa$pairs$p.value[w] - 1e-12))
  expect_error(predictmeans(fm, "nitro", adj = "westfall", plot = FALSE, prtnum = FALSE), "needs permutation")
  expect_error(predictmeans(fm, "nitro", permlist = r, adj = "scheffe", plot = FALSE, prtnum = FALSE), "scheffe")
})

test_that("contrastmeans(): max-|t| from the refits of a permlist, and other adjustments", {
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  fm <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = Oats)
  r <- permmodels(fm, nperm = 49, seed = 1, prt = FALSE, ncore = 1, method = "raw")
  cm <- rbind(c(-1, 1, 0, 0), c(-1, 0, 1, 0), c(-1, 0, 0, 1))
  out <- contrastmeans(fm, "nitro", ctrmatrix = cm, permlist = r, adj = "westfall")[[1]]
  # by hand: |t| of each refit for the three contrasts
  rK <- cm %*% Kmatrix(fm, "nitro")$K
  T <- pm_perm_tstats(r$permlist, rK)
  expect_equal(unname(out[, "Permuted Pr(>|t|)"]),
               round(pm_perm_maxT(abs(out[, "t value"]), T, stepdown = TRUE), 4), tolerance = 1e-3)
  none <- contrastmeans(fm, "nitro", ctrmatrix = cm, permlist = r)[[1]]
  holm <- contrastmeans(fm, "nitro", ctrmatrix = cm, permlist = r, adj = "holm")[[1]]
  expect_equal(unname(holm[, 4]), round(p.adjust(none[, 4], "holm"), 4), tolerance = 1e-4)
})
