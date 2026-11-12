## contrastmeans(joint = TRUE): all the contrasts tested together

test_that("the joint F test of contrasts", {
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  cn <- rbind(c(-1, 1, 0, 0), c(-1, 0, 1, 0), c(-1, 0, 0, 1))
  # lm, balanced: the contrasts span the nitrogen differences, so the joint F
  # is the nitro F of the ANOVA table
  fm <- lm(yield ~ Block + nitro * Variety, data = Oats)
  jt <- contrastmeans(fm, "nitro", ctrmatrix = cn, joint = TRUE)[["Joint F test of the contrasts"]]
  a <- anova(fm)
  expect_equal(jt$F.ratio, round(a["nitro", "F value"], 4))
  expect_equal(c(jt$df1, jt$df2), c(3, a["Residuals", "Df"]))
  # lmer: the same as pm_joint() (Kenward-Roger); "join" works by partial matching
  fl <- lme4::lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
  j2 <- contrastmeans(fl, "nitro", ctrmatrix = cn, join = TRUE)[[3]]
  pj <- pm_joint(fl, "nitro", ctrmatrix = cn)
  expect_equal(j2$F.ratio, round(pj$F.ratio, 4))
  expect_equal(j2$df2, round(pj$df2, 4))
  # one contrast: F = t^2; dependent rows are dropped (df1 = rank)
  one <- contrastmeans(fm, "nitro", ctrmatrix = cn[1, , drop = FALSE], joint = TRUE)
  expect_equal(one[[3]]$F.ratio, unname(round(one[[1]][1, "t value"]^2, 4)), tolerance = 1e-3)
  dep <- contrastmeans(fm, "nitro", ctrmatrix = rbind(cn, cn[1, ] + cn[2, ]), joint = TRUE)[[3]]
  expect_identical(dep$df1, 3)
  # without joint the result is as before (two components)
  expect_length(contrastmeans(fm, "nitro", ctrmatrix = cn), 2)
})

test_that("the joint F test by permutation", {
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  cn <- rbind(c(-1, 1, 0, 0), c(-1, 0, 1, 0), c(-1, 0, 0, 1))
  fl <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = Oats)
  r <- permmodels(fl, nperm = 29, seed = 1, prt = FALSE, ncore = 1, method = "raw")
  jt <- contrastmeans(fl, "nitro", ctrmatrix = cn, permlist = r, joint = TRUE)[[3]]
  # by hand: the Wald F of the refits against the model's
  L <- cn %*% Kmatrix(fl, "nitro")$K
  wald <- function(b, V) { d <- L %*% b; as.numeric(t(d) %*% solve(L %*% V %*% t(L), d)) / 3 }
  Fs <- vapply(r$permlist, function(x) wald(x$coef, as.matrix(x$vcov)), 0)
  Fo <- wald(lme4::fixef(fl), as.matrix(vcov(fl)))
  expect_equal(jt$Perm_p_value, round((sum(round(Fs, 6) >= round(Fo, 6)) + 1) / 30, 4))
  # Freedman-Lane (lm, V fixed) and with refits (lmer)
  fm <- lm(yield ~ Block + nitro + Variety, data = Oats)
  rf <- permmodels(fm, nperm = 99, seed = 2, prt = FALSE)
  expect_true(contrastmeans(fm, "nitro", ctrmatrix = cn, permlist = rf, joint = TRUE)[[3]]$Perm_p_value <= 0.02)
  rr <- permmodels(fl, nperm = 19, seed = 3, prt = FALSE, ncore = 1)
  expect_true(contrastmeans(fl, "nitro", ctrmatrix = cn, permlist = rr, joint = TRUE)[[3]]$Perm_p_value <= 0.1)
})
