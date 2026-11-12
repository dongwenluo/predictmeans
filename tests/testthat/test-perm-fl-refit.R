## Freedman-Lane with refits for lmer/lme/gls (permmodels(method = "freedman-lane"), refit = TRUE)

test_that("Freedman-Lane with refits equals a direct computation (lmer, type III)", {
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  m <- lme4::lmer(yield ~ nitro + Variety + (1 | Block), data = Oats)
  r <- permmodels(m, nperm = 19, type = "III", seed = 4, prt = FALSE, ncore = 1)
  # by hand: the model without Variety, its residuals decorrelated with its V
  m0 <- lme4::lmer(yield ~ nitro + (1 | Block), data = Oats)
  W0 <- pm_whiten(lme4::getME(m0, "Z") %*% lme4::getME(m0, "Lambda"))
  f0 <- fitted(m0) - as.vector(lme4::getME(m0, "Z") %*% lme4::getME(m0, "b"))
  e0 <- W0$solve(Oats$yield - f0)
  idx <- pm_perm_stream(nrow(Oats), 4)(19)
  Fv <- function(fit) { a <- pm_anova(fit, "III", "Chisq"); a["Variety", "Chisq"] / a["Variety", "Df"] }
  Fs <- apply(idx, 2, function(i) Fv(lme4::refit(m, f0 + W0$mult(e0[i]))))
  expect_equal(r$ANOVA["Variety", "Perm_p_value"],
               round((sum(round(Fs, 6) >= round(Fv(m), 6)) + 1) / 20, 4))
  # pairs of Variety: the family null is the term's null, so pair 1 v 2 is the coefficient's test
  pm <- predictmeans(m, "Variety", permlist = r, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_equal(pm$pairs$p.value[1], unname(r$COEFFICENTS["VarietyMarvellous", "Perm_p_value"]), tolerance = 1e-4)
  expect_match(attr(r$ANOVA, "heading"), "refitted")
})

test_that("Freedman-Lane with refits for lme and gls; refit = FALSE holds V", {
  skip_on_cran()
  Oats <- as.data.frame(nlme::Oats); Oats$nitro <- factor(Oats$nitro)
  m <- nlme::lme(yield ~ nitro + Variety, random = ~ 1 | Block, data = Oats)
  r <- permmodels(m, nperm = 19, type = "II", seed = 1, prt = FALSE, ncore = 1)
  expect_equal(r$perm$scheme, "freedman-lane-refit")
  expect_true(all(r$ANOVA$Perm_p_value > 0))
  r0 <- permmodels(m, nperm = 19, type = "II", seed = 1, prt = FALSE, refit = FALSE)
  expect_equal(r0$perm$scheme, "freedman-lane")
  expect_equal(r$ANOVA[, "F value"], r0$ANOVA[, "F value"])
  g <- nlme::gls(yield ~ nitro + Variety, correlation = nlme::corCompSymm(form = ~ 1 | Block), data = Oats)
  rg <- permmodels(g, nperm = 19, type = "III", seed = 1, prt = FALSE, ncore = 1)
  expect_identical(rownames(rg$ANOVA), c("nitro", "Variety"))
  pmg <- predictmeans(g, "nitro", permlist = rg, plot = FALSE, prtnum = FALSE, ncore = 1)
  expect_true(all(pmg$pairs$p.value > 0 & pmg$pairs$p.value <= 1))
})
