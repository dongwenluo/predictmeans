## df_method = "KR-SAS": Kenward-Roger with the observed information of the
## variance parameters, as SAS PROC MIXED DDFM=KR. Reference values from the
## SAS run of dev/sas/pm-reference.sas (dev/sas/sas-results.csv, run 15:
## the split plot without Block II x Victory x 0.6).

oats71 <- function() {
  d <- as.data.frame(nlme::Oats)
  d$nitro <- factor(d$nitro)
  droplevels(subset(d, !(Block == "II" & Variety == "Victory" & nitro == "0.6")))
}

test_that("the observed information equals the Hessian of the REML log-likelihood", {
  d <- oats71()
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block) + (1 | Block:Variety), data = d)
  X <- as.matrix(lme4::getME(m, "X")); y <- lme4::getME(m, "y"); Z <- lme4::getME(m, "Z"); Gp <- lme4::getME(m, "Gp")
  Zk <- lapply(1:2, function(k) as.matrix(Z[, (Gp[k] + 1):Gp[k + 1]]))
  reml2 <- function(s) {                 # -2 REML log-likelihood (up to a constant)
    V <- diag(s[3], length(y)) + s[1] * tcrossprod(Zk[[1]]) + s[2] * tcrossprod(Zk[[2]])
    Vi <- solve(V)
    XVX <- crossprod(X, Vi %*% X)
    r <- y - X %*% solve(XVX, crossprod(X, Vi %*% y))
    as.numeric(determinant(V)$modulus + determinant(XVX)$modulus + crossprod(r, Vi %*% r))
  }
  vc <- as.data.frame(lme4::VarCorr(m))
  s <- c(vc$vcov[vc$grp == "Block:Variety"], vc$vcov[vc$grp == "Block"], vc$vcov[vc$grp == "Residual"])
  H <- numDeriv::hessian(reml2, s)
  expect_equal(attr(pm_vcov_kr(m, "observed"), "W"), 2 * solve(H), tolerance = 1e-5)
  expect_identical(attr(pm_vcov_kr(m, "observed"), "information"), "observed")
})

test_that("predictmeans(), contrastmeans() and pm_joint() with df_method = \"KR-SAS\" equal SAS DDFM=KR", {
  d <- oats71()
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d)
  pv <- predictmeans(m, "Variety", df_method = "KR-SAS", plot = FALSE, prtnum = FALSE)
  expect_equal(pv$means$df, c(8.58, 8.58, 8.7), tolerance = 2e-3)
  expect_equal(pv$pairs$df, c(9.79, 9.97, 9.97), tolerance = 2e-3)
  expect_equal(pv$pairs$std.error, c(6.8298, 6.8647, 6.8647), tolerance = 2e-4)
  expect_identical(pv$settings$df_method, "KR-SAS")
  expect_true(pv$settings$vcov_adjust)
  # pbkrtest's version differs here by about 1%
  pk <- predictmeans(m, "Variety", plot = FALSE, prtnum = FALSE)
  expect_true(all(abs(pk$pairs$df / pv$pairs$df - 1) > 0.005))
  cn <- rbind(c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
  ct <- contrastmeans(m, "Variety:nitro", ctrmatrix = cn, df_method = "KR-SAS")[[1]]
  expect_equal(unname(ct[, "df"]), c(31.3, 34.1), tolerance = 2e-3)
  expect_equal(unname(abs(ct[, "Std. Error"])), c(9.5577, 9.9495), tolerance = 2e-4)
  jt <- pm_joint(m, "Variety", df_method = "KR-SAS")
  expect_equal(c(jt$df2, jt$F.ratio), c(9.91, 1.38), tolerance = 5e-3)
  expect_identical(attr(jt, "df_method"), "KR-SAS")
  # balanced data: the same as pbkrtest's
  dd <- as.data.frame(nlme::Oats); dd$nitro <- factor(dd$nitro)
  mb <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = dd)
  expect_equal(predictmeans(mb, "nitro", df_method = "KR-SAS", plot = FALSE, prtnum = FALSE)$pairs[, c("std.error", "df")],
               predictmeans(mb, "nitro", plot = FALSE, prtnum = FALSE)$pairs[, c("std.error", "df")], tolerance = 1e-5)
  L <- rbind(c(0, -1, 0, rep(0, 9)))
  expect_equal(df_term(m, ctrmatrix = L, type = "KR-SAS"),
               df_term(m, ctrmatrix = L, vcov_adj = pm_vcov_kr(m, "observed")))
})
