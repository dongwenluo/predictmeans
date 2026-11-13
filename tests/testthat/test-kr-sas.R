## Kenward-Roger (df_method = "Kenward-Roger", the default for REML lmer fits)
## uses the observed information of the variance parameters, as SAS PROC MIXED
## DDFM=KR. Reference values from the
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

test_that("predictmeans(), contrastmeans() and pm_joint() with df_method = \"Kenward-Roger\" equal SAS DDFM=KR", {
  d <- oats71()
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = d)
  pv <- predictmeans(m, "Variety", df_method = "Kenward-Roger", plot = FALSE, prtnum = FALSE)
  expect_equal(pv$means$df, c(8.58, 8.58, 8.7), tolerance = 2e-3)
  expect_equal(pv$pairs$df, c(9.79, 9.97, 9.97), tolerance = 2e-3)
  expect_equal(pv$pairs$std.error, c(6.8298, 6.8647, 6.8647), tolerance = 2e-4)
  expect_identical(pv$settings$df_method, "Kenward-Roger")
  expect_true(pv$settings$vcov_adjust)
  # it is the default; pbkrtest's version ("KR-expected") differs here by about 1%
  expect_identical(predictmeans(m, "Variety", plot = FALSE, prtnum = FALSE)$pairs, pv$pairs)
  pk <- predictmeans(m, "Variety", df_method = "KR-expected", plot = FALSE, prtnum = FALSE)
  expect_true(all(abs(pk$pairs$df / pv$pairs$df - 1) > 0.005))
  cn <- rbind(c(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0), c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1))
  ct <- contrastmeans(m, "Variety:nitro", ctrmatrix = cn, df_method = "Kenward-Roger")[[1]]
  expect_equal(unname(ct[, "df"]), c(31.3, 34.1), tolerance = 2e-3)
  expect_equal(unname(abs(ct[, "Std. Error"])), c(9.5577, 9.9495), tolerance = 2e-4)
  jt <- pm_joint(m, "Variety", df_method = "Kenward-Roger")
  expect_equal(c(jt$df2, jt$F.ratio), c(9.91, 1.38), tolerance = 5e-3)
  expect_identical(attr(jt, "df_method"), "Kenward-Roger")
  # balanced data: the same as pbkrtest's
  dd <- as.data.frame(nlme::Oats); dd$nitro <- factor(dd$nitro)
  mb <- lme4::lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = dd)
  expect_equal(predictmeans(mb, "nitro", df_method = "Kenward-Roger", plot = FALSE, prtnum = FALSE)$pairs[, c("std.error", "df")],
               predictmeans(mb, "nitro", df_method = "KR-expected", plot = FALSE, prtnum = FALSE)$pairs[, c("std.error", "df")],
               tolerance = 1e-5)
  L <- rbind(c(0, -1, 0, rep(0, 9)))
  expect_equal(df_term(m, ctrmatrix = L, type = "Kenward-Roger"),
               df_term(m, ctrmatrix = L, vcov_adj = pm_vcov_kr(m, "observed")))
})

## "KR-average": Kenward-Roger with the average information of the REML
## iterations of Genstat and ASReml, A_ij = y'P G_i P G_j P y / 2 (Gilmour,
## Thompson and Cullis 1995). The Genstat runs of dev/genstat/README.md gave
## 8.7 df for Variety (oats, a whole plot missing) and 25.8 for Group x Time
## (AR(1), dev/sas/ar1.txt), where this gives 8.75 and 25.76 (notes only:
## that output was pasted, not uploaded).

test_that("KR-average: the average information is the mean of the observed and expected ones for variance components", {
  d <- oats71()
  m <- lme4::lmer(yield ~ Variety * nitro + (1 | Block) + (1 | Block:Variety), data = d)
  Wo <- attr(pm_vcov_kr(m, "observed"), "W")
  We <- attr(pm_vcov_kr(m, "expected"), "W")
  Wa <- pm_vcov_kr(m, "average")
  expect_identical(attr(Wa, "information"), "average")
  expect_equal(solve(attr(Wa, "W")), (solve(Wo) + solve(We)) / 2, tolerance = 1e-8)
  e <- nlme::lme(yield ~ Variety * nitro, random = ~ 1 | Block/Variety, data = d)
  expect_equal(pm_joint(e, "Variety", df_method = "KR-average")$df2,
               pm_joint(m, "Variety", df_method = "KR-average")$df2, tolerance = 1e-5)
})

test_that("KR-average for gls with AR(1): the average information y'P G_i P G_j P y / 2, computed densely", {
  d <- read.table(test_path("_data", "sas", "ar1.txt"), header = TRUE)
  d$Subj <- factor(d$Subj); d$Group <- factor(d$Group); d$t <- d$Time; d$Time <- factor(d$Time)
  g <- nlme::gls(Y ~ Group * Time, data = d, correlation = nlme::corAR1(form = ~ t | Subj))
  R <- pm_nlme_reml(g)
  S0 <- as.matrix(R$Sig(R$th0))
  G <- lapply(seq_len(R$k), function(i) {
    if (i == R$k) return(2 * S0)
    h <- 1e-4 * max(1, abs(R$th0[i]))
    e <- replace(numeric(R$k), i, h)
    as.matrix(R$Sig(R$th0 + e) - R$Sig(R$th0 - e)) / (2 * h)
  })
  Si <- solve(S0)
  P <- Si - Si %*% R$X %*% solve(crossprod(R$X, Si %*% R$X), crossprod(R$X, Si))
  Py <- P %*% R$y
  AI <- outer(seq_len(R$k), seq_len(R$k), Vectorize(function(i, j) 0.5 * drop(crossprod(Py, G[[i]] %*% P %*% G[[j]] %*% Py))))
  expect_equal(attr(pm_vcov_kr_nlme(g, "average"), "W"), solve(AI), tolerance = 1e-6)
  last <- rbind(c(1, 0, 0, -1), c(0, 1, 0, -1), c(0, 0, 1, -1))
  C <- kronecker(rbind(c(1, -1)), last)
  ddf <- vapply(c("KR-expected", "KR-average", "Kenward-Roger"), function(dm)
    pm_joint(g, "Group:Time", df_method = dm, ctrmatrix = C)$df2, numeric(1))
  expect_true(ddf[1] < ddf[2] && ddf[2] < ddf[3])
  p <- predictmeans(g, "Group", df_method = "KR-average", vcov_adjust = FALSE, plot = FALSE, prtnum = FALSE)
  expect_identical(p$settings$df_method, "KR-average")
  expect_false(p$settings$vcov_adjust)
})
