## Bands of semipred() for generalized fits, against independent calculations.
##   2. wageTUnion (binary): M. P. Wand's script (original/wageTUnionLMM.R)
##      fits by PQL (glmmPQL) and draws mu +- 2 sd, sd = sqrt(diag(Cg (C'WC +
##      lambda D)^-1 Cg')) mu'(eta). semireg() fits by Laplace (glmer), so
##      (a) the two fits are compared, (b) Wand's formula is applied to the
##      semireg() fit and compared with semipred()'s link-scale SE and its
##      delta-method SE of the probability.
##   3. A difference of two curves of a GLMM (Poisson counts, split plot with
##      curves over time): semipred(contr = ) against the same quantities
##      computed from lme4's own pieces of the fit (X, Z, Lambda, working
##      weights): V = (C'WC + blockdiag(0, (Lambda Lambda')^-1))^-1, the
##      contrast rows of the prediction design; on the log scale, as a ratio,
##      and as a difference of means (delta method); and roughly against
##      mgcv::gam() with the same model (its own smoothing parameters).
## Result (2026-11-17; semipred(unconditional = FALSE), as the references are
## given the variance parameters; the default for GLMMs adds the
## smoothing-parameter correction):
##   2. fits: spline SD PQL 0.0788, Laplace 0.0799; probability curves within
##      0.003 (probabilities up to 0.32). Wand's formula on the semireg() fit
##      equals semipred(): probability, link-scale SE and SE of the
##      probability to 1e-14. Limits plogis(eta -+ 1.96 SE) (inside 0-1),
##      where the script draws mu -+ 2 SE_mu clipped at 0.001.
##   3. semipred(contr = ) equals the calculation from lme4's pieces to 1e-15:
##      log-scale difference and SE, ratio and limits = exp(), difference of
##      means and its delta-method SE; limits t on the 10 whole-plot df.
##      mgcv::gam (own basis and smoothing): difference within 0.011, SE 1-4%.
##   Rscript dev/semireg-golden/check-bands-glmm.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages({library(nlme); library(MASS)})
q <- function(e) suppressWarnings(suppressMessages(e))
orig <- new.env(); sys.source("dev/semireg-golden/original/ZOSull.r", orig); ZOSull <- orig$ZOSull
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

## ---- 2. wageTUnion ---------------------------------------------------------
wu <- read.table("inst/extdata/semireg/wageTUnion.txt", header = TRUE)
x <- wu[, 1]; y <- wu[, 2]
X <- cbind(1, x)
kn <- quantile(unique(x), seq(0, 1, length = 17))[-c(1, 17)]
range.x <- range(x)
Z <- ZOSull(x, range.x, kn)
group <- rep(1, length(x))
data.fr <- groupedData(y ~ x | group, data = data.frame(x, y))
fitP <- q(glmmPQL(y ~ -1 + X, random = list(group = pdIdent(~ -1 + Z)), data = data.fr, family = binomial,
                  verbose = FALSE))
xg <- seq(min(x), max(x), length = 101)
## Wand's band for a fit given by (beta, u, sd_u) on its own basis (Zd at the
## data, Zg at the grid)
wand <- function(beta, u, sd_u, Zd, Zg) {
  eta <- as.vector(X %*% beta + Zd %*% u)
  w <- exp(eta) / (1 + exp(eta))^2
  C <- cbind(X, Zd)
  M <- crossprod(C * w, C) + diag(c(0, 0, rep(1 / sd_u^2, ncol(Zd))))
  Cg <- cbind(1, xg, Zg)
  etag <- as.vector(Cg %*% c(beta, u))
  se_eta <- sqrt(rowSums((Cg %*% solve(M)) * Cg))
  list(mu = plogis(etag), se_eta = se_eta, se_mu = se_eta * exp(etag) / (1 + exp(etag))^2)
}
sdP <- as.numeric(fitP$sigma * exp(unlist(fitP$modelStruct)))
wP <- wand(fitP$coef$fixed, unlist(fitP$coef$random), sdP, Z, ZOSull(xg, range.x, kn))

d <- data.frame(y = y, x = x)
fm <- q(semireg(y ~ s(x, bs = "os", k = 15, xt = list(range = range.x)), family = binomial, data = d, prt = FALSE))
blk <- names(fm$u_lst)[1]
Zs <- function(v) do.call(cbind, pm_sm_predict(fm$smooths, data.frame(x = v), blk))
sdS <- as.numeric(attr(lme4::VarCorr(fm$semer)[[blk]], "stddev"))
sc <- unname(coef(stats::lm(Zs(x)[, 1] ~ Z[, 1] - 1)))          # semireg's basis = Wand's x sc
wS <- wand(lme4::fixef(fm$semer), fm$u_lst[[1]], sdS, Zs(x), Zs(xg))
sp <- q(semipred(fm, "x", "x", covariateV = xg, unconditional = FALSE, prt = FALSE))$pred_df
se_link <- sp$SE                                                # the link scale
stopifnot(isTRUE(all.equal(sp$LL, plogis(qlogis(sp$y) - qnorm(0.975) * sp$SE), tolerance = 1e-10)))
cat("2. wageTUnion (binary, no random term besides the smooth: normal limits)\n")
cat(sprintf("   (a) fits: spline SD PQL %.4f, Laplace %.4f (on Wand's basis; semireg's basis is his x %.4f);\n",
            sdP, sdS * sc, sc),
    sprintf("      probability curve max |Laplace - PQL| = %.4f (probabilities up to %.2f)\n",
            max(abs(wS$mu - wP$mu)), max(wP$mu)))
cat(sprintf("   (b) Wand's formula on the semireg() fit vs semipred(): probability %.1e, link-scale SE %.1e, SE of the probability %.1e\n",
            rel(sp$y, wS$mu), rel(se_link, wS$se_eta), rel(sp$SE_mu, wS$se_mu)))
cat(sprintf("       band SE of the probability, Laplace vs PQL fit: ratio %.2f to %.2f\n",
            min(wS$se_mu / wP$se_mu), max(wS$se_mu / wP$se_mu)))
cat("       limits: semipred plogis(eta -+ 1.960 SE) (inside 0-1, asymmetric); Wand mu -+ 2 SE_mu (clipped at 0.001)\n")

## ---- 3. a difference of two curves of a Poisson GLMM ---------------------------
set.seed(31)
dd <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
eta <- 1 + sin(dd$t / 2) / 2 + 0.04 * (dd$A == "a2") * dd$t + rnorm(6, 0, 0.3)[dd$Block] +
  rnorm(18, 0, 0.2)[interaction(dd$Block, dd$A)]
dd$y <- rpois(nrow(dd), exp(eta))
fp <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), family = poisson, data = dd,
                prt = FALSE))
m <- fp$semer
tg <- seq(0, 10, by = 0.25)
cs <- function(scale) q(semipred(fp, "A", "t", covariateV = tg, contr = c("a2", "a1"), unconditional = FALSE,
                                 contr_scale = scale, prt = FALSE))
pr_model <- cs("model")
pr_resp <- cs("response")
## the same from lme4's pieces: all coefficients (beta, b), V = (C'WC + P)^-1
Xm <- lme4::getME(m, "X"); Zm <- lme4::getME(m, "Z")
Lam <- lme4::getME(m, "Lambda")
b <- lme4::getME(m, "b")
etaf <- as.vector(Xm %*% lme4::fixef(m) + Zm %*% b)
fam <- poisson()
w <- fam$mu.eta(etaf)^2 / fam$variance(fam$linkinv(etaf))
C <- cbind(Xm, as.matrix(Zm))
LL <- as.matrix(Matrix::tcrossprod(Lam))
P <- matrix(0, ncol(C), ncol(C)); kx <- ncol(Xm)
P[kx + seq_len(ncol(Zm)), kx + seq_len(ncol(Zm))] <- solve(LL)
V <- solve(crossprod(C * w, C) + P)
## the contrast rows: the prediction design of a2 and a1 (fixed part and the
## smooth blocks, mapped to the columns of Z by the block names)
des <- q(pm_semi_design(fp, "A", "t", tg, NULL, fp$data, unconditional = FALSE))
D <- as.matrix(des$D_mx); g <- des$pred_df$A
gp <- lme4::getME(m, "Gp"); fl <- names(lme4::getME(m, "flist")); cn <- names(lme4::getME(m, "cnms"))
cols <- unlist(lapply(names(fp$u_lst), function(bk) { i <- match(bk, cn); kx + (gp[i] + 1):gp[i + 1] }))
Lfull <- matrix(0, sum(g == "a2"), ncol(C))
Lfull[, c(seq_len(kx), cols)] <- D[g == "a2", ] - D[g == "a1", ]
est <- as.vector(Lfull %*% c(lme4::fixef(m), as.vector(b)))
se <- sqrt(rowSums((Lfull %*% V) * Lfull))
## ---- 3. a difference of two curves of a Poisson GLMM ---------------------------
set.seed(31)
dd <- expand.grid(t = 0:10, A = factor(c("a1", "a2", "a3")), Block = factor(1:6))
eta <- 1 + sin(dd$t / 2) / 2 + 0.04 * (dd$A == "a2") * dd$t + rnorm(6, 0, 0.3)[dd$Block] +
  rnorm(18, 0, 0.2)[interaction(dd$Block, dd$A)]
dd$y <- rpois(nrow(dd), exp(eta))
fp <- q(semireg(y ~ s(t, by = A, bs = "os", k = 6) + (1 | Block) + (1 | Block:A), family = poisson, data = dd,
                prt = FALSE))
m <- fp$semer
tg <- seq(0, 10, by = 0.25)
pr <- q(semipred(fp, "A", "t", covariateV = tg, contr = c("a2", "a1"), unconditional = FALSE, prt = FALSE))$pred_df
## the same from lme4's pieces: all coefficients (beta, b), V = (C'WC + P)^-1,
## P = (Lambda Lambda')^-1 for the random part, W the working weights at the fit
Xm <- lme4::getME(m, "X"); Zm <- lme4::getME(m, "Z")
b <- lme4::getME(m, "b")
etaf <- as.vector(Xm %*% lme4::fixef(m) + Zm %*% b)
fam <- poisson()
w <- fam$mu.eta(etaf)^2 / fam$variance(fam$linkinv(etaf))
C <- cbind(Xm, as.matrix(Zm))
kx <- ncol(Xm)
P <- matrix(0, ncol(C), ncol(C))
P[kx + seq_len(ncol(Zm)), kx + seq_len(ncol(Zm))] <- solve(as.matrix(Matrix::tcrossprod(lme4::getME(m, "Lambda"))))
V <- solve(crossprod(C * w, C) + P)
## the rows of the prediction design for a2 and a1 (fixed part and the smooth
## blocks, mapped to the columns of Z by the blocks' names)
des <- q(pm_semi_design(fp, "A", "t", tg, NULL, fp$data, unconditional = FALSE))
D <- as.matrix(des$D_mx); g <- des$pred_df$A
gp <- lme4::getME(m, "Gp"); cn <- names(lme4::getME(m, "cnms"))
cols <- unlist(lapply(names(fp$u_lst), function(bk) { i <- match(bk, cn); kx + (gp[i] + 1):gp[i + 1] }))
full <- function(Drows) { L <- matrix(0, nrow(Drows), ncol(C)); L[, c(seq_len(kx), cols)] <- Drows; L }
La <- full(D[g == "a2", , drop = FALSE]); Lb <- full(D[g == "a1", , drop = FALSE])
coefs <- c(lme4::fixef(m), as.vector(b))
ea <- as.vector(La %*% coefs); eb <- as.vector(Lb %*% coefs)
se_link <- sqrt(rowSums(((La - Lb) %*% V) * (La - Lb)))
G <- exp(ea) * La - exp(eb) * Lb                       # delta method for exp(eta_a) - exp(eta_b)
se_resp <- sqrt(rowSums((G %*% V) * G))
cat("3. Poisson split plot, difference a2 - a1 of the curves over time (given the variance parameters)\n")
cat(sprintf("   log scale: difference %.1e, SE %.1e (semipred vs lme4's pieces)\n", rel(pr$diff, ea - eb), rel(pr$ses, se_link)))
cat(sprintf("   ratio exp(diff) and its limits exp(LL), exp(UL): %s\n",
            isTRUE(all.equal(c(pr$ratio, pr$ratio.LL, pr$ratio.UL), exp(c(pr$diff, pr$LL, pr$UL)), tolerance = 1e-10))))
cat(sprintf("   response scale: difference of means %.1e, its SE (delta method) %.1e\n",
            rel(pr$response.diff, exp(ea) - exp(eb)), rel(pr$response.ses, se_resp)))
crit <- (pr$UL - pr$diff) / pr$ses
cat(sprintf("   limits: diff -+ %.3f SE (t with %.1f df, the whole-plot stratum; z would be %.3f); response limits the same multiple: %s\n",
            crit[1], uniroot(function(v) qt(0.975, v) - crit[1], c(1, 1e4))$root, qnorm(0.975),
            isTRUE(all.equal(pr$response.UL - pr$response.diff, crit * pr$response.ses, tolerance = 1e-8))))
## roughly: mgcv::gam() with the same model, its own smoothing parameters and
## estimation (PQL-type REML), Bayesian covariance Vp
g1 <- mgcv::gam(y ~ A + s(t, by = A, bs = "cr", k = 8) + s(Block, bs = "re") + s(Block, A, bs = "re"),
                family = poisson, method = "REML", data = dd)
nd <- function(a) data.frame(t = tg, A = factor(a, levels = levels(dd$A)), Block = dd$Block[1])
Xa <- predict(g1, nd("a2"), type = "lpmatrix", exclude = c("s(Block)", "s(Block,A)"))
Xb <- predict(g1, nd("a1"), type = "lpmatrix", exclude = c("s(Block)", "s(Block,A)"))
re <- grepl("s\\(Block", colnames(Xa)); Xa[, re] <- 0; Xb[, re] <- 0
dg <- as.vector((Xa - Xb) %*% coef(g1)); sg <- sqrt(rowSums(((Xa - Xb) %*% g1$Vp) * (Xa - Xb)))
cat(sprintf("   mgcv::gam (other basis and smoothing parameters): difference within %.3f, SE ratio semipred / gam %.2f to %.2f\n",
            max(abs(pr$diff - dg)), min(pr$ses / sg), max(pr$ses / sg)))
