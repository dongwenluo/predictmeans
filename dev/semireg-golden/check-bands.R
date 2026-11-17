## Bands of semipred() against M. P. Wand's own calculation
## (dev/semireg-golden/original/growthIndianaCont.R): the contrast curve
## (black - white mean height of the Indiana boys) with its band
## qnorm(0.975) * sigma * sqrt(diag(Cg (C'C + Lambda)^-1 Cg')), the covariance
## of the fixed and spline coefficients given the variance parameters.
##   A. Wand's model (base curve + black x contrast curve + subject curves),
##      fitted by lme() as in the script and by semireg() (s(age, by = black)
##      with black numeric); contrast and band from each fit.
##   B. semipred(contr = ) itself, which needs a factor: a curve per group
##      (s(age, by = blackF)), fitted by lme() (one pdIdent block per group,
##      the band by Wand's formula) and by semireg(); the difference of the
##      two curves and its SE from semipred(contr = c("1", "0")).
##   Rscript dev/semireg-golden/check-bands.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(nlme))
q <- function(e) suppressWarnings(suppressMessages(e))
gi <- read.table("inst/extdata/semireg/growthIndiana.txt", header = TRUE)
gm <- gi[gi$male == 1, ]
idnum <- match(gm$idnum, unique(gm$idnum))
age <- (gm$age - min(gm$age)) / diff(range(gm$age))
height <- (gm$height - min(gm$height)) / diff(range(gm$height))
black <- gm$black
n <- length(height); nsub <- max(idnum)
range.age <- c(1.01 * min(age) - 0.01 * max(age), 1.01 * max(age) - 0.01 * min(age))
kn20 <- quantile(unique(age), seq(0, 1, length = 22))[-c(1, 22)]
kn10 <- quantile(unique(age), seq(0, 1, length = 12))[-c(1, 12)]
Z <- ZOSull(age, range.age, kn20)
ZSubj <- ZOSull(age, range.age, kn10)
ageg <- seq(range.age[1], range.age[2], length = 101)
Xg <- cbind(1, ageg)
Zg <- ZOSull(ageg, range.age, kn20)
## the subject columns of C (intercept and slope, then the spline) and the
## penalty of the subject part, as in the script
subj_cols <- function() {
  L <- matrix(0, n, 2 * nsub); S <- matrix(0, n, ncol(ZSubj) * nsub)
  for (i in 1:nsub) {
    r <- idnum == i
    L[r, 2 * i - 1:0] <- cbind(1, age[r])
    S[r, (i - 1) * ncol(ZSubj) + seq_len(ncol(ZSubj))] <- ZSubj[r, ]
  }
  cbind(L, S)
}
Csub <- subj_cols()
## Wand's band from an lme fit: blocks = penalty blocks for the spline columns
## (in the order of 'cols'), the subject part from the last two variance
## parameters; L = rows of the contrast for (fixed, splines)
wand_band <- function(fit, Cfix, Zcols, lam_spl, L) {
  s2 <- fit$sigma^2
  sig.u <- as.numeric(fit$sigma * exp(unlist(fit$modelStruct)))
  Sig <- crossprod(rbind(c(sig.u[2], sig.u[3]), c(0, sig.u[4])))    # pdSymm(~age)
  pen <- c(rep(0, ncol(Cfix)), unlist(lapply(seq_along(Zcols), function(j) rep(s2 / lam_spl[j]^2, ncol(Zcols[[j]])))))
  Lam <- diag(c(pen, rep(0, ncol(Csub))))
  k1 <- length(pen)
  Lam[k1 + 1:(2 * nsub), k1 + 1:(2 * nsub)] <- s2 * kronecker(diag(nsub), solve(Sig))
  ks <- k1 + 2 * nsub
  diag(Lam)[ks + seq_len(ncol(ZSubj) * nsub)] <- s2 / sig.u[1]^2
  C <- cbind(Cfix, do.call(cbind, Zcols), Csub)
  V <- solve(crossprod(C) + Lam)[1:k1, 1:k1]
  sqrt(s2) * sqrt(rowSums((L %*% V) * L))
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

## ---- A. Wand's model -------------------------------------------------------
X <- cbind(1, age, black, age * black)
Zbase <- Z; Zcont <- black * Z
dummyId <- factor(rep(1, n))
data.fr <- groupedData(height ~ X[, -1] | rep(1, length = n), data = data.frame(height, X, idnum))
t0 <- proc.time()[[3]]
fitA <- lme(height ~ -1 + X, data = data.fr,
            random = list(dummyId = pdIdent(~ -1 + Zbase), dummyId = pdIdent(~ -1 + Zcont),
                          idnum = pdSymm(~ age), idnum = pdIdent(~ -1 + ZSubj)))
cat(sprintf("A. lme fit %.0f s\n", proc.time()[[3]] - t0))
sig.uA <- as.numeric(fitA$sigma * exp(unlist(fitA$modelStruct)))   # (subject spline, pdSymm x3, contrast, base)
ContA <- Xg %*% fitA$coef$fixed[3:4] + Zg %*% as.vector(fitA$coef$random[[2]])
LA <- cbind(0, 0, Xg, matrix(0, 101, ncol(Z)), Zg)
sdA <- wand_band(fitA, X, list(Zbase, Zcont), c(sig.uA[6], sig.uA[5]), LA)

d <- data.frame(height = height, age = age, black = black, idnum = factor(idnum))
fmA <- q(semireg(height ~ s(age, bs = "os", k = 20, xt = list(range = range.age)) +
                   s(age, by = black, bs = "os", k = 20, xt = list(range = range.age)) +
                   (s(age, bs = "os", k = 10, xt = list(range = range.age)) | idnum), data = d, prt = FALSE))
cat(sprintf("   logLik lme %.4f, semireg %.4f\n", as.numeric(logLik(fitA)), as.numeric(logLik(fmA$semer))))
b <- lme4::fixef(fmA$semer)
blk <- grep("black", names(fmA$u_lst), value = TRUE)
Zc <- do.call(cbind, pm_sm_predict(fmA$smooths, data.frame(age = ageg, black = 1), blk))
Lf <- matrix(0, 101, length(b), dimnames = list(NULL, names(b)))
Lf[, grep("^black$", names(b))] <- 1
Lf[, grep("age:black|black:age", names(b))] <- ageg
ci <- fmA$Cov_ind
idx <- c(seq_along(b), (ci[match(blk, names(ci)) - 1] + 1):ci[match(blk, names(ci))])
Lc <- cbind(Lf, Zc)
ContS <- Lc %*% c(b, fmA$u_lst[[blk]])
sdS <- sigma(fmA$semer) * sqrt(rowSums((Lc %*% as.matrix(fmA$CovMat[idx, idx])) * Lc))
cat(sprintf("   contrast curve: max |semireg - lme| / max |lme| = %.1e\n", rel(ContS, ContA)))
cat(sprintf("   band SE:        max |semireg - lme| / max |lme| = %.1e\n", rel(sdS, sdA)))

## ---- B. a curve per group: semipred(contr = ) ------------------------------
blackF <- factor(black)
XB <- cbind(1, age, black, age * black)              # = (white curve, black - white) fixed part
Zw <- (1 - black) * Z; Zb <- black * Z
data.fr$Zw <- NULL
t0 <- proc.time()[[3]]
fitB <- lme(height ~ -1 + X, data = data.fr,
            random = list(dummyId = pdIdent(~ -1 + Zw), dummyId = pdIdent(~ -1 + Zb),
                          idnum = pdSymm(~ age), idnum = pdIdent(~ -1 + ZSubj)))
cat(sprintf("B. lme fit %.0f s\n", proc.time()[[3]] - t0))
sig.uB <- as.numeric(fitB$sigma * exp(unlist(fitB$modelStruct)))   # (subject spline, pdSymm x3, black, white)
uW <- as.vector(fitB$coef$random[[1]]); uB <- as.vector(fitB$coef$random[[2]])
ContB <- Xg %*% fitB$coef$fixed[3:4] + Zg %*% (uB - uW)
LB <- cbind(0, 0, Xg, -Zg, Zg)
sdB <- wand_band(fitB, XB, list(Zw, Zb), c(sig.uB[6], sig.uB[5]), LB)

d$blackF <- blackF
fmB <- q(semireg(height ~ s(age, by = blackF, bs = "os", k = 20, xt = list(range = range.age)) +
                   (s(age, bs = "os", k = 10, xt = list(range = range.age)) | idnum), data = d, prt = FALSE))
cat(sprintf("   logLik lme %.4f, semireg %.4f\n", as.numeric(logLik(fitB)), as.numeric(logLik(fmB$semer))))
sp <- q(semipred(fmB, "blackF", "age", contr = c("1", "0"), covariateV = ageg, prt = FALSE))$pred_df
cat(sprintf("   contrast curve: max |semipred - lme| / max |lme| = %.1e\n", rel(sp$height, ContB)))
cat(sprintf("   band SE:        max |semipred - lme| / max |lme| = %.1e\n", rel(sp$SE_contr, sdB)))
cat(sprintf("   limits: semipred t with df %.0f-%.0f (qt %.4f-%.4f) vs Wand's qnorm %.4f;",
            min(sp$df), max(sp$df), qt(0.975, max(sp$df)), qt(0.975, min(sp$df)), qnorm(0.975)),
    sprintf("UL - estimate = qt(0.975, df) x SE: %s\n",
            isTRUE(all.equal(sp$UL - sp$height, qt(0.975, sp$df) * sp$SE_contr, tolerance = 1e-8))))
