## Reference values of the golden semi-parametric examples (the maintainer's
## lme examples, in original/; data in inst/extdata/semireg/). Each model is
## fitted exactly as in its original script (the scripts' own ZOSull.r and
## Ztps.r, nlme::lme, MASS::glmmPQL, mgcv::gamm / gam), without the plots.
## Saves log-likelihood, sigma, fitted values, key coefficients and the time
## of each fit to tests/testthat/fixtures/semireg-golden.rds.
##
##   Rscript dev/semireg-golden/make-reference.R
suppressMessages({library(nlme); library(mgcv); library(MASS)})
here <- "dev/semireg-golden"
dat <- function(f, ...) read.table(file.path("inst/extdata/semireg", f), ...)
orig <- new.env()
sys.source(file.path(here, "original/ZOSull.r"), orig)
sys.source(file.path(here, "original/Ztps.r"), orig)
ZOSull <- orig$ZOSull
Ztps <- orig$Ztps
out <- list()
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "time") <- proc.time()[["elapsed"]] - t0; v }
keep <- function(fit, extra = list()) {
  c(list(logLik = as.numeric(logLik(fit)), sigma = fit$sigma, fitted = unname(fitted(fit)),
         fixed = fixef(fit), time = attr(fit, "time")), extra)
}

## 1. addMix: additive mixed model, simulated (Smith and Wand, 2008)
set.seed(39402); m <- 250; nVals <- sample(1:4, m, replace = TRUE)
betaVal <- 0.3; sigU <- 0.5; sigEps <- 0.2
f <- function(x) -sin(2*pi*x)
U <- rnorm(m, 0, sigU)
sVals <- NULL; xVals <- NULL; Uvals <- NULL; idNum <- NULL
for (i in 1:m) {
  idNum <- c(idNum, rep(i, nVals[i]))
  stt <- runif(1, 0, 1 - 0.05*(nVals[i] - 1))
  sVals <- c(sVals, seq(stt, by = 0.05, length = nVals[i]))
  xCurr <- sample(c(0, 1), 1)
  xVals <- c(xVals, rep(xCurr, nVals[i]))
  Uvals <- c(Uvals, rep(U[i], nVals[i]))
}
yVals <- f(sVals) + betaVal*xVals + Uvals + rnorm(sum(nVals), 0, sigEps)
addMix <- data.frame(y = yVals, s = sVals, x = xVals, id = idNum)
intKnots <- quantile(unique(sVals), seq(0, 1, length = 17)[-c(1, 17)])
X <- cbind(1, sVals, xVals)
ZSpline <- ZOSull(sVals, c(0, 1), intKnots)
groupVec <- factor(rep(1, length(yVals)))
dataFr_lme <- groupedData(yVals ~ xVals | groupVec, data = data.frame(yVals, X, ZSpline, idNum))
fit <- tm(lme(yVals ~ -1 + X, data = dataFr_lme, random = list(groupVec = pdIdent(~ ZSpline - 1), idNum = pdIdent(~ 1))))
# the script's streamlined covariance of (beta, u)
p <- 1; numKnots <- ncol(ZSpline)
sig.eps.hat <- fit$sigma
sig.u.hat <- intervals(fit)$reStruct$groupVec$est
sig.U.hat <- intervals(fit)$reStruct$idNum$est
M11inv <- diag(c(rep(0, p + 2), rep(sig.eps.hat^2/sig.u.hat^2, numKnots))) + crossprod(cbind(X, ZSpline))
M.21 <- rowsum(cbind(X, ZSpline), idNum)
for (i in 1:m) M11inv <- M11inv - sig.U.hat^2*tcrossprod(M.21[i, ])/(sig.U.hat^2*M.21[i, 1] + sig.eps.hat^2)
covMat <- sig.eps.hat^2*solve(M11inv)
out$addMix <- keep(fit, list(data = addMix, x_effect = unname(fixef(fit)[3]), x_se = sqrt(covMat[3, 3])))

## 2. copper: geoadditive model, thin-plate (xcoord, ycoord) + O'Sullivan zcoord
copper <- dat("copper.txt", header = TRUE)
copperKnots <- as.matrix(dat("copperKnots.txt"))
n <- nrow(copper)
Zgeog <- Ztps(cbind(copper$xcoord, copper$ycoord), copperKnots)
Zz <- ZOSull(copper$zcoord, range(copper$zcoord), seq(80, 120, by = 5))
X <- cbind(1, copper$xcoord, copper$ycoord, copper$zcoord)
Z <- cbind(Zgeog, Zz)
blk <- list(1:ncol(Zgeog), ncol(Zgeog) + 1:ncol(Zz))
Z.block <- lapply(blk, function(b) as.formula(paste0("~Z[,c(", paste(b, collapse = ","), ")]-1")))
group <- rep(1, n); grade <- copper$grade
dataFr <- groupedData(grade ~ -1 + X | group, data = data.frame(grade, X, Z))
fit <- tm(lme(grade ~ -1 + X, random = list(group = pdBlocked(Z.block, pdClass = rep("pdIdent", 2))), data = dataFr))
out$copper <- keep(fit)

## 3. growthIndiana: black males, mean curve + subject curves ((age | id) + spline)
gi <- dat("growthIndiana.txt", header = TRUE)
gbm <- gi[gi$male == 1 & gi$black == 1, ]
gbm$idnum <- match(gbm$idnum, unique(gbm$idnum))
age <- gbm$age; height <- gbm$height; idnum <- gbm$idnum
numObs <- length(height)
X <- model.matrix(height ~ age)
range.age <- c(5.5, 20)
Z <- ZOSull(age, range.age, quantile(unique(age), seq(0, 1, length = 22))[-c(1, 22)])
ZSubj <- ZOSull(age, range.age, quantile(unique(age), seq(0, 1, length = 12))[-c(1, 12)])
dummyId <- factor(rep(1, numObs))
Z.block <- list(dummyId = pdIdent(~ -1 + Z), idnum = pdSymm(~ age), idnum = pdIdent(~ -1 + ZSubj))
data.fr <- groupedData(height ~ X[, -1] | rep(1, length = numObs), data = data.frame(height, X, Z, ZSubj, idnum))
fit <- tm(lme(height ~ -1 + X, data = data.fr, random = Z.block))
out$growthIndiana <- keep(fit)

## 4. growthIndianaCont: all males, mean curve + black-white contrast curve +
## subject curves (as HRW's demo maleGrowthIndiananlme, scaled by the range)
gm <- gi[gi$male == 1, ]
gm$idnum <- match(gm$idnum, unique(gm$idnum))
gm$age01 <- (gm$age - min(gm$age))/diff(range(gm$age))
gm$height01 <- (gm$height - min(gm$height))/diff(range(gm$height))
age <- gm$age01; height <- gm$height01; black <- gm$black; idnum <- gm$idnum
numObs <- length(height)
X <- cbind(rep(1, numObs), age, black, age*black)
range.age <- c(1.01*min(age) - 0.01*max(age), 1.01*max(age) - 0.01*min(age))
Zbase <- ZOSull(age, range.age, quantile(unique(age), seq(0, 1, length = 22))[-c(1, 22)])
Zcont <- black*Zbase
ZSubj <- ZOSull(age, range.age, quantile(unique(age), seq(0, 1, length = 12))[-c(1, 12)])
dummyId <- factor(rep(1, numObs))
Z.block <- list(dummyId = pdIdent(~ -1 + Zbase), dummyId = pdIdent(~ -1 + Zcont),
                idnum = pdSymm(~ age), idnum = pdIdent(~ -1 + ZSubj))
data.fr <- groupedData(height ~ X[, -1] | rep(1, length = numObs), data = data.frame(height, X, idnum))
fit <- tm(lme(height ~ -1 + X, data = data.fr, random = Z.block))
out$growthIndianaCont <- keep(fit)

## 5. ozoneCalifLMM: additive model, three O'Sullivan smooths
oz <- dat("ozoneCalif.txt", header = TRUE)
y <- oz$ozone
xs <- oz[c("inversion.base.height", "daggett.pressure.gradient", "inversion.base.temp")]
X <- cbind(1, as.matrix(xs))
Zs <- lapply(xs, function(x) ZOSull(x, range(x), quantile(unique(x), seq(0, 1, length = 17))[-c(1, 17)]))
Z <- do.call(cbind, Zs)
K <- cumsum(c(0, sapply(Zs, ncol)))
Z.block <- lapply(1:3, function(i) as.formula(paste0("~Z[,c(", paste((K[i] + 1):K[i + 1], collapse = ","), ")]-1")))
group <- rep(1, length(y))
data.fr <- groupedData(y ~ -1 + X | group, data = data.frame(y, X, Z))
fit <- tm(lme(y ~ -1 + X, random = list(group = pdBlocked(Z.block, pdClass = rep("pdIdent", 3))), data = data.fr))
out$ozoneCalifLMM <- keep(fit)

## 6. scallop: thin-plate surface of asinh(catch)
sc <- dat("scallop.txt", header = TRUE)
scallopKnots <- as.matrix(dat("scallopKnots.txt"))
x <- sc$longitude; y <- sc$latitude; z <- asinh(sc$tot.catch)
n <- length(z)
X <- cbind(rep(1, n), x, y)
Z <- Ztps(cbind(x, y), scallopKnots)
group <- rep(1, n)
data.fr <- groupedData(y ~ x | group, data = data.frame(x, y))
fit <- tm(lme(z ~ -1 + X, random = pdIdent(~ -1 + Z), data = data.fr))
out$scallop <- keep(fit)

## 7. wageTUnionLMM: binary response, penalised quasi-likelihood (glmmPQL)
wu <- dat("wageTUnion.txt", header = TRUE)
x <- wu[, 1]; y <- wu[, 2]
X <- cbind(rep(1, length(x)), x)
Z <- ZOSull(x, range(x), quantile(unique(x), seq(0, 1, length = 17))[-c(1, 17)])
group <- rep(1, length(x))
data.fr <- groupedData(y ~ x | group, data = data.frame(x, y))
fit <- tm(suppressMessages(glmmPQL(y ~ -1 + X, random = list(group = pdIdent(~ -1 + Z)), data = data.fr,
                                   family = binomial, verbose = FALSE)))
out$wageTUnion <- list(fitted_prob = plogis(unname(fitted(fit))),   # lme fitted(): the link scale
                        fixed = fixef(fit), time = attr(fit, "time"),
                       sd_u = as.numeric(fit$sigma*exp(unlist(fit$modelStruct))))

## 8. SBMDfemGCV: gamm with a random intercept per girl; the script's
## likelihood ratio test of ethnicity
SBMDfem <- dat("SBMDfem.txt", header = TRUE)
fit <- tm(gamm(spnbmd ~ s(age) + black + hispanic + white, random = list(idnum = ~ 1), data = SBMDfem))
fitReduced <- gamm(spnbmd ~ s(age), random = list(idnum = ~ 1), data = SBMDfem)
lr <- 2*(fit$lme$logLik - fitReduced$lme$logLik)
out$SBMDfem <- list(logLik = fit$lme$logLik, sigma = fit$lme$sigma, fitted = unname(fitted(fit$lme)),
                    fixed = coef(fit$gam)[c("black", "hispanic", "white")], time = attr(fit, "time"),
                    lrt = lr, lrt_p = 1 - pchisq(lr, df = 3))

## 9. ozoneCalifGAM: gam with the default criterion (GCV) and with REML
fit <- tm(gam(ozone ~ s(inversion.base.height) + s(daggett.pressure.gradient) + s(inversion.base.temp), data = oz))
fitR <- gam(ozone ~ s(inversion.base.height) + s(daggett.pressure.gradient) + s(inversion.base.temp), data = oz,
            method = "REML")
out$ozoneCalifGAM <- list(fitted_gcv = unname(fitted(fit)), fitted_reml = unname(fitted(fitR)),
                          edf_gcv = sum(fit$edf), edf_reml = sum(fitR$edf), time = attr(fit, "time"))

dir.create("tests/testthat/fixtures", showWarnings = FALSE)
saveRDS(out, "tests/testthat/fixtures/semireg-golden.rds", version = 2)
print(sapply(out, function(o) c(logLik = if (is.null(o$logLik)) NA else o$logLik, time = o$time)))
