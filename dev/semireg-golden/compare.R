## The golden examples refitted with semireg() (smooths in the formula, lme4
## engine) and compared with the reference fits of make-reference.R:
## REML log-likelihood, residual SD, fitted values (all random effects
## included), key coefficients, time.
##
##   Rscript dev/semireg-golden/compare.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
ref <- readRDS("tests/testthat/fixtures/semireg-golden.rds")
dat <- function(f, ...) read.table(file.path("inst/extdata/semireg", f), ...)
tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; attr(v, "time") <- proc.time()[["elapsed"]] - t0; v }
## a gaussian fit, with semireg()'s default optimizer (nloptwrap since 2.0)
## and again with 1.x's bobyqa (time and log-likelihood kept)
tm2 <- function(expr) {
  cl <- substitute(expr)
  fit <- tm(eval(cl, parent.frame()))
  cl$control <- quote(lme4::lmerControl(optimizer = "bobyqa"))
  fit2 <- tm(eval(cl, parent.frame()))
  attr(fit, "time_bobyqa") <- attr(fit2, "time")
  attr(fit, "logLik_bobyqa") <- as.numeric(logLik(fit2$semer))
  fit
}
row <- function(case, fm, r, y, note = "") {
  m <- fm$semer
  data.frame(case = case, logLik = as.numeric(logLik(m)), logLik_ref = if (is.null(r$logLik)) NA else r$logLik,
             sigma_ratio = if (is.null(r$sigma)) NA else sigma(m)/r$sigma,
             fitted_diff = max(abs(fitted(m) - r$fitted))/sd(y), time = attr(fm, "time"),
             time_bobyqa = if (is.null(attr(fm, "time_bobyqa"))) NA else attr(fm, "time_bobyqa"),
             dlogLik_bobyqa = if (is.null(attr(fm, "logLik_bobyqa"))) NA else attr(fm, "logLik_bobyqa") - as.numeric(logLik(m)),
             time_ref = r$time, note = note)
}
res <- list()

## 1. addMix
d <- ref$addMix$data
d$id <- factor(d$id); d$t <- d$s
fm <- tm2(semireg(y ~ x + s(t, bs = "os", k = 15, xt = list(range = c(0, 1))) + (1 | id), data = d, prt = FALSE))
res$addMix <- row("addMix", fm, ref$addMix, d$y,
                  sprintf("x effect %.4f (SE %.4f); script %.4f (%.4f)", lme4::fixef(fm$semer)[["x"]],
                          sqrt(diag(as.matrix(vcov(fm$semer))))[["x"]], ref$addMix$x_effect, ref$addMix$x_se))

## 2. copper
copper <- dat("copper.txt", header = TRUE)
kn <- as.matrix(dat("copperKnots.txt")); zk <- seq(80, 120, by = 5); rz <- range(copper$zcoord)
fm <- tm2(semireg(grade ~ s(xcoord, ycoord, bs = "ztps", xt = list(knots = kn)) +
                   s(zcoord, bs = "os", xt = list(knots = zk, range = rz)), data = copper, prt = FALSE))
res$copper <- row("copper", fm, ref$copper, copper$grade)

## 3. growthIndiana (black males)
gi <- dat("growthIndiana.txt", header = TRUE)
d <- gi[gi$male == 1 & gi$black == 1, ]; d$idnum <- factor(match(d$idnum, unique(d$idnum)))
fm <- tm2(semireg(height ~ s(age, bs = "os", k = 20, xt = list(range = c(5.5, 20))) +
                   (s(age, k = 10, xt = list(range = c(5.5, 20))) | idnum), data = d, prt = FALSE))
res$growthIndiana <- row("growthIndiana", fm, ref$growthIndiana, d$height)

## 4. growthIndianaCont (all males, contrast curve)
d <- gi[gi$male == 1, ]; d$idnum <- factor(match(d$idnum, unique(d$idnum)))
d$age01 <- (d$age - min(d$age))/diff(range(d$age)); d$height01 <- (d$height - min(d$height))/diff(range(d$height))
fm <- tm2(semireg(height01 ~ s(age01, bs = "os", k = 20) + s(age01, by = black, bs = "os", k = 20) +
                   (s(age01, k = 10) | idnum), data = d, prt = FALSE))
res$growthIndianaCont <- row("growthIndianaCont", fm, ref$growthIndianaCont, d$height01)

## 5. ozoneCalifLMM
oz <- dat("ozoneCalif.txt", header = TRUE)
r1 <- range(oz$inversion.base.height); r2 <- range(oz$daggett.pressure.gradient); r3 <- range(oz$inversion.base.temp)
fm <- tm2(semireg(ozone ~ s(inversion.base.height, bs = "os", k = 15, xt = list(range = r1)) +
                   s(daggett.pressure.gradient, bs = "os", k = 15, xt = list(range = r2)) +
                   s(inversion.base.temp, bs = "os", k = 15, xt = list(range = r3)), data = oz, prt = FALSE))
res$ozoneCalifLMM <- row("ozoneCalifLMM", fm, ref$ozoneCalifLMM, oz$ozone)

## 6. scallop
sc <- dat("scallop.txt", header = TRUE); sc$z <- asinh(sc$tot.catch)
sk <- as.matrix(dat("scallopKnots.txt"))
fm <- tm2(semireg(z ~ s(longitude, latitude, bs = "ztps", xt = list(knots = sk)), data = sc, prt = FALSE))
res$scallop <- row("scallop", fm, ref$scallop, sc$z)

## 7. wageTUnion: Laplace (glmer) here, penalised quasi-likelihood in the script
wu <- dat("wageTUnion.txt", header = TRUE)
rw <- range(wu$wage)
fm <- tm(semireg(unionMember ~ s(wage, bs = "os", k = 15, xt = list(range = rw)), family = binomial, data = wu,
                 prt = FALSE))
res$wageTUnion <- data.frame(case = "wageTUnion", logLik = as.numeric(logLik(fm$semer)), logLik_ref = NA,
                             sigma_ratio = NA, fitted_diff = max(abs(fitted(fm$semer) - ref$wageTUnion$fitted_prob)),
                             time = attr(fm, "time"), time_bobyqa = NA, dlogLik_bobyqa = NA,
                             time_ref = ref$wageTUnion$time,
                             note = "Laplace vs PQL; fitted_diff on the probability scale")

## 8. SBMDfem: gamm(s(age) + ethnicity, random intercept)
sb <- dat("SBMDfem.txt", header = TRUE); sb$idnum <- factor(sb$idnum)
sb$ethnicity <- factor(sb$ethnicity)
fm <- tm2(semireg(spnbmd ~ s(age) + black + hispanic + white + (1 | idnum), data = sb, prt = FALSE))
fe <- suppressMessages(semireg(spnbmd ~ s(age) + ethnicity + (1 | idnum), data = sb, prt = FALSE))
a <- anova(fe, smooth = FALSE)$fixed
res$SBMDfem <- row("SBMDfem", fm, ref$SBMDfem, sb$spnbmd,
                   sprintf("ethnicity: KR F p = %.3g; script's REML LRT p = %.3g", a["ethnicity", "Pr(>F)"],
                           ref$SBMDfem$lrt_p))

## 9. ozoneCalifGAM: semireg (REML) vs gam REML and gam GCV
fm <- tm2(semireg(ozone ~ s(inversion.base.height) + s(daggett.pressure.gradient) + s(inversion.base.temp),
                 data = oz, prt = FALSE))
res$ozoneCalifGAM <- data.frame(case = "ozoneCalifGAM", logLik = as.numeric(logLik(fm$semer)), logLik_ref = NA,
                                sigma_ratio = NA,
                                fitted_diff = max(abs(fitted(fm$semer) - ref$ozoneCalifGAM$fitted_reml))/sd(oz$ozone),
                                time = attr(fm, "time"), time_bobyqa = attr(fm, "time_bobyqa"),
                                dlogLik_bobyqa = attr(fm, "logLik_bobyqa") - as.numeric(logLik(fm$semer)),
                                time_ref = ref$ozoneCalifGAM$time,
                                note = sprintf("vs gam REML; vs gam GCV: %.3f",
                                               max(abs(fitted(fm$semer) - ref$ozoneCalifGAM$fitted_gcv))/sd(oz$ozone)))

tab <- do.call(rbind, res)
rownames(tab) <- NULL
print(format(tab, digits = 6), right = FALSE)
saveRDS(tab, file.path(tempdir(), "semireg-golden-compare.rds"))
