## The maintainer's golden semi-parametric examples (lme fits of penalised
## splines; scripts in dev/semireg-golden/original/, data in
## inst/extdata/semireg/): semireg() with the smooths in the formula gives the
## same fit. Reference values from real runs of the original model code
## (dev/semireg-golden/make-reference.R): REML log-likelihood, residual SD,
## fitted values. The fits are exact reproductions (lme and lmer agree to
## 1e-7), so the tolerances are tight.

golden <- function() readRDS(test_path("fixtures", "semireg-golden.rds"))
gdata <- function(f, ...) {
  p <- system.file("extdata", "semireg", f, package = "predictmeans")
  if (!nzchar(p)) p <- test_path("..", "..", "inst", "extdata", "semireg", f)
  utils::read.table(p, ...)
}
q <- function(e) suppressWarnings(suppressMessages(e))
expect_same_fit <- function(fit, r, y, label) {
  m <- fit$semer
  expect_equal(as.numeric(logLik(m)), r$logLik, tolerance = 1e-5, label = paste(label, "logLik"))
  expect_equal(sigma(m), r$sigma, tolerance = 1e-4, label = paste(label, "sigma"))
  expect_lt(max(abs(fitted(m) - r$fitted))/stats::sd(y), 1e-4)
}

test_that("golden examples: additive models, thin-plate surfaces, additive mixed model (lme)", {
  ref <- golden()
  # addMix: O'Sullivan spline + binary covariate + random intercept (Smith and Wand, 2008)
  d <- ref$addMix$data
  d$id <- factor(d$id); d$t <- d$s
  fm <- q(semireg(y ~ x + s(t, bs = "os", k = 15, xt = list(range = c(0, 1))) + (1 | id), data = d, prt = FALSE))
  expect_same_fit(fm, ref$addMix, d$y, "addMix")
  expect_equal(lme4::fixef(fm$semer)[["x"]], ref$addMix$x_effect, tolerance = 1e-5)
  expect_equal(sqrt(diag(as.matrix(vcov(fm$semer))))[["x"]], ref$addMix$x_se, tolerance = 1e-4)

  # copper: thin-plate surface with given knots + O'Sullivan spline whose covariate has fewer distinct
  # values than basis columns
  copper <- gdata("copper.txt", header = TRUE)
  kn <- as.matrix(gdata("copperKnots.txt"))
  zk <- seq(80, 120, by = 5)
  rz <- range(copper$zcoord)
  fm <- q(semireg(grade ~ s(xcoord, ycoord, bs = "ztps", xt = list(knots = kn)) +
                    s(zcoord, bs = "os", xt = list(knots = zk, range = rz)), data = copper, prt = FALSE))
  expect_same_fit(fm, ref$copper, copper$grade, "copper")

  # ozoneCalifLMM: three O'Sullivan splines, range of the data
  oz <- gdata("ozoneCalif.txt", header = TRUE)
  r1 <- range(oz$inversion.base.height); r2 <- range(oz$daggett.pressure.gradient)
  r3 <- range(oz$inversion.base.temp)
  fm <- q(semireg(ozone ~ s(inversion.base.height, bs = "os", k = 15, xt = list(range = r1)) +
                    s(daggett.pressure.gradient, bs = "os", k = 15, xt = list(range = r2)) +
                    s(inversion.base.temp, bs = "os", k = 15, xt = list(range = r3)), data = oz, prt = FALSE))
  expect_same_fit(fm, ref$ozoneCalifLMM, oz$ozone, "ozoneCalifLMM")

  # scallop: thin-plate surface of asinh(catch)
  sc <- gdata("scallop.txt", header = TRUE)
  sc$z <- asinh(sc$tot.catch)
  sk <- as.matrix(gdata("scallopKnots.txt"))
  fm <- q(semireg(z ~ s(longitude, latitude, bs = "ztps", xt = list(knots = sk)), data = sc, prt = FALSE))
  expect_same_fit(fm, ref$scallop, sc$z, "scallop")
})

test_that("golden example growthIndiana: mean curve + subject curves (s(age) | idnum), as lme", {
  skip_on_cran()
  ref <- golden()
  gi <- gdata("growthIndiana.txt", header = TRUE)
  d <- gi[gi$male == 1 & gi$black == 1, ]
  d$idnum <- factor(match(d$idnum, unique(d$idnum)))
  fm <- q(semireg(height ~ s(age, bs = "os", k = 20, xt = list(range = c(5.5, 20))) +
                    (s(age, k = 10, xt = list(range = c(5.5, 20))) | idnum), data = d, prt = FALSE))
  expect_same_fit(fm, ref$growthIndiana, d$height, "growthIndiana")
  # (age | idnum) is a correlated intercept and slope per boy, as pdSymm(~age)
  vc <- as.data.frame(lme4::VarCorr(fm$semer))
  expect_true(any(vc$grp == "idnum" & vc$var1 == "(Intercept)" & vc$var2 %in% "age"))
  # the subject curves are random effects: the population curve leaves them out
  sp <- q(semipred(fm, "age", "age", prt = FALSE))
  expect_true(all(is.finite(sp$pred_df$SE)))
})

test_that("golden examples fitted by gamm() and gam(): semireg() equals gam(method = 'REML')", {
  ref <- golden()
  # SBMDfem: gamm(s(age) + ethnicity, random intercept); gamm's lme stops short of the REML optimum
  sb <- gdata("SBMDfem.txt", header = TRUE)
  sb$idnum <- factor(sb$idnum)
  fm <- q(semireg(spnbmd ~ s(age) + black + hispanic + white + (1 | idnum), data = sb, prt = FALSE))
  g <- mgcv::gam(spnbmd ~ s(age) + black + hispanic + white + s(idnum, bs = "re"), data = sb, method = "REML")
  expect_lt(max(abs(fitted(fm$semer) - fitted(g)))/stats::sd(sb$spnbmd), 1e-5)
  expect_lt(max(abs(fitted(fm$semer) - ref$SBMDfem$fitted))/stats::sd(sb$spnbmd), 5e-3)
  # ozoneCalifGAM: the script's gam uses GCV; semireg is REML
  oz <- gdata("ozoneCalif.txt", header = TRUE)
  fm <- q(semireg(ozone ~ s(inversion.base.height) + s(daggett.pressure.gradient) + s(inversion.base.temp),
                  data = oz, prt = FALSE))
  expect_lt(max(abs(fitted(fm$semer) - ref$ozoneCalifGAM$fitted_reml))/stats::sd(oz$ozone), 1e-5)
})

test_that("the package data sets of the golden examples equal their files (dev/semireg-golden/make-data.R)", {
  get <- function(nm) { e <- new.env(); utils::data(list = nm, package = "predictmeans", envir = e); e[[nm]] }
  same <- function(x, y) expect_equal(unname(as.list(x)), unname(as.list(y)), tolerance = 1e-12)
  same(get("growthIndiana"), gdata("growthIndiana.txt", header = TRUE))
  sb <- gdata("SBMDfem.txt", header = TRUE)
  sb$ethnicity <- factor(sb$ethnicity)
  same(get("SBMDfem"), sb)
  same(get("ozoneCalif"), gdata("ozoneCalif.txt", header = TRUE))
  same(get("scallop"), gdata("scallop.txt", header = TRUE))
  same(get("scallopKnots"), gdata("scallopKnots.txt"))
  same(get("scallopBdry"), gdata("scallopBdry.txt"))
  same(get("copperKnots"), gdata("copperKnots.txt"))
  same(get("wageTUnion"), gdata("wageTUnion.txt", header = TRUE))
  same(get("copper"), gdata("copper.txt", header = TRUE))
})
