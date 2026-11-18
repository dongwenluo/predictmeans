# The 1.x smoothZ = list(smZ(...)) forms of semireg()/semireg_tmb() against
# their 2.0 formula forms (vignette("migrating")): the REML log-likelihood,
# the fitted values and the semipred() curves (estimates and SEs).
# Run on a build that still has smZ() (before its removal, 2026-11-18):
#   Rscript --no-init-file dev/check-smz-migration.R <library>
lib <- commandArgs(TRUE)[1]
suppressPackageStartupMessages(library(predictmeans, lib.loc = lib))
q <- function(e) suppressMessages(suppressWarnings(e))
set.seed(11)
n <- 300
d <- data.frame(x = runif(n, 0, 3), z = runif(n, 0, 2), f = gl(2, n/2, labels = c("a", "b")), id = gl(25, 1, n))
d$y <- sin(2*d$x) + (d$f == "b")*d$x/2 + 0.5*d$z^2 + rnorm(25, 0, 0.3)[d$id] + rnorm(n, sd = 0.3)
d$cnt <- rpois(n, exp(0.5 + sin(d$x)/2 + 0.3*(d$f == "b")))
kn <- c(0.5, 1, 1.5, 2, 2.5); r <- c(-0.05, 3.05)
bs_kn <- attr(splines::bs(d$x, df = 6), "knots"); bs_r <- range(d$x)
xs <- c(0.3, 1, 2, 2.7)
cases <- list(
  default = list(quote(semireg(y ~ x, smoothZ = list(sm = smZ(x)), data = d)),
                 quote(semireg(y ~ s(x, bs = "os"), data = d)), NULL, "x"),
  k10 = list(quote(semireg(y ~ x, smoothZ = list(sm = smZ(x, k = 10)), data = d)),
             quote(semireg(y ~ s(x, bs = "os", k = 10), data = d)), NULL, "x"),
  by_group = list(quote(semireg(y ~ x*f, smoothZ = list(sm = smZ(x, by = f, group = TRUE)), data = d)),
                  quote(semireg(y ~ s(x, by = f, bs = "os"), data = d)), "f", "x"),
  by_one = list(quote(semireg(y ~ x*f, smoothZ = list(sm = smZ(x, by = f)), data = d)),
                quote(semireg(y ~ s(x, by = f, bs = "os", id = 1), data = d)), "f", "x"),
  two_smooths = list(quote(semireg(y ~ x + z, smoothZ = list(sx = smZ(x), sz = smZ(z)), data = d)),
                     quote(semireg(y ~ s(x, bs = "os") + s(z, bs = "os"), data = d)), NULL, "x"),
  knots = list(quote(semireg(y ~ x, smoothZ = list(sm = smZ(x, intKnots = kn, range.x = r)), data = d)),
               quote(semireg(y ~ s(x, bs = "os", xt = list(knots = kn, range = r)), data = d)), NULL, "x"),
  random = list(quote(semireg(y ~ x + (1 | id), smoothZ = list(sm = smZ(x)), data = d)),
                quote(semireg(y ~ s(x, bs = "os") + (1 | id), data = d)), NULL, "x"),
  ztps = list(quote(semireg(y ~ x + z, smoothZ = list(sm = smZ(cbind(x, z), type = "Ztps")), data = d)),
              quote(semireg(y ~ s(x, z, bs = "ztps", k = 6), data = d)), NULL, NULL),
  poisson = list(quote(semireg(cnt ~ x*f, smoothZ = list(sm = smZ(x, by = f, group = TRUE)), family = poisson,
                               data = d)),
                 quote(semireg(cnt ~ s(x, by = f, bs = "os"), family = poisson, data = d)), "f", "x"),
  bs = list(quote(semireg(y ~ x, smoothZ = list(sm = smZ(x, type = "bs", k = 6)), data = d)),
            quote(semireg(y ~ s(x, bs = "user", xt = list(basis = function(x)
              splines::bs(x, knots = bs_kn, Boundary.knots = bs_r))), data = d)), NULL, "x"),
  tmb = list(quote(semireg_tmb(cnt ~ x*f, smoothZ = list(sm = smZ(x, by = f, group = TRUE)), family = poisson,
                               data = d)),
             quote(semireg_tmb(cnt ~ s(x, by = f, bs = "os"), family = poisson, data = d)), "f", "x"))
cmp <- function(a, b, tol = 1e-6) isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
out <- data.frame()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  m1 <- tryCatch(q(eval(cs[[1]])), error = function(e) conditionMessage(e))
  m2 <- tryCatch(q(eval(cs[[2]])), error = function(e) conditionMessage(e))
  if (is.character(m1) || is.character(m2)) {
    out <- rbind(out, data.frame(case = nm, logLik_1x = NA, logLik_2 = NA, fitted = NA, curves = NA,
                                 note = paste(c(if (is.character(m1)) paste("1.x:", m1), if (is.character(m2)) paste("2.0:", m2)), collapse = " ")))
    next
  }
  l1 <- as.numeric(logLik(m1$semer)); l2 <- as.numeric(logLik(m2$semer))
  fv <- max(abs(fitted(m1$semer) - fitted(m2$semer)))
  cv <- if (!is.null(cs[[4]])) {
    sp <- function(m) if ("at" %in% names(formals(semipred))) semipred(m, cs[[3]], cs[[4]], at = xs) else
      semipred(m, cs[[3]], cs[[4]], covariateV = xs, prt = FALSE)
    p1 <- q(sp(m1))$pred_df
    p2 <- q(sp(m2))$pred_df
    resp <- setdiff(names(p1), c(cs[[3]], cs[[4]], "SE", "LL", "UL", "df", "SE_mu", "Df"))[1]
    max(abs(p1[[resp]] - p2[[resp]]), abs(p1$SE - p2$SE))
  } else NA
  out <- rbind(out, data.frame(case = nm, logLik_1x = l1, logLik_2 = l2, fitted = fv, curves = cv, note = ""))
}
print(out, digits = 10, row.names = FALSE)

# ---- The models of vignettes/articles/overview.Rmd: the 1.1.1 article's code
# (git 075eb43) against the 2.0 article's code ----
if (requireNamespace("HRW", quietly = TRUE)) {
  data(WarsawApts, package = "HRW"); data(indonRespir, package = "HRW"); data(ragweed, package = "HRW")
  data(copper, package = "predictmeans")
  WarsawApts$district <- factor(WarsawApts$district)
  names(WarsawApts)[6] <- "apartment.price"
  indonRespir[, c(1, 5)] <- lapply(indonRespir[, c(1, 5)], factor)
  indonRespir$sex <- indonRespir$female
  levels(indonRespir$sex) <- c("male", "female")
  ragweed$year <- factor(ragweed$year)
  nbinom1 <- glmmTMB::nbinom1
  # Caschool (in predictmeans since 2.0; for an older library give the file in CASCHOOL_RDATA)
  if (nzchar(Sys.getenv("CASCHOOL_RDATA"))) load(Sys.getenv("CASCHOOL_RDATA")) else data(Caschool, package = "predictmeans")
  Caschool$log.avginc <- log(Caschool$avginc)
  art <- list(
    warsaw1 = list(quote(semireg(apartment.price ~ construction.date, smoothZ = list(grp = smZ(construction.date, k = 25)),
                                 data = WarsawApts)),
                   quote(semireg(apartment.price ~ s(construction.date, bs = "os", k = 25), data = WarsawApts))),
    warsaw2 = list(quote(semireg(apartment.price ~ construction.date*district,
                                 smoothZ = list(grp = smZ(construction.date, k = 10, by = district, group = TRUE)),
                                 data = WarsawApts)),
                   quote(semireg(apartment.price ~ s(construction.date, by = district, bs = "os", k = 10),
                                 data = WarsawApts))),
    caschool = list(quote(semireg_tmb(mathscr ~ calwpct + log.avginc + compstu + (1|county) + (1|grspan),
                                      smoothZ = list(cal_grp = smZ(calwpct, k = 12), avg_grp = smZ(log.avginc, k = 5),
                                                     comp_grp = smZ(compstu, k = 12)), data = Caschool)),
                    quote(semireg_tmb(mathscr ~ s(calwpct, bs = "os", k = 12) + s(log.avginc, bs = "os", k = 5) +
                                        s(compstu, bs = "os", k = 12) + (1 | county) + (1 | grspan),
                                      data = Caschool))),
    indon1 = list(quote(semireg_tmb(respirInfec ~ age + vitAdefic + sex + height + stunted + visit2 + visit3 + visit4 +
                                      visit5 + visit6 + (1|idnum), smoothZ = list(grp = smZ(age, k = 25)),
                                    family = binomial, data = indonRespir)),
                  quote(semireg_tmb(respirInfec ~ s(age, bs = "os", k = 25) + vitAdefic + sex + height + stunted +
                                      visit2 + visit3 + visit4 + visit5 + visit6 + (1 | idnum),
                                    family = binomial, data = indonRespir))),
    indon2 = list(quote(semireg_tmb(respirInfec ~ age + vitAdefic + sex + height + stunted + visit2 + visit3 + visit4 +
                                      visit5 + visit6 + (1|idnum),
                                    smoothZ = list(grp = smZ(age, k = 25, by = sex, group = TRUE)),
                                    family = binomial, data = indonRespir)),
                  quote(semireg_tmb(respirInfec ~ s(age, by = sex, bs = "os", k = 25) + vitAdefic + height + stunted +
                                      visit2 + visit3 + visit4 + visit5 + visit6 + (1 | idnum),
                                    family = binomial, data = indonRespir))),
    pollen = list(quote(semireg_tmb(pollenCount ~ rain + year + temperatureResidual + dayInSeason + windSpeed,
                                    smoothZ = list(s_day = smZ(dayInSeason, by = year, group = TRUE),
                                                   s_temp = smZ(temperatureResidual), s_wind = smZ(windSpeed)),
                                    family = nbinom1, data = ragweed)),
                  quote(semireg_tmb(pollenCount ~ rain + s(dayInSeason, by = year, bs = "os") +
                                      s(temperatureResidual, bs = "os") + s(windSpeed, bs = "os"),
                                    family = nbinom1, data = ragweed))),
    copper = list(quote(semireg(log(grade) ~ (xcoord + ycoord + zcoord)^2,
                                smoothZ = list(grp_z = smZ(zcoord), grp_geo = smZ(cbind(xcoord, ycoord), type = "Ztps")),
                                data = copper)),
                  quote(semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6) +
                                  xcoord:ycoord + xcoord:zcoord + ycoord:zcoord, data = copper))))
  aout <- data.frame()
  for (nm in names(art)) {
    m1 <- tryCatch(q(eval(art[[nm]][[1]])), error = function(e) conditionMessage(e))
    m2 <- tryCatch(q(eval(art[[nm]][[2]])), error = function(e) conditionMessage(e))
    if (is.character(m1) || is.character(m2)) {
      aout <- rbind(aout, data.frame(model = nm, logLik_1x = NA, logLik_2 = NA, fitted = NA,
                                     note = paste(m1[is.character(m1)], m2[is.character(m2)])))
      next
    }
    aout <- rbind(aout, data.frame(model = nm, logLik_1x = as.numeric(logLik(m1$semer)),
                                   logLik_2 = as.numeric(logLik(m2$semer)),
                                   fitted = max(abs(fitted(m1$semer) - fitted(m2$semer))), note = ""))
  }
  cat("\nThe overview article's models\n")
  print(aout, digits = 10, row.names = FALSE)
}
