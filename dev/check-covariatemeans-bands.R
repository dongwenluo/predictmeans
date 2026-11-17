## Bands of covariatemeans() against emmeans and by hand: curves (means
## along a covariate per level of a factor) and the difference of two curves
## (contr), pointwise and simultaneous (Scheffe), for lmer (Kenward-Roger),
## glmer (Poisson) and glm (binomial).
## Result (2026-11-17):
##   lmer: curves equal to emmeans (KR): means exact, SE 3e-6; contrast
##     difference exact, SE 1e-6; df 1% apart (observed-information KR, as SAS,
##     against emmeans' expected information: the package's decision);
##     simultaneous band = sqrt(2 F(0.95; 2, df)) row by row (Scheffe, rank 2).
##   glmer Poisson: link SE and means exact; ratio = exp(log difference),
##     difference of means and its SE equal to emmeans' regrid(); limits t on
##     the inner-outer df (250 here) where emmeans uses z.
##   glm binomial: link SE, probability and limits plogis(eta -+ z SE) exact.
##   Rscript dev/check-covariatemeans-bands.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(lme4))
q <- function(e) suppressWarnings(suppressMessages(e))
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))
line <- function(...) cat(sprintf(...), "\n")

## ---- lmer, Kenward-Roger --------------------------------------------------
cw <- as.data.frame(ChickWeight)
cw$Chick <- factor(cw$Chick, ordered = FALSE)
m1 <- lmer(weight ~ Diet * Time + (1 | Chick), data = cw)
xs <- c(4, 10, 16)
cm <- q(covariatemeans(m1, "Diet", "Time", covariateV = xs))$pltdf
em <- summary(emmeans::emmeans(m1, ~ Diet | Time, at = list(Time = xs), lmer.df = "kenward-roger"))
o <- match(paste(em$Diet, em$Time), paste(cm$Diet, cm$xvar))
crit <- (cm$UL - cm$Mean) / cm$ses
line("lmer curves:  mean %.1e, SE %.1e, limits %.1e (emmeans KR); t quantile of %s df",
     rel(cm$Mean[o], em$emmean), rel(cm$ses[o], em$SE), rel(cm$UL[o], em$upper.CL),
     paste(range(round(em$df)), collapse = "-"))
cc <- q(covariatemeans(m1, "Diet", "Time", covariateV = xs, contr = c("2", "1")))$pltdf
emc <- summary(emmeans::contrast(emmeans::emmeans(m1, ~ Diet | Time, at = list(Time = xs), lmer.df = "kenward-roger"),
                                 list(d21 = c(-1, 1, 0, 0))), infer = TRUE)
line("lmer contrast 2 - 1: difference %.1e, SE %.1e, df %.1e, limits %.1e (emmeans KR)",
     rel(cc$diff, emc$estimate), rel(cc$ses, emc$SE), rel(cc$Df, emc$df), rel(cc$UL, emc$upper.CL))
## simultaneous: Scheffe over the covariate, the rank of the contrast rows
## (a straight line: rank 2)
xg <- seq(0, 21, by = 1)
cs <- q(covariatemeans(m1, "Diet", "Time", covariateV = xg, contr = c("2", "1"), band = "simultaneous"))$pltdf
crit_s <- (cs$UL - cs$diff) / cs$ses
line("lmer simultaneous: critical value %s; sqrt(2 F(0.95; 2, df)) with the df of the rows: %s",
     paste(format(range(crit_s), digits = 5), collapse = "-"),
     paste(format(range(sqrt(2 * qf(0.95, 2, cs$Df))), digits = 5), collapse = "-"))

## ---- glmer, Poisson -----------------------------------------------------------
set.seed(3)
d <- expand.grid(x = seq(0, 10, length = 11), A = factor(c("a1", "a2", "a3")), site = factor(1:8))
d$y <- rpois(nrow(d), exp(0.5 + 0.1 * d$x + 0.2 * (d$A == "a2") + rnorm(8, 0, 0.3)[d$site]))
m2 <- glmer(y ~ A * x + (1 | site), family = poisson, data = d)
xs2 <- c(2, 5, 8)
cg <- q(covariatemeans(m2, "A", "x", covariateV = xs2))$pltdf
eg <- summary(emmeans::emmeans(m2, ~ A | x, at = list(x = xs2)))                     # link scale
eg_r <- summary(emmeans::emmeans(m2, ~ A | x, at = list(x = xs2), type = "response"))
o2 <- match(paste(eg$A, eg$x), paste(cg$A, cg$xvar))
crit2 <- unique(round((log(cg$UL) - log(cg$Mean)) / cg$ses, 4))
line("glmer curves: link SE %.1e (emmeans); mean exp(eta) %.1e; limits exp(eta -+ crit SE), crit %s = qt(0.975, %s) (emmeans z 1.960)",
     rel(cg$ses[o2], eg$SE), rel(cg$Mean[o2], eg_r$rate), paste(crit2, collapse = ","),
     paste(round(uniroot(function(v) qt(0.975, v) - crit2[1], c(1, 1e6))$root), collapse = ","))
cgc <- q(covariatemeans(m2, "A", "x", covariateV = xs2, contr = c("a2", "a1")))$pltdf
egc <- summary(emmeans::contrast(emmeans::emmeans(m2, ~ A | x, at = list(x = xs2)), list(r = c(-1, 1, 0))))
egr <- summary(emmeans::contrast(emmeans::regrid(emmeans::emmeans(m2, ~ A | x, at = list(x = xs2))), list(r = c(-1, 1, 0))))
line("glmer contrast a2 / a1: log ratio %.1e, SE %.1e; ratio = exp %s; difference of means %.1e, its SE %.1e (emmeans regrid)",
     rel(cgc$diff, egc$estimate), rel(cgc$ses, egc$SE), isTRUE(all.equal(cgc$ratio, exp(cgc$diff))),
     rel(cgc$response.diff, egr$estimate), rel(cgc$response.ses, egr$SE))

## ---- glm, binomial -------------------------------------------------------------
data(Drug)
g <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
xs3 <- c(0, 0.5, 1)
cb <- q(covariatemeans(g, "Drug", "log10Dose", covariateV = xs3))$pltdf
eb <- summary(emmeans::emmeans(g, ~ Drug | log10Dose, at = list(log10Dose = xs3)))
eb_r <- summary(emmeans::emmeans(g, ~ Drug | log10Dose, at = list(log10Dose = xs3), type = "response"))
o3 <- match(paste(eb$Drug, eb$log10Dose), paste(cb$Drug, cb$xvar))
line("glm curves:   link SE %.1e, probability %.1e, limits plogis(eta -+ z SE) %.1e (emmeans: z, back-transformed)",
     rel(cb$ses[o3], eb$SE), rel(cb$Mean[o3], eb_r$prob), rel(cb$UL[o3], eb_r$asymp.UCL))
