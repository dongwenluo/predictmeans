## Coverage of the band of the difference between two smooth curves,
## semipred(contr = ): pointwise (each covariate value), simultaneous (the
## whole curve) and, for generalized models, the response-scale difference.
## Run: Rscript dev/sim-curve-difference.R   (about 4 minutes)
##
## Result of 2026-11-17 (n = 400, 150 data sets each, nominal 0.95), without
## and with the correction for the uncertainty of the smoothing parameters
## (semipred(unconditional = TRUE), R/semireg-engine.R):
##   family    pointwise     (range over x)   whole curve   response-scale
##   gaussian  0.968 / 0.976  0.92-0.99        0.973 / 0.980  -
##   poisson   0.933 / 0.957  0.88-0.99 / 0.91-1   0.907 / 0.940  0.938
##   binomial  0.936 / 0.955  0.89-0.97 / 0.92-0.98   0.893 / 0.947  0.926
## Generalized models were below the nominal level because the covariance was
## conditional on the estimated smoothing parameters; with the correction they
## are at it. Gaussian bands were already at or above it. Hence the default
## of semipred() and predictmeans() on semireg fits: the correction for
## generalized models, not for gaussian ones. (Result of 2026-11-16, before
## the semireg engine fixes: gaussian 0.958/0.959, poisson 0.922/0.865,
## binomial 0.936/0.876, conditional only.)
args <- commandArgs(TRUE)       # optional: the package source, for pkgload::load_all()
if (length(args)) suppressMessages(pkgload::load_all(args[1], quiet = TRUE)) else suppressMessages(library(predictmeans))
grDevices::pdf(NULL)
n <- 400; xs <- seq(0.2, 2.8, 0.2); nsim <- 150
fa <- function(x) 0.4*sin(2*x); fb <- function(x) 0.4*sin(2*x) + 0.5*x/3
q <- function(e) suppressWarnings(suppressMessages(e))
run <- function(fam) {
  hit <- hits <- matrix(NA, nsim, length(xs)); hr <- hitu <- hitsu <- hit
  for (i in seq_len(nsim)) {
    d <- data.frame(x = runif(n, 0, 3), g = gl(2, n/2, labels = c("a", "b")))
    eta <- ifelse(d$g == "a", fa(d$x), fb(d$x))
    d$y <- switch(fam, gaussian = eta + rnorm(n, 0, .5), poisson = rpois(n, exp(0.5 + eta)), binomial = rbinom(n, 1, plogis(-0.5 + 2*eta)))
    tr <- switch(fam, gaussian = fb(xs) - fa(xs), poisson = fb(xs) - fa(xs), binomial = 2*(fb(xs) - fa(xs)))
    trr <- switch(fam, gaussian = tr, poisson = exp(.5 + fb(xs)) - exp(.5 + fa(xs)), binomial = plogis(-.5 + 2*fb(xs)) - plogis(-.5 + 2*fa(xs)))
    sm <- tryCatch(q(if (fam == "gaussian") semireg(y ~ x*g, smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d, prt = FALSE) else
      if (fam == "poisson") semireg(y ~ x*g, smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d, family = poisson, prt = FALSE) else
      semireg(y ~ x*g, smoothZ = list(grp = smZ(x, k = 8, by = g, group = TRUE)), data = d, family = binomial, prt = FALSE)), error = function(e) NULL)
    if (is.null(sm)) next
    p <- tryCatch(q(semipred(sm, "g", "x", contr = c("b", "a"), covariateV = xs, prt = FALSE))$pred_df, error = function(e) NULL)
    ps <- tryCatch(q(semipred(sm, "g", "x", contr = c("b", "a"), covariateV = xs, prt = FALSE, simultaneous = TRUE))$pred_df, error = function(e) NULL)
    # with the uncertainty of the smoothing parameters (unconditional = TRUE)
    pu <- tryCatch(q(semipred(sm, "g", "x", contr = c("b", "a"), covariateV = xs, prt = FALSE, unconditional = TRUE))$pred_df, error = function(e) NULL)
    psu <- tryCatch(q(semipred(sm, "g", "x", contr = c("b", "a"), covariateV = xs, prt = FALSE, simultaneous = TRUE, unconditional = TRUE))$pred_df, error = function(e) NULL)
    if (is.null(p) || is.null(ps) || is.null(pu) || is.null(psu)) next
    hitu[i, ] <- pu$LL <= tr & tr <= pu$UL; hitsu[i, ] <- psu$LL <= tr & tr <= psu$UL
    if (fam == "gaussian") { hit[i, ] <- p$LL <= tr & tr <= p$UL; hits[i, ] <- ps$LL <= tr & tr <= ps$UL }
    else { hit[i, ] <- p$LL <= tr & tr <= p$UL; hits[i, ] <- ps$LL <= tr & tr <= ps$UL; hr[i, ] <- p$response.LL <= trr & trr <= p$response.UL }
  }
  ok <- !is.na(hit[, 1])
  cat(fam, ": fits", sum(ok), "| pointwise coverage (mean over x)", round(mean(hit[ok, ]), 3), "| range over x", round(range(colMeans(hit[ok, ])), 3),
      "| simultaneous: whole curve covered", round(mean(apply(hits[ok, , drop = FALSE], 1, all)), 3),
      if (fam != "gaussian") paste("| response-scale pointwise", round(mean(hr[ok, ]), 3)),
      "\n   unconditional: pointwise", round(mean(hitu[ok, ]), 3), "| range over x", round(range(colMeans(hitu[ok, ])), 3),
      "| simultaneous", round(mean(apply(hitsu[ok, , drop = FALSE], 1, all)), 3), "\n")
}
set.seed(11); for (f in c("gaussian", "poisson", "binomial")) run(f)
