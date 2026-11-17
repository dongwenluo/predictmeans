## The smoothing-parameter correction of semipred(unconditional = TRUE)
## against mgcv's (predict.gam(unconditional = TRUE), Vc; Wood, Pya and
## Safken 2016), the same model and O'Sullivan basis fitted by semireg() and
## gam(method = "REML"); gaussian and Poisson curves.
## Result (2026-11-17):
##   gaussian: conditional SE equal to 0.1%; with the correction, SE equal to
##     0.1%, the added variance 0.98-1.00 of mgcv's where it adds > 1% to the SE.
##   Poisson (glmer Laplace vs gam REML, estimates 0.014 apart on the log
##     scale): SE with the correction equal to 0.5%; added variance 1.00-1.12
##     of mgcv's where it matters (mgcv has a second-order term besides the
##     first-order J V J' of semireg). Pointwise ratios of the added variance
##     are unstable only where it is near 0.
##   Rscript dev/check-uncond-mgcv.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(11)
n <- 120
d <- data.frame(x = c(0, 10, runif(n - 2, 0, 10)))
xs <- seq(0.5, 9.5, by = 0.5)
for (fam in c("gaussian", "poisson")) {
  if (fam == "gaussian") { d$y <- sin(d$x / 1.5) + rnorm(n, 0, 0.4); F <- gaussian() } else
    { d$y <- rpois(n, exp(1 + sin(d$x / 1.5) / 2)); F <- poisson() }
  f <- if (fam == "gaussian") q(semireg(y ~ s(x, bs = "os", k = 10), data = d, prt = FALSE)) else
    q(semireg(y ~ s(x, bs = "os", k = 10), family = poisson, data = d, prt = FALSE))
  g <- mgcv::gam(y ~ s(x, bs = "os", k = 10), family = F, method = "REML", data = d)
  pc <- q(semipred(f, "x", "x", covariateV = xs, unconditional = FALSE, prt = FALSE))$pred_df
  pu <- q(semipred(f, "x", "x", covariateV = xs, unconditional = TRUE, prt = FALSE))$pred_df
  se_c <- if (fam == "gaussian") pc$SE else pc$SE                      # link scale
  se_u <- pu$SE
  gc <- predict(g, data.frame(x = xs), se.fit = TRUE)
  gu <- predict(g, data.frame(x = xs), se.fit = TRUE, unconditional = TRUE)
  ## the first-order part of mgcv's correction alone: Vp + J V_rho J' (Vc1 = Vc without the second-order term)
  X <- predict(g, data.frame(x = xs), type = "lpmatrix")

  cat(sprintf("%-8s fits: estimate %.1e apart | conditional SE: semipred / gam %.4f-%.4f\n", fam,
              max(abs((if (fam == "gaussian") pc$y else log(pc$y)) - gc$fit)), min(se_c / gc$se.fit), max(se_c / gc$se.fit)))
  inc_s <- se_u^2 - se_c^2; inc_g <- gu$se.fit^2 - gc$se.fit^2
  cat(sprintf("         added variance: semipred / gam %.2f-%.2f (median %.2f); SE increase semipred %.1f-%.1f%%, gam %.1f-%.1f%%\n",
              min(inc_s / inc_g), max(inc_s / inc_g), median(inc_s / inc_g),
              100 * min(se_u / se_c - 1), 100 * max(se_u / se_c - 1), 100 * min(gu$se.fit / gc$se.fit - 1), 100 * max(gu$se.fit / gc$se.fit - 1)))
  cat(sprintf("         corrected SE: semipred / gam %.4f-%.4f; where the correction adds > 1%% to the SE: added variance ratio %.2f-%.2f\n",
              min(se_u / gu$se.fit), max(se_u / gu$se.fit),
              min((inc_s / inc_g)[gu$se.fit / gc$se.fit > 1.01]), max((inc_s / inc_g)[gu$se.fit / gc$se.fit > 1.01])))
  if (!is.null(g$Vc1)) {
    s1 <- sqrt(rowSums((X %*% g$Vc1) * X))
    inc_1 <- s1^2 - gc$se.fit^2
    cat(sprintf("         mgcv's first-order term alone (Vc1): semipred / gam %.2f-%.2f\n", min(inc_s / inc_1), max(inc_s / inc_1)))
  }
  cat(sprintf("         smoothing parameter: mgcv sp %.4g; semireg sigma^2/sigma_u^2 %.4g\n", g$sp,
              sigma(f$semer)^2 / as.numeric(attr(lme4::VarCorr(f$semer)[[names(f$u_lst)[1]]], "stddev"))^2))
}
