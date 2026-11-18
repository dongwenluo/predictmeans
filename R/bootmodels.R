#' Parametric Bootstrap Tests
#'
#' @description
#' Parametric bootstrap p-values for the terms and coefficients of a model. Responses are simulated from the model
#' fitted without its fixed effects (only an intercept, any offset, and the same random effects and covariance
#' structure: the global null hypothesis), the full model is refitted to each, and a p-value is the proportion of
#' simulated statistics at least as large as the observed one. Useful for GLMs and GLMMs and small samples, where
#' Wald or F tests can be poor.
#' @details
#' The result has the same form as that of \code{\link{permmodels}}: its \code{permlist} can be given to
#' \code{\link{predictmeans}} or \code{\link{contrastmeans}} (argument \code{permlist}) for bootstrap p-values of
#' pairwise comparisons or contrasts (\code{p.method} \code{"bootstrap"} in the \code{pairs} table).
#'
#' In a simulation study (500 data sets per case, 5\% level; \code{dev/sim-boot.R} in the package sources) the
#' bootstrap corrected the liberal Wald chi-square test of a whole-plot factor in a split plot with 6 blocks (10\%
#' to 6-7\%, as Kenward-Roger) and kept the level in Poisson and binomial GLMMs, also when other terms had effects.
#' Like the likelihood-ratio test it relies on the model: with overdispersed counts fitted as Poisson it rejected
#' 37\% of true null hypotheses, where \code{\link{permmodels}} (Freedman-Lane on the working scale) gave 7\%.
#' @usage
#' bootmodels(model, nsim = 999, type = c("I", "II", "III", "1", "2", "3"), seed = NULL,
#'            ncore = 1, ...)
#' @param model A model fitted by \code{lm}, \code{aov} (without \code{Error()}), \code{glm} (a family with a
#'   distribution: gaussian, binomial, Poisson, Gamma; also \code{MASS::glm.nb}), \code{gls}, \code{lme},
#'   \code{lmer}, \code{glmer} or \code{glmmTMB}. Responses of \code{lme} and \code{gls} models are simulated from
#'   the null model's covariance of the data (random effects, correlation and variance structures).
#' @param nsim Number of simulated data sets, the default is 999.
#' @param type The type of ANOVA table (as in \code{\link{permmodels}}): Type I gives sequential F statistics (for a
#'   \code{glm} with a known dispersion, the deviance; \code{glmmTMB} has none and uses Type II), Types II and III
#'   Wald chi-square statistics (\code{lm}: F; \code{glm}: likelihood ratio).
#' @param seed Optional random number seed, for reproducible results; the user's random number stream is left as it was.
#' @param ncore Number of cores for refitting, the default is 1.
#' @param ... Not used: \code{prt} (before 2.0) stops with its 2.0 spelling; \code{print()} of the result shows the tables.
#' @return
#' Invisibly, a list with \code{permlist} (the coefficients and covariance matrices of the refits),
#' \code{COEFFICENTS} and \code{ANOVA} (the model's tables with a \code{Boot_p_value} column). Attribute
#' \code{"resampling"} is \code{"parametric bootstrap"}.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   \donttest{
#'   bt <- bootmodels(fm, nsim = 99, seed = 1)   # use 999 or more in practice
#'   predictmeans(fm, "nitro", permlist = bt, plot = FALSE, ncore = 1)$pairs
#'   }
#' @name bootmodels
NULL

bootmodels <- function(model, nsim = 999, type = c("I", "II", "III", "1", "2", "3"), seed = NULL, ncore = 1, ...) {
  pm_defunct("bootmodels", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  type <- match.arg(as.character(type)[1], c("I", "II", "III", "1", "2", "3"))
  if (inherits(model, "aovlist")) {
    stop("bootmodels() does not handle aov models with Error(): fit the error strata as random effects (lmer or ",
         "lme) instead.", call. = FALSE)
  }
  if (inherits(model, "gam")) {
    stop("bootmodels() does not handle gam models (a null model for the smooths is not defined); use permmodels().",
         call. = FALSE)
  }
  if (!inherits(model, c("lm", "lme", "gls", "lmerMod", "glmerMod", "glmmTMB"))) {
    stop("bootmodels() handles lm, aov, glm, gls, lme, lmer, glmer and glmmTMB models.", call. = FALSE)
  }
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  plain <- function(m) if (inherits(m, "lmerModLmerTest")) methods::as(m, "lmerMod") else m
  model <- plain(model)
  if (inherits(model, "glmmTMB") && type %in% c("I", "1")) {
    message("glmmTMB models have no sequential (Type I) tests; Type II Wald tests are used.")
    type <- "II"
  }
  # the null model, simulation and refits of this class (R/boot-parts.R)
  bp <- pm_boot_parts(model)
  ysim <- if (is.null(seed)) bp$simulate(nsim) else pm_with_seed(seed, bp$simulate(nsim))

  stats_of <- function(m) {
    mp <- mymodelparm(m, df = 0)[1:2]
    list(mp, pm_boot_anova(m, type))
  }
  one <- function(y) {
    rf <- tryCatch(bp$refit(y), error = function(e) NULL)
    if (is.null(rf)) NULL else stats_of(rf)
  }
  res <- pm_perm_drop(pm_perm_map(ysim, one, ncore, packages = bp$packages))
  failures <- res$failures
  res <- res$res
  B <- length(res)
  if (B == 0) stop("No refit to simulated data succeeded.")

  obs <- stats_of(model)
  # coefficients: |t|
  t_obs <- abs(obs[[1]]$coef/sqrt(diag(obs[[1]]$vcov)))
  t_sim <- vapply(res, function(r) abs(r[[1]]$coef/sqrt(diag(r[[1]]$vcov))), numeric(length(t_obs)))
  ctab <- pm_coef_table(model)
  bp <- pm_perm_p(t_obs, t_sim)
  if ("(Intercept)" %in% rownames(ctab)) bp[rownames(ctab) == "(Intercept)"] <- NA
  COEF <- round(cbind(ctab, Boot_p_value = bp), 4)
  # terms
  aT <- obs[[2]]
  scol <- pm_boot_stat_col(aT)
  s_sim <- vapply(res, function(r) r[[2]][, scol], numeric(nrow(aT)))
  ap <- pm_perm_p(aT[, scol], s_sim)
  ap[rownames(aT) %in% c("(Intercept)", "Residuals")] <- NA     # the null model keeps the intercept
  ANOVA <- round(cbind(aT, Boot_p_value = ap), 4)
  attr(ANOVA, "heading") <- attr(aT, "heading")
  out <- pm_permmodels_result(list(permlist = lapply(res, `[[`, 1), COEFFICENTS = COEF, ANOVA = ANOVA), NULL,
                              paste0("Note: Boot_p_value is from ", sQuote(B),
                                     " parametric bootstrap samples under the global null model",
                                     pm_perm_failnote(failures), "."))
  attr(out, "resampling") <- "parametric bootstrap"
  if (failures > 0) attr(out, "failures") <- failures
  out
}
