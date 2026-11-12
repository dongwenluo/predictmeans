#' Parametric Bootstrap Tests for Mixed Models
#'
#' @description
#' Parametric bootstrap p-values for the terms and coefficients of an \code{lmer} or \code{glmer} model. Responses are
#' simulated from the model fitted without its fixed effects (only an intercept and the same random effects: the
#' global null hypothesis), the full model is refitted to each, and a p-value is the proportion of simulated statistics
#' at least as large as the observed one. Useful for GLMMs and small samples, where Wald or F tests can be poor.
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
#'            ncore = 1, prt = TRUE)
#' @param model A model fitted by \code{lmer} or \code{glmer}.
#' @param nsim Number of simulated data sets, the default is 999.
#' @param type The type of ANOVA table (as in \code{\link{permmodels}}): Type I gives sequential F statistics, Types II
#'   and III Wald chi-square statistics.
#' @param seed Optional random number seed, for reproducible results; the user's random number stream is left as it was.
#' @param ncore Number of cores for refitting, the default is 1.
#' @param prt Logical: print the tables? The default is TRUE.
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

bootmodels <- function(model, nsim = 999, type = c("I", "II", "III", "1", "2", "3"), seed = NULL, ncore = 1,
                       prt = TRUE) {
  type <- match.arg(as.character(type)[1], c("I", "II", "III", "1", "2", "3"))
  if (!inherits(model, c("lmerMod", "glmerMod"))) stop("bootmodels() needs a model fitted by lmer or glmer.")
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  plain <- function(m) if (inherits(m, "lmerModLmerTest")) methods::as(m, "lmerMod") else m
  model <- plain(model)
  ff <- stats::formula(model)
  bars <- pm_findbars(ff)
  if (length(bars) == 0) stop("The model has no random effects.")
  null_f <- stats::as.formula(paste(deparse(ff[[2]]), "~ 1 +",
                                    paste0("(", vapply(bars, function(b) paste(deparse(b), collapse = ""), ""), ")",
                                           collapse = " + ")),
                              env = environment(ff))
  null_fit <- suppressMessages(suppressWarnings(pm_refit(model, formula = null_f)))
  ysim <- if (is.null(seed)) stats::simulate(null_fit, nsim = nsim) else
    pm_with_seed(seed, stats::simulate(null_fit, nsim = nsim))

  stats_of <- function(m) {
    mp <- mymodelparm(m, df = 0)[1:2]
    aT <- if (type %in% c("I", "1")) stats::anova(m) else pm_anova(m, type = type, test.statistic = "Chisq")
    list(mp, aT)
  }
  one <- function(y) {
    rf <- tryCatch(suppressMessages(suppressWarnings(lme4::refit(model, newresp = y))), error = function(e) NULL)
    if (is.null(rf)) NULL else stats_of(rf)
  }
  res <- pm_perm_drop(pm_perm_map(ysim, one, ncore, packages = "lme4"))
  failures <- res$failures
  res <- res$res
  B <- length(res)
  if (B == 0) stop("No refit to simulated data succeeded.")

  obs <- stats_of(model)
  # coefficients: |t|
  t_obs <- abs(obs[[1]]$coef/sqrt(diag(obs[[1]]$vcov)))
  t_sim <- vapply(res, function(r) abs(r[[1]]$coef/sqrt(diag(r[[1]]$vcov))), numeric(length(t_obs)))
  ctab <- stats::coef(summary(model))
  bp <- pm_perm_p(t_obs, t_sim)
  if ("(Intercept)" %in% rownames(ctab)) bp[rownames(ctab) == "(Intercept)"] <- NA
  COEF <- round(cbind(ctab, Boot_p_value = bp), 4)
  # terms
  aT <- obs[[2]]
  scol <- if (type %in% c("I", "1")) "F value" else "Chisq"
  s_sim <- vapply(res, function(r) r[[2]][, scol], numeric(nrow(aT)))
  ap <- pm_perm_p(aT[, scol], s_sim)
  ap[rownames(aT) == "(Intercept)"] <- NA            # the null model keeps the intercept
  ANOVA <- round(cbind(aT, Boot_p_value = ap), 4)
  attr(ANOVA, "heading") <- attr(aT, "heading")
  if (prt) {
    cat("\nCoefficients of (fixed) effects:\n")
    print(COEF)
    cat("\nANOVA:\n")
    print(ANOVA)
    cat("\nNote: Boot_p_value is from ", sQuote(B), " parametric bootstrap samples under the global null model",
        pm_perm_failnote(failures), ".\n\n", sep = "")
  }
  out <- list(permlist = lapply(res, `[[`, 1), COEFFICENTS = COEF, ANOVA = ANOVA)
  attr(out, "resampling") <- "parametric bootstrap"
  if (failures > 0) attr(out, "failures") <- failures
  invisible(out)
}
