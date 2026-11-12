#' Permutation Test of random or fixed effects for \code{lmer}, \code{glmer} and \code{glmmTMB} models.
#'
#' @description
#' This function provides permutation tests for the terms in a linear mixed model of \code{lmer}, and in a generalized
#' linear mixed model of \code{glmer} or \code{glmmTMB} (also as \code{permglmer()} and \code{permglmmTMB()}).
#' @usage
#' permlmer(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
#'          method=c("auto", "working", "quantile"))
#' permglmer(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
#'           method=c("auto", "working", "quantile"))
#' permglmmTMB(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
#'             method=c("auto", "working", "quantile"))
#' @param lmer0 \code{lmer} (or \code{glmer}, \code{glmmTMB}) model under H0, note that \code{lmer0} model must nest within
#'   \code{lmer1} model.
#' @param lmer1 \code{lmer} (or \code{glmer}, \code{glmmTMB}) model under H1, note that \code{lmer0} model must nest within
#'   \code{lmer1} model.
#' @param nperm Number of permutation, the default value is 999.
#' @param ncore Number of cores for parallel computing, the default value is 1 (before predictmeans 2.0: 3).
#' @param plot Plot permutation distribution or not, the default value is FALSE.
#' @param seed Specify a random number generator seed, for reproducible results. The user's random number stream is
#'   left as it was.
#' @param perms Optional matrix of permutations to use instead of free ones, one row per observation and one column
#'   per permutation (e.g. from \code{\link{permindex}}); \code{nperm} is then its number of columns.
#' @param method For GLMMs only: \code{"working"} permutes on the working (link) scale, \code{"quantile"} permutes
#'   quantile residuals of a 0/1 response and refits the GLMMs; \code{"auto"} (the default) uses \code{"quantile"} for
#'   0/1 responses and \code{"working"} otherwise. See Details.
#' @return
#' The \code{anova()} table of the two models with the permutation p-value in column \code{Perm-p}; for GLMMs,
#' attribute \code{"method"} says which method was used, and attribute \code{"failures"} gives the number of refits
#' that failed and were left out (if any).
#' @details
#' The residuals of \code{lmer0} are decorrelated with the symmetric square root of their covariance under
#' \code{lmer0}, permuted, correlated again and added to the fitted values of \code{lmer0} (Lee and Braun, 2012); both models
#' are refitted to each permuted data set. The likelihood-ratio statistics are compared after rounding to 4 decimals, so an
#' observed statistic of 0 (the tested variance estimated at 0) gives a p-value of 1. Since predictmeans 2.0 the square root
#' is symmetric, so the results do not depend on the order of the data (1.x used the Cholesky factor).
#'
#' For GLMMs (\code{glmer} or \code{glmmTMB}; since predictmeans 2.0), \code{method = "working"} does the same on the
#' working (link) scale: \code{lmer0} is linearised at its fit (working response \eqn{z = \eta + (y - \mu) d\eta/d\mu},
#' working variance \eqn{\phi V(\mu)/(a (d\mu/d\eta)^2)}; the covariance of z also has the random effects, as in
#' \code{\link{R2_glmm}}), the decorrelated working residuals are permuted, and the statistic is the likelihood-ratio
#' statistic of weighted linear mixed models fitted to the working responses (REML when the two models have the same
#' fixed effects, ML otherwise). In simulations it kept the 5\% level for counts and for binomial data with several
#' trials, and was much more robust than the usual likelihood-ratio test when counts were overdispersed. For 0/1
#' responses the working scale fails when events are rare or clusters tiny (36\% rejections at the 5\% level), so
#' \code{method = "quantile"} is used: randomised quantile residuals (Dunn and Smyth, 1996) of the response given the
#' conditional means of \code{lmer0} are permuted and mapped back to 0/1 responses, and both models are refitted
#' (\code{glmer} fits with \code{nAGQ = 0}, for the observed data too, so the statistic can differ from the
#' \code{Chisq} column); it kept the level (5.6-6.2\%). A random effect with one level per observation (often used
#' for overdispersion) cannot be used on the working scale; \code{\link{bootmodels}} can test it. \code{glmmTMB} models must have no zero inflation and the
#' default dispersion model; families poisson, binomial, nbinom1, nbinom2, Gamma, gaussian and beta are handled.
#' @references
#' Oliver E. Lee and Thomas M. Braun (2012),
#'     \emph{Permutation Tests for Random Effects in Linear Mixed Models. Biometrics}, Journal 68(2).
#'
#' Peter K. Dunn and Gordon K. Smyth (1996),
#'     \emph{Randomized Quantile Residuals}. Journal of Computational and Graphical Statistics, 5(3), 236-244.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' ## Test random effects (99 permutations to keep the example quick;
#' ## use 999 or more in practice)
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
#' fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
#' anova(fm1, fm2, fm3)
#' \donttest{
#' permlmer(fm3, fm2, nperm=99, ncore=1)
#' permlmer(fm2, fm1, nperm=99, ncore=1)
#' }
#'
#' ## Test fixed effects
#' Oats$nitro <- factor(Oats$nitro)
#' fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' \donttest{
#' permlmer(fm0, fm, nperm=99, ncore=1)
#'
#' ## A GLMM: do tick counts vary between locations, beyond broods?
#' gm1 <- glmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD), family = poisson, data = lme4::grouseticks)
#' gm2 <- glmer(TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | LOCATION), family = poisson,
#'              data = lme4::grouseticks)
#' permglmer(gm1, gm2, nperm = 99, ncore = 1, seed = 1)
#' }
#' @aliases permglmer permglmmTMB
#' @name permlmer
NULL

permlmer <- function(lmer0, lmer1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
                     method=c("auto", "working", "quantile")){

  c_deparse <- function (...)
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  if (inherits(lmer0, c("glmerMod", "glmmTMB")) || inherits(lmer1, c("glmerMod", "glmmTMB"))) {
    return(pm_permglmm(lmer0, lmer1, nperm, ncore, if (missing(seed)) NULL else seed, perms,
                       match.arg(method), plot, lmernames[1:2]))
  }
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) {
    stop("The models must be lmer, glmer or glmmTMB fits!")
  }
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) {
    stop("Please check the response in your model!")
  }

  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1)
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1)

  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }

  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))

  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) {
    stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  }

  if(setequal(fixef0name, fixef1name)) {
    ref <- TRUE
  } else {
    ref <- FALSE
  }

  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  Lambda1n <- getME(lmer1, "Lambda")
  Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # V = Z Lambda Lambda' Z' + I (H0 variances in H1's structure); W is its
  # order-invariant square root (pm_whiten(), R/whiten.R).
  W <- pm_whiten(getME(lmer1, "Z") %*% Lambda1n)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta

  wterrors <- W$solve(errors)

  # Random effects (same fixed effects): the decorrelated residuals are
  # permuted. Fixed effects: fitted values and decorrelated residuals are
  # permuted together.
  if (!is.null(perms)) nperm <- ncol(as.matrix(perms))
  idx <- pm_perm_source(length(errors), if (missing(seed)) NULL else seed, perms)(nperm)
  permy <- function(k) {
    r <- idx[, k]
    if (ref) xbeta + W$mult(wterrors[r]) else xbeta[r] + W$mult(wterrors[r])
  }

  # Likelihood-ratio statistics. They are compared after rounding to 4
  # decimals (pm_perm_p()): when the tested variance is estimated at 0 the
  # LRT is 0 up to optimizer noise, and without rounding permuted values of 0
  # could count as smaller than an observed 1e-9.
  lrtest1 <- max(as.numeric(2*(logLik(lmer1, REML=ref)-logLik(lmer0, REML=ref))), 0)
  lrtest2 <- pm_perm_map(seq_len(nperm), function(k) {
    y <- permy(k)
    as.numeric(2*(logLik(suppressMessages(refit(lmer1, y)), REML=ref) -
                  logLik(suppressMessages(refit(lmer0, y)), REML=ref)))
  }, ncore, packages = "lme4")
  lrtest2 <- pm_perm_drop(lrtest2)
  lrtest <- pmax(unlist(lrtest2$res), 0)
  ok <- !is.na(lrtest)
  failures <- lrtest2$failures + sum(!ok)
  lrtest <- round(lrtest[ok], 4)
  lrtest1 <- round(lrtest1, 4)
  perm_p <- pm_perm_p(lrtest1, lrtest, digits = 4)
  aod <- anova(lmer0, lmer1, refit=!ref)
  aod$'Perm-p' <- c(NA, perm_p)
  if (failures > 0) attr(aod, "failures") <- failures

  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

permglmer <- permlmer

permglmmTMB <- permlmer
