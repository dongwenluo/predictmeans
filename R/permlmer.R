#' Permutation Test of Random or Fixed Effects of Mixed Models
#'
#' @description
#' Permutation likelihood ratio test of two nested \code{lmer}, \code{glmer} or \code{glmmTMB} models. The p-value
#' is (b + 1)/(B + 1), b of the B permuted statistics being at least as large as the observed one; results are the
#' same for any \code{ncore}.
#' @usage
#' permlmer(model0, model1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
#'          method=c("auto", "working", "quantile"), ...)
#' @param model0 \code{lmer} (or \code{glmer}, \code{glmmTMB}) model under H0 (\code{lmer0} before 2.0); it must
#'   be nested within \code{model1}.
#' @param model1 \code{lmer} (or \code{glmer}, \code{glmmTMB}) model under H1 (\code{lmer1} before 2.0).
#' @param ... Not used: \code{lmer0}, \code{lmer1} (before 2.0) stop with their 2.0 spelling. \code{permglmer()} and
#'   \code{permglmmTMB()}, the same function under other names before 2.0, were removed: use \code{permlmer()}.
#' @param nperm Number of permutations; default 999.
#' @param ncore Number of cores (default 1; before predictmeans 2.0: 3).
#' @param plot \code{TRUE}: plot the permutation distribution of the statistic. Default \code{FALSE}.
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
#' The residuals of \code{model0} are decorrelated with the symmetric square root of their covariance under
#' \code{model0}, permuted, correlated again and added to the fitted values of \code{model0} (Lee and Braun, 2012); both models
#' are refitted to each permuted data set. The likelihood-ratio statistics are compared after rounding to 4 decimals, so an
#' observed statistic of 0 (the tested variance estimated at 0) gives a p-value of 1. Since predictmeans 2.0 the square root
#' is symmetric, so the results do not depend on the order of the data (1.x used the Cholesky factor).
#' Models with the same fixed effects (tests of random effects) are compared by the REML likelihood ratio; models
#' with different fixed effects are first refitted by ML, so the permuted statistic is the ML likelihood ratio
#' that \code{anova()} prints in column \code{Chisq} (before 2.0 it was the ML likelihood at the REML estimates).
#' Check: the identity permutation gives \code{Perm-p} = 1 (\file{tests/testthat/test-permlmer.R}).
#'
#' For GLMMs (\code{glmer} or \code{glmmTMB}; since predictmeans 2.0), \code{method = "working"} does the same on the
#' working (link) scale: \code{model0} is linearised at its fit (working response \eqn{z = \eta + (y - \mu) d\eta/d\mu},
#' working variance \eqn{\phi V(\mu)/(a (d\mu/d\eta)^2)}; the covariance of z also has the random effects, as in
#' \code{\link{R2_model}}), the decorrelated working residuals are permuted, and the statistic is the likelihood-ratio
#' statistic of weighted linear mixed models fitted to the working responses (REML when the two models have the same
#' fixed effects, ML otherwise). In simulations it kept the 5\% level for counts and for binomial data with several
#' trials, and was much more robust than the usual likelihood-ratio test when counts were overdispersed. For 0/1
#' responses the working scale fails when events are rare or clusters tiny (36\% rejections at the 5\% level), so
#' \code{method = "quantile"} is used: randomised quantile residuals (Dunn and Smyth, 1996) of the response given the
#' conditional means of \code{model0} are permuted and mapped back to 0/1 responses, and both models are refitted
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
#' @seealso \code{\link{permmodels}}, \code{\link{varcomp}} (\code{test = TRUE}), \code{\link{permindex}}
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
#' permlmer(gm1, gm2, nperm = 99, ncore = 1, seed = 1)
#' }
#' @name permlmer
NULL

#' @export
permlmer <- function(model0, model1, nperm = 999, ncore=1, plot=FALSE, seed, perms=NULL,
                     method=c("auto", "working", "quantile"), ...){
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("permlmer", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  lmer0 <- model0
  lmer1 <- model1

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

  # Different fixed effects: the LRT needs ML fits (REML likelihoods of models
  # with different fixed effects are not comparable). Both models are refitted
  # by ML once, so the observed statistic, the null model that generates the
  # permuted data and every permuted refit are ML, and the observed statistic
  # equals the printed anova() chi-square. (Before 2.0.0 the REML fits were
  # kept and logLik(REML = FALSE) taken at their REML estimates: a different
  # statistic from the printed one.)
  if (!ref) {
    if (lme4::isREML(lmer0)) lmer0 <- suppressMessages(pm_refit(lmer0, REML = FALSE))
    if (lme4::isREML(lmer1)) lmer1 <- suppressMessages(pm_refit(lmer1, REML = FALSE))
    theta0 <- getME(lmer0, "theta")
    theta1 <- getME(lmer1, "theta")
    fixef0 <- fixef(lmer0)
  }

  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  Lambda1n <- getME(lmer1, "Lambda")
  Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # V = Z Lambda Lambda' Z' + diag(1/w) (H0 variances in H1's structure; w the
  # prior weights, 1 without); W is its order-invariant square root
  # (R/whiten.R). The offset is part of the fitted values: refit() keeps the
  # model's weights and offset, so only the response changes.
  pp <- pm_mer_perm_parts(lmer0, W = pm_mer_whiten(lmer1, Lambda1n), beta = fixef0)
  W <- pp$W
  xbeta <- pp$xbeta
  errors <- pp$errors
  wterrors <- pp$wterrors

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
  # lme4's convergence warnings of single permuted refits are counted and
  # reported once (attr "convergence"), not printed for each refit
  lrtest2 <- pm_perm_map(seq_len(nperm), function(k) {
    y <- permy(k)
    nconv <- 0L
    quiet_refit <- function(m) withCallingHandlers(suppressMessages(pm_mer_refit(m, y)), warning = function(w) {
      if (grepl("converge", conditionMessage(w))) { nconv <<- nconv + 1L; invokeRestart("muffleWarning") }
    })
    lr <- as.numeric(2*(logLik(quiet_refit(lmer1), REML=ref) - logLik(quiet_refit(lmer0), REML=ref)))
    structure(lr, conv = nconv > 0)
  }, ncore, packages = "lme4")
  nconv <- sum(vapply(lrtest2, function(r) isTRUE(attr(r, "conv")), logical(1)))
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
  if (nconv > 0) {
    attr(aod, "convergence") <- nconv
    message(nconv, " of ", nperm, " permuted refits had an lme4 convergence warning (usually a gradient just above ",
            "lme4's tolerance; the statistic is kept).")
  }

  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

## Removed in 2.0 (they were permlmer() itself): permlmer().
#' @export
permglmer <- function(...) pm_defunct_fun("permglmer")
#' @export
permglmmTMB <- function(...) pm_defunct_fun("permglmmTMB")
