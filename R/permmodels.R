#' Permutation Test of Linear Model
#'
#' @description
#' This function provides permutation t-tests for coefficients of (fixed) effects and permutation F-tests 
#' for the terms in a linear model such as \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @usage
#' permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
#'   test.statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, 
#'   fo=NULL, prt=TRUE, ncore=3, seed, method=c("raw", "freedman-lane"), perms=NULL)
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @param nperm The number of permutations. The default is 4999.
#' @param type type of ANOVA test, "I", "II", "III", 1, 2, or 3. Roman numerals are equivalent to the corresponding Arabic numerals.
#' @param test.statistic For type I ANOVA, F test is applied to all models, while for type II and III ANOVA, F test is performed for \code{lm}, 
#'   Chisq test for \code{lme} and \code{gls} model, Chisq or F tests for \code{lmer} model and Likelihood ratio, Wald or F tests for \code{glm} model.
#' @param exact A logical variable to indicate whether or not exact no. of permutations will be used (applicable only to free the permutation case). The default is FALSE.
#' @param data In some cases, you need to provide the data set used in model fitting, especially when you have applied some variable trnasformation in the model.
#' @param fo A model formula used in the \code{model}; \code{fo!=NULL} when the formula is specified by function \code{formula}.
#' @param prt A logical variable to indicate whether or not to print output on the screen. The default is TRUE.
#' @param ncore Number of core for parallel computing, the default value is 3.
#' @param seed Specify a random number generator seed, for reproducible results.
#' @param method How the data are permuted: \code{"raw"} (the default) or \code{"freedman-lane"}. See Details.
#' @param perms Optional matrix of permutations to use instead of free ones, one row per observation used in the
#'   model and one column per permutation (e.g. from \code{\link{permindex}}, for permutations within blocks or of
#'   whole plots); \code{nperm} is then its number of columns. Use \code{exact = TRUE} when its columns are all the
#'   permutations the design allows.
#' @return
#' The function produces permutation t-test table for coefficients of (fixed) effects, permutation ANOVA table for model terms and a model parameter list \code{permlist}, a list containing \code{nsim=4999} times permutation refitted \code{model} parameters which are used in functions \code{predictmeans} and \code{contrastmeans}.
#' With \code{method = "freedman-lane"}, \code{permlist} is \code{NULL} and the result has a fourth component,
#' \code{perm}: what is needed to test any comparison with the same permutations. Give the whole result to
#' \code{\link{predictmeans}} or \code{\link{contrastmeans}} (argument \code{permlist}) for Freedman-Lane p-values of
#' pairwise comparisons and contrasts. If refits fail, the number left out is in attribute \code{"failures"}.
#' @details
#' With \code{method = "raw"}, the response is permuted for \code{lm}, \code{aov} and \code{glm} models; for
#' \code{gls}, \code{lme} and \code{lmer} models the fitted values and the decorrelated residuals of the model are
#' permuted together. Each permutation refits the model, and the refits are kept in \code{permlist}. The residuals are
#' decorrelated with the symmetric square root of their covariance, so the results do not depend on the order of the
#' data (predictmeans 1.x used the Cholesky factor, which made tests of terms that vary within blocks reject too often).
#'
#' With \code{method = "freedman-lane"} (Freedman and Lane, 1983), each test has its own permutations: the
#' residuals of the model without the tested term (for type I and II tests also without the terms that contain it;
#' for a t test, without that coefficient) are permuted and added back to that model's fitted values, and the
#' statistic is computed again. The statistics come from QR projections, so no model is refitted and it is fast;
#' weights and offsets are allowed. Use it when the model has a covariate with extreme values (skewed, or with
#' outliers) that is related to the tested term: in a simulation study, permuting the response then rejected
#' 6 to 7\% of true null hypotheses at the 5\% level, and Freedman-Lane 4 to 5\% (Anderson and Legendre (1999) found
#' the same). In other designs the two methods agreed.
#'
#' For mixed models with a covariate with extreme values related to the tested term, the same study found
#' the \code{"raw"} scheme slightly liberal (6.3\% at the 5\% level), while the Satterthwaite F test kept the level
#' (4.6\%); prefer the Kenward-Roger or Satterthwaite F test there (e.g. \code{anova()} of an \pkg{lmerTest} fit, or
#' \code{\link{pm_joint}}), or \code{method = "freedman-lane"} (4.6\%).
#'
#' Since predictmeans 2.0, \code{method = "freedman-lane"} also handles \code{lmer}, \code{lme} and \code{gls} models,
#' and \code{glm}, \code{glmer} and \code{glmmTMB} models. For mixed models the covariance of the data, V, is fixed at
#' the fitted model: the data and the design are decorrelated with the symmetric square root of V, so least squares
#' on them gives the model's estimates and its Wald t and F statistics (unadjusted covariance), and Freedman-Lane
#' permutations are done as for \code{lm}, without refits. GLMs and GLMMs are first linearised at the fit (working
#' response and working variances on the link scale, as in \code{\link{permlmer}}); their t and F statistics are Wald
#' statistics on the working scale with the dispersion estimated from the working residuals, so overdispersion
#' affects them much less than likelihood-ratio tests (in a pilot study with overdispersed counts: 8-10\% against 35\%
#' at the 5\% level). For these classes the tables show these statistics and the permutation
#' p-values; type I tests are sequential in the GLS metric (as \code{anova()} of \pkg{lme4} and \pkg{nlme}). The
#' result can be given to \code{\link{predictmeans}} and \code{\link{contrastmeans}} as for \code{lm}. For \code{glmer} and
#' \code{glmmTMB} models only \code{method = "freedman-lane"} is available.
#' @references
#' Anderson, M. J. and Legendre, P. (1999). An empirical comparison of permutation methods for tests of partial
#' regression coefficients in a linear model. \emph{Journal of Statistical Computation and Simulation}, 62,
#' 271-303.
#'
#' Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of reported significance levels.
#' \emph{Journal of Business & Economic Statistics}, 1, 292-298.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'
#' \donttest{
#' ## Permutation test for model terms (99 permutations to keep the example
#' ## quick; use the default or more in practice)
#' permlme <- permmodels(model=fm, nperm=99, ncore=1)
#'
#' ## Permutation test for multiple comparisons
#' predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH",
#'   permlist=permlme, plot=FALSE, ncore=1)
#'
#' ## Permutation test for specified contrasts
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
#' }
#'
#' ## Freedman-Lane permutations for a mixed model (no refits: fast)
#' fm3 <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#' permmodels(fm3, nperm = 999, type = "II", seed = 1, method = "freedman-lane")
#'
#' ## Freedman-Lane permutations for a linear model with a covariate
#' mtcars$cyl <- factor(mtcars$cyl)
#' fm2 <- lm(mpg ~ wt + cyl, data = mtcars)
#' pfl <- permmodels(fm2, nperm = 999, type = "II", seed = 1, method = "freedman-lane")
#' ## ... and for the pairwise comparisons of the cyl means
#' predictmeans(fm2, "cyl", permlist = pfl, plot = FALSE)$pairs
#' @name permmodels
NULL

## Statistics kept from one permutation refit: coefficients and covariance
## (for permlist), the ANOVA table and the t table. For lmer refits the t
## values are computed from that covariance (as summary() does) instead of
## calling summary(), which would compute it again.
pm_perm_stats <- function(rfitmodel, type, test.statistic, Ls = NULL) {
  mp <- mymodelparm(rfitmodel, df = 0)[1:2]
  if (type %in% c("I", "1")) {
    # lmer: lmerTest's type I F (Wald tests of the hypotheses Ls), as in the observed table
    aT <- if (!is.null(Ls)) pm_type1_F(Ls, mp$coef, mp$vcov) else anova(rfitmodel)
  } else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
  if (inherits(rfitmodel, "lmerMod")) {
    se <- sqrt(diag(mp$vcov))
    coefT <- cbind("Estimate" = mp$coef, "Std. Error" = se, "t value" = mp$coef/se)
  } else coefT <- coef(summary(rfitmodel))
  list(mp, aT, coefT)
}

permmodels <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                       test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
					   data=NULL, fo=NULL, prt=TRUE, ncore=3, seed, method=c("raw", "freedman-lane"), perms=NULL) {
  
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  method <- match.arg(method)
  seed <- if (missing(seed)) NULL else seed
  if (!is.null(perms)) nperm <- ncol(as.matrix(perms))
  if (method == "freedman-lane") {
    if (class(model)[1] == "aovlist") stop("Plese use model 'lme' instead of 'aov'!")
    return(invisible(pm_perm_fl(model, nperm, type, exact, seed, prt, perms)))
  }
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) test.statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test.statistic <- "Chisq"
  if (type %in% c("I", "1")) test.statistic <- "F"	  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, c("glmerMod", "glmmTMB"))) {
    stop("For glmer and glmmTMB models use method = \"freedman-lane\" (permutations on the working scale).",
         call. = FALSE)
  }
  # 'fo': the model formula is an object; refit in this session only
  if (!is.null(fo)) ncore <- 1

  # Each branch sets 'idx' (permutation indices, one column per permutation),
  # 'one' (refit to the k-th permuted data set and return its statistics, or
  # NULL if the refit failed) and the packages a Windows worker needs.
  if (any(inherits(model, "lm"), inherits(model, "aov"), inherits(model, "glm"))) {
    # the response is permuted
    if (!is.null(data)) mod_df <- data else {
    if (inherits(model, "glm")) mod_df <- as.data.frame(model$data) else mod_df <- as.data.frame(model.frame(model))
	}
    if (inherits(model, "mlm")) stop("permmodels() handles models with one response.", call. = FALSE)
    # the variables in the response, permuted together (log(y): y;
    # cbind(R, N - R) or cbind(R, F): both columns)
    yname <- all.vars(attr(terms(model), "variables")[[2]])
    if (!all(yname %in% names(mod_df))) {
      stop("Can't find the response variable(s) ", paste(sQuote(yname), collapse = ", "),
           " in the data; please give the data used to fit the model in 'data'.", call. = FALSE)
    }
    idx <- pm_perm_source(nrow(mod_df), seed, perms)(nperm)
    one <- function(k) {
      d <- mod_df
      d[yname] <- mod_df[idx[, k], yname, drop = FALSE]
      rfitmodel <- pm_refit(model, data = d)
      if (class(rfitmodel)[1] == class(model)[1]) pm_perm_stats(rfitmodel, type, test.statistic)
    }
    packages <- character()
  }
  
  if (any(inherits(model, "gls"), inherits(model, "lme"))) { 
    # fitted values and decorrelated residuals are permuted together
    model <- pm_refit(model, data=as.data.frame(pm_get_data(model)), na.action=na.omit, check=TRUE)
    mod_df <- as.data.frame(pm_get_data(model))
    X <- model.matrix(model)
    if (inherits(model, "gls")) beta <- coef(model) else beta <- nlme::fixef(model)
    # check for columns dropped from model
    col_names <- names(beta)
    if (ncol(X) != length(col_names)) X <- X[,col_names,drop=FALSE]    
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]
    y <- nlme::getResponse(model)
    X_y <- cbind(X, y)
    V_list <- build_Sigma_mats(model, sigma_scale =FALSE)
    ord_X_y <- get_order(A=V_list, B=X_y)
    X_y <- ord_X_y[[1]]
    mod_df <- mod_df[ord_X_y[[2]], ]
    xbeta <- as.vector(X_y[, col_names]%*%beta)
    errors <- X_y[, "y"] - xbeta
    # order-invariant square root of the block-diagonal V (R/whiten.R)
    W <- pm_whiten_blocks(V_list)
    wterrors <- W$solve(errors)    
    idx <- pm_perm_source(length(errors), seed, perms)(nperm)
    if (!is.null(perms)) {
      # the data were reordered by group: map the user's permutations
      ord <- ord_X_y[[2]]
      idx <- matrix(order(ord)[idx[ord, , drop = FALSE]], nrow = length(ord))
    }
    one <- function(k) {
      r <- idx[, k]
      d <- mod_df
      d[, yname] <- xbeta[r] + W$mult(wterrors[r])
      rfitmodel <- pm_refit(model, data = d)
      if (class(rfitmodel)[1] == class(model)[1]) pm_perm_stats(rfitmodel, type, test.statistic)
    }
    packages <- "nlme"
  }
  
  if (inherits(model, "lmerMod")) { 
    # fitted values and decorrelated residuals are permuted together
    lmer_var_names <- all.vars(stats::formula(model))   # (the call may hold a formula object)
	lmer_var_names <- lmer_var_names[lmer_var_names!="pi"]
	mod_df <- as.data.frame(nlme::getData(model))
	mod_df <- na.omit(mod_df[, lmer_var_names])	
	model <- pm_refit(model, data=mod_df, check=TRUE)
	# lmerTest tables (Satterthwaite df and p-values) however the model was fitted
	if (!inherits(model, "lmerModLmerTest")) model <- as_lmerModLmerTest(model)	
    fixef <- fixef(model) 
    Lambda <- getME(model, "Lambda")
    # V = Z Lambda Lambda' Z' + I; W is its order-invariant square root (R/whiten.R)
    W <- pm_whiten(getME(model, "Z") %*% Lambda)
    xbeta <- as.vector(getME(model, "X")%*%fixef)
    errors <- getME(model, "y") - xbeta
    wterrors <- W$solve(errors)
    idx <- pm_perm_source(length(errors), seed, perms)(nperm)
    Ls <- if (type %in% c("I", "1")) pm_type1_L(getME(model, "X"), attr(terms(model), "term.labels")) else NULL
    one <- function(k) {
      r <- idx[, k]
      rfitmodel <- refit(model, xbeta[r] + W$mult(wterrors[r]))
      if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) pm_perm_stats(rfitmodel, type, test.statistic, Ls)
    }
    packages <- "lme4"
  }

  permmod <- pm_perm_drop(pm_perm_map(seq_len(nperm), one, ncore, packages))
  failures <- permmod$failures
  permmod <- permmod$res
  nperm <- length(permmod)
  if (nperm==0) stop("permutation can't produce useful data!")
  failnote <- pm_perm_failnote(failures)

  if (!(inherits(model, "aov"))) {
	model_tTable <- coef(summary(model))
    Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    Tper.p <- pm_perm_p(abs(model_tTable[, Tvalue]), sapply(permmod, function(x) abs(x[[3]][, Tvalue])), exact)
	if ("(Intercept)" %in% rownames(model_tTable)) Tper.p[1] <- NA
	COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p),4)
    if (prt) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEFFICENTS)
      cat("\nNote: Perm_p_value of t test is obtained using ", sQuote(nperm), " permutations", failnote, ".\n", sep = "")
    }
  }else COEFFICENTS <- NULL
  
  if (type %in% c("I", "1")) {
    # (lmerTest's anova() gives type III unless asked for type I)
    key_anova <- if (inherits(model, "lmerModLmerTest")) anova(model, type = "I") else anova(model)
    stats_value <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
    if (class(model)[1] %in% c("glm")) stats_value <- "Deviance"
  }else{
    key_anova <- pm_anova(model, type=type, test.statistic=test.statistic)
    if (class(model)[1] %in% c("lm", "aov")) stats_value <- "F value"  else stats_value <-  test.statistic 
    if (class(model)[1] %in% c("glm")) stats_value <- switch(test.statistic, LR = "LR Chisq", Wald = "Chisq", F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
  }
  stats.p <- pm_perm_p(key_anova[, stats_value], sapply(permmod, function(x) x[[2]][, stats_value]), exact)
  if ("(Intercept)" %in% rownames(key_anova)) stats.p[1] <- NA
  ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4)
  if (!(type %in% c("I", "1"))) attr(ANOVA, "heading") <- attr(key_anova, "heading")
  
  ANOVA[is.na(ANOVA)] <-""
  if (prt) {
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep=""))
    print(ANOVA)
    cat("\nNote: Perm_p_value of ", test.statistic, " test is obtained using ", sQuote(nperm), " permutations",
        failnote, ".\n\n", sep = "")
  }
  permlist <- lapply(permmod, `[[`, 1)
  out <- list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA)
  if (failures > 0) attr(out, "failures") <- failures
  return(invisible(out))
} 
