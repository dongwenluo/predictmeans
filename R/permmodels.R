#' Permutation Test of Linear Model
#'
#' @description
#' This function provides permutation t-tests for coefficients of (fixed) effects and permutation F-tests 
#' for the terms in a linear model such as \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @usage
#' permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
#'   test.statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, 
#'   fo=NULL, prt=TRUE, ncore=3, seed, method=c("raw", "freedman-lane"))
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
#' @param method How the data are permuted: \code{"raw"} (the default) or \code{"freedman-lane"} (\code{lm} and
#'   \code{aov} models only). See Details.
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
#' \code{\link{pm_joint}}).
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
pm_perm_stats <- function(rfitmodel, type, test.statistic) {
  mp <- mymodelparm(rfitmodel, df = 0)[1:2]
  if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- pm_anova(rfitmodel, type=type, test.statistic=test.statistic)
  if (inherits(rfitmodel, "lmerMod")) {
    se <- sqrt(diag(mp$vcov))
    coefT <- cbind("Estimate" = mp$coef, "Std. Error" = se, "t value" = mp$coef/se)
  } else coefT <- coef(summary(rfitmodel))
  list(mp, aT, coefT)
}

permmodels <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                       test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
					   data=NULL, fo=NULL, prt=TRUE, ncore=3, seed, method=c("raw", "freedman-lane")) {
  
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  method <- match.arg(method)
  seed <- if (missing(seed)) NULL else seed
  if (method == "freedman-lane") {
    if (!inherits(model, "lm") || inherits(model, "glm")) {
      stop("method = \"freedman-lane\" is available for lm and aov models.", call. = FALSE)
    }
    return(invisible(pm_perm_fl(model, nperm, type, exact, seed, prt)))
  }
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) test.statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test.statistic <- "Chisq"
  if (type %in% c("I", "1")) test.statistic <- "F"	  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "glmerMod")) stop("This function is not applied to 'glmer' yet!")
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
    idx <- pm_perm_stream(nrow(mod_df), seed)(nperm)
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
    idx <- pm_perm_stream(length(errors), seed)(nperm)
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
    lmer_var_names <- all.vars(model@call$formula)
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
    idx <- pm_perm_stream(length(errors), seed)(nperm)
    one <- function(k) {
      r <- idx[, k]
      rfitmodel <- refit(model, xbeta[r] + W$mult(wterrors[r]))
      if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) pm_perm_stats(rfitmodel, type, test.statistic)
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
    key_anova <- anova(model)
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

## Freedman-Lane permutation tests for lm/aov models (method = "freedman-lane").
## For each test (a term of the ANOVA table, or one coefficient), fit the
## reduced model without it (for type I and II tests also without the terms
## that contain it), permute the reduced model's residuals, add them back to
## its fitted values and recompute the statistic. Weighted fits are handled on
## the sqrt(weights) scale. The statistics come from QR projections of whole
## blocks of permuted responses, so no model is refitted. The same
## permutations are used for every test. Returns the tables of permmodels()
## with permlist = NULL, and 'perm', a permutation object with which
## predictmeans() and contrastmeans() test pairs and contrasts by
## Freedman-Lane permutations (pm_perm_fl_t()).
pm_perm_fl <- function(model, nperm, type, exact, seed, prt) {
  mf <- stats::model.frame(model)
  y <- stats::model.response(mf)
  if (!is.null(dim(y))) stop("method = \"freedman-lane\" needs a single response.", call. = FALSE)
  X <- stats::model.matrix(model)
  asgn <- attr(X, "assign")
  b <- stats::coef(model)
  ok <- !is.na(b)
  X <- X[, ok, drop = FALSE]
  asgn <- asgn[ok]
  off <- stats::model.offset(mf)
  z <- if (is.null(off)) y else y - off
  w <- stats::weights(model)
  if (!is.null(w)) {
    keep <- w > 0
    X <- X[keep, , drop = FALSE] * sqrt(w[keep])
    z <- z[keep] * sqrt(w[keep])
  }
  n <- nrow(X)
  df_res <- n - ncol(X)
  if (df_res < 1) stop("residual df = 0", call. = FALSE)

  Terms <- stats::terms(model)
  labs <- attr(Terms, "term.labels")
  fac <- attr(Terms, "factors")
  all_cols <- seq_len(ncol(X))
  cols_of <- function(terms) which(asgn %in% terms)
  qr_cache <- list()
  get_qr <- function(cols) {
    key <- paste(cols, collapse = ",")
    if (is.null(qr_cache[[key]])) qr_cache[[key]] <<- qr(X[, cols, drop = FALSE])
    qr_cache[[key]]
  }
  rss <- function(cols, Y) {
    if (!length(cols)) return(colSums(Y^2))
    colSums(qr.resid(get_qr(cols), Y)^2)
  }
  rank_of <- function(cols) if (length(cols)) get_qr(cols)$rank else 0L

  ## one test: columns of the reduced model (red) and the two nested models
  ## whose difference in RSS is the test's sum of squares (A within B)
  tests <- list()
  for (j in seq_along(labs)) {
    if (!any(asgn == j)) next
    rel <- pm_relatives(labs[j], labs, fac)
    red <- setdiff(all_cols, cols_of(c(j, rel)))
    AB <- switch(type,
      I = , "1" = list(cols_of(0:(j - 1)), cols_of(0:j)),
      II = , "2" = list(red, setdiff(all_cols, cols_of(rel))),
      list(setdiff(all_cols, cols_of(j)), all_cols))
    if (type %in% c("III", "3")) red <- AB[[1]]
    tests[[labs[j]]] <- list(red = red, A = AB[[1]], B = AB[[2]])
  }
  for (k in all_cols) {
    if (colnames(X)[k] == "(Intercept)") next
    tests[[paste0("coef:", colnames(X)[k])]] <- list(red = setdiff(all_cols, k), A = setdiff(all_cols, k), B = all_cols)
  }
  stat <- function(tst, Y) {
    dfk <- rank_of(tst$B) - rank_of(tst$A)
    rB <- rss(tst$B, Y)
    rF <- if (length(tst$B) == length(all_cols)) rB else rss(all_cols, Y)
    ((rss(tst$A, Y) - rB) / dfk) / (rF / df_res)
  }
  for (k in names(tests)) {
    r <- tests[[k]]$red
    f0 <- if (length(r)) qr.fitted(get_qr(r), z) else rep(0, n)
    tests[[k]]$f0 <- f0
    tests[[k]]$e0 <- z - f0
    tests[[k]]$obs <- round(stat(tests[[k]], matrix(z)), 6)
    tests[[k]]$count <- 0
  }

  perms <- pm_perm_stream(n, seed)
  chunk <- max(1L, min(nperm, floor(2e6 / n)))
  done <- 0L
  while (done < nperm) {
    m <- min(chunk, nperm - done)
    idx <- perms(m)
    for (k in names(tests)) {
      tk <- tests[[k]]
      Y <- tk$f0 + matrix(tk$e0[idx], n)
      tests[[k]]$count <- tk$count + sum(round(stat(tk, Y), 6) >= tk$obs)
    }
    done <- done + m
  }
  pval <- vapply(tests, function(tk) (tk$count + !exact) / (nperm + !exact), numeric(1))

  COEF <- NULL
  if (!inherits(model, "aov")) {        # as with method = "raw"
    ctab <- stats::coef(summary(model))
    COEF <- round(cbind(ctab, Perm_p_value = unname(pval[paste0("coef:", rownames(ctab))])), 4)
  }
  aT <- if (type %in% c("I", "1")) stats::anova(model) else pm_anova(model, type = type, test.statistic = "F")
  ANOVA <- round(cbind(aT, Perm_p_value = unname(pval[rownames(aT)])), 4)
  attr(ANOVA, "heading") <- attr(aT, "heading")
  ANOVA[is.na(ANOVA)] <- ""
  if (prt) {
    if (!is.null(COEF)) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEF)
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "Freedman-Lane permutations.\n")
    }
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep = ""))
    print(ANOVA)
    cat("\nNote: Perm_p_value of F test is obtained using", sQuote(nperm), "Freedman-Lane permutations.\n\n")
  }
  list(permlist = NULL, COEFFICENTS = COEF, ANOVA = ANOVA,
       perm = pm_perm_fl_object(X, z, nperm, exact, attr(perms, "start"), pm_perm_coef_names(model)))
}
