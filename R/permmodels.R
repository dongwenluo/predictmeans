#' Permutation Tests of the Terms and Coefficients of a Model
#'
#' @description
#' Permutation t tests of the coefficients of the fixed effects and permutation F (or chi-square) tests of the
#' terms of a model fitted by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, \code{lmer}, \code{glmer} or
#' \code{glmmTMB}. p-values are (b + 1)/(B + 1): b of the B permuted statistics at least as large as the observed
#' one. Results are the same for any \code{ncore}.
#' @usage
#' permmodels(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3), 
#'   test_statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, ncore=1, seed,
#'   method=c("auto", "raw", "freedman-lane"), perms=NULL, refit=TRUE, robust=FALSE, ...)
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, \code{lmer},
#'   \code{glmer}, \code{glmmTMB} or \code{ordinal::clmm}. For a \code{clmm}, each main-effect term of a variable
#'   applied to whole clusters (constant within the clusters of a random effect) is tested by re-randomising it:
#'   its values are relabelled among the clusters within strata (the other cluster-level variables, such as fixed
#'   blocks, and coarser random effects the clusters are nested in), the model refitted and its Wald chi-square
#'   compared. This is exact for the null hypothesis that the variable has no effect at all (main effect and
#'   interactions), whatever the within-cluster effects and cluster sizes. The clusters are those of the coarsest
#'   random effect within which the variable is constant; numeric covariates stay with their clusters. \code{perms}
#'   is refused and \code{exact} not used. When the allocations (one per relabelling of the values, which give the
#'   same statistic, for equal group sizes) are at most \code{nperm} all are used (exact p-values); refits that do
#'   not converge count as at least as extreme. A note gives, for each variable, the strata, the number of
#'   allocations and the smallest attainable p-value (6 clusters split 3/3: 20 allocations, p at least 0.1);
#'   other terms get none. The result's attribute \code{"resampling"} is \code{"cluster randomisation"}.
#'   \code{\link{bootmodels}} tests every term, by a parametric bootstrap under the global null.
#' @param nperm The number of permutations. The default is 4999.
#' @param type The type of the tests: "I", "II" or "III" (or 1, 2, 3).
#'   Type III tests are those of sum-to-zero contrasts whatever contrasts the model was fitted with (as SAS Type 3);
#'   \code{car::Anova()} gives them only for models fitted with such contrasts (with R's default treatment contrasts
#'   it tests each main effect at the first level of the other factors).
#' @param test_statistic The statistic of Type II and III tests with \code{method = "raw"}: F for \code{lm};
#'   chi-square for \code{lme} and \code{gls}; \code{"Chisq"} or \code{"F"} for \code{lmer}; \code{"LR"},
#'   \code{"Wald"} or \code{"F"} for \code{glm}. Type I tests use F. A statistic the model does not have is
#'   replaced by the one it has, with a message; with \code{method = "freedman-lane"} the F (or t) statistics are
#'   permuted and a different \code{test_statistic} is ignored with a message.
#'   The F tests of an \code{lmer} model are Kenward-Roger F tests with the observed information of the variance
#'   parameters, as in the rest of predictmeans (\code{car::Anova()} uses the expected information; the two agree
#'   for balanced data).
#' @param exact \code{TRUE}: \code{perms} holds every permutation the design allows (a complete enumeration), so
#'   p-values are the proportion of them with a statistic at least as large as the observed one. Otherwise (the
#'   default \code{FALSE}) the observed data count as one of the permutations: p = (b + 1)/(B + 1), never 0.
#'   \code{exact = TRUE} without \code{perms}, or with repeated columns, is ignored with a warning.
#' @param data The data the model was fitted to, when they cannot be found from the model's call (e.g. after
#'   transforming variables outside the data).
#' @param ... Not used: \code{prt} and \code{fo} (before 2.0) stop with their 2.0 spelling. Nothing is printed during
#'   the call: \code{print()} of the result shows the tables. A model whose formula is held in a variable
#'   (\code{lm(f, data)}; 1.x's \code{fo}) is detected and refitted in this session only.
#' @param ncore Number of cores (default 1; before predictmeans 2.0: 3). Each permutation has its own seed, so
#'   the results do not depend on \code{ncore}.
#' @param seed A random number seed, for reproducible results; the user's random number stream is left as it was.
#' @param method How the data are permuted: \code{"freedman-lane"} or \code{"raw"}; \code{"auto"} (the default)
#'   chooses the more robust one for the model class, which is \code{"freedman-lane"} for every class at present.
#'   \code{"raw"} is the scheme of predictmeans 1.x. See Details.
#' @param perms Optional matrix of permutations to use instead of free ones, one row per observation used in the
#'   model and one column per permutation (e.g. from \code{\link{permindex}}, for permutations within blocks or of
#'   whole plots); \code{nperm} is then its number of columns. Use \code{exact = TRUE} when its columns are all the
#'   permutations the design allows.
#' @param refit For \code{method = "freedman-lane"} with mixed models (\code{lmer}, \code{lme}, \code{gls}, and
#'   \code{glmer}, \code{glmmTMB} through their working linear mixed model): refit the model to each permuted data set
#'   (the default), so the variance parameters are estimated again each time.
#'   \code{FALSE} holds the covariance of the data at the fit (no refits, much faster), which rejects too often
#'   when the variance parameters rest on few units. In a simulation (test of the whole-plot factor of a split plot,
#'   true null, 5\% level) \code{refit = FALSE} rejected 8.6\%, 7.2\%, 6.0\% and 5.6\% with 6, 12, 24 and 48 blocks
#'   (whole-plot error df 15, 33, 69 and 141), and refits 5.0\%, 6.2\%, 5.4\% and 4.6\% (Kenward-Roger 4.6\% to
#'   6.4\%); at the 1\% level \code{refit = FALSE} rejected up to 2.8\%. So \code{refit = FALSE} suits a quick look,
#'   or a final test when the error df of the tested terms exceed about 100; otherwise keep the refits (for 199
#'   permutations of the split plot's three terms on one core: 8 s with 6 blocks, 25 s with 48).
#' @param robust For \code{lm}, \code{aov} and \code{glm} models with \code{method = "freedman-lane"}: tests for
#'   error variances that are unequal while the model assumes them equal. The statistics are
#'   heteroscedasticity-robust (HC3 sandwich) t and Wald F statistics (DiCiccio and Romano, 2017), and the signs of
#'   the reduced model's residuals are flipped at random instead of permuted (a wild bootstrap; Davidson and
#'   Flachaire, 2008), so each residual stays with its own unit and keeps its variance. In a simulation with
#'   group SDs 3, 1, 1 (2000 data sets per case) the usual and permutation F tests rejected 19-24\% of true null
#'   hypotheses at the 5\% level, these tests 5-7\% (normal errors) or about 9\% (strongly skewed errors, against
#'   11-14\% for the HC3 F test), and 5.4\% with equal variances. The tables show the robust statistics; pairs
#'   and contrasts from the
#'   result use them too (robust t, and confidence limits from robust standard errors). Type II tests match
#'   \code{car::Anova(model, type = "II", white.adjust = "hc3")}. Not with \code{perms}. The default is
#'   \code{FALSE}.
#' @return
#' A list of class \code{"pm_permmodels"}: \code{COEFFICENTS} (permutation t tests of the coefficients),
#' \code{ANOVA} (permutation tests of the terms, of the \code{type} asked for) and \code{permlist} (with
#' \code{method = "raw"}: the coefficients and covariance of each refit, for \code{\link{predictmeans}} and
#' \code{\link{contrastmeans}}). Values are stored rounded to 4 decimals; the notes give the smallest attainable
#' p-value, 1/(B + 1). Advice on the number of permutations is kept in \code{attr(, "advice")} (\code{\link{advice}}).
#' With \code{method = "freedman-lane"}, \code{permlist} is \code{NULL} and the result has a fourth component,
#' \code{perm}: what is needed to test any comparison with the same permutations. Give the whole result to
#' \code{\link{predictmeans}} or \code{\link{contrastmeans}} (argument \code{permlist}) for Freedman-Lane p-values of
#' pairwise comparisons and contrasts. If refits fail, the number left out is in attribute \code{"failures"}.
#' @details
#' With \code{method = "auto"} (the default since predictmeans 2.0) the method is chosen by model class from a
#' simulation study (\code{dev/permutation-plan.md} in the package sources): at present Freedman-Lane for every
#' class. Use \code{method = "raw"} for the results of predictmeans 1.x.
#'
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
#' (4.6\%); prefer the Kenward-Roger or Satterthwaite F test there (\code{\link{anova.lmerMod}} or
#' \code{\link{pm_joint}}), or \code{method = "freedman-lane"} (4.6\%).
#'
#' Since predictmeans 2.0, \code{method = "freedman-lane"} also handles \code{lmer}, \code{lme} and \code{gls} models,
#' and \code{glm}, \code{glmer} and \code{glmmTMB} models. For \code{lmer}, \code{lme} and \code{gls} models (with
#' \code{refit = TRUE}, the default) the model without the tested term is fitted, its residuals are decorrelated with
#' the symmetric square root of its covariance, permuted, correlated again and added to its fitted values, and the
#' model is refitted to each permuted data set; the statistics are the model's Wald t and F statistics (unadjusted
#' covariance). One set of refits serves a term and its coefficients; pairs and contrasts in \code{predictmeans()}
#' get one set of refits for their family (the model with all the compared means equal). With \code{refit = FALSE}
#' the covariance of the data, V, is fixed at the fitted model: the data and the design are decorrelated with the
#' symmetric square root of V, so least squares on them gives the same statistics, and Freedman-Lane permutations
#' are done as for \code{lm}, without refits. That is much faster, but the permuted statistics then behave as if V
#' were known, so comparisons whose error rests on few units rejected too many true null hypotheses in a simulation
#' study (whole-plot factor of a split plot with 6 blocks: 8.6\% at the 5\% level, against 5.0\% with refits; pairs of
#' whole-plot means 6.2\%). GLMs and GLMMs are first linearised at the fit (working
#' response and working variances on the link scale, as in \code{\link{permlmer}}); GLMMs are then treated as the
#' working linear mixed model (weights from the working variances), refitted to each permuted working response when
#' \code{refit = TRUE}. Their t and F statistics are Wald statistics on the working scale with the dispersion
#' estimated from the working residuals, so overdispersion
#' affects them much less than likelihood-ratio tests (in a pilot study with overdispersed counts: 8-10\% against 35\%
#' at the 5\% level). For these classes the tables show these statistics and the permutation
#' p-values; type I tests are sequential in the GLS metric (as \code{anova()} of \pkg{lme4} and \pkg{nlme}). The
#' result can be given to \code{\link{predictmeans}} and \code{\link{contrastmeans}} as for \code{lm}. For \code{glmer} and
#' \code{glmmTMB} models only \code{method = "freedman-lane"} is available.
#' @references
#' Davidson, R. and Flachaire, E. (2008). The wild bootstrap, tamed at last. \emph{Journal of Econometrics},
#' 146, 162-169.
#'
#' DiCiccio, C. J. and Romano, J. P. (2017). Robust permutation tests for correlation and regression
#' coefficients. \emph{Journal of the American Statistical Association}, 112, 1211-1220.
#'
#' Anderson, M. J. and Legendre, P. (1999). An empirical comparison of permutation methods for tests of partial
#' regression coefficients in a linear model. \emph{Journal of Statistical Computation and Simulation}, 62,
#' 271-303.
#'
#' Freedman, D. and Lane, D. (1983). A nonstochastic interpretation of reported significance levels.
#' \emph{Journal of Business & Economic Statistics}, 1, 292-298.
#' @seealso \code{\link{permlmer}}, \code{\link{permindex}}, \code{\link{permstats}}, \code{\link{permcombine}}, \code{\link{bootmodels}}, \code{\link{advice}}
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
#'   permlist=permlme)
#'
#' ## Permutation test for specified contrasts
#' cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#'             c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#' contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)
#'
#' ## The same by lmer, type II tests; refit = FALSE is much faster (no refits)
#' fm3 <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#' permmodels(fm3, nperm = 99, type = "II", seed = 1, ncore = 1)
#' permmodels(fm3, nperm = 999, type = "II", seed = 1, refit = FALSE)
#' }
#'
#' ## Freedman-Lane permutations for a linear model with a covariate
#' mtcars$cyl <- factor(mtcars$cyl)
#' fm2 <- lm(mpg ~ wt + cyl, data = mtcars)
#' pfl <- permmodels(fm2, nperm = 999, type = "II", seed = 1, method = "freedman-lane")
#' ## ... and for the pairwise comparisons of the cyl means
#' predictmeans(fm2, "cyl", permlist = pfl)$pairs
#' ## Robust (HC3) statistics, for error variances that differ between the cyl groups
#' permmodels(fm2, nperm = 999, type = "II", seed = 1, robust = TRUE)
#' @name permmodels
NULL

## Statistics kept from one permutation refit: coefficients and covariance
## (for permlist), the ANOVA table and the t table. For lmer refits the t
## values are computed from that covariance (as summary() does) instead of
## calling summary(), which would compute it again. 'map3': the Type III map
## of the observed fit (R/anova.R), the same for every refit.
pm_perm_stats <- function(rfitmodel, type, test_statistic, Ls = NULL, map3 = NULL) {
  mp <- pm_model_parm(rfitmodel, df = 0)[1:2]
  if (type %in% c("I", "1")) {
    # lmer: lmerTest's type I F (Wald tests of the hypotheses Ls), as in the observed table
    aT <- if (!is.null(Ls)) pm_type1_F(Ls, mp$coef, mp$vcov) else pm_anova_seq(rfitmodel)
  } else aT <- pm_anova(rfitmodel, type=type, test_statistic=test_statistic, type3_map=map3)
  if (inherits(rfitmodel, "lmerMod")) {
    se <- sqrt(diag(mp$vcov))
    coefT <- cbind("Estimate" = mp$coef, "Std. Error" = se, "t value" = mp$coef/se)
  } else coefT <- pm_coef_table(rfitmodel)                     # (R/boot-parts.R; clmm: latent coefficients)
  list(mp, aT, coefT)
}

#' @export
permmodels <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                       test_statistic=c("Chisq", "F", "LR", "Wald"), exact=FALSE, data=NULL, ncore=1, seed,
                       method=c("auto", "raw", "freedman-lane"), perms=NULL, refit=TRUE, robust=FALSE, ...) {
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("permmodels", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  # a model whose formula is held in a variable (lm(f, data)) is refitted in
  # this session only: the parallel workers would not find the variable
  # (1.x asked for it in 'fo')
  if (pm_formula_is_object(model)) ncore <- 1
  
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)
  type <- pm_anova_type(type[1])
  ts_given <- !missing(test_statistic)
  test_statistic <- match.arg(test_statistic)
  method <- match.arg(method)
  if (method == "auto") method <- pm_perm_auto_method(model)
  seed <- if (missing(seed)) NULL else seed
  if (inherits(model, "clmm")) {                # randomisation of whole clusters (R/perm-clmm.R)
    if (isTRUE(robust)) stop("robust = TRUE needs method = \"freedman-lane\", which clmm models do not have.",
                             call. = FALSE)
    if (!is.null(perms)) stop("permmodels() of a clmm permutes the treatments of whole clusters itself; leave ",
                              "'perms' out.", call. = FALSE)
    if (method != "raw") stop("clmm models are permuted as whole clusters (method = \"raw\").", call. = FALSE)
    if (type %in% c("I", "1")) {
      message("clmm models have no sequential (Type I) tests; Type II Wald tests are used.")
      type <- "II"
    }
    if (ts_given && test_statistic != "Chisq") message("clmm tests are Wald chi-square; test_statistic = \"",
                                                        test_statistic, "\" is not used.")
    if (!missing(exact)) message("For clmm the allocations are enumerated (exact) whenever there are at most ",
                                 "'nperm'; 'exact' is not used.")
    return(pm_perm_clmm(model, nperm, type, data, seed, ncore))
  }
  if (!is.null(perms)) nperm <- ncol(as.matrix(perms))
  exact <- pm_perm_exact(exact, perms)
  if (isTRUE(robust)) {
    if (method != "freedman-lane") stop("robust = TRUE needs method = \"freedman-lane\".", call. = FALSE)
    return(pm_perm_fl_robust(model, nperm, type, exact, seed, perms))
  }
  if (method == "freedman-lane") {
    if (ts_given && test_statistic != "F")
      message("test_statistic = \"", test_statistic, "\" is used only with method = \"raw\"; ",
              "Freedman-Lane permutes F (or t) statistics.")
    if (class(model)[1] == "aovlist") stop("permmodels() does not handle aov models with Error(): fit the error strata as random effects (lmer or ",
                                                 "lme) instead.", call. = FALSE)
    return(pm_perm_fl(model, nperm, type, exact, seed, perms, refit, ncore))
  }
  asked <- test_statistic
  if (inherits(model, "glm") && !(test_statistic %in% c("LR", "Wald", "F"))) test_statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test_statistic <- "Chisq"
  if (type %in% c("I", "1")) test_statistic <- "F"
  # tell only a user who chose the statistic, not one who left the default
  if (ts_given && test_statistic != asked)
    message("test_statistic = \"", asked, "\" is not available here; \"", test_statistic, "\" is used (",
            if (type %in% c("I", "1")) "Type I tests are F tests" else
              if (inherits(model, c("gls", "lme"))) "gls/lme tests are Chisq" else
                "glm tests are LR, Wald or F", ").")
  if (class(model)[1]=="aovlist") stop("permmodels() does not handle aov models with Error(): fit the error strata as random effects (lmer or ",
                                                 "lme) instead.", call. = FALSE)
  if (inherits(model, c("glmerMod", "glmmTMB"))) {
    stop("For glmer and glmmTMB models use method = \"freedman-lane\" (permutations on the working scale).",
         call. = FALSE)
  }

  # each class sets the permutation indices, a function that refits to the
  # k-th permuted data set and returns its statistics (NULL if the refit
  # failed), and the packages a Windows worker needs (pm_perm_raw())
  parts <- pm_perm_raw(model, data, seed, perms, nperm, type, test_statistic)
  model <- parts$model

  # the Type III map depends only on the design, which the refits share
  map3 <- if (type %in% c("III", "3")) tryCatch({
    b0 <- pm_coef(model)
    pm_type3_map(model, names(b0)[!is.na(b0)])
  }, error = function(e) NULL)
  one <- function(k) parts$one(k, map3)
  permmod <- pm_perm_drop(pm_perm_map(seq_len(nperm), one, ncore, parts$packages))
  failures <- permmod$failures
  permmod <- permmod$res
  nperm <- length(permmod)
  if (nperm==0) stop("permutation can't produce useful data!")
  failnote <- pm_perm_failnote(failures)

  # lmer: the tables with Satterthwaite df and p-values (as lmerTest's), however
  # the model was fitted (R/anova.R)
  lt <- if (inherits(model, "lmerMod")) tryCatch(pm_lmer_tables(model), error = function(e) NULL)
  if (!(inherits(model, "aov"))) {
	model_tTable <- if (!is.null(lt)) lt$coef else pm_coef_table(model)
    Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    Tper.p <- pm_perm_p(abs(model_tTable[, Tvalue]), sapply(permmod, function(x) abs(x[[3]][, Tvalue])), exact)
	if ("(Intercept)" %in% rownames(model_tTable)) Tper.p[1] <- NA
	COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p),4)
  }else COEFFICENTS <- NULL
  
  if (type %in% c("I", "1")) {
    key_anova <- if (!is.null(lt)) lt$anova else pm_anova_seq(model)
    stats_value <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
    if (class(model)[1] %in% c("glm")) stats_value <- "Deviance"
  }else{
    key_anova <- pm_anova(model, type=type, test_statistic=test_statistic)
    if (class(model)[1] %in% c("lm", "aov")) stats_value <- "F value"  else stats_value <-  test_statistic 
    if (class(model)[1] %in% c("glm")) stats_value <- switch(test_statistic, LR = "LR Chisq", Wald = "Chisq", F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
  }
  stats.p <- pm_perm_p(key_anova[, stats_value], sapply(permmod, function(x) x[[2]][, stats_value]), exact)
  if ("(Intercept)" %in% rownames(key_anova)) stats.p[1] <- NA
  ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4)
  if (!(type %in% c("I", "1"))) attr(ANOVA, "heading") <- attr(key_anova, "heading")
  
  ANOVA[is.na(ANOVA)] <-""
  permlist <- lapply(permmod, `[[`, 1)
  out <- pm_permmodels_result(
    list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA),
    if (!is.null(COEFFICENTS)) paste0("Note: Perm_p_value of t test is obtained using ", sQuote(nperm), " permutations",
                                      failnote, "."),
    paste0("Note: Perm_p_value of ", test_statistic, " test is obtained using ", sQuote(nperm), " permutations",
           failnote, "."), nperm = nperm, exact = exact)
  if (failures > 0) attr(out, "failures") <- failures
  out
}

## TRUE when the model's formula is held in a variable (lm(f, data)) rather
## than written in the call.
pm_formula_is_object <- function(model) {
  cl <- tryCatch(stats::getCall(model), error = function(e) NULL)
  if (is.null(cl)) return(FALSE)
  f <- if (!is.null(cl$formula)) cl$formula else if (!is.null(cl$fixed)) cl$fixed else if (length(cl) > 1) cl[[2]]
  # only a formula written out (y ~ x) or stored as a formula object is safe
  # for the workers; anything else (a name, a call) runs on one core
  !(inherits(f, "formula") || (is.call(f) && identical(f[[1]], as.name("~"))))
}

## The raw permutation scheme of permmodels() for each class: the model (nlme
## fits are refitted on their complete rows first), the permutation indices
## 'idx' (one column per permutation), one(k, map3) refitting the model to the
## k-th permuted data set and returning its statistics (NULL if the refit
## failed), and the packages a Windows worker needs.
pm_perm_raw <- function(model, data, seed, perms, nperm, type, test_statistic) UseMethod("pm_perm_raw")

#' @exportS3Method NULL
pm_perm_raw.default <- function(model, data, seed, perms, nperm, type, test_statistic) {
  stop("permmodels() handles lm, aov, glm, gls, lme and lmer models with method = \"raw\"; use method = ",
       "\"freedman-lane\" for glmer and glmmTMB.", call. = FALSE)
}

#' @exportS3Method NULL
pm_perm_raw.lm <- function(model, data, seed, perms, nperm, type, test_statistic) {
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
  one <- function(k, map3) {
    d <- mod_df
    d[yname] <- mod_df[idx[, k], yname, drop = FALSE]
    rfitmodel <- pm_refit(model, data = d)
    if (class(rfitmodel)[1] == class(model)[1]) pm_perm_stats(rfitmodel, type, test_statistic, map3 = map3)
  }
  packages <- character()
  list(model = model, idx = idx, one = one, packages = packages)
}

#' @exportS3Method NULL
pm_perm_raw.lme <- function(model, data, seed, perms, nperm, type, test_statistic) {
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
  # rows by group, and the order-invariant square root of the block-diagonal V (R/whiten.R)
  nw <- pm_nlme_whiten(model, cbind(X, y))
  X_y <- nw$B
  mod_df <- mod_df[nw$ord, ]
  xbeta <- as.vector(X_y[, col_names]%*%beta)
  errors <- X_y[, "y"] - xbeta
  W <- nw$W
  wterrors <- W$solve(errors)    
  idx <- pm_perm_source(length(errors), seed, perms)(nperm)
  if (!is.null(perms)) {
    # the data were reordered by group: map the user's permutations
    ord <- nw$ord
    idx <- matrix(order(ord)[idx[ord, , drop = FALSE]], nrow = length(ord))
  }
  one <- function(k, map3) {
    r <- idx[, k]
    d <- mod_df
    d[, yname] <- xbeta[r] + W$mult(wterrors[r])
    rfitmodel <- pm_refit(model, data = d)
    if (class(rfitmodel)[1] == class(model)[1]) pm_perm_stats(rfitmodel, type, test_statistic, map3 = map3)
  }
  packages <- "nlme"
  list(model = model, idx = idx, one = one, packages = packages)
}

#' @exportS3Method NULL
pm_perm_raw.gls <- pm_perm_raw.lme

## The distinct orderings of a multiset (integers 1..k with repeats), one per
## column.
pm_multiset_perms <- function(v) {
  rec <- function(left) {
    if (!length(left)) return(list(integer(0)))
    out <- list()
    for (u in unique(left)) {
      rest <- left[-match(u, left)]
      out <- c(out, lapply(rec(rest), function(r) c(u, r)))
    }
    out
  }
  do.call(cbind, rec(sort(v)))
}

## A clmm refit that converged (a gradient near 0): refits that did not are
## failures (left out), not statistics.
pm_clmm_ok <- function(fit) {
  # (the gradient only: optimisers code success differently, e.g. nlminb 0, ucminf 1 or 2)
  inherits(fit, "clmm") && all(is.finite(fit$gradient)) && max(abs(fit$gradient)) < 1e-2
}

## The rows a clmm was fitted to (its data without the rows it dropped for NA).
pm_clmm_rows <- function(model) {
  d <- as.data.frame(pm_get_data(model))
  n <- length(model$gfList[[1]])
  if (nrow(d) != n && !is.null(model$na.action)) d <- d[-as.integer(model$na.action), , drop = FALSE]
  if (nrow(d) != n) stop("predictmeans could not find the rows the clmm was fitted to; give its data without ",
                         "missing values.", call. = FALSE)
  d
}

#' @exportS3Method NULL
pm_perm_raw.lmerMod <- function(model, data, seed, perms, nperm, type, test_statistic) {
  # fitted values and decorrelated residuals are permuted together. refit()
  # replaces only the response of the rows used, so prior weights, offsets and
  # NA handling stay as fitted (2.0 refitted to the formula variables: weights lost)
  # V = Z Lambda Lambda' Z' + diag(1/w); W its order-invariant square root (R/whiten.R)
  pp <- pm_mer_perm_parts(model)
  W <- pp$W
  xbeta <- pp$xbeta
  wterrors <- pp$wterrors
  idx <- pm_perm_source(length(xbeta), seed, perms)(nperm)
  Ls <- if (type %in% c("I", "1")) pm_type1_L(getME(model, "X"), attr(terms(model), "term.labels")) else NULL
  one <- function(k, map3) {
    r <- idx[, k]
    rfitmodel <- pm_mer_refit(model, xbeta[r] + W$mult(wterrors[r]))
    if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) pm_perm_stats(rfitmodel, type, test_statistic, Ls, map3)
  }
  packages <- "lme4"
  list(model = model, idx = idx, one = one, packages = packages)
}
