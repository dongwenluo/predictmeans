#' Predicted Means of a Linear Model
#'
#' @description
#' This function obtains predicted means, SE of means, SED of means, LSDs and plots of means
#'   with SE bar or LSD bar for parametric models such as \code{aov}, \code{lm}, 
#'   \code{glm}, \code{gls}, \code{lme}, and \code{lmer}. The function also perfomrs pairwise comparisons
#'   and permutation tests.
#' @usage
#' predictmeans(model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL,
#'   lsd_bar=TRUE, level=NULL, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = I,
#'   transOff = 0, responsen=NULL, count=FALSE, plotord=NULL, lineplot=TRUE, plottitle=NULL, 
#'   plotxlab=NULL, plotylab=NULL, which.plots = c("mean", "pval", "back"),
#'   mplot=NULL, barplot=NULL, pplot=NULL, bkplot=NULL, plot=TRUE, jitterv=0.2, basesz=12, prtnum=TRUE, prtplt=TRUE, newwd=FALSE, permlist=NULL, 
#'   ncore=1, ndecimal=4, cov_reduce=mean, exclude=NULL, vcov_adjust=TRUE,
#'   weights=c("equal", "proportional", "outer", "cells", "flat"),
#'   df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-SAS"), vcov=NULL)
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, \code{lmer}, \code{glmer},
#'   \code{glmmTMB} and \code{gam} (package \pkg{mgcv}), or a \code{\link{semireg}} fit. For a \code{semireg} fit the
#'   means include the smooth terms (as \code{\link{semipred}}), at the covariates' means or at the values of
#'   \code{covariate} (a named vector, e.g. \code{c(construction.date = 1980)}); the result holds only the
#'   \code{means}, \code{pairs} and \code{settings} tables, with normal-based tests and any \code{adj}.
#' @param modelterm Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
#'   The \code{modelterm} must be factors and given exactly as it appears in the printed model, e.g. "A" or "A:B".
#' @param data In some cases, you need to provide the data set used in model fitting, especially when you have applied some 
#'   variable trnasformation in the model.
#' @param pairwise An option for showing pair-wise LSDs and p-values, or not. The default is FALSE.
#' @param atvar When \code{pairwise = TRUE}, a quoted name indicating within levels of which variable in 
#'   \code{modelterm} the multiple comparison will be performed.
#' @param adj Name (in "quote") for indicating a method for adjusting p-values of pairwise comparisons. 
#'   The choices are "none", "tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "sidak", "scheffe"
#'   and "mvt". "sidak", "scheffe" (rank of the family: number of means - 1) and "mvt" (single-step multivariate t, the exact
#'   counterpart of SAS \code{ADJUST=SIMULATE}; needs package \pkg{mvtnorm}) are as in \pkg{emmeans}. With
#'   \code{permlist}, "tukey" and "mvt" give the single-step max-|t| adjustment of the permutation p-values (the
#'   permutation version of Tukey's test: each pair's |t| is compared with the largest |t| of the family in each
#'   permuted data set), and "westfall" the step-down max-|t| of Westfall and Young (1993), which is uniformly more
#'   powerful; both use the correlation between the comparisons. "westfall" needs \code{permlist}; "scheffe" is not
#'   available with it. 
#'   The default method is "none"; letter case does not matter ("Tukey" works). Note that LSD can't be adjusted
#'   except for "bonferroni" method.
#' @param Df A degree of freedom for calculating LSD. For the above models, Df is obtained from the function automatically.
#' @param lsd_bar A logical variable to indicate to print an average LSD or SED bar on the means plot. The default is TRUE.
#' @param level A significant level for calculating LSD, CI etc. The default value is 0.05.
#' @param covariate A numerical vector to specify values of covariates for calculating predicted means. The default 
#'   values are the means of the associated covariates.
#' @param meandecr A logical variable to indicate whether to print letters for multiple comparisons by decreasing order of means in the mean_table. 
#'   The default is NULL which indicates the mean order follows the associated factor levels.
#' @param letterCI A logical variable to indicate printed letters for multiple comparisons by whether or not CI overlap in the mean_table. 
#'   The default is FALSE. With \code{level = NULL} the CIs are 83.4\% CIs: for two uncorrelated means with equal SEs
#'   they fail to overlap when the means differ at about the 5\% level. Otherwise overlap only approximates the
#'   pairwise test: it finds fewer differences when the means are positively correlated, as in most mixed models,
#'   and more when their SEs differ (Maghsoodloo and Huang, 2010). With \code{letterCI = TRUE}, \code{level} sets the CIs
#'   (e.g. \code{level = 0.05} gives 95\% CIs), \code{adj} is "none", and with \code{atvar} the letters compare
#'   the means within each \code{atvar} group. The note under the mean table gives the CI level and the
#'   approximate test level.
#' @param trans A function object for calculating the back transformed means, e.g. \code{trans=exp}.
#' @param transOff When you use \code{trans=exp(x+1)}, then \code{transOff=1}, the default is 0.
#' @param responsen Name (in "quotes") of the back transformed response variable in the \code{model}.
#' @param count An option for indicating the back transformed mean values are counts or not. The default is FALSE.
#' @param plotord A numeric (or character) vector specifying the order of plotting for two or three way interaction (e.g. 
#'   \code{plotord = c(2, 1, 3)} or \code{plotord = c("B", "A", "C")} will put the second variable 'B' in 
#'   \code{modelterm} on the \code{X} axis, the first variable 'A' as the grouping variable, and the third one 
#'   'C' as the panel variable). The defaults are \code{c(1, 2)} and \code{c(1, 2, 3)} for two and three way 
#'   interactions.
#' @param lineplot An option for drawing a line chart, or dot chart. The default is TRUE.
#' @param plottitle A character vector specifying the main title for plot(s). The default is NULL.
#' @param plotxlab A character vector specifying the x label for plot(s). The default is NULL.
#' @param plotylab A character vector specifying the y label for plot(s). The default is NULL.
#' @param which.plots Which plots to draw: any of \code{"mean"} (means plot), \code{"bar"} (bar chart), \code{"pval"}
#'   (p-value plot, when there are more than six p-values) and \code{"back"} (back-transformed plot).
#'   The default is \code{c("mean", "pval", "back")}.
#' @param mplot Deprecated, use \code{which.plots}. An option for drawing a means plot, or not.
#' @param barplot Deprecated, use \code{which.plots}. An option for drawing a bar chart, or not.
#' @param pplot Deprecated, use \code{which.plots}. An option for drawing a p-values plot, or not.
#' @param bkplot Deprecated, use \code{which.plots}. An option for drawing back transformed plot, or not.
#' @param plot An option for drawing plots, or not. The default is TRUE.
#' @param jitterv A degree of jitter in x and y direction in the back transformed means graph. The default is zero.
#' @param basesz The base font size. The default is 12.
#' @param prtnum An option for printing covariate information on the screen, or not. The default is TRUE.
#' @param prtplt An option for printing plots on the screen, or not. The default is TRUE.
#' @param newwd A logical variable to indicate whether to print graph in a new window. The default is FALSE.
#' @param permlist A model parameter list produced by the function \code{permmodels}. When \code{permlist != NULL}, 
#'   the option \code{Df} will be non-functional. This is a key option for pairwise comparisons via 
#'   permutation tests. With the result of \code{permmodels(method = "freedman-lane")} (\code{lm} and \code{aov}
#'   models), each pair is tested by Freedman-Lane permutations of the model in which that pair's means are equal.
#'   With a \code{permlist}, \code{vcov_adjust} is not used: the observed t statistics are compared with those of the
#'   permuted refits, which use the unadjusted covariance.
#' @param ncore Number of cores for parallel computing when \code{permlist != NULL}, the default value is 1 (before predictmeans 2.0: 3).
#' @param ndecimal An option for specifying number of decimal point to be print at predicted means table. The default is 4.
#' @param cov_reduce How covariates that are not varied are handled when building the predicted means: a function such as
#'   \code{mean} (the default) or \code{median}, or \code{"average"} to average the predictions over the
#'   observed covariate values (e.g. means adjusted for a spatial smooth in a \code{gam}). See \code{\link{Kmatrix}}.
#' @param exclude For \code{gam} models: smooth terms to leave out of the predictions. The default excludes random-effect
#'   smooths. See \code{\link{Kmatrix}}.
#' @param vcov_adjust For \code{lmer} models fitted by REML: if \code{TRUE} (the default), standard errors use the Kenward-Roger
#'   adjusted covariance matrix of the fixed effects, as SAS \code{DDFM=KR} and \pkg{emmeans} do, to go with the
#'   Kenward-Roger degrees of freedom. \code{FALSE} uses the unadjusted covariance and reproduces results of
#'   predictmeans 1.x. It has no effect for other models.
#' @param weights Passed to \code{\link{Kmatrix}}: how the cells of the factors that are not in \code{modelterm} are weighted when averaging over them:
#'   \code{"equal"} (the default; equal weights), \code{"proportional"} (weights proportional to the frequencies of
#'   those cells in the data, the same for every mean, as SAS \code{LSMEANS / OM}), \code{"outer"} (products of each
#'   factor's own frequencies), \code{"cells"} (each mean weights the cells by their own frequencies) or
#'   \code{"flat"} (equal weights over the non-empty cells). The definitions are those of \pkg{emmeans}.
#' @param df_method Denominator degrees of freedom for \code{lmer} models. \code{"auto"} (the default) uses Kenward-Roger for
#'   REML fits and Satterthwaite for ML fits. \code{"containment"} gives the SAS \code{DDFM=CONTAIN} df: the
#'   smallest rank contribution to \eqn{[X Z]} of the random effects that contain \code{modelterm}, or the
#'   residual df when none does; all means and comparisons of the term then share that df, and the unadjusted
#'   covariance is used unless \code{vcov_adjust = TRUE} is given. \code{"between-within"} gives the SAS \code{DDFM=BETWITHIN} df, for \code{lmer}, \code{lme} and
#'   \code{gls} models: the subject is the coarsest random-effect grouping (\code{lmer}, \code{lme}) or the grouping of
#'   the correlation structure (\code{gls}); an effect constant within subjects gets the between-subject df, and the
#'   others the df of their subject-by-effect interaction (for example Block x Variety in a split plot); with an
#'   unstructured covariance (\code{gls} with \code{corSymm}) and no random effects every effect gets the
#'   between-subject df, as in SAS. Checked against SAS PROC MIXED. As for containment, the unadjusted covariance is
#'   used unless \code{vcov_adjust = TRUE} is given.
#'   \code{"KR-SAS"} is Kenward-Roger as SAS PROC MIXED \code{DDFM=KR} computes it: with the observed
#'   information of the variance parameters where \code{"Kenward-Roger"} (as \pkg{pbkrtest} and \pkg{emmeans}) uses
#'   the expected information; the two agree in balanced designs and differ slightly otherwise (about 1\% in the df
#'   of a split plot with one plot missing).
#'   Other choices are ignored for other model classes.
#' @param vcov Optional covariance of the coefficients to use instead of the model's. For \code{lm}, \code{aov}
#'   and \code{glm} models: \code{"HC0"}, \code{"HC1"}, \code{"HC2"} or \code{"HC3"} (heteroscedasticity-consistent,
#'   as \code{sandwich::vcovHC()}), or a one-sided formula naming a cluster variable of the data, e.g. \code{~ Block}
#'   (cluster-robust CR1, as \code{sandwich::vcovCL()}). For any model: a covariance matrix of the coefficients, or
#'   a function that returns one from the model (for mixed models e.g. from the \pkg{clubSandwich} package). The
#'   degrees of freedom stay those of the model (the residual df; for \code{lmer} Kenward-Roger or Satterthwaite df
#'   computed from the model), and \code{vcov_adjust} is not used. With a \code{permlist}, \code{vcov} is not used
#'   (for robust permutation tests see \code{permmodels(robust = TRUE)}). The default \code{NULL} uses the model's
#'   covariance.
#' @return
#' \item{Predicted Means}{A table of predicted means.}
#'  \item{Standard Error of Means}{A table of standard errors of predicted means.}
#'  \item{Standard Error of Differences}{Standard errors of differences between predicted means.}
#'  \item{LSD}{Least significant differences between predicted means.}
#'  \item{Pairwise p-value}{A matrix with t-values above the diagonal and p-values below the diagonal, or 
#'  matrix of pairwise comparison p-values for each level of \code{atvar}.}
#'  \item{mean_table}{A summary of predicted means result including 'Predicted means', 'Standard error', 'Df' 
#'  and 'CIs'. When \code{trans!=NULL} or \code{trans!=I}, a table of back transformed means with CIs are also shown.}
#'  \item{predictmeansPlot}{ggplot of predicted means.}
#'  \item{predictmeansBKPlot}{ggplot of back transformed means.}
#'  \item{predictmeansBarPlot}{gg bar plot of predicted means.}
#'  \item{p_valueMatrix}{p_value matrix for pairwise comparison.}
#'  \item{means}{A data frame with one row per predicted mean and fixed column names: the factor (and covariate)
#'  columns, \code{estimate}, \code{std.error}, \code{df}, \code{conf.low}, \code{conf.high}; when the means are
#'  back-transformed also \code{response}, \code{response.std.error} (delta method, as \pkg{emmeans} with
#'  \code{type = "response"}), \code{response.conf.low} and \code{response.conf.high}; \code{group}
#'  (letters, when pairwise comparisons are made), \code{n} (number of observations) and \code{estimable}. Values are
#'  not rounded. Also returned by \code{tidy()} and \code{as.data.frame()}.}
#'  \item{pairs}{A data frame with one row per pairwise comparison: \code{level1}, \code{level2}, \code{by} (the
#'  level of \code{atvar}, when given; \code{NA} for pairs across levels), \code{estimate} (mean of \code{level1}
#'  minus mean of \code{level2}), \code{std.error}, \code{df}, \code{statistic} (t), \code{p.value}, \code{p.adj}
#'  (adjusted by \code{adj}, within each level of \code{atvar} when given), \code{conf.low} and \code{conf.high} (limits of
#'  the difference at \code{level}), \code{lsd} (the half-width of those limits) and \code{p.method}; with \code{permlist}
#'  the limits are permutation-t limits (estimate +- c SE, c the critical value of the permutation test; simultaneous
#'  within each family with \code{adj} "tukey", "mvt" or "westfall"); for a log
#'  link or \code{trans = exp} also \code{ratio} (of the back-transformed means, \code{exp(estimate)}),
#'  \code{ratio.conf.low} and \code{ratio.conf.high}, and for a logit link \code{odds.ratio} with its limits.}
#'  \item{settings}{A list with the options used: model term, significance level, adjustment method, \code{atvar},
#'  degrees-of-freedom method, \code{vcov_adjust}, \code{vcov} (\code{"model"}, \code{"HC3"}, ...), number of permutations, covariate values, model class and
#'  package version.}
#'  The result has class \code{c("predictmeans", "pdmlist")}. Printing shows the same output as earlier versions;
#'  the \code{means}, \code{pairs} and \code{settings} slots are not printed.
#' @note
#' The \code{predictmeans} function becomes confused if a factor or covariate is changed to the other
#'   in a model formula. Consequently, formulae that include calls \code{as.factor}, \code{factor}, or \code{numeric}
#'   (e.g. \code{as.factor(income)}) will cause errors. Instead, create the modified variables outside of the model
#'   formula (e.g., \code{fincome <- as.factor(income)}) and then use them in the model formula.
#'   
#'   Factors cannot have colons in level names (e.g., \code{"level:A"}) for the legacy tables; the \code{predictmeans}
#'   function will confuse the colons with interactions; rename levels to avoid colons. The \code{means} and
#'   \code{pairs} tables are not affected.
#'   
#'   For \code{predictmeans} function, it is assumed that methods \code{coef}, \code{vcov}, \code{model.matrix}, \code{model.frame} and \code{terms} are available for \code{model}.
#' @references
#' Torsten Hothorn, Frank Bretz and Peter Westfall (2008),
#'     \emph{Simultaneous Inference in General Parametric Models. Biometrical}, Journal 50(3), 346--363.
#'   
#' 	Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004),
#' 	\emph{Prediction in linear mixed models}, Australian and New Zealand Journal of Statistics, 46(3), 325-347. 
#' 	
#' 	Vanessa, C. (05 October 2022),
#' 	\emph{Confidence tricks: the 83.4\% confidence interval for comparing means}, https://vsni.co.uk/blogs/confidence_trick.
#' 
#' 	Maghsoodloo, S. and Huang, C.-Y. (2010), \emph{Comparing the overlapping of two independent confidence intervals
#' 	with a single confidence interval for two normal population parameters}, Journal of Statistical Planning and
#' 	Inference, 140(11), 3295-3305.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   predictmeans(fm, "nitro", adj="BH")
#'   predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", line=FALSE)
#'   predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", 
#'     barplot=TRUE, line=FALSE)
#'   names(predictout)
#'   print(predictout$predictmeansPlot)
#'   print(predictout$predictmeansBarPlot)
#' @name predictmeans
NULL

predictmeans <- function (model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL, lsd_bar=TRUE,
                          level=NULL, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = I, transOff=0,
						  responsen=NULL, count=FALSE, plotord=NULL, lineplot=TRUE, plottitle=NULL, plotxlab=NULL,
						  plotylab=NULL,
						  which.plots = c("mean", "pval", "back"),
						  mplot=NULL, barplot=NULL, pplot=NULL, bkplot=NULL, plot=TRUE, jitterv=0.2,
						  basesz=12, prtnum=TRUE, prtplt=TRUE, newwd=FALSE, permlist=NULL, ncore=1, ndecimal=4,
                          cov_reduce=mean, exclude=NULL, vcov_adjust=TRUE,
                          weights=c("equal", "proportional", "outer", "cells", "flat"),
                          df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-SAS"), vcov=NULL) {
  vcov_adjust_default <- missing(vcov_adjust)

  ## Bias output towards fixed notation. See ?options for more info
  op <- options(scipen = 6)
  on.exit(options(op), add = TRUE)

  if(any(missing(model), missing(modelterm))){
    stop("The arguments 'model', and 'modelterm' must be provided!")
  }

  # semireg() fits: means with the smooth terms included (see semireg-means.R)
  if (inherits(model, "semireg")) {
    return(pm_semireg_means(model, modelterm, covariate = covariate, level = if (is.null(level)) 0.05 else level,
                            adj = adj, atvar = atvar))
  }

  weights <- match.arg(weights)
  df_method <- match.arg(df_method)
  # SAS DDFM=CONTAIN uses the unadjusted covariance
  if (df_method %in% c("containment", "between-within") && vcov_adjust_default) vcov_adjust <- FALSE
  vcov_adjust <- pm_perm_vcov_adjust(model, permlist, vcov_adjust, vcov_adjust_default)
  vcov <- pm_vcov_with_perm(vcov, permlist)
  if (!is.null(vcov)) vcov_adjust <- FALSE     # the given covariance replaces the adjusted one

  # The calculation runs in stages that share this call's state (see
  # predictmeans-stages.R).
  st <- environment()
  pm_stage_setup(st)
  pm_stage_estimate(st)
  pm_stage_compare(st)
  pm_stage_mean_table(st)
  pm_stage_backtransform(st)
  if (plot) {
    pm_stage_plots(st)
  }
  pm_stage_assemble(st)
  get("pm_out", envir = st)
}
