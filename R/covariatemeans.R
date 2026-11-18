#' Predicted Means of a Linear Model with Covariate Variable(s)
#'
#' @description
#' Predicted means along a covariate for each level of a factor term (curves) with their confidence limits, or
#' the difference between two of the curves with its band. The result holds the tables in the layout of
#' \code{\link{predictmeans}} and the plot, which \code{autoplot()} returns.
#' @usage
#' covariatemeans(model, modelterm=NULL, covariate, as.is=FALSE, at=NULL, data=NULL,
#'   conf.level=0.95, trans=NULL, cov_reduce=mean, exclude=NULL,
#'   weights=c("equal", "proportional", "outer", "cells", "flat"),
#'   df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
#'   vcov_adjust=TRUE, contr=NULL, contr_scale=c("model", "response"),
#'   band=c("pointwise", "simultaneous"), bias_adjust=FALSE, sigma=NULL,
#'   mean_type=c("conditional", "overall"), plot=getOption("predictmeans.plot", FALSE),
#'   plot_args=list(), ...)
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
#'   The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B".
#' @param covariate Name (in "quotes") of one the covariate variables in the \code{model}.
#' @param as.is \code{TRUE}: each line only over the range of the covariate observed in its group (rows outside it are
#'   left out of the result); the default \code{FALSE} uses the values of \code{at} for every group.
#' @param at The values of \code{covariate} at which the lines are predicted (\code{covariateV} before 2.0); by
#'   default 50 values over the observed range. See \code{\link{Kmatrix}}.
#' @param data In some cases, you need to provide the data set used in model fitting, especially when you have applied some variable trnasformation in the model.
#' @param conf.level The confidence level of the limits, the default is 0.95 (\code{level = 0.05} before 2.0). The limits
#'   are for each line on its own: whether the bands of two groups overlap is not a test of their difference
#'   (the lines are estimated from the same model and are correlated; see \code{letterCI} in
#'   \code{\link{predictmeans}}). Use \code{contr} for the band of the difference between two groups, or compare
#'   groups at chosen covariate values with \code{\link{predictmeans}} (\code{covariate = }) and their slopes
#'   with \code{\link{pm_slopes}}.
#' @param trans A function object for calculating the back transformed means, e.g. \code{trans=exp}. When it is
#'   not given, the inverse link of the model is used, or the inverse of a transformation of the response in the
#'   model formula (\code{log(y)}, \code{sqrt(y)}, ...; see \code{\link{predictmeans}}), with a message;
#'   \code{trans = I} keeps the model scale. The back-transformed line is the back-transformed fitted value.
#'   (Before 2.0 half the standard error of the fitted value was added inside \code{exp()} for a log.)
#' @param plot_args A named list of plot settings (names may be abbreviated):
#'   \describe{
#'     \item{\code{trellis}}{\code{TRUE} (default): one panel per level.}
#'     \item{\code{plotord}}{The order of the factors for two or three way terms (e.g. \code{c(2, 1)} or
#'       \code{c("B", "A")}): the first gives the colours, the second (and third) the panels. Giving it draws the plot.}
#'     \item{\code{mtitle}}{The main title.}
#'     \item{\code{ci}}{Draw the confidence band (default \code{TRUE}).}
#'     \item{\code{point}}{Draw the data points (default \code{TRUE}).}
#'     \item{\code{jitterv}}{Jitter of the points in x and y (default 0).}
#'   }
#'   Before 2.0 these were arguments of their own; given so they stop with this spelling.
#' @param ... Not used: arguments removed in 2.0 (\code{level}, \code{Df}, \code{covariateV}, \code{transOff},
#'   \code{responsen}, \code{newwd}, the plot settings) stop with their 2.0 spelling (see \code{vignette("migrating")}).
#' @param cov_reduce How covariates that are not varied are handled when building the predicted means: a function such as
#'   \code{mean} (the default) or \code{median}, or \code{"average"} to average the predictions over the
#'   observed covariate values (e.g. means adjusted for a spatial smooth in a \code{gam}). See \code{\link{Kmatrix}}.
#' @param exclude For \code{gam} models: smooth terms to leave out of the predictions. The default excludes random-effect
#'   smooths. See \code{\link{Kmatrix}}.
#' @param weights Passed to \code{\link{Kmatrix}}: how the cells of the factors that are not in \code{modelterm} are weighted when averaging over them:
#'   \code{"equal"} (the default; equal weights), \code{"proportional"} (weights proportional to the frequencies of
#'   those cells in the data, the same for every mean, as SAS \code{LSMEANS / OM}), \code{"outer"} (products of each
#'   factor's own frequencies), \code{"cells"} (each mean weights the cells by their own frequencies) or
#'   \code{"flat"} (equal weights over the non-empty cells). The definitions are those of \pkg{emmeans}.
#' @param df_method,vcov_adjust As in \code{\link{predictmeans}} (a number in \code{df_method} is a fixed df, \code{Df} before 2.0): by default Kenward-Roger df (each point of the curves
#'   has its own) and adjusted standard errors for models fitted by REML (\code{lmer}, \code{lme}, \code{gls}).
#'   predictmeans 1.x used the unadjusted covariance and, for \code{lmer}, one df per term; \code{vcov_adjust = FALSE}
#'   and, for \code{lme} and \code{gls}, \code{df_method = "containment"} give its results.
#' @param contr Two levels of \code{modelterm}, by their names as in the data (e.g.
#'   \code{c("Golden Rain", "Victory")}) or by their numbers (\code{c(1, 3)}; for a term with several factors the combined labels such as
#'   \code{"a:x"}): the difference between the first and the second along the covariate, with its confidence
#'   band, is returned and plotted instead of the lines. See \sQuote{Comparing two lines}.
#' @param contr_scale With \code{contr}: what is plotted. \code{"model"} (the default) is the difference on
#'   the scale of the model, where the test is made; for a log link or a log response it is shown as the ratio
#'   of the two lines and for a logit link as the odds ratio (exact: exp of the difference and of its limits).
#'   \code{"response"} is the difference of the back-transformed lines, in the units of the response, with a
#'   delta-method band; it describes the size of the difference, the test stays on the model scale.
#' @param band With \code{contr}: \code{"pointwise"} (the default) gives limits at the confidence level \code{conf.level} for
#'   each covariate value on its own; \code{"simultaneous"} gives Scheffe limits that hold for all covariate
#'   values together, for saying over which range of the covariate the lines differ.
#' @param mean_type Only for \code{glmmTMB} fits with zero inflation (\code{ziformula}; ignored, with a message,
#'   for other models): \code{"conditional"} (default), curves of the means of the count part, \eqn{\mu};
#'   \code{"overall"}, curves of the overall means \eqn{(1 - \pi)\mu} (\eqn{\pi} the probability of an extra zero),
#'   their bands \eqn{\exp(\log m \mp crit\,SE)} with the SE of \eqn{\log m} from the joint covariance of the
#'   count and zero-inflation coefficients (delta method), and with \code{contr} the ratio and the difference of the
#'   overall means. Log-link count families only.
#' @param plot Draw the plot now? The default is \code{FALSE} (since 2.0; 1.x drew it): the plot is in the result
#'   (\code{autoplot()} returns it). \code{plot = TRUE}, or \code{options(predictmeans.plot = TRUE)} for a whole
#'   session, draws it as 1.x did; \code{plot_args = list(plotord = )} also draws it.
#' @param bias_adjust,sigma For back-transformed lines (a link, a transformed response or \code{trans}), as in
#'   \code{\link{predictmeans}}. By default a line is the back-transformed fitted value: the median of the response
#'   for a transformed response (the geometric mean for a log), and the response of a typical cluster (random effects
#'   at zero) for a generalized linear mixed model. \code{bias_adjust = TRUE} gives the mean of the response
#'   instead (for a mixed model: averaged over the clusters), computed over the variation on the model scale:
#'   the residual variance for a transformed response plus the variances of the random intercepts, or
#'   \code{sigma^2} when \code{sigma} is given (needed with random slopes). It is exact (for a log the factor
#'   \eqn{exp(\sigma^2/2)}, for a square root \eqn{+ \sigma^2}); \code{"taylor"} is the second-order formula of
#'   emmeans. Use it for totals and expected amounts; ratios and tests do not need it. The limits treat the
#'   variance as known. With \code{contr} the response-scale difference is that of the adjusted lines.
#' @section Comparing two lines:
#' Whether the bands of two lines overlap does not say whether the lines differ: the lines are estimated from
#' the same model and are correlated. The band of their difference does, and it is computed from the contrast
#' of the two lines, not from their limits.
#'
#' \emph{The scale of the test.} The difference is tested on the scale of the model (the link, or the
#' transformed response): there the model is additive, the estimates are closest to normal and the standard
#' error does not depend on the mean. \eqn{d = (K_a - K_b)\hat\beta}, \eqn{SE^2 = (K_a - K_b) V (K_a - K_b)'},
#' limits \eqn{d \pm c\,SE} (columns \code{diff}, \code{ses}, \code{Df}, \code{LL}, \code{UL}). With a log
#' the difference back-transforms exactly to a ratio and with a logit to an odds ratio (columns \code{ratio}
#' or \code{odds.ratio} with their limits): the limits leave out 1 exactly when those of \code{diff} leave out
#' 0, so nothing is approximated. For other transformations (square root, probit, a power) the
#' back-transformed difference has no meaning and the model scale is shown.
#'
#' \emph{Why back-transform.} Only to read the size of the difference in the units of the response; it is not
#' needed to decide whether the lines differ. The response-scale difference \eqn{h(\eta_a) - h(\eta_b)} is
#' another quantity: its standard error is an approximation (the delta method, gradient
#' \eqn{h'(\eta_a) K_a - h'(\eta_b) K_b}; columns \code{response.diff}, \code{response.ses},
#' \code{response.LL}, \code{response.UL}), it depends on the values of the other variables of the model, and
#' lines that are parallel on the model scale are not parallel on the response scale. Its band can disagree
#' with the test near the boundary, so with \code{contr_scale = "response"} the covariate values where the
#' lines differ are still those of the model scale (marked along the axis). For a generalized linear mixed
#' model it is the difference for a typical cluster (random effects at zero), as the lines are.
#'
#' \emph{Significance at a level.} Three questions have three answers:
#' \enumerate{
#' \item At one covariate value chosen beforehand: the pointwise band (\code{band = "pointwise"}); the lines
#'   differ at the level \code{1 - conf.level} when it leaves out zero (one for a ratio).
#' \item Over which range of the covariate: looking at many values needs wider limits. \code{band =
#'   "simultaneous"} uses Scheffe's critical value \eqn{c = \sqrt{r F_{1-\alpha;\,r,\,df}}}, \eqn{r} the rank
#'   of the contrasts along the covariate (2 for straight lines: the Working-Hotelling band). It holds for
#'   every covariate value at once and is slightly conservative over a limited range.
#' \item Whether the lines differ anywhere: the joint F test that all the differences along the covariate
#'   are zero (\code{test} in the result).
#' }
#' Each set of limits in \code{pltdf} has a column that says whether it includes the value of no difference:
#' \code{includes0} for \code{diff}, \code{ratio.includes1} or \code{odds.ratio.includes1} (always equal to
#' \code{includes0}) and \code{response.includes0} for the response-scale difference (which can differ from
#' \code{includes0} near the boundary). Column \code{significant} (\code{!includes0}) and \code{significant} in
#' the result (the ranges of the covariate) follow the band asked for, on the model scale.
#' @section Zero inflation (glmmTMB):
#' For a \code{glmmTMB} fit with zero inflation (\code{ziformula}) the means are those of the conditional (count)
#' part, \eqn{\mu}, as \code{emmeans} and \code{predict(type = "conditional")} give, not the overall mean
#' \eqn{(1 - \pi)\mu} (\eqn{\pi} the probability of an extra zero); ratios are of these means, and equal the
#' ratios of overall means only when \eqn{\pi} is the same. A message says so and \code{settings$component} is
#' \code{"conditional"}. The two can order treatments differently when \eqn{\pi} differs between them.
#' \code{mean_type = "overall"} (\code{predictmeans()}, \code{covariatemeans()}, \code{semipred()}) gives the overall
#' means and their ratios.
#' @return
#' An object of class \code{c("covariatemeans", "pm_curves")}, a list. Its tables have the layout of
#' \code{\link{predictmeans}} (columns \code{estimate}, \code{std.error}, \code{df}, \code{conf.low},
#' \code{conf.high}, ..., the covariate under its own name; see \code{\link{tidy.pm_curves}}); \code{tidy()} and
#' \code{as.data.frame()} return them, \code{autoplot()} the plot, and \code{print()} shows a summary.
#' \item{means}{The curves: one row per level of \code{modelterm} and covariate value (\code{NULL} with
#'   \code{contr}).}
#' \item{pairs}{With \code{contr}: the difference of the two curves at each covariate value.}
#' \item{settings}{The options that produced them (df method, level, band, back-transformation, ...).}
#' \item{plt}{The plot (a \code{ggplot}); drawn only with \code{plot = TRUE}.}
#' \item{pltdf}{The data of the plot as in 1.x: the fitted values with their limits, or with \code{contr} the
#'   differences on the model scale and, when there is a back-transformation, on the response scale.}
#' \item{test, significant}{With \code{contr}: the joint F test that the two lines are the same (numerator and
#'   denominator df, F, p-value), and the ranges of the covariate where the band leaves out zero.}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   data(Oats, package="nlme")
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # library(lme4)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   cm <- covariatemeans(fm, "Variety", covariate="nitro")
#'   ggplot2::autoplot(cm)
#'   covariatemeans(fm, "Variety", covariate="nitro", at=seq(0, 0.6, 0.3))
#' # the difference between two varieties along nitro, with a simultaneous band
#' # (the two levels by name, or by number: contr=c(1, 3))
#'   cd <- covariatemeans(fm, "Variety", covariate="nitro", contr=c("Golden Rain", "Victory"),
#'                        band="simultaneous")
#'   cd$settings$test
#'   head(as.data.frame(cd))
#' @name covariatemeans
NULL

covariatemeans <- function (model, modelterm=NULL, covariate, as.is=FALSE, at=NULL, data=NULL, conf.level=0.95, trans=NULL,
                            cov_reduce=mean, exclude=NULL,
                            weights=c("equal", "proportional", "outer", "cells", "flat"),
                            df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
                            vcov_adjust=TRUE, contr=NULL, contr_scale=c("model", "response"),
                            band=c("pointwise", "simultaneous"), bias_adjust=FALSE, sigma=NULL,
                            mean_type=c("conditional", "overall"), plot=getOption("predictmeans.plot", FALSE),
                            plot_args=list(), ...) {
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("covariatemeans", match.call(), strict = TRUE, formals_of = formals(sys.function()))
  level <- pm_conf_alpha(conf.level)
  vcov_adjust_default <- missing(vcov_adjust)
  dfa <- pm_df_arg(df_method, eval(formals(sys.function())$df_method))
  df_method <- dfa$method
  Df <- dfa$Df
  pa <- pm_plot_args(plot_args, pm_cov_plot_defaults)
  trellis <- pa$trellis; plotord <- pa$plotord; mtitle <- pa$mtitle
  ci <- pa$ci; point <- pa$point; jitterv <- pa$jitterv
  if (inherits(model, "aovlist")) model <- pm_aov_carrier(model)     # by strata (R/aovlist.R)
  mean_type <- pm_zi_check(model, match.arg(mean_type))               # R/zero-inflation.R
  if (identical(mean_type, "overall") && !isFALSE(bias_adjust)) {
    stop("mean_type = \"overall\" cannot be combined with 'bias_adjust'.", call. = FALSE)
  }
  pm_zi_message(model, mean_type)

  if (is.null(modelterm) || all(modelterm%in%c("NULL", ""))) {
    modelterm <- covariate
    trellis=FALSE
  }

  if (!is.null(at) && length(covariate) > 1){
      stopifnot(is.matrix(at), dim(at)[2] == length(covariate))
      colnames(at) <- covariate
  }

  vars <- unlist(strsplit(modelterm, "\\:"))
  if (!is.null(plotord) && !unique(plotord%in%c("NULL", ""))) {
    if(length(plotord) != length(vars)){
      stop(paste("plotord must be a character vector of length", length(vars)))
    }
    if(is.character(plotord)){
      i = match(plotord, vars)
      if(any(is.na(i))){
        stop(paste("The entries", paste0(plotord[is.na(i)], sep = ", "), "do not match any term in ", modelterm))
      }else{
        plotord = i
      }
    }
  }
  # draw the plot now? (as predictmeans(): only when asked; a plot order
  # asks for it); the plot is in the result either way
  draw <- isTRUE(plot) || (!is.null(plotord) && !all(plotord %in% c("NULL", "")))
  weights <- match.arg(weights)
  ctr.matrix <- Kmatrix(model, modelterm, covariate, at, data, cov_reduce=cov_reduce, exclude=exclude, weights=weights)
  KK <- ctr.matrix$K
  pltdf <- ctr.matrix$fctnames
  response <- ctr.matrix$response
  preddf <- ctr.matrix$preddf
  inf <- pm_inference(model, df_method, vcov_adjust, vcov_adjust_default)     # R/inference.R
  mp <- inf$mp
  if (mean_type == "overall") {
    # the curves of overall means (1 - pi) mu, as log m (R/zero-inflation.R);
    # df from the count model's columns
    bz <- glmmTMB::fixef(model)$zi
    b <- mp$coef[mp$estimable]
    Vf <- as.matrix(stats::vcov(model, full = TRUE))[seq_len(length(b) + length(bz)), seq_len(length(b) + length(bz))]
    a <- pm_zi_augment(KK[, mp$estimable, drop = FALSE], b, pm_zi_X(model, pltdf, names(bz)), bz, Vf)
    KK <- a$K
    mp$coef <- a$coef
    mp$vcov <- a$vcov
    mp$estimable <- rep(TRUE, length(a$coef))
    df0 <- inf$df
    inf$df <- function(L, terms) df0(as.matrix(L)[, seq_len(a$ncond), drop = FALSE], terms)
    inf$mp <- mp
  }
  bhat <- mp$coef

  # the difference between two levels along the covariate, with its own band
  if (!is.null(contr)) {
    return(pm_cov_difference(contr, KK[, mp$estimable, drop = FALSE], pltdf, vars, covariate, modelterm, mp, inf,
                             level, Df, model, trans, mtitle, match.arg(contr_scale), match.arg(band),
                             bias_adjust, sigma, draw = draw, mean_type = mean_type))
  }

  if (!setequal(modelterm, covariate) & as.is) {
    agg_formula <- formula(paste(covariate, "~", paste(vars, collapse="+")))
    min_df <- aggregate(agg_formula, preddf, min, na.rm=TRUE)
    names(min_df)[ncol(min_df)] <- "min_value"
    max_df <- aggregate(agg_formula, preddf, max, na.rm=TRUE)
    names(max_df)[ncol(max_df)] <- "max_value"
    pltdf <- merge(merge(pltdf, min_df, by=vars, sort=FALSE), max_df, by=vars, sort=FALSE)
    pltdf <- pltdf[pltdf[[covariate]] >= pltdf[["min_value"]] & pltdf[[covariate]] <= pltdf[["max_value"]], ]
    KK <- KK[as.numeric(rownames(pltdf)),]
  }

  # We'll work only with the non-NA elements of bhat
  KK <- KK[, mp$estimable, drop=FALSE]
  pltdf$yhat <- as.numeric(KK%*%bhat)
  pltdf$ses <- sqrt(base::diag(KK %*% tcrossprod(mp$vcov, KK)))

  Df_arg <- Df
  if (is.null(Df) || all(Df%in%c("NULL", ""))) {
    Df <- inf$df(KK, modelterm)          # each point its own df (Kenward-Roger, Satterthwaite)
    if (any(Df==0)) {
      stop("You need provide Df for this model!")
    }
  }

  trans_arg <- trans
  if (is.null(trans)) {
    trans <- pm_link_trans(model)   # NULL for an identity link
    # a transformed response in the formula (log(y), ...; R/back-transform.R)
    if (is.null(trans)) {
      resp_trans <- pm_response_trans(model)
      if (!is.null(resp_trans)) {
        trans <- resp_trans$fun
        message("The response is ", resp_trans$label, ": the fitted lines are back-transformed to the original scale ",
                "(trans = I for the model scale).")
      }
    }
  }

  Mean <- LL <- UL <- xvar <- factors <- bky <- NULL
  bt <- NULL
  if (is.null(trans)) {
    if (!isFALSE(bias_adjust)) pm_cov_back(model, NULL, bias_adjust, sigma)     # warns: nothing to adjust
    pltdf$Mean <- pltdf$yhat
    pltdf$LL <- pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses
    pltdf$UL <- pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses
  } else {
    # (1.x added half the SE of the fitted value inside exp() for a log:
    # not a bias adjustment, and the line then depended on its own precision)
    # bias_adjust: the mean of the back-transformed response, not its median
    # (or the typical cluster of a GLMM); the same code as predictmeans()
    bt <- pm_cov_back(model, trans_arg, bias_adjust, sigma)
    hb <- bt$h
    if (is.null(hb)) hb <- trans
    lim <- cbind(hb(pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses), hb(pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses))
    pltdf$Mean <- hb(pltdf$yhat)
    pltdf$LL <- pmin(lim[, 1], lim[, 2])
    pltdf$UL <- pmax(lim[, 1], lim[, 2])
  }

  # the stable table of the curves (R/curves.R), the rows of pltdf
  lab_cols <- c(setdiff(intersect(vars, names(pltdf)), covariate), intersect(covariate, names(pltdf)))
  means_tab <- pm_curve_means(pltdf[, lab_cols, drop = FALSE], pltdf$yhat, pltdf$ses, Df, qt(1 - level/2, df = Df),
                              h = bt$h, d = bt$d)

  # pltdf$yhat <- pltdf$ses <- NULL
  pltdf$yhat <- NULL

  if (setequal(modelterm, covariate)) {
    pltdf$factors <- factor(1)
  } else {
    pltdf$factors <- factor(do.call("paste", c(pltdf[, vars, drop=FALSE], sep=":")))
  }
  colnames(pltdf)[colnames(pltdf)==covariate] <- "xvar"

  # delete empty factor combinations
  if (!is.null(data)) {
    mdf <- data
  } else {
    mdf <- tryCatch(model.frame(model), error = function(e) ctr.matrix$data)
  }
  if (length(setdiff(response, names(mdf)))!=0) {
    mdf[,response] <- eval(parse(text=response), mdf)
  }
  mdf <- cbind(mdf, preddf[, !names(preddf)%in%names(mdf), drop=FALSE])

  ndf <- data.frame(table(mdf[, vars, drop = FALSE]))
  if (any(ndf$Freq==0)) {
    ndf0 <- ndf[ndf$Freq==0, , drop=FALSE]
    ndf0$factors <- factor(do.call("paste", c(ndf0[, vars, drop=FALSE], sep=":")))
    keep <- !pltdf$factors%in%ndf0$factors
    pltdf <- pltdf[keep, ]
    means_tab <- means_tab[keep, , drop = FALSE]
    rownames(means_tab) <- NULL
  }


  if (is.null(mtitle) || mtitle%in%c("NULL", "")) {
    mtitle <- paste("Fitted and observed relationship with", (1-level)*100, "% CI")
  }
  if (is.null(trans)) {
    mdf$bky <- mdf[, response]
  } else {
    if (length(setdiff(response, names(mdf)))==0) {    ## Transformed y before modelling
      if (!is.null(pm_link_trans(model))) {    # models with a link: data on the response scale
        if (inherits(mdf[, response], "factor")) {
          mdf$bky <- as.numeric(mdf[, response])-1
        }else if (!is.null(dim(mdf[, response]))) {
          mdf$bky <- mdf[, response][,1]/rowSums(mdf[, response])
          response <- "Probability"
        } else {
          mdf$bky <- mdf[, response]
        }
        if (any( identical(trans, function(x) x, ignore.environment=TRUE), identical(trans, I, ignore.environment=TRUE))) {
          if (inherits(model, "glm")) {
            mdf$bky <- ifelse(mdf$bky >=1 | mdf$bky <= 0, NA, model$family$linkfun(mdf$bky))
          }
          if (inherits(model, "glmerMod")) {
            mdf$bky <- ifelse(mdf$bky >=1 | mdf$bky <= 0, NA, slot(model, "resp")$family$linkfun(mdf$bky))
          }
          if (inherits(model, "glmmTMB")) {
            mdf$bky <- ifelse(mdf$bky >=1 | mdf$bky <= 0, NA, model$modelInfo$family$linkfun(mdf$bky))
          }
          response <- "Response"
        }
      } else {
        mdf$bky <- trans(mdf[, response])
        response <- paste("Transformed", response)
      }
    } else {       ## Transformed y within modelling
      response <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
      if (length(setdiff(response, names(mdf)))!=0) {
        # the variable of the response in the formula (e.g. log(y + 1))
        response <- intersect(all.vars(stats::formula(model)[[2]]), names(mdf))[1]
        if (is.na(response)) {
          stop("The plot needs the response variable in the data; use the 'means' table instead.", call. = FALSE)
        }
      }
      mdf$bky <- mdf[, response]
    }
  }

  if (all(is.na(mdf$bky))) {
    point <- FALSE
  }
  if (setequal(modelterm, covariate)) {
    mdf$factors <- factor(1)
  } else {
    mdf$factors <- do.call("paste", c(mdf[, vars, drop=FALSE], sep=":"))
  }
  names(mdf)[names(mdf)==covariate] <- "xvar"
  pltdf[, c("Mean", "LL", "UL")] <- lapply(pltdf[, c("Mean", "LL", "UL")], as.numeric)

  if (!trellis) {
    if (setequal(modelterm, covariate))  {
      plt <- ggplot(pltdf, aes(xvar, Mean))+
        labs(title=paste(pm_wrap_title(mtitle, 11), "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
        geom_line(linewidth=0.5)+
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.2, data=pltdf, stat="identity")
      }
      if (point) {
        plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      }
    }else{
      plt <- ggplot(pltdf, aes(xvar, Mean, colour=factors))+
        labs(title=paste(pm_wrap_title(mtitle, 11), "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
        geom_line(linewidth=0.5)+
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(aes(ymin = LL, ymax = UL, fill=factors), alpha = 0.2, data=pltdf, stat="identity")
      }
      if (point) {
        plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      }
      plt <- plt+guides(col = guide_legend(modelterm), fill=guide_legend(modelterm))
    }
    if (draw) print(plt)
  }else{
    if (length(vars)==1) {
      plt <- ggplot(pltdf, aes(xvar, Mean, colour=factors))+
        labs(title=paste(pm_wrap_title(mtitle, 11), "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
        geom_line(linewidth=0.5)+
        facet_wrap(~  factors)+
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(aes(ymin = LL, ymax = UL, fill=factors), alpha = 0.2, data=pltdf, stat="identity")
      }
      if (point) {
        plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      }
      plt <- plt+guides(col = guide_legend(modelterm), fill=guide_legend(modelterm))
      if (setequal(modelterm, covariate))  {
        plt <- plt+ theme(legend.position="none")
      }
      if (draw) print(plt)
    }
    if (length(vars)==2) {
      if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
        plotord <- 1:2
      }
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      plt <- ggplot(pltdf, aes(xvar, Mean, colour=factor(eval(parse(text = fact1)))))+
        labs(title=paste(pm_wrap_title(mtitle, 11), "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
        geom_line(linewidth=0.5)+
        facet_grid(eval(parse(text = paste("~",fact2, sep=""))))+
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(aes(ymin = LL, ymax = UL, fill=factor(eval(parse(text = fact1)))), alpha = 0.2, data=pltdf, stat="identity")
      }
      if (point) {
        plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      }
      plt <- plt+guides(col = guide_legend(fact1), fill=guide_legend(fact1))
      if (draw) print(plt)
    }
    if (length(vars)==3) {
      if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
        plotord <- 1:3
      }
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      fact3 <- (vars[plotord])[3]
      plt <- ggplot(pltdf, aes(xvar, Mean, colour=factor(eval(parse(text = fact1)))))+
        labs(title=paste(pm_wrap_title(mtitle, 11), "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
        geom_line(linewidth=0.5)+
        facet_grid(eval(parse(text = paste(fact2, "~",fact3, sep=""))))+
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(aes(ymin = LL, ymax = UL, fill=factor(eval(parse(text = fact1)))), alpha = 0.2, data=pltdf, stat="identity")
      }
      if (point) {
        plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      }
      plt <- plt+guides(col = guide_legend(fact1), fill=guide_legend(fact1))
      if (draw) print(plt)
    }
  }

  settings <- pm_curve_settings(model, modelterm, covariate, level, pm_cov_df_label(inf, Df_arg), mean_type,
                                pm_curve_back_label(model, trans_arg, bt), vcov_adjust = isTRUE(inf$vcov_adjust),
                                vcov = if (is.null(inf$vcov_label)) "model" else inf$vcov_label)
  pm_curve_result(means = means_tab, settings = settings, legacy = list(plt = plt, pltdf = pltdf))
}

## The difference between two levels of the term along the covariate, with
## its confidence band (covariatemeans(contr = )). KK: the rows of K of the
## grid 'pltdf' (the levels of 'vars' by the covariate values).
##
## The band of a difference is not made from the bands of the two lines: the
## lines come from the same model and are correlated. Two bands are right:
##  * model (link) scale: d = (Ka - Kb) b, SE^2 = (Ka - Kb) V (Ka - Kb)',
##    limits d +- t SE with the df of the contrast. This is the scale of the
##    test. Back-transformed it is a ratio for a log link (exp of d and of its
##    limits) and an odds ratio for a logit link. For other links the
##    back-transformed limits of d mean nothing.
##  * response scale: h(eta_a) - h(eta_b), h the inverse link or the back-
##    transformation. Not h of the limits of d, and not a difference of the
##    lines' limits: by the delta method, with gradient
##    g = h'(eta_a) Ka - h'(eta_b) Kb, SE^2 = g V g', limits +- t SE. In a
##    GLMM this is the difference for a typical cluster (random effects at
##    zero), as the lines of covariatemeans() are.
pm_cov_difference <- function(contr, KK, pltdf, vars, covariate, modelterm, mp, inf, level, Df, model, trans,
                              mtitle, contr_scale = "model", band = "pointwise", bias_adjust = FALSE,
                              sigma = NULL, draw = FALSE, mean_type = "conditional") {
  if (length(covariate) != 1 || setequal(modelterm, covariate)) {
    stop("'contr' compares two levels of 'modelterm' along one covariate.", call. = FALSE)
  }
  lab <- do.call(paste, c(lapply(pltdf[vars], as.character), sep = ":"))
  levs <- unique(lab)
  if (length(contr) != 2) stop("'contr' must name two levels of '", modelterm, "'.", call. = FALSE)
  k <- if (is.numeric(contr)) contr else match(contr, levs)
  if (anyNA(k) || any(!k %in% seq_along(levs)) || k[1] == k[2]) {
    stop("'contr' must be two different levels of '", modelterm, "': ", paste(levs, collapse = ", "), ".", call. = FALSE)
  }
  x <- pltdf[[covariate]]
  a <- which(lab == levs[k[1]])
  b <- which(lab == levs[k[2]])
  b <- b[match(x[a], x[b])]
  Ka <- KK[a, , drop = FALSE]
  Kb <- KK[b, , drop = FALSE]
  L <- Ka - Kb
  df <- if (is.null(Df) || all(Df %in% c("NULL", ""))) inf$df(L, modelterm) else rep_len(Df, nrow(L))
  # the lines differ anywhere? the joint F test of the differences (a basis of L)
  qL <- qr(t(L))
  r <- qL$rank
  test <- tryCatch(stats::setNames(inf$ftest(L[qL$pivot[seq_len(r)], , drop = FALSE], modelterm),
                                   c("numDF", "denDF", "F", "p.value")), error = function(e) NULL)
  # pointwise: t; simultaneous over the covariate: Scheffe with the rank of L
  # (straight lines: Working-Hotelling), from pm_adjust_crit() (R/p-adjust.R)
  crit <- if (band == "simultaneous") {
    pm_adjust_crit("scheffe", level, df, list(list(idx = seq_len(nrow(L)), ok = seq_len(nrow(L)),
                                                     m = nrow(L), k = 2, r = r)))
  } else stats::qt(1 - level/2, df = df)

  # the back-transformation: the inverse link, the response transformation
  # in the formula, or 'trans' (trans = I: the model scale only)
  bt <- pm_cov_back(model, trans, bias_adjust, sigma)
  res <- pm_diff_table(x[a], Ka, Kb, mp$coef, mp$vcov, df, crit, h = bt$h, d1 = bt$d, kind = bt$kind)
  plt <- pm_diff_plot(res, covariate, levs[k], modelterm, level, band, contr_scale, test, mtitle,
                      typical = inherits(model, c("merMod", "glmmTMB", "lme")) && !bt$adjusted, adjusted = bt$adjusted)
  if (draw) print(plt)
  # the stable table (R/curves.R): pointwise limits, and with a simultaneous
  # band the Scheffe p-values and limits
  simult <- band == "simultaneous"
  tt <- res$table$diff/res$table$ses
  p_adj <- if (simult) pm_adjust_p("scheffe", tt, 2*stats::pt(abs(tt), df, lower.tail = FALSE), df,
                                   list(list(idx = seq_along(tt), ok = which(is.finite(tt)), m = length(tt), k = 2, r = r)))
  pairs <- pm_diff_pairs(res$table, covariate, levs[k], level, kind = res$kind, back = res$back,
                         crit_s = if (simult) crit, p_adj = p_adj)
  out <- res$table
  names(out)[1] <- covariate
  settings <- pm_curve_settings(model, modelterm, covariate, level, pm_cov_df_label(inf, Df), mean_type,
                                pm_curve_back_label(model, trans, bt), contr = levs[k], band = band,
                                test = if (!is.null(test)) data.frame(df1 = test[["numDF"]], df2 = test[["denDF"]],
                                                                      F.ratio = test[["F"]], p.value = test[["p.value"]]),
                                vcov_adjust = isTRUE(inf$vcov_adjust),
                                vcov = if (is.null(inf$vcov_label)) "model" else inf$vcov_label)
  pm_curve_result(pairs = pairs, settings = settings,
                  legacy = list(plt = plt, pltdf = out, test = test, significant = res$significant))
}

## The plot settings of covariatemeans() ('plot_args').
pm_cov_plot_defaults <- list(trellis = TRUE, plotord = NULL, mtitle = NULL, ci = TRUE, point = TRUE, jitterv = 0)

## The df label of covariatemeans() settings: "user" when Df is given.
pm_cov_df_label <- function(inf, Df) {
  if (!is.null(Df) && !all(Df %in% c("NULL", ""))) "user" else inf$label
}

## The back-transformation of the fitted lines of a model: the inverse link,
## the response transformation in the formula, or 'trans' (I: none), with
## the mean of the back-transformed response when bias_adjust is asked
## (pm_bias_funs(); the variance of pm_bias_sigma2() or 'sigma').
## Returns list(h, d (NULL: no back-transformation), kind ("ratio",
## "odds.ratio" or "none": what exp of a difference is), link_scale, adjusted).
pm_cov_back <- function(model, trans, bias_adjust = FALSE, sigma = NULL) {
  ba <- isTRUE(bias_adjust) || identical(bias_adjust, "taylor")
  if (!ba && !isFALSE(bias_adjust)) stop("'bias_adjust' must be TRUE, FALSE or \"taylor\".", call. = FALSE)
  link <- pm_link_trans(model)
  if (is.null(trans)) trans <- link
  if (is.null(trans)) trans <- pm_response_trans(model)$fun
  if (is.null(trans) || identical(trans, I)) {
    if (ba) warning("'bias_adjust' applies to back-transformed lines (a link or a transformed response); ignored.",
                    call. = FALSE)
    return(list(h = NULL, d = NULL, kind = "none", link_scale = FALSE, adjusted = FALSE))
  }
  link_scale <- !is.null(link) && identical(trans, link)
  link_name <- if (link_scale) tryCatch(stats::family(model)$link, error = function(e) NA_character_) else NA_character_
  kind <- if (identical(trans, exp) || identical(link_name, "log")) "ratio" else
    if (identical(trans, stats::plogis) || identical(link_name, "logit")) "odds.ratio" else "none"
  bf <- pm_back_funs(trans, trans, model, link_scale)
  h <- bf$h
  d <- bf$d1
  if (ba) {
    s2 <- if (!is.null(sigma)) sigma^2 else if (inherits(model, "pm_aovlist")) {
      # aov with Error(): the stratum variance components (negative estimates as zero)
      sum(pmax(suppressMessages(varcomp(attr(model, "aovlist")))$vcov, 0))
    } else pm_bias_sigma2(model, link_scale)
    bias <- pm_bias_funs(bf, s2, if (identical(bias_adjust, "taylor")) "taylor" else "exact")
    h <- bias$h
    d <- bias$d
  }
  list(h = h, d = d, kind = kind, link_scale = link_scale, adjusted = ba)
}

## The difference between two curves a and b at the values x, the one
## calculation of covariatemeans(contr = ) and semipred(contr = ).
## Ka, Kb: the rows that give the curves from the coefficients (covariance V);
## df, crit: of each difference; h, d1: the back-transformation and its
## derivative (NULL: none); kind: "ratio" / "odds.ratio" when exp of the
## difference is one.
##  * model scale: d = (Ka - Kb) b, SE^2 = (Ka - Kb) V (Ka - Kb)', limits
##    d +- crit SE: the scale of the test. Exp of d and of its limits is a
##    ratio (log) or an odds ratio (logit).
##  * response scale: h(eta_a) - h(eta_b) with the delta-method SE from the
##    gradient h'(eta_a) Ka - h'(eta_b) Kb: not h of the limits of d, and not
##    a difference of the curves' limits (the curves are correlated).
## Each set of limits has a column that says whether it includes no
## difference; 'significant' is that of the model scale.
## Returns list(table, significant (the ranges of x), kind, back).
pm_diff_table <- function(x, Ka, Kb, coef, V, df, crit, h = NULL, d1 = NULL, kind = "none") {
  Ka <- as.matrix(Ka)
  Kb <- as.matrix(Kb)
  V <- as.matrix(V)
  coef <- as.numeric(coef)
  L <- Ka - Kb
  out <- data.frame(xvar = x, diff = as.numeric(L %*% coef), ses = sqrt(rowSums((L %*% V) * L)))
  out$Df <- df
  out$LL <- out$diff - crit*out$ses
  out$UL <- out$diff + crit*out$ses
  out$includes0 <- out$LL <= 0 & out$UL >= 0
  back <- !is.null(h)
  if (back) {
    if (kind %in% c("ratio", "odds.ratio")) {
      out[[kind]] <- exp(out$diff)
      out[[paste0(kind, ".LL")]] <- exp(out$LL)
      out[[paste0(kind, ".UL")]] <- exp(out$UL)
      out[[paste0(kind, ".includes1")]] <- out$includes0
    }
    eta_a <- as.numeric(Ka %*% coef)
    eta_b <- as.numeric(Kb %*% coef)
    G <- d1(eta_a) * Ka - d1(eta_b) * Kb
    out$response.diff <- as.numeric(h(eta_a)) - as.numeric(h(eta_b))
    out$response.ses <- sqrt(rowSums((G %*% V) * G))
    out$response.LL <- out$response.diff - crit*out$response.ses
    out$response.UL <- out$response.diff + crit*out$response.ses
    out$response.includes0 <- out$response.LL <= 0 & out$response.UL >= 0
  }
  out$significant <- !out$includes0
  runs <- rle(out$significant)
  end <- cumsum(runs$lengths)
  start <- end - runs$lengths + 1L
  sig <- data.frame(from = out$xvar[start[runs$values]], to = out$xvar[end[runs$values]],
                    direction = ifelse(out$diff[start[runs$values]] > 0, ">", "<"))
  list(table = out, significant = sig, kind = if (back) kind else "none", back = back)
}

## The plot of pm_diff_table(): contr_scale "model" draws the difference, as
## the ratio or odds ratio when there is one (log axis, line at 1);
## "response" the response-scale difference, with marks where the curves
## differ on the model scale.
pm_diff_plot <- function(res, covariate, levs, modelterm, level, band, contr_scale = "model", test = NULL,
                         mtitle = NULL, typical = FALSE, adjusted = FALSE, style = pm_curve_style()) {
  out <- res$table
  ratio <- if (res$kind %in% c("ratio", "odds.ratio")) res$kind else character(0)
  what <- if (res$back && contr_scale == "response") "response" else if (length(ratio)) ratio else "diff"
  pd <- switch(what,
    response = data.frame(y = out$response.diff, lo = out$response.LL, hi = out$response.UL),
    diff = data.frame(y = out$diff, lo = out$LL, hi = out$UL),
    data.frame(y = out[[what]], lo = out[[paste0(what, ".LL")]], hi = out[[paste0(what, ".UL")]]))
  pd$x <- out$xvar
  pair <- paste0(levs[1], if (what %in% c("ratio", "odds.ratio")) " / " else " - ", levs[2])
  ylab <- switch(what, response = paste0("Difference (", pair, "), response scale"),
                 diff = paste0("Difference (", pair, ")", if (res$back) ", model scale" else ""),
                 ratio = paste0("Ratio (", pair, ")"), odds.ratio = paste0("Odds ratio (", pair, ")"))
  if (is.null(mtitle) || mtitle %in% c("NULL", "")) {
    mtitle <- paste0(switch(what, ratio = "Ratio of '", odds.ratio = "Odds ratio of '", "Difference between '"),
                     modelterm, "' ", levs[1], if (what %in% c("ratio", "odds.ratio")) " to " else " and ", levs[2],
                     " along ", covariate, " with ", (1 - level)*100, "% ",
                     if (band == "simultaneous") "simultaneous " else "", "CI")
  }
  sub <- c(if (!is.null(test) && is.finite(test[["p.value"]])) paste0(
             "Lines equal: F(", test[["numDF"]], ", ", format(round(test[["denDF"]], 1)), ") = ",
             format(round(test[["F"]], 2)), ", p ", if (test[["p.value"]] < 0.001) "< 0.001" else
               paste("=", format(round(test[["p.value"]], 3)))),
           if (what == "response") paste0(
             "delta-method band", if (adjusted) " (bias-adjusted means)" else if (typical) " (random effects at zero)" else "",
             "; marks: where the lines differ on the model scale"),
           if (what != "response" && any(out$significant, na.rm = TRUE))
             paste0("marks: where the band leaves out ", if (what %in% c("ratio", "odds.ratio")) 1 else 0))
  ref <- if (what %in% c("ratio", "odds.ratio")) 1 else 0
  # the band, the reference line, then the difference on top (R/plot-style.R)
  plt <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), fill = style$diff_band, alpha = style$alpha,
                         colour = NA) +
    ggplot2::geom_hline(yintercept = ref, linetype = 2, colour = style$ref) +
    ggplot2::geom_line(colour = style$diff_line, linewidth = style$linewidth) +
    ggplot2::labs(title = paste0(pm_wrap_title(mtitle, 11), "\n"), x = covariate, y = ylab,
                  subtitle = if (length(sub)) pm_wrap_title(paste(sub, collapse = "; "), 11, rel = 1)) +
    ggplot2::theme_bw()
  if (ref == 1) plt <- plt + ggplot2::scale_y_log10()            # ratios: a log axis
  # marks under the curve where the band leaves out the reference (on the
  # model scale for response-scale plots)
  plt <- plt + pm_sig_bars(res$significant, pd$x, style)
  plt
}
