#' Predicted Curves and Surfaces of a Semiparametric Model
#'
#' @description
#' Predicted curves (or surfaces) of a semiparametric model along a covariate, with pointwise or simultaneous
#' bands; their slopes (\code{deriv}), or the difference between two curves (\code{contr}). The result holds the
#' tables in the layout of \code{\link{predictmeans}} (see \code{\link{tidy.pm_curves}}) and the plot, which
#' \code{autoplot()} returns.
#' @usage
#' semipred(model, modelterm=NULL, covariate, sm_term=NULL, contr=NULL, at=NULL,
#'          boundary=NULL, conf_level=0.95, trans=NULL, contr_scale=c("model", "response"),
#'          band=c("pointwise", "simultaneous"),
#'          unconditional=getOption("predictmeans.semi_unconditional"), deriv=0,
#'          plot_args=list(), mean_type=c("conditional", "overall"),
#'          bias_adjust=FALSE, sigma=NULL, plot=getOption("predictmeans.plot", FALSE), ...)
#' @param model A fit of \code{\link{semireg}} or \code{\link{semireg_tmb}} (\code{semireg} before 2.0).
#' @param modelterm Name (in "quotes") of the factor term whose curves are predicted, as in the model formula, e.g. "A" or "A:B"
#'   (any number of factors: every combination is computed; \code{autoplot(x, group = , panel = )} chooses the layout).
#'   With \code{NULL} or the name of the covariate, one curve without grouping. Factors that are not in \code{modelterm}, including \code{by} factors of the smooth terms, are averaged over with equal weights, in both the parametric and the smooth parts of the model.
#' @param covariate Name (in "quotes") of the covariate, or two names for a surface, e.g. \code{c("x1", "x2")}.
#' @param sm_term The smooths used in the prediction; the default \code{NULL} uses all of them. \code{print(model)}
#'   (or \code{tidy(model)}) lists the smooths of the fit with the values to use here. Each entry may be a
#'   term as written in the formula (\code{"s(pressure, by = QB)"}, all its pieces), a smooth as labelled
#'   (\code{"s(pressure):QB200"}, column \code{smooth}), a covariate (\code{"pressure"}: every smooth of it, column
#'   \code{covariate}) or a random-effect term (\code{"s_pressure_QB200"}, column \code{term}). Smooths of other covariates enter the curves as
#'   a constant; the advice \code{smooth_other_covariate} says when the prediction includes such smooths.
#'   \code{settings$sm_term} records the blocks used.
#' @param contr Two levels of \code{modelterm} (numbers or names, e.g. \code{c(4, 1)} or \code{c("d", "a")}): the difference of their curves (level 4 minus level 1), with limits, along \code{covariate}.
#' @param at The values of \code{covariate} to predict at (\code{covariateV} before 2.0): a numeric vector (one covariate) or a matrix with one named column per covariate (two covariates, e.g. for a \code{Ztps} smooth); each row of the matrix is one prediction point.
#' @param boundary For two covariates (a surface): the region where predictions are shown. A matrix or data frame of two
#'   columns (the covariates in the order of \code{covariate}), the vertices of a polygon, or of several rings separated
#'   by a row of \code{NA} (a point is inside when it is inside an odd number of them: separate pieces and holes), such
#'   as \code{\link{semiboundary}()} returns; or \code{TRUE} for \code{semiboundary()} of the model's data (the region
#'   the data cover, following their shape). \code{NULL} (default): the whole grid.
#' @param conf_level Confidence level of the two-sided limits, default 0.95 (\code{level = 0.05} before 2.0;
#'   before 2.0 the limits used \code{qnorm(1 - level)}, i.e. 90\%). For gaussian fits the limits are t limits
#'   with Satterthwaite df at each point (column \code{df}; see Details); for generalized fits t limits with the
#'   df of the design's strata (inner-outer, as for GLMMs), normal limits when the model has no random term
#'   besides the smooths.
#' @param trans A function that back-transforms the curves, e.g. \code{trans = exp}. Without it the curves are
#'   back-transformed by the inverse link of a generalized fit, or by the inverse of a response transformed in the
#'   formula (\code{log(y)}, \code{sqrt(y)}, \code{log(y + 1)}, ...), as in \code{\link{covariatemeans}} (one rule
#'   for both); \code{trans = I} keeps the model scale.
#' @param bias_adjust,sigma For back-transformed curves, as in \code{\link{covariatemeans}}: by default a curve is the
#'   back-transformed fitted value (the median for a transformed response, a typical cluster for a GLMM);
#'   \code{bias_adjust = TRUE} gives the mean of the response, over the residual variance (transformed response) and
#'   the variances of the random intercepts other than the smooths (the smooths' random effects are the curves
#'   themselves), or \code{sigma^2}.
#' @param mean_type Only for \code{glmmTMB} fits with zero inflation (\code{ziformula}; ignored, with a message,
#'   for other models): \code{"conditional"} (default), curves of the means of the count part, \eqn{\mu};
#'   \code{"overall"}, curves of the overall means \eqn{(1 - \pi)\mu} (\eqn{\pi} the probability of an extra zero),
#'   their bands \eqn{\exp(\log m \mp crit\,SE)} with the SE of \eqn{\log m} from the joint covariance of the
#'   count and zero-inflation coefficients (delta method), and with \code{contr} the ratio and the difference of the
#'   overall means. Log-link count families only.
#' @param plot_args Superseded by the arguments of \code{\link{autoplot.pm_curves}}, kept for 1.x scripts (before
#'   2.0 arguments of their own; given so they stop with this spelling). A named list: \code{trellis} (\code{TRUE}, default: one panel per level), \code{scales} (\code{"fixed"},
#'   default, \code{"free"}, \code{"free_x"}, \code{"free_y"}: the axes of the panels), \code{plotord} (the order
#'   of the two factors of a two-way term, e.g. \code{c(2, 1)}; giving it draws the plot), \code{ci} (draw the
#'   bands, \code{TRUE}), \code{point} (draw the data, \code{TRUE}), \code{jitterv} (jitter of the points, 0),
#'   \code{threeD} (a 3-D surface with \pkg{plotly}, \code{FALSE}); and colours and sizes, any of: \code{line} and \code{band} (one
#'   curve; dark and light blue), \code{diff_line} and \code{diff_band} (a difference of two curves, \code{contr};
#'   dark and light green), \code{ref} (the dashed reference line at 0, or 1 for ratios, and the marks under a
#'   difference where its band leaves out the reference; red), \code{alpha} (of the bands, 0.35), \code{linewidth}
#'   (1), \code{points}, \code{point_alpha}, \code{point_size} (the data: grey, 0.6, 1.3), \code{palette} (curves
#'   by group: eight colours that colour-blind readers can tell apart; ggplot2's own beyond eight groups),
#'   \code{surface_points} (the data on surfaces: red). Bands are drawn first and their lines on top. For
#'   example \code{plot_args = list(line = "black", band = "grey70")}; the plot is a ggplot, so anything else can
#'   be added to it.
#' @param plot Draw the plot now? The default is \code{FALSE} (since 2.0; 1.x drew it): the plot is in the result
#'   (\code{autoplot()} returns it). \code{plot = TRUE}, or \code{options(predictmeans.plot = TRUE)} for a whole
#'   session, draws it; so does \code{plot_args = list(plotord = )}.
#' @param ... Not used: arguments removed in 2.0 (\code{semireg}, \code{covariateV}, \code{level}, \code{prt},
#'   \code{simultaneous}, the plot settings) stop with their 2.0 spelling (see \code{vignette("migrating")}).
#' @param band \code{"pointwise"} (the default) or \code{"simultaneous"}: simultaneous confidence bands over the
#'   whole curve(s) instead of pointwise limits; the critical value is the \code{1 - level} quantile of the largest
#'   standardised deviation over the predicted points, from 10000 draws of the estimates' distribution (Ruppert,
#'   Wand and Carroll, 2003); for gaussian fits a multivariate t with the df of each point (Details). The same argument as in \code{\link{covariatemeans}}, where lines of a parametric
#'   model use Scheffe's critical value instead; for smooths that would be far too wide.
#' @param unconditional If \code{TRUE}, the covariance of the coefficients includes the uncertainty of the
#'   smoothing parameters (the variances of the smooth terms), as mgcv's \code{unconditional = TRUE}: the first-order
#'   correction \eqn{J V_\theta J'} is added, \eqn{J} the derivative of the estimated coefficients in the variance
#'   parameters and \eqn{V_\theta} the covariance of their estimates (from the Hessian of the deviance; for
#'   \code{semireg_tmb()} fits TMB's joint precision of all parameters). The default (\code{NULL}) is \code{TRUE}
#'   for generalized models and \code{FALSE} for gaussian ones: in a simulation of the difference between two curves
#'   (\code{dev/sim-curve-difference.R}, 150 data sets, nominal 0.95) the bands without it covered 0.93 pointwise and
#'   0.89-0.91 over the whole curve for Poisson and binary data, and 0.955-0.957 and 0.94-0.95 with it; for gaussian
#'   data 0.97 without it (0.98 with it). \code{options(predictmeans.semi_unconditional = TRUE)} (or \code{FALSE}) sets
#'   it for \code{semipred()} and \code{predictmeans()} on \code{semireg} fits.
#' @param deriv \code{0} (the default): the curves; \code{1}: their slopes, \code{2}: their curvatures (first and
#'   second derivatives along \code{covariate}, one covariate), with bands as for the curves, and with \code{contr}
#'   the difference of the slopes (curvatures) of two curves. See \sQuote{Slopes and curvature}.
#' @param contr_scale With \code{contr} and a generalized model: \code{"model"} (the default) plots the difference
#'   on the link scale, as the ratio (log link) or the odds ratio (logit link) of the two curves;
#'   \code{"response"} the difference of the back-transformed curves with a delta-method band. See
#'   \sQuote{Comparing two lines} in \code{\link{covariatemeans}}: the test is on the model scale.
#' @section Slopes and curvature:
#' With \code{deriv = 1} (\code{2}) each row of the prediction design is replaced by its first (second) derivative
#' along the covariate, by central differences of the whole design (the fixed part and every smooth basis, on the
#' same grid), so the slopes are those of the fitted curves for any basis; for the O'Sullivan basis (\code{bs = "os"},
#' 1.x \code{smZ()}) they equal the exact derivatives of the basis (\code{ZOSull(drv = )}),
#' to 1e-8. Standard errors, df and pointwise or simultaneous bands are computed as for the curves; generalized
#' models give the slopes of the linear predictor (link scale). \code{pred_df} holds the slopes in the column of the
#' response, with \code{includes0} and \code{significant}; the ranges of the covariate where the band lies above or
#' below zero (the curve rises or falls) are printed with the result and kept in \code{attr(pred_df, "significant")}. A covariate
#' scaled by \code{semireg(nscale = TRUE)} is differentiated in its own units. The table \code{slopes} of the
#' result has the slopes with \code{statistic} and \code{p.value} of slope = 0 (and \code{p.adj} with a
#' simultaneous band), and \code{significant} the ranges. In a simulation
#' (\file{dev/sim-semipred-deriv.R}; \eqn{y = \sin(x/1.5)} plus noise, 150 points, 500 data sets) the bands covered
#' the true slope 0.97 pointwise and over the whole range with the default 6 knots, 0.99-1.00 with 15; the
#' curvature 0.98-1.00.
#' @details
#' Degrees of freedom (gaussian fits, both engines): a point of a curve, or of a difference of curves, is a
#' combination of fixed and smooth (random) coefficients. Its df are Satterthwaite's for that combination (as
#' \code{\link{predictmeans}} on the fit and SAS for estimates that involve random effects), counting the
#' variance components of the design (blocks, plots, ...) and the residual; the smoothing parameters are held at
#' their estimates (their uncertainty is what \code{unconditional} adds), and components estimated at zero are left
#' out. So the df are the residual df when there is no other random term, and few where a design stratum dominates
#' (curves of treatments applied to whole plots). Pointwise limits use \code{qt(1 - level/2, df)} at each point.
#' Simultaneous bands use a multivariate t: the simulated deviations are divided by variance estimates that move
#' together over the points (each point with its own df), the points are put on one scale by their t
#' probabilities, and each point gets the t quantile of its df at the simulated \code{1 - level} quantile (with one
#' df for all points: the exact multivariate t).
#'
#' Generalized fits (both engines) take the df that \code{\link{predictmeans}} gives GLMMs: the inner-outer df of
#' the stratum of the design's random effects in which the fixed coefficients of a point vary (a curve of a
#' treatment applied to whole plots: the whole-plot df), the smooth blocks left out (they are not strata of the
#' design); with no random term besides the smooths, normal limits (as for a \code{glm} with a known dispersion).
#' In a simulation (\file{dev/sim-semireg-glmm-df.R}; Poisson counts in the split plot above, 1000 data sets per
#' setting) normal limits rejected 5.8-9.2\% of true equalities of two treatments at a time and covered a
#' treatment's curve 0.91-0.93, the t limits (10 df) 3.1-5.9\% and 0.955-0.970; bands of differences covered
#' 0.92-0.96 (normal) and 0.95-0.996 (t).
#'
#' In a simulation (\file{dev/sim-semipred-df.R}; split plot, 6 blocks, 3 treatments on whole plots, curves over 11
#' times; 500 data sets per setting) normal limits of a treatment's curve covered 0.91 pointwise and 0.92 over the
#' whole curve, these t limits (about 8-11 df) 0.94-0.95 and 0.96-0.97; for differences of curves (15-31 df) normal
#' 0.94-0.96, t 0.96-0.97. With curves by group and no design strata both covered 0.95-0.98.
#'
#' With \code{contr} the difference between two curves is computed as in \code{\link{covariatemeans}}, from the
#' contrast of the two curves and the covariance of the fixed and smooth coefficients. For a linear model
#' \code{pred_df} has the difference (in the column of the response), \code{LL}, \code{UL}, \code{SE_contr} and
#' the indicators \code{includes0} (the limits include zero) and \code{significant}. For a generalized model it
#' has the columns of \code{covariatemeans(contr = )}: \code{diff}, \code{ses}, \code{LL}, \code{UL} and
#' \code{includes0} on the link scale; the \code{ratio} or \code{odds.ratio} with its limits and
#' \code{..includes1}; \code{response.diff}, \code{response.ses}, \code{response.LL}, \code{response.UL} and
#' \code{response.includes0} on the response scale; and \code{significant} (the link scale). Before 2.0 the band of
#' a generalized model's difference used the standard error of the link-scale difference multiplied by a weight
#' evaluated at the difference of the means, which is not a standard error of any difference. The same difference
#' is in the table \code{pairs} in the layout of \code{predictmeans()}. The ranges of the
#' covariate where the curves differ are printed with the result and kept in \code{significant} and
#' \code{attr(pred_df, "significant")}; with
#' \code{band = "simultaneous"} they hold for the whole curve, and the curves differ somewhere at level \code{1 - conf_level}
#' when there is such a range. (A joint F test, as \code{covariatemeans()} gives, is not given: Wald tests of
#' penalised smooth coefficients are not reliable.)
#' @section Zero inflation (glmmTMB):
#' For a \code{glmmTMB} fit with zero inflation (\code{ziformula}) the means are those of the conditional (count)
#' part, \eqn{\mu}, as \code{emmeans} and \code{predict(type = "conditional")} give, not the overall mean
#' \eqn{(1 - \pi)\mu} (\eqn{\pi} the probability of an extra zero); ratios are of these means, and equal the
#' ratios of overall means only when \eqn{\pi} is the same. A message says so and \code{settings$component} is
#' \code{"conditional"}. The two can order treatments differently when \eqn{\pi} differs between them.
#' \code{mean_type = "overall"} (\code{predictmeans()}, \code{covariatemeans()}, \code{semipred()}) gives the overall
#' means and their ratios.
#' @return
#' An object of class \code{c("semipred", "pm_curves")}, a list. \code{tidy()} and \code{as.data.frame()} return its
#' table, \code{autoplot()} the plot, and \code{print()} shows a summary (see \code{\link{tidy.pm_curves}}); notes in \code{attr(, "advice")} (\code{\link{advice}}).
#' \item{means}{The curves: the factor columns of \code{modelterm}, the covariate(s), \code{estimate},
#'   \code{std.error}, \code{df}, \code{conf.low}, \code{conf.high} on the model scale (pointwise limits), with a
#'   simultaneous band \code{adj.conf.low}, \code{adj.conf.high}; for generalized fits and with \code{trans} also
#'   \code{response}, \code{response.std.error}, \code{response.conf.low}, \code{response.conf.high} (and
#'   \code{response.adj.conf.*}). \code{NULL} with \code{contr} or \code{deriv}.}
#' \item{slopes}{With \code{deriv = 1} (\code{2}): the slopes (curvatures) in the same columns, with
#'   \code{statistic}, \code{p.value} and, with a simultaneous band, \code{p.adj}.}
#' \item{pairs}{With \code{contr}: the difference of the two curves, as \code{covariatemeans()$pairs}; with a
#'   simultaneous band \code{p.adj} is the share of the simulated maxima at least as large (below 0.05 exactly where
#'   the band leaves out zero; the smallest value is 1/10000).}
#' \item{settings}{The options: \code{df_method}, \code{band}, \code{deriv}, \code{unconditional},
#'   \code{back_transform}, ...; with \code{contr}, \code{test}: whether the curves differ anywhere along the
#'   covariate, the max-|t| test (\code{statistic} the largest |t| over the \code{points}, \code{p.value} the share of
#'   simulated maxima at least as large, from the covariance of the fixed and smooth coefficients: the smallest level at
#'   which the simultaneous band leaves out zero somewhere). A Wald F test of the smooth coefficients is not used:
#'   they are shrunken, and the statistic would not have an F distribution. In a simulation of equal curves
#'   (200 data sets) it rejected 3.0\% at the 5\% level.}
#' \item{significant}{With \code{contr} or \code{deriv}: the ranges of the covariate where the band asked for
#'   leaves out zero.}
#' \item{plt}{\code{autoplot()} of the result with \code{plot_args}: a ggplot (one plotly widget for \code{threeD});
#'   drawn only when \code{plot = TRUE}.}
#' \item{pred_df}{The 1.x table: a data.frame with the predicted values, their limits \code{LL} and \code{UL} (the band asked for), \code{SE} and
#'   \code{df}; for generalized fits the means and limits are on the response scale (limits \code{h(eta -+ crit SE)}),
#'   \code{SE} on the link scale and \code{SE_mu} on the response scale (delta method), \code{df} when the limits
#'   are t; with \code{contr} the columns described in Details.}
#' @references
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003), \emph{Semiparametric Regression}, Cambridge University Press.
#' @seealso \code{\link{semireg}}, \code{\link{tidy.pm_curves}}, \code{\link{semiboundary}}, \code{\link{covariatemeans}}, \code{\link{advice}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' data(Dialyzer, package = "nlme")
#' library(ggplot2)
#' ggplot(Dialyzer, aes(x = pressure, y = rate, col = QB)) +
#'   geom_line() +
#'   facet_wrap(vars(Subject))
#'
#' # a curve per blood flow rate, a random intercept and slope per dialyzer
#' fm <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
#'               data = Dialyzer)
#' residplot(fm$semer, group = "QB")
#' anova(fm)   # fixed effects (Kenward-Roger) and a test of each smooth term
#'
#' # the curves with pointwise bands (t limits with Satterthwaite df), then simultaneous ones
#' out <- semipred(fm, "QB", "pressure")
#' out
#' autoplot(out)
#' \donttest{
#' autoplot(semipred(fm, "QB", "pressure", band = "simultaneous"), panel = character(0))
#' # the smooths of one term only, as written in the formula
#' semipred(fm, "QB", "pressure", sm_term = "s(pressure, by = QB)")$settings$sm_term
#' # at given pressures
#' semipred(fm, "QB", "pressure", at = c(0.5, 1, 1.5, 2))
#'
#' # the difference of the two curves, and where they differ
#' dif <- semipred(fm, "QB", "pressure", contr = c("300", "200"), band = "simultaneous")
#' dif$significant
#' head(as.data.frame(dif))
#'
#' # slopes of the curves, and where they differ from zero over the whole range
#' slp <- semipred(fm, "QB", "pressure", deriv = 1, band = "simultaneous")
#' slp$significant
#'
#' # a surface: copper grade over the plane and a smooth in depth (data(copper))
#' data(copper)
#' fcu <- semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6),
#'                data = copper)
#' autoplot(semipred(fcu, covariate = "zcoord", trans = exp))
#' autoplot(semipred(fcu, covariate = c("xcoord", "ycoord"), trans = exp))
#'
#' # a binary response: the probability scale, and the odds ratio of two curves
#' set.seed(4)
#' db <- expand.grid(x = seq(0, 1, length = 25), g = factor(c("a", "b")), cl = factor(1:12))
#' db$y <- rbinom(nrow(db), 1, plogis(-1 + 2*sin(pi*db$x) + 0.8*(db$g == "b")*db$x +
#'                                      rnorm(12, 0, 0.5)[db$cl]))
#' fbin <- semireg(y ~ s(x, by = g, bs = "os", k = 6) + (1 | cl), family = binomial, data = db)
#' autoplot(semipred(fbin, "g", "x"))
#' head(as.data.frame(semipred(fbin, "g", "x", contr = c("b", "a"))))
#' }
#' @name semipred
NULL

#' @export
semipred <- function(model, modelterm=NULL, covariate, sm_term=NULL, contr=NULL, at=NULL, boundary=NULL,
                     conf_level=0.95, trans=NULL, contr_scale=c("model", "response"), band=c("pointwise", "simultaneous"),
                     unconditional=getOption("predictmeans.semi_unconditional"), deriv=0, plot_args=list(),
                     mean_type=c("conditional", "overall"), bias_adjust=FALSE, sigma=NULL,
                     plot=getOption("predictmeans.plot", FALSE), ...) {
  # arguments removed in 2.0, conf_level, plot settings, zero inflation (one
  # preamble with covariatemeans(), R/covariatemeans.R)
  if (!deriv %in% 0:2) stop("'deriv' must be 0, 1 or 2.", call. = FALSE)
  pre <- pm_curve_preamble("semipred", match.call(), formals(sys.function()), conf_level, plot_args, plot,
                           model, mean_type, deriv)
  semireg <- model
  covariateV <- at
  level <- pre$level; mean_type <- pre$mean_type; draw <- pre$draw
  if (deriv > 0) {
    if (length(covariate) != 1) stop("'deriv' needs one covariate.", call. = FALSE)
    if (!is.null(trans)) stop("'trans' is not used with 'deriv': derivatives are on the model scale.", call. = FALSE)
  }
  trans_arg <- trans
  contr_scale <- match.arg(contr_scale)
  band <- match.arg(band)
  simultaneous <- band == "simultaneous"

  stopifnot(inherits(semireg, "semireg"), !is.null(semireg$CovMat))

  if (is.null(modelterm) || all(modelterm%in%c("NULL", ""))) {
    modelterm <- covariate
  }
  vars_tab <- if (!setequal(modelterm, covariate)) unlist(strsplit(modelterm, ":"))
  # the plot settings as autoplot() arguments (also checks plotord)
  pargs <- pm_curve_plot_args(pre$pa, if (is.null(vars_tab)) character(0) else vars_tab,
                              if (!is.null(contr)) contr_scale)

  mod_df <- semireg$data
  if (isFALSE(boundary)) boundary <- NULL
  if (isTRUE(boundary)) {
    # the region covered by the data (semiboundary())
    if (length(covariate) != 2) stop("'boundary = TRUE' needs two covariates.", call. = FALSE)
    boundary <- semiboundary(mod_df[, covariate])
  }
  if (!is.null(contr)) {
    trans <- NULL
    if (length(unlist(strsplit(modelterm, "\\:"))) != 1 || setequal(modelterm, covariate)) {
      stop("'contr' compares two levels of one factor: 'modelterm' must be a single factor.", call. = FALSE)
    }
    if (length(contr) != 2) stop("'contr' must name two levels of '", modelterm, "'.", call. = FALSE)
    if (is.character(contr)) {
      grp_fct_levl <- levels(factor(mod_df[[modelterm]]))
      i <- match(contr, grp_fct_levl)
      if (anyNA(i)) stop("'contr' must be two levels of '", modelterm, "': ", paste(grp_fct_levl, collapse = ", "),
                         "; not ", paste(contr[is.na(i)], collapse = ", "), ".", call. = FALSE)
      contr <- i
    }
  }

  if (!is.null(covariateV) && length(covariate) > 1) {
    if (!(is.matrix(covariateV) && ncol(covariateV) == length(covariate))) {
      stop("With ", length(covariate), " covariates 'at' must be a matrix with one column per covariate ",
           "(each row a point).", call. = FALSE)
    }
    colnames(covariateV) <- covariate
  }

  sm_term <- pm_sm_resolve(semireg, sm_term)
  des <- pm_semi_design(semireg, modelterm, covariate, covariateV, sm_term, mod_df, unconditional = unconditional)
  semer <- des$semer; pred_df <- des$pred_df
  response <- des$response
  if (deriv > 0) des <- pm_semi_deriv(semireg, des, modelterm, covariate, sm_term, mod_df, deriv)
  D_mx <- des$D_mx; beta_u <- des$beta_u; sm_Cov <- des$sm_Cov; sig.epsHat <- des$sig.epsHat
  ncond <- ncol(D_mx)
  if (mean_type == "overall") {
    # overall means (1 - pi) mu of a zero-inflated fit: the rows of log m,
    # with the joint covariance of the count (fixed + smooth) and
    # zero-inflation coefficients from TMB's joint precision (unconditional:
    # its whole inverse) (R/zero-inflation.R)
    nzc <- ncol(lme4::getME(semer, "Z"))
    nb <- length(glmmTMB::fixef(semer)$cond)
    bz <- glmmTMB::fixef(semer)$zi
    Vj <- pm_tmb_cond_cov(semer, nzc, zi = TRUE, full = if (is.null(unconditional)) TRUE else isTRUE(unconditional))
    idx <- c(des$cov_idx, nb + nzc + seq_along(bz))
    a <- pm_zi_augment(D_mx, beta_u, pm_zi_X(semer, pred_df, names(bz)), bz, Vj[idx, idx, drop = FALSE])
    D_mx <- a$K; beta_u <- a$coef; sm_Cov <- a$vcov / sig.epsHat^2
  }

  # two-sided (1 - level) limits, pointwise or simultaneous over the curve(s),
  # with the df of pm_sr_ddf() (as predictmeans() on the fit): gaussian fits
  # Satterthwaite at each point, generalized fits the inner-outer df of the
  # design's strata (normal without strata)
  sr_terms <- if (!is.null(modelterm) && !setequal(modelterm, covariate)) modelterm
  df_label <- "normal (semireg)"
  df_of <- function(D) {
    D <- D[, seq_len(ncond), drop = FALSE]         # the count model's columns (mean_type = "overall")
    df <- if (is.null(des$cov_idx)) rep(Inf, nrow(D)) else pm_sr_ddf(semireg, D, des$cov_idx, sr_terms)
    if (!is.null(attr(df, "method"))) df_label <<- attr(df, "method")
    df <- as.numeric(df)
    df[is.na(df)] <- Inf
    df
  }
  # the critical values: pointwise t, or the simultaneous band simulated by
  # pm_sim_crit() (R/p-adjust.R), whose draws also give the adjusted p-values
  sim_of <- function(D, df) {
    if (simultaneous) pm_sim_crit(D, sm_Cov, level, df = df, draws = TRUE) else list(crit = stats::qt(1 - level/2, df))
  }
  # (with contr only the difference of two curves is reported: its own df
  # and critical value below)
  df_mx <- if (is.null(contr)) df_of(D_mx) else rep(Inf, nrow(D_mx))
  sim_m <- if (is.null(contr)) sim_of(D_mx, df_mx)
  crit <- if (is.null(contr)) sim_m$crit else stats::qnorm(1 - level/2)
  stderr <- sig.epsHat*sqrt(rowSums(as.matrix(D_mx %*% sm_Cov) * as.matrix(D_mx)))
  y_hat <- as.numeric(D_mx%*%beta_u)
  eta_hat <- y_hat                                 # the model scale, for the stable tables
  unc <- if (is.null(unconditional)) {
    if (inherits(semer, "glmmTMB")) semer$modelInfo$family$family != "gaussian" else !lme4::isLMM(semer)
  } else isTRUE(unconditional)
  LL <- y_hat - crit*stderr
  UL <- y_hat + crit*stderr
  glmm <- !isLMM(semer) && deriv == 0              # a generalized fit's curves: back-transformed by the inverse link
  fam <- if (!isLMM(semer)) (if (inherits(semer, "glmmTMB")) semer$modelInfo$family else slot(semer, "resp")$family)

  # the back-transformation of the curves: one rule with covariatemeans()
  # (pm_cov_back()): 'trans', else the inverse link, else the inverse of a
  # response transformed in the formula (log(y), sqrt(y), ...); bias_adjust
  # gives the mean of the response instead of the back-transformed fitted value
  # (the smooths' random effects are the curves, not variation to average over)
  bt <- if (deriv == 0) pm_cov_back(semer, trans, bias_adjust, sigma, exclude = names(semireg$smooths$maps)) else list(h = NULL, d = NULL, kind = "none",
                                                                             link_scale = FALSE, adjusted = FALSE)
  if (deriv > 0 && !isFALSE(bias_adjust)) warning("'bias_adjust' is not used with 'deriv' (model scale).", call. = FALSE)
  if (is.null(trans) && !is.null(bt$h) && !glmm) trans <- bt$h      # a transformed response: its inverse
  if (is.null(trans) || (glmm && isTRUE(bt$link_scale) && !isTRUE(bt$adjusted))) {
    if (glmm) {
      # the means and limits back-transformed with the inverse link, as
      # predictmeans() and covariatemeans() do: h(eta), h(eta -+ crit SE), and
      # the delta-method SE |h'(eta)| SE. (Before 2.0: limits h(eta) -+ crit
      # SE_mu, clipped to [0, 1] for binomial, with SE_mu from the weight
      # h'(eta)^2 / V(mu), which is |h'(eta)| only for canonical links.)
      stderr_mu <- abs(fam$mu.eta(y_hat))*stderr
      lim <- cbind(fam$linkinv(LL), fam$linkinv(UL))
      y_hat <- fam$linkinv(y_hat)
      LL <- pmin(lim[, 1], lim[, 2])
      UL <- pmax(lim[, 1], lim[, 2])
      pred_df$SE <- stderr                         # the link scale (the limits are h(eta -+ crit SE))
      pred_df$SE_mu <- stderr_mu
      if (any(is.finite(df_mx))) pred_df$df <- df_mx
    } else {
      pred_df$SE <- stderr
      pred_df$df <- df_mx
    }
  } else {
    h1 <- if (!is.null(bt$h)) bt$h else trans
    y_hat <- h1(y_hat)
    lim <- cbind(h1(LL), h1(UL))
    LL <- pmin(lim[, 1], lim[, 2])
    UL <- pmax(lim[, 1], lim[, 2])
    response <- "Trans_y"                          # the 1.x column of pred_df
  }
  # the back-transformation of the stable tables: 'trans', or the inverse
  # link of a generalized fit (curves only; slopes stay on the link scale)
  glm_back <- isTRUE(bt$link_scale)
  h_tab <- bt$h
  d_tab <- bt$d
  back_info <- if (!is.null(h_tab)) {
    if (isTRUE(bt$adjusted)) list(by = if (glm_back) "the inverse link" else if (!is.null(trans_arg)) "trans" else
                                    "the inverse of the response transformation",
                                  estimate = "mean (bias adjusted)", average = "link")
    else list(by = if (!is.null(trans_arg)) "trans" else if (glm_back) "the inverse link" else
                "the inverse of the response transformation",
              estimate = if (glm_back) (if (pm_sr_design_re(semireg)) "typical cluster (random effects at zero)" else
                                          "mean") else "median of the original scale (bias_adjust = TRUE gives the mean)",
              average = "link")
  }

  pred_df[[response]] <- y_hat
  pred_df$LL <- LL
  pred_df$UL <- UL
  for (i in names(pred_df)) {
    if (is(pred_df[, i], "AsIs")) pred_df[, i] <- as.numeric(pred_df[, i])
  }

  # the stable table of the curves or slopes (R/curves.R): model-scale
  # values, pointwise limits, the simultaneous band in adj.conf.*, the
  # covariates in their own units
  curve_tab <- pairs_tab <- significant <- contr_lev <- curve_test <- NULL
  inside <- function(d) if (length(covariate) == 2 && !is.null(boundary)) pm_points_in_poly(d[, covariate], boundary) else
    rep(TRUE, nrow(d))
  if (is.null(contr)) {
    curve_tab <- pm_curve_means(pred_df[, c(vars_tab, covariate), drop = FALSE], eta_hat, stderr, df_mx,
                                stats::qt(1 - level/2, df_mx), h = h_tab, d = d_tab,
                                crit_s = if (simultaneous) sim_m$crit, test = deriv > 0,
                                p_adj = if (simultaneous && deriv > 0) pm_sim_padj(sim_m, eta_hat/stderr, df_mx))
    if (length(covariate) == 2 && !is.null(boundary)) {
      curve_tab <- curve_tab[inside(curve_tab), , drop = FALSE]
      rownames(curve_tab) <- NULL
      pred_df <- droplevels(pred_df[inside(pred_df), , drop = FALSE])
    }
  } else {
    # the difference between the two curves (surfaces): the one calculation of
    # covariatemeans(contr = ) (pm_diff_table(), R/covariatemeans.R), with the
    # covariance of the fixed and smooth coefficients
    grp_fct <- pred_df[[modelterm]]
    grp_fct_levl <- levels(grp_fct)
    plt_df_lst <- lapply(split(pred_df, grp_fct), function(x) x[, !names(x) %in% modelterm])
    plt_df_contr <- plt_df_lst[[grp_fct_levl[contr[1]]]] - plt_df_lst[[grp_fct_levl[contr[2]]]]
    diff_covarn <- setdiff(names(plt_df_contr), c("SE", response, "LL", "UL"))
    plt_df_contr[, diff_covarn] <- plt_df_lst[[grp_fct_levl[contr[1]]]][, diff_covarn]
    D_mx_lst <- lapply(split(data.frame(as.matrix(D_mx)), grp_fct), as.matrix)
    D_a <- D_mx_lst[[grp_fct_levl[contr[1]]]]
    D_b <- D_mx_lst[[grp_fct_levl[contr[2]]]]
    df_c <- df_of(D_a - D_b)
    sim_c <- sim_of(D_a - D_b, df_c)
    crit_c <- sim_c$crit
    one <- length(covariate) == 1
    res <- pm_diff_table(if (one) plt_df_contr[[covariate]] else seq_len(nrow(D_a)), D_a, D_b, beta_u,
                         sig.epsHat^2*as.matrix(sm_Cov), df_c, crit_c,
                         h = bt$h, d1 = bt$d, kind = bt$kind)
    tab <- res$table
    pairs_tab <- pm_diff_pairs(tab, covariate, grp_fct_levl[contr], level, kind = res$kind, back = res$back,
                               crit_s = if (simultaneous) crit_c,
                               p_adj = if (simultaneous) pm_sim_padj(sim_c, tab$diff/tab$ses, df_c),
                               x = plt_df_contr[, covariate, drop = FALSE])
    contr_lev <- grp_fct_levl[contr]
    # do the curves differ anywhere? The max-|t| test over the covariate
    # values: the p-value is the smallest level at which the simultaneous band
    # of the difference leaves out zero somewhere (simulated with the
    # covariance of the fixed and smooth coefficients, as the band). A Wald F
    # test of the shrunken smooth coefficients would not have an F
    # distribution; this one needs none.
    sim_t <- if (simultaneous) sim_c else pm_sim_crit(D_a - D_b, sm_Cov, level, df = df_c, draws = TRUE)
    t_c <- tab$diff/tab$ses
    p_any <- suppressWarnings(min(pm_sim_padj(sim_t, t_c, df_c), na.rm = TRUE))
    curve_test <- if (is.finite(p_any)) data.frame(statistic = max(abs(t_c), na.rm = TRUE), method = "max-|t| (simulated)",
                                                     points = sum(is.finite(t_c)), p.value = p_any)
    if (one) {
      significant <- res$significant
    } else if (!is.null(boundary)) {
      pairs_tab <- pairs_tab[inside(pairs_tab), , drop = FALSE]
      rownames(pairs_tab) <- NULL
    }
    if (glmm) {
      # a GLMM: the columns of covariatemeans(contr = ); SE_contr on the model scale
      if (!one) {
        resp_sc <- contr_scale == "response"
        plt_df_contr[[response]] <- if (resp_sc) tab$response.diff else tab$diff
        plt_df_contr$LL <- if (resp_sc) tab$response.LL else tab$LL
        plt_df_contr$UL <- if (resp_sc) tab$response.UL else tab$UL
      }
      pred_df <- cbind(plt_df_contr[setdiff(diff_covarn, c("SE_mu", "significant"))],
                       tab[setdiff(names(tab), c("xvar", "Df"))])
      pred_df$SE_contr <- tab$ses
    } else {
      # a linear model: the 1.x columns, then the indicator columns
      plt_df_contr[[response]] <- tab$diff
      plt_df_contr$LL <- tab$LL
      plt_df_contr$UL <- tab$UL
      if (!one) plt_df_contr$significant <- tab$significant
      pred_df <- plt_df_contr
      pred_df$SE_contr <- tab$ses
      pred_df$df <- tab$Df
      pred_df$SE <- NULL
      pred_df$SE_mu <- NULL
      pred_df$includes0 <- tab$includes0
      pred_df$significant <- tab$significant
    }
    if (one) attr(pred_df, "significant") <- res$significant
  }

  if (deriv > 0 && is.null(contr)) {
    pred_df$includes0 <- pred_df$LL <= 0 & pred_df$UL >= 0
    pred_df$significant <- !pred_df$includes0
    sig <- pm_deriv_ranges(pred_df, vars_tab, covariate)
    attr(pred_df, "significant") <- sig
    significant <- if (is.null(sig)) data.frame(level = character(0), from = numeric(0), to = numeric(0),
                                                direction = character(0), stringsAsFactors = FALSE) else sig
  }

  rownames(pred_df) <- NULL
  settings <- pm_curve_settings(semireg, modelterm, covariate, level, df_label, mean_type, back_info,
                                contr = contr_lev, band = band, test = curve_test,
                                vcov = if (unc) "unconditional (smoothing-parameter uncertainty included)" else
                                  "conditional on the smoothing parameters",
                                deriv = deriv, unconditional = unc, boundary = !is.null(boundary),
                                sm_term = des$sm_term,
                                response = pm_response_label_info(semireg, trans_arg, glm_back, !is.null(h_tab), from_formula = is.null(trans_arg)),
                                contr_scale = if (!is.null(contr)) contr_scale)
  out <- pm_curve_result(means = if (deriv == 0) curve_tab, slopes = if (deriv > 0) curve_tab, pairs = pairs_tab,
                         settings = settings, legacy = list(plt = FALSE, pred_df = pred_df, significant = significant),
                         fun = "semipred")
  attr(out, "plot_data") <- pm_curve_plot_data(semireg, c(vars_tab, covariate), trans = trans_arg)
  out <- pm_add_advice(out, tryCatch(pm_semipred_advice(out, semireg), error = function(e) NULL),
                       pm_sm_covariate_advice(semireg, des$sm_term, covariate))
  pm_curve_add_plot(out, pargs, draw)
}

## The smooth blocks of 'sm_term' (NULL: all). Each entry may be a block name
## (names(fit$u_lst), e.g. "s_pressure_QB200"), a smooth as mgcv labels it
## ("s(pressure):QB200"), a term as written in the formula
## ("s(pressure, by = QB, k = 4)": all its blocks, matched by its covariates
## and its 'by' or factor term) or a covariate ("pressure": every smooth of
## it). Returns block names in the order of the fit.
pm_sm_resolve <- function(fit, sm_term = NULL) {
  if (is.null(sm_term)) return(NULL)
  sm <- fit$smooths
  if (is.null(sm)) return(sm_term)                 # (1.x smoothZ fits: pm_semi_design() stops)
  blk <- intersect(names(fit$u_lst), names(sm$maps))
  # per block: the mgcv labels, covariates and 'by' variables of its smooths
  info <- lapply(blk, function(b) {
    s <- sm$sms[vapply(sm$maps[[b]], `[[`, 1L, "sm")]
    by <- unique(vapply(s, `[[`, "", "by"))
    by <- unique(c(by, unlist(sm$derived[intersect(by, names(sm$derived))])))
    list(labels = vapply(s, `[[`, "", "label"), vars = unique(unlist(lapply(s, `[[`, "term"))), by = by,
         pm = any(startsWith(vapply(s, function(z) as.character(if (is.null(z$id)) "" else z$id), ""), "pm:")))
  })
  names(info) <- blk
  norm <- function(z) gsub("\\s+", "", z)
  as_term <- function(t) {
    e <- tryCatch(str2lang(t), error = function(err) NULL)
    if (!is.call(e) || !deparse(e[[1]]) %in% c("s", "te", "t2", "ti")) return(NULL)
    a <- as.list(e)[-1]
    nm <- if (is.null(names(a))) rep("", length(a)) else names(a)
    un <- lapply(a[nm == ""], all.vars)
    v <- unlist(un)
    fac <- if ("by" %in% nm) all.vars(a[["by"]]) else character(0)
    # s(x, f1*f2): the factor term is an unnamed argument; its variables are not covariates
    cov_vars <- unique(unlist(lapply(info, `[[`, "vars")))
    fac <- c(fac, setdiff(v, cov_vars))
    # derived factors (f1_f2 of s(x, f1*f2)) count when all their factors are given
    fac <- c(fac, names(Filter(function(d) all(d %in% fac), sm$derived)))
    list(vars = intersect(v, cov_vars), fac = fac)
  }
  pick <- lapply(sm_term, function(t) {
    if (t %in% blk) return(t)
    hit <- blk[vapply(info, function(i) norm(t) %in% norm(i$labels), logical(1))]
    if (length(hit)) return(hit)
    tm <- as_term(t)
    if (!is.null(tm) && length(tm$vars)) {
      # a factor term: the blocks of its factors, and the common curve of an
      # expanded s(x, f1*f2); none: the plain smooth only
      hit <- blk[vapply(info, function(i) setequal(i$vars, tm$vars) &&
                          if (length(tm$fac)) all(i$by %in% c(tm$fac, "NA")) && (any(i$by %in% tm$fac) || i$pm) else
                            all(i$by == "NA") && !i$pm, logical(1))]
      if (length(hit)) return(hit)
    }
    hit <- blk[vapply(info, function(i) t %in% i$vars, logical(1))]
    if (length(hit)) return(hit)
    stop("'sm_term' \"", t, "\" is not a smooth of this fit. Give smooths as labelled (",
         paste(unique(unlist(lapply(info, `[[`, "labels"))), collapse = ", "), "), terms as written in the formula, ",
         "block names (", paste(blk, collapse = ", "), ") or covariates (",
         paste(unique(unlist(lapply(info, `[[`, "vars"))), collapse = ", "), ").", call. = FALSE)
  })
  blk[blk %in% unlist(pick)]
}

## Advice: smooths in the prediction that do not involve 'covariate' (they add
## a constant, their average over the grid, to the curve).
pm_sm_covariate_advice <- function(fit, blocks, covariate) {
  sm <- fit$smooths
  if (is.null(sm) || !length(blocks)) return(NULL)
  other <- blocks[vapply(blocks, function(b) {
    v <- unique(unlist(lapply(sm$sms[vapply(sm$maps[[b]], `[[`, 1L, "sm")], `[[`, "term")))
    !any(covariate %in% v)
  }, logical(1))]
  if (!length(other)) return(NULL)
  pm_advice("smooth_other_covariate", paste0(
    "The smooth(s) ", paste(other, collapse = ", "), " do not involve ", paste(covariate, collapse = ", "),
    ": they enter the curves as a constant (at the reference values of their own covariates). ",
    "sm_term = \"", covariate[1], "\" uses the smooths of ", covariate[1], " only."), "note")
}

## Has a semireg fit random terms besides its smooths (blocks, plots, ...)?
pm_sr_design_re <- function(semireg) {
  semer <- semireg$semer
  cn <- tryCatch(if (inherits(semer, "glmmTMB")) semer$modelInfo$reTrms$cond$cnms else lme4::getME(semer, "cnms"),
                 error = function(e) NULL)
  length(setdiff(names(cn), names(semireg$refit$smoothZ))) > 0
}




## The design of semipred(): rows D = [K, Z] of the fixed part and the smooth
## basis on the reference grid of 'modelterm' (and covariate values), the
## estimates (beta, u) and their covariance up to sigma^2. Also used by
## predictmeans() for semireg fits.
pm_semi_design <- function(semireg, modelterm, covariate, covariateV, sm_term, mod_df, unconditional = NULL,
                           cov = TRUE) {
  semer <- semireg$semer 
  vcov_indN <- semireg$Cov_indN
  cov_lst <- semireg$cov_lst
  fullCovMat <- semireg$CovMat
  
  if (inherits(semer, "glmmTMB")) betaHat <- fixef(semer)$cond else betaHat <- fixef(semer)
  sig.epsHat <- sigma(semer) 
  
  # Kmatrix() understands 'semireg' objects directly: it uses the semireg
  # data and carries smooth-only covariates in its reference grid.
  KK <- Kmatrix(semireg, modelterm, covariate, covariateV)

  response <- response0 <- KK$response
  K <- KK$K[, colnames(model.matrix(semer)), drop=FALSE]
  pred_df <- KK$fctnames
  if (!(response %in% names(mod_df))) {   # all.vars(str2expression(response0))]
    mod_df[["response_y"]] <- eval(parse(text=response), mod_df)
    names(pred_df)[names(pred_df)==response] <- "response_y" 
    response <- "response_y"
  }
  
  # the rows of CovMat of the fixed effects and of the given smooth blocks, in
  # the order of the fit (that of u_lst and of Z); with unconditional = TRUE
  # the covariance includes the uncertainty of the smoothing parameters
  # (pm_sr_uncond(), R/semireg-engine.R; computed once per fit)
  ci <- semireg$Cov_ind
  # default: with the correction for generalized models (not gaussian ones),
  # as the coverage simulation showed (dev/sim-curve-difference.R)
  if (is.null(unconditional)) {
    unconditional <- if (inherits(semer, "glmmTMB")) semer$modelInfo$family$family != "gaussian" else !lme4::isLMM(semer)
  }
  cov_idx <- function(blocks) {
    c(seq_len(length(betaHat)), unlist(lapply(blocks, function(b) (ci[match(b, names(ci)) - 1] + 1):ci[match(b, names(ci))])))
  }
  cov_of <- function(blocks) {
    if (!cov) return(NULL)                     # (the design only: derivatives)
    idx <- cov_idx(blocks)
    V <- fullCovMat[idx, idx]
    if (isTRUE(unconditional)) {
      corr <- pm_sr_uncond_cached(semireg)
      if (!is.null(corr)) V <- as.matrix(V) + corr[idx, idx]
    }
    V
  }
  if (is.null(semireg$smooths)) {
    stop("This fit has its smooths in 'smoothZ' (smZ(), before 2.0): refit it with the smooths in the formula ",
         "(see vignette(\"migrating\")).", call. = FALSE)
  }
  # the blocks by name, their design for the grid from the stored mgcv smooth
  # objects (R/semireg-smooth.R); random curves (s(x) | g) are random
  # effects, not part of the population curves
  blk <- intersect(names(semireg$u_lst), names(semireg$smooths$maps))
  if (is.null(sm_term)) sm_term <- blk
  bad <- setdiff(sm_term, blk)
  if (length(bad)) stop("'sm_term' must be among the smooth blocks: ", paste(blk, collapse = ", "), call. = FALSE)
  blocks <- blk[blk %in% sm_term]
  sm_Cov <- cov_of(blocks)
  smooth_basis <- function(nd) do.call(cbind, pm_sm_predict(semireg$smooths, nd, blocks))
  ZN <- pm_refgrid_average(KK$refgrid, smooth_basis)
  list(semer = semer, KK = KK, K = K, ZN = ZN, mod_df = mod_df, pred_df = pred_df, response = response,
       response0 = response0, D_mx = cbind(K, ZN), beta_u = c(betaHat, do.call("c", unname(semireg$u_lst[blocks]))),
       sm_Cov = sm_Cov, sig.epsHat = sig.epsHat, cov_idx = cov_idx(blocks), sm_term = blocks)
}

## Derivatives of the curves of semipred(): the rows of the design D(x) on
## the grid replaced by dD/dx (deriv = 1) or d2D/dx2 (deriv = 2), by central
## differences of the whole design (fixed part and every smooth basis; for
## the O'Sullivan basis they equal ZOSull(drv = ) to 1e-8, tested), the same
## grid rebuilt at x + h and x - h, per unit of the covariate.
pm_semi_deriv <- function(semireg, des, modelterm, covariate, sm_term, mod_df, deriv) {
  xg <- unique(des$pred_df[[covariate]])
  r <- diff(range(xg))
  if (!is.finite(r) || r <= 0) stop("'deriv' needs more than one value of the covariate.", call. = FALSE)
  h <- pm_fd_h(0, r, order = deriv)            # R/deriv.R: 1e-5 or 1e-4 of the range
  at <- function(v) pm_semi_design(semireg, modelterm, covariate, v, sm_term, mod_df, cov = FALSE)
  dp <- at(xg + h)
  dm <- at(xg - h)
  same <- function(a) identical(a$pred_df[names(a$pred_df) != covariate], des$pred_df[names(des$pred_df) != covariate])
  if (!same(dp) || !same(dm)) stop("The grid could not be rebuilt for the derivative.", call. = FALSE)
  D <- if (deriv == 1) (as.matrix(dp$D_mx) - as.matrix(dm$D_mx))/(2*h) else
    (as.matrix(dp$D_mx) - 2*as.matrix(des$D_mx) + as.matrix(dm$D_mx))/h^2
  dimnames(D) <- dimnames(as.matrix(des$D_mx))
  des$D_mx <- D
  des
}

## The ranges of the covariate where a derivative band leaves out zero, per
## level of the factor (pred_df in the units of the data).
pm_deriv_ranges <- function(pred_df, by, covariate) {
  grp <- if (length(by)) interaction(pred_df[by], sep = ":", drop = TRUE) else factor(rep("", nrow(pred_df)))
  out <- lapply(levels(grp), function(g) {
    d <- pred_df[grp == g, , drop = FALSE]
    runs <- rle(d$significant)
    end <- cumsum(runs$lengths)
    start <- end - runs$lengths + 1L
    k <- runs$values
    if (!any(k)) return(NULL)
    data.frame(level = g, from = d[[covariate]][start[k]], to = d[[covariate]][end[k]],
               direction = ifelse(d$LL[start[k]] > 0, "> 0", "< 0"), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
