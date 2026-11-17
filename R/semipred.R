#' Predicted Means of a Semi Paramatric Model with Covariate Variable(s)
#'
#' @description
#' This function produces predicted means with graph for a semi paramatric model with new set of covariate values.
#' @usage
#' semipred(semireg, modelterm=NULL, covariate, sm_term=NULL, contr=NULL,
#'            covariateV=NULL, boundary=NULL, level=0.05, trans=NULL, trellis=TRUE, 
#' 		   scales=c("fixed", "free", "free_x", "free_y"),
#' 		   plotord=NULL, ci=TRUE, point=TRUE, jitterv=0, threeD=FALSE, prt=TRUE,
#' 		   simultaneous=FALSE, contr_scale=c("model", "response"),
#' 		   band=c("pointwise", "simultaneous"),
#' 		   unconditional=getOption("predictmeans.semi_unconditional"), deriv=0)
#' @param semireg A list object returned by \code{semireg}.
#' @param modelterm Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
#'   The \code{modelterm} must be given exactly as it appears in \code{semireg} model, e.g. "A" or "A:B". In case \code{modelterm} is the same as \code{covariate} or \code{NULL}, them \code{semipred} will product predictmeans with CI based on \code{covariate} only with out any grouping. Factors that are not in \code{modelterm}, including \code{by} factors of the smooth terms, are averaged over with equal weights, in both the parametric and the smooth parts of the model.
#' @param covariate Name (in "quotes") of one or two (for \code{Ztps} smooth) the covariate variables in the \code{semireg} shuch as "x1" or c("x1", "x2").
#' @param sm_term Names (in "quotes") of smooth terms (from \code{smoothZ} list in \code{semireg} model) used in the prediction such as "sm1_grp" or c("sm1_grp", "sm2_grp"). The default is using all smooth terms which is \code{sm_term=NULL}.
#' @param contr A numeric (or character) vector with length of two (e.g. c(4, 1) or c("d", "a")) which indicates to produce predicted mean with CI for difference between \code{modelterm} level 4 ("d") vs level 1 ("a") along \code{covariate}.
#' @param covariateV A numeric vector (one covariate) or a matrix with one named column per covariate (two covariates, e.g. for a \code{Ztps} smooth), then semipred will produce the result for \code{covariate} at value of \code{covariateV}. Each row of the matrix is one prediction point.
#' @param boundary For two covariates (a surface): the region where predictions are shown. A matrix or data frame of two
#'   columns (the covariates in the order of \code{covariate}), the vertices of a polygon, or of several rings separated
#'   by a row of \code{NA} (a point is inside when it is inside an odd number of them: separate pieces and holes), such
#'   as \code{\link{semiboundary}()} returns; or \code{TRUE} for \code{semiboundary()} of the model's data (the region
#'   the data cover, following their shape). \code{NULL} (default): the whole grid.
#' @param level Significance level of the two-sided confidence limits: the default 0.05 gives 95\% limits (before
#'   version 2.0 the limits used \code{qnorm(1 - level)}, i.e. 90\%). For gaussian fits the limits are t limits
#'   with Satterthwaite df at each point (column \code{df}; see Details); for generalized fits t limits with the
#'   df of the design's strata (inner-outer, as for GLMMs), normal limits when the model has no random term
#'   besides the smooths.
#' @param trans A function object for calculating the back transformed means, e.g. \code{trans=exp}.
#' @param trellis A logical scalar. If set to TRUE (default), a trellis plots of predicted means with CI will be drawn.
#' @param scales Should scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x", "free_y") in a trellis graph?
#' @param plotord A numeric vector specifying the order of plotting for two way interaction (e.g. 
#'   \code{plotord = c(2, 1)} will put the second variable in \code{modelterm} on the \code{X} 
#'   axis, the first variable as the grouping variable, and the third one as the panel variable). 
#'   The defaults are \code{c(1, 2)} for two way interactions.
#' @param ci A logical scalar to indicate whether to print confidence interval. The default value is TRUE.
#' @param point A logical scalar to indicate whether to print raw data points. The default value is TRUE.
#' @param jitterv A degree of jitter in x and y direction in the graph. The default is zero.
#' @param threeD A logical scalar to indicate whether to produce a 3-D plot or not. The default value is FALSE. 3-D plots need package \pkg{plotly}.
#' @param prt A logical scalar to indicate whether to produce plots on the screen or not. The default value is TRUE.
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
#' @param simultaneous \code{TRUE} is the same as \code{band = "simultaneous"} (the name before 2.0).
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
#' \code{smZ()}) they equal the exact derivatives of the basis, \code{ZOSull(drv = )} (as \code{smZ(drv = )}),
#' to 1e-8. Standard errors, df and pointwise or simultaneous bands are computed as for the curves; generalized
#' models give the slopes of the linear predictor (link scale). \code{pred_df} holds the slopes in the column of the
#' response, with \code{includes0} and \code{significant}; the ranges of the covariate where the band lies above or
#' below zero (the curve rises or falls) are printed and kept in \code{attr(pred_df, "significant")}. A covariate
#' scaled by \code{semireg(nscale = TRUE)} is differentiated in its own units. In a simulation
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
#' evaluated at the difference of the means, which is not a standard error of any difference. The ranges of the
#' covariate where the curves differ are printed and kept in \code{attr(pred_df, "significant")}; with
#' \code{simultaneous = TRUE} they hold for the whole curve, and the curves differ somewhere at level \code{level}
#' when there is such a range. (A joint F test, as \code{covariatemeans()} gives, is not given: Wald tests of
#' penalised smooth coefficients are not reliable.)
#' @return
#' \item{plt}{A ggplot object.}
#' \item{pred_df}{A data.frame with the predicted values, their limits \code{LL} and \code{UL} and (gaussian fits)
#'   \code{SE} and \code{df}; with \code{contr} the columns described in Details.}
#' @references
#' Ruppert, D., Wand, M. P. and Carroll, R. J. (2003), \emph{Semiparametric Regression}, Cambridge University Press.
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
#'               data = Dialyzer, prt = FALSE)
#' residplot(fm$semer, group = "QB")
#' anova(fm)   # fixed effects (Kenward-Roger) and a test of each smooth term
#'
#' # the curves with pointwise bands (t limits with Satterthwaite df), then simultaneous ones
#' out <- semipred(fm, "QB", "pressure")
#' head(out$pred_df)
#' semipred(fm, "QB", "pressure", band = "simultaneous")
#' # at given pressures
#' semipred(fm, "QB", "pressure", covariateV = c(0.5, 1, 1.5, 2), prt = FALSE)$pred_df
#'
#' # the difference of the two curves, and where they differ
#' dif <- semipred(fm, "QB", "pressure", contr = c("300", "200"), band = "simultaneous")
#' attr(dif$pred_df, "significant")
#'
#' # slopes of the curves, and where they differ from zero over the whole range
#' slp <- semipred(fm, "QB", "pressure", deriv = 1, band = "simultaneous")
#' attr(slp$pred_df, "significant")
#' \donttest{
#' # a surface: copper grade over the plane and a smooth in depth (data(copper))
#' data(copper)
#' fcu <- semireg(log(grade) ~ s(zcoord, bs = "os") + s(xcoord, ycoord, bs = "ztps", k = 6),
#'                data = copper, prt = FALSE)
#' semipred(fcu, covariate = "zcoord", trans = exp, prt = FALSE)$plt
#' semipred(fcu, covariate = c("xcoord", "ycoord"), trans = exp, prt = FALSE)$plt
#'
#' # a binary response: the probability scale, and the odds ratio of two curves
#' set.seed(4)
#' db <- expand.grid(x = seq(0, 1, length = 25), g = factor(c("a", "b")), cl = factor(1:12))
#' db$y <- rbinom(nrow(db), 1, plogis(-1 + 2*sin(pi*db$x) + 0.8*(db$g == "b")*db$x +
#'                                      rnorm(12, 0, 0.5)[db$cl]))
#' fbin <- semireg(y ~ s(x, by = g, bs = "os", k = 6) + (1 | cl), family = binomial, data = db,
#'                 prt = FALSE)
#' semipred(fbin, "g", "x", prt = FALSE)$plt
#' head(semipred(fbin, "g", "x", contr = c("b", "a"), prt = FALSE)$pred_df)
#' }
#' @name semipred
NULL

semipred <- function(semireg, modelterm=NULL, covariate, sm_term=NULL, contr=NULL,
                     covariateV=NULL, boundary=NULL, level=0.05, trans=NULL, trellis=TRUE, 
                     scales=c("fixed", "free", "free_x", "free_y"), 
                     plotord=NULL, ci=TRUE, point=TRUE, jitterv=0, threeD=FALSE, prt=TRUE, simultaneous=FALSE,
                     contr_scale=c("model", "response"), band=c("pointwise", "simultaneous"),
                     unconditional=getOption("predictmeans.semi_unconditional"), deriv=0) {
  if (!deriv %in% 0:2) stop("'deriv' must be 0, 1 or 2.", call. = FALSE)
  if (deriv > 0) {
    if (length(covariate) != 1) stop("'deriv' needs one covariate.", call. = FALSE)
    if (!is.null(trans)) stop("'trans' is not used with 'deriv': derivatives are on the model scale.", call. = FALSE)
    point <- FALSE
  }
  contr_scale <- match.arg(contr_scale)
  # band = "simultaneous" as in covariatemeans(); simultaneous = TRUE (1.x) still works
  band <- match.arg(band)
  if (isTRUE(simultaneous)) band <- "simultaneous"
  simultaneous <- band == "simultaneous"
  
  stopifnot(inherits(semireg, "semireg"), !is.null(semireg$CovMat))  
  scales <- as.character(scales)
  scales <- match.arg(scales)
  
  if (is.null(modelterm) || all(modelterm%in%c("NULL", ""))) {  
    modelterm <- covariate
  }
  
  mod_df <- semireg$data
  if (isFALSE(boundary)) boundary <- NULL
  if (isTRUE(boundary)) {
    # the region covered by the data (semiboundary()), in the covariates' units
    if (length(covariate) != 2) stop("'boundary = TRUE' needs two covariates.", call. = FALSE)
    bx <- mod_df[, covariate]
    for (i in covariate) {
      s <- attr(mod_df, paste(i, "sd", sep = "_"))
      m <- attr(mod_df, paste(i, "mean", sep = "_"))
      if (!is.null(s)) bx[[i]] <- bx[[i]] * s
      if (!is.null(m)) bx[[i]] <- bx[[i]] + m
    }
    boundary <- semiboundary(bx)
  }
  if (!is.null(contr)) {
    trans <- NULL
    stopifnot("'modelterm' must be a single term"=length(unlist(strsplit(modelterm, "\\:")))==1) 
    
    if(length(contr) != 2){
      stop("contr must be a character vector of length 2")
    }
    
    if(is.character(contr)){
      grp_fct <- mod_df[[modelterm]]
      grp_fct_levl <- levels(grp_fct)
      
      i = match(contr, grp_fct_levl)
      if(any(is.na(i))){
        stop(paste("The entries", paste0(contr[is.na(i)], sep = ", "), "do not match any level in ", modelterm))
      }else{
        contr = i
      }
    }
  }
  
  if (!is.null(covariateV)){
    if (length(covariate) > 1) {
      stopifnot(is.matrix(covariateV), dim(covariateV)[2] == length(covariate))
      colnames(covariateV) <- covariate
    }
    
    for (i in covariate) {
      if (!is.null(attr(mod_df, paste(i, "mean", sep="_")))) {
        if (!is.matrix(covariateV)) covariateV <- covariateV-attr(mod_df, paste(i, "mean", sep="_"))
        else covariateV[,i] <- covariateV[,i]-attr(mod_df, paste(i, "mean", sep="_"))
      }
      if (!is.null(attr(mod_df, paste(i, "sd", sep="_")))) {
        if (length(covariate) ==1) covariateV <- covariateV/attr(mod_df, paste(i, "sd", sep="_"))
        else covariateV[,i] <- covariateV[,i]/attr(mod_df, paste(i, "sd", sep="_"))
      } 
    }
  }
  
  des <- pm_semi_design(semireg, modelterm, covariate, covariateV, sm_term, mod_df, unconditional = unconditional)
  semer <- des$semer; mod_df <- des$mod_df; pred_df <- des$pred_df
  response <- des$response; response0 <- des$response0
  if (deriv > 0) des <- pm_semi_deriv(semireg, des, modelterm, covariate, sm_term, mod_df, deriv)
  D_mx <- des$D_mx; beta_u <- des$beta_u; sm_Cov <- des$sm_Cov; sig.epsHat <- des$sig.epsHat

  # two-sided (1 - level) limits, pointwise or simultaneous over the curve(s),
  # with the df of pm_sr_ddf() (as predictmeans() on the fit): gaussian fits
  # Satterthwaite at each point, generalized fits the inner-outer df of the
  # design's strata (normal without strata)
  sr_terms <- if (!is.null(modelterm) && !setequal(modelterm, covariate)) modelterm
  df_of <- function(D) {
    df <- if (is.null(des$cov_idx)) rep(Inf, nrow(D)) else pm_sr_ddf(semireg, D, des$cov_idx, sr_terms)
    df[is.na(df)] <- Inf
    as.numeric(df)
  }
  crit_of <- function(D, df) {
    if (simultaneous) pm_sim_crit(D, sm_Cov, level, df = df) else stats::qt(1 - level/2, df)
  }
  # (with contr only the difference of two curves is reported: its own df
  # and critical value below)
  df_mx <- if (is.null(contr)) df_of(D_mx) else rep(Inf, nrow(D_mx))
  crit <- if (is.null(contr)) crit_of(D_mx, df_mx) else stats::qnorm(1 - level/2)
  stderr <- sig.epsHat*sqrt(diag(D_mx%*%sm_Cov%*%t(D_mx)))
  y_hat <- as.numeric(D_mx%*%beta_u)
  LL <- y_hat - crit*stderr
  UL <- y_hat + crit*stderr
  
  if (is.null(trans)) {
    if (!isLMM(semer) && deriv == 0) {
      # the means and limits back-transformed with the inverse link, as
      # predictmeans() and covariatemeans() do: h(eta), h(eta -+ crit SE), and
      # the delta-method SE |h'(eta)| SE. (Before 2.0: limits h(eta) -+ crit
      # SE_mu, clipped to [0, 1] for binomial, with SE_mu from the weight
      # h'(eta)^2 / V(mu), which is |h'(eta)| only for canonical links.)
      fam <- if (inherits(semer, "glmmTMB")) semer$modelInfo$family else slot(semer, "resp")$family
      family_n <- fam$family
      stderr_mu <- abs(fam$mu.eta(y_hat))*stderr
      lim <- cbind(fam$linkinv(LL), fam$linkinv(UL))
      y_hat <- fam$linkinv(y_hat)
      LL <- pmin(lim[, 1], lim[, 2])
      UL <- pmax(lim[, 1], lim[, 2])
      if (family_n == "binomial"){
        if (is.factor(mod_df[[response]])) {
          mod_df[[response]] <- as.numeric(mod_df[[response]]) - 1
        }else if (!is.null(dim(mod_df[, response])) && ncol(mod_df[, response])==2) {
          mod_df[["Proportion"]] <- mod_df[, response][,1]/rowSums(mod_df[, response])
          response <- "Proportion"
        }		
      } 
      pred_df$SE_mu <- stderr_mu
      if (any(is.finite(df_mx))) pred_df$df <- df_mx
    }else {
      pred_df$SE <- stderr
      pred_df$df <- df_mx
    }
  }else{
    y_hat <- trans(y_hat)
    LL <- trans(LL)
    UL <- trans(UL)
    mod_df[["Trans_y"]] <- trans(mod_df[[response]])
    response <- "Trans_y"
  }
  
  pred_df[[response]] <- y_hat
  pred_df$LL <- LL
  pred_df$UL <- UL
  
  # Undo semireg()'s centring and scaling. Numeric covariates that are not in
  # the prediction grid (averaged over) are only in mod_df.
  if (!is.null(attr(mod_df, "numeric_var"))) {
    for (i in attr(mod_df, "numeric_var")) {
      in_pred <- i %in% names(pred_df)
      if (!is.null(attr(mod_df, paste(i, "sd", sep="_")))) {
        mod_df[[i]] <- mod_df[[i]]*attr(mod_df, paste(i, "sd", sep="_"))
        if (in_pred) pred_df[[i]] <- pred_df[[i]]*attr(mod_df, paste(i, "sd", sep="_"))
      }

      if (!is.null(attr(mod_df, paste(i, "mean", sep="_")))) {
        mod_df[[i]] <- mod_df[[i]]+attr(mod_df, paste(i, "mean", sep="_"))
        if (in_pred) pred_df[[i]] <- pred_df[[i]]+attr(mod_df, paste(i, "mean", sep="_"))
      }
    }
  }
  
  if (!is.null(attr(mod_df, paste(response, "sd", sep="_")))) {
    pred_df$LL <- pred_df$LL*attr(mod_df, paste(response, "sd", sep="_")) 
    pred_df$UL <- pred_df$UL*attr(mod_df, paste(response, "sd", sep="_")) 
    if (deriv > 0) {
      pred_df[[response]] <- pred_df[[response]]*attr(mod_df, paste(response, "sd", sep="_"))
      pred_df$SE <- pred_df$SE*attr(mod_df, paste(response, "sd", sep="_"))
    }
  }
  
  if (deriv == 0 && !is.null(attr(mod_df, paste(response, "mean", sep="_")))) {
    pred_df$LL <- pred_df$LL+attr(mod_df, paste(response, "mean", sep="_")) 
    pred_df$UL <- pred_df$UL+attr(mod_df, paste(response, "mean", sep="_")) 
  }
  
  for (i in names(pred_df)){
    if (is(pred_df[,i], "AsIs")) pred_df[,i] <- as.numeric(pred_df[,i])
  }
  
  if (setequal(modelterm, covariate)) {
    if (length(covariate)==1) {
      plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response))))+  
        geom_line(linewidth=1, col="red") +
        geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.6, stat="identity", fill="#00FFFF", col="#00FFFF")+
        labs(x=covariate, y=response)+
        theme_bw()
      if (point) plt <- plt + geom_point(data=mod_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response))), col="DeepPink", position = position_jitter(width = jitterv, height = jitterv))
      if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))
    }else if(length(covariate)==2) {
      if (!is.null(boundary)){
        inBdry <- pm_points_in_poly(pred_df[, covariate], boundary)
        pred_df <- droplevels(pred_df[inBdry,])
      }
      
      if (threeD) {
        pm_need("plotly", "3-D plots (threeD = TRUE)")
        plt <- plotly::plot_ly(pred_df, x = formula(paste("~", covariate[1])), 
                               y = formula(paste("~", covariate[2])),
                               z = formula(paste("~", response)), 
                               intensity = formula(paste("~", response)), 
                               type = 'mesh3d',
                               opacity = 0.7) 
        
        if (point) plt <- plotly::add_trace(plt, x = mod_df[[covariate[1]]], y = mod_df[[covariate[2]]], 
                                            z = mod_df[[response]], mode = "markers", type = "scatter3d", 
                                            marker = list(size = 2.5, color = "red", symbol = 104))						
      }else{		
        plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate[1])), y=eval(parse(text = covariate[2])), z=eval(parse(text = response)))) +
          geom_contour_filled() +
          labs(x=covariate[1], y=covariate[2], z=response)
        if (point) plt <- plt + geom_point(data=mod_df, col="red", size=1.5)   
        if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))		
      }
      
    }
  }else{
    vars <- unlist(strsplit(modelterm, "\\:"))
    if (length(vars)==1) {
      if (length(covariate)==1) {
        if (is.null(contr)) {
          plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = modelterm))))+
            geom_line(aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = modelterm))), linewidth=1) +
            labs(x=covariate, y=response, col=modelterm) +
            theme_bw()	
          if (ci) plt <- plt + 
              geom_ribbon(aes(ymin = LL, ymax = UL, fill=eval(parse(text = modelterm))), alpha = 0.3, stat="identity")+
              labs(fill=modelterm)
          if (point) plt <- plt + 
              geom_point(data=mod_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = modelterm))), position = position_jitter(width = jitterv, height = jitterv))+
              labs(col=modelterm)
          if (trellis)  plt <- plt + 
              facet_wrap(formula(paste("~", modelterm)), scales=scales) 
          if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))  
        }else{
          grp_fct <- pred_df[[modelterm]]
          grp_fct_levl <- levels(grp_fct)
          plt_df_lst <- lapply(split(pred_df, grp_fct), function(x) x[, !names(x) %in% modelterm])
          plt_df_contr <- plt_df_lst[[grp_fct_levl[contr[1]]]] - plt_df_lst[[grp_fct_levl[contr[2]]]]		  
          diff_covarn <- setdiff(names(plt_df_contr), c("SE", response, "LL", "UL"))
          plt_df_contr[,diff_covarn] <- plt_df_lst[[grp_fct_levl[contr[1]]]][, diff_covarn]
          
          # the difference between the two curves: the one calculation of
          # covariatemeans(contr = ) (pm_diff_table(), R/covariatemeans.R),
          # with the covariance of the fixed and smooth coefficients
          D_mx_lst <- lapply(split(data.frame(as.matrix(D_mx)), grp_fct), as.matrix)
          D_a <- D_mx_lst[[grp_fct_levl[contr[1]]]]
          D_b <- D_mx_lst[[grp_fct_levl[contr[2]]]]
          df_c <- df_of(D_a - D_b)
          crit_c <- crit_of(D_a - D_b, df_c)
          glmm <- !isLMM(semer) && deriv == 0
          fam <- if (glmm) (if (inherits(semer, "glmmTMB")) semer$modelInfo$family else slot(semer, "resp")$family)
          res <- pm_diff_table(plt_df_contr[[covariate]], D_a, D_b, beta_u, sig.epsHat^2*as.matrix(sm_Cov), df_c, crit_c,
                               h = if (glmm) fam$linkinv, d1 = if (glmm) fam$mu.eta,
                               kind = if (glmm) switch(fam$link, log = "ratio", logit = "odds.ratio", "none") else "none")
          # semireg() may have standardised the response of a linear model
          y_sd <- attr(mod_df, paste(response, "sd", sep = "_"))
          if (!glmm && !is.null(y_sd)) {
            for (cn in c("diff", "ses", "LL", "UL")) res$table[[cn]] <- res$table[[cn]]*y_sd
          }
          band_c <- if (simultaneous) "simultaneous" else "pointwise"
          plt <- pm_diff_plot(res, covariate, grp_fct_levl[contr], modelterm, level, band_c, contr_scale,
                              typical = glmm)
          if (response == "response_y") plt <- plt + labs(y = paste("The difference of ", response0, "\n", sep = ""))
          lv_c <- grp_fct_levl[contr]
          if (deriv > 0) lv_c <- c(paste("The", if (deriv == 1) "slope" else "curvature", "of", lv_c[1]),
                                   paste("that of", lv_c[2]))
          if (prt) pm_diff_report(res, covariate, lv_c, level, band_c)
          tab <- res$table
          if (glmm) {
            # a GLMM: the columns of covariatemeans(contr = ); SE_contr on the model scale
            pred_df <- cbind(plt_df_contr[setdiff(diff_covarn, "SE_mu")], tab[setdiff(names(tab), c("xvar", "Df"))])
            pred_df$SE_contr <- tab$ses
          } else {
            # a linear model: the 1.x columns, then the indicator columns
            plt_df_contr[[response]] <- tab$diff
            plt_df_contr$LL <- tab$LL
            plt_df_contr$UL <- tab$UL
            pred_df <- plt_df_contr
            pred_df$SE_contr <- tab$ses
            pred_df$df <- tab$Df
            pred_df$SE <- NULL
            pred_df$SE_mu <- NULL
            pred_df$includes0 <- tab$includes0
            pred_df$significant <- tab$significant
          }
          attr(pred_df, "significant") <- res$significant
        }
      }else if(length(covariate)==2) {
        if (is.null(contr)) {
          if (!is.null(boundary)){
            inBdry <- pm_points_in_poly(pred_df[, covariate], boundary)
            pred_df <- droplevels(pred_df[inBdry,])
          }
          
          if (threeD) {
            pm_need("plotly", "3-D plots (threeD = TRUE)")
            grp_fct <- pred_df[[modelterm]]
            plt <- lapply(split(pred_df, grp_fct), function(x) {
              plotly::plot_ly(droplevels(x), x = formula(paste("~", covariate[1])), 
                              y = formula(paste("~", covariate[2])),
                              z = formula(paste("~", response)), 
                              intensity = formula(paste("~", response)), 
                              type = 'mesh3d',
                              opacity = 0.7
              ) 
            })
            
            if (point) {
              mod_df_lst <- split(mod_df, mod_df[[modelterm]])
              for (i in names(mod_df_lst)) {
                plt[[i]] <- plotly::add_trace(plt[[i]], x = mod_df_lst[[i]][[covariate[1]]], y = mod_df_lst[[i]][[covariate[2]]], 
                                              z = mod_df_lst[[i]][[response]], mode = "markers", type = "scatter3d", 
                                              marker = list(size = 2.5, color = "red", symbol = 104)) 
                
              }
            }
          }else{		
            plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate[1])), y=eval(parse(text = covariate[2])), z=eval(parse(text = response)))) +
              geom_contour_filled()+
              facet_wrap(formula(paste("~", modelterm)), scales=scales)+
              labs(x=covariate[1], y=covariate[2], z=response) 
            if (point) plt <- plt + geom_point(data=mod_df, col="red", size=1.5)
            if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))			
          }          
        }else{
          grp_fct <- pred_df[[modelterm]]
          grp_fct_levl <- levels(grp_fct)
          plt_df_lst <- lapply(split(pred_df, grp_fct), function(x) x[, !names(x) %in% modelterm])
          plt_df_contr <- plt_df_lst[[grp_fct_levl[contr[1]]]] - plt_df_lst[[grp_fct_levl[contr[2]]]]
          diff_covarn <- setdiff(names(plt_df_contr), c("SE", response, "LL", "UL"))
          plt_df_contr[,diff_covarn] <- plt_df_lst[[grp_fct_levl[contr[1]]]][, diff_covarn]
          
          # the difference of the two surfaces: pm_diff_table() as for one
          # covariate (R/covariatemeans.R); the plot shows the model scale (or
          # the response scale, contr_scale = "response") and marks the points
          # where the model-scale band leaves out zero
          D_mx_lst <- lapply(split(data.frame(as.matrix(D_mx)), grp_fct), as.matrix)
          D_a <- D_mx_lst[[grp_fct_levl[contr[1]]]]
          D_b <- D_mx_lst[[grp_fct_levl[contr[2]]]]
          df_c <- df_of(D_a - D_b)
          crit_c <- crit_of(D_a - D_b, df_c)
          glmm <- !isLMM(semer)
          fam <- if (glmm) (if (inherits(semer, "glmmTMB")) semer$modelInfo$family else slot(semer, "resp")$family)
          res <- pm_diff_table(seq_len(nrow(D_a)), D_a, D_b, beta_u, sig.epsHat^2*as.matrix(sm_Cov), df_c, crit_c,
                               h = if (glmm) fam$linkinv, d1 = if (glmm) fam$mu.eta,
                               kind = if (glmm) switch(fam$link, log = "ratio", logit = "odds.ratio", "none") else "none")
          tab <- res$table
          y_sd <- attr(mod_df, paste(response, "sd", sep = "_"))
          if (!glmm && !is.null(y_sd)) {
            for (cn in c("diff", "ses", "LL", "UL")) tab[[cn]] <- tab[[cn]]*y_sd
          }
          resp_sc <- glmm && contr_scale == "response"
          plt_df_contr[[response]] <- if (resp_sc) tab$response.diff else tab$diff
          plt_df_contr$LL <- if (resp_sc) tab$response.LL else tab$LL
          plt_df_contr$UL <- if (resp_sc) tab$response.UL else tab$UL
          plt_df_contr$significant <- tab$significant
          stderr_contr <- tab$ses

          if (!is.null(boundary)){
            inBdry <- pm_points_in_poly(plt_df_contr[, covariate], boundary)
            plt_df_contrN <- droplevels(plt_df_contr[inBdry,])
          }else plt_df_contrN <- plt_df_contr
          plt_df_contrSub <- plt_df_contrN[plt_df_contrN$significant, , drop = FALSE]
          
          plt <- ggplot(plt_df_contrN, aes(x=eval(parse(text = covariate[1])), y=eval(parse(text = covariate[2])), z=eval(parse(text = response)))) + 
            geom_contour_filled() +
            geom_point(data=plt_df_contrSub, col="red", size=2.5, shape=22) +
            ggtitle(pm_wrap_title(paste("Factor", modelterm, "level", grp_fct_levl[contr[1]], "vs level", grp_fct_levl[contr[2]], "with Points at the Area out of", (1-level) * 100, if (simultaneous) "% simultaneous band" else "% CI"), 11)) +
            labs(x=covariate[1], y=paste("The difference of ", response, "\n", sep=""), z=response) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5)) 			
          if (response == "response_y") plt <- plt + labs(y=paste("The difference of ", response0, "\n", sep=""))
          if (glmm) {
            pred_df <- cbind(plt_df_contr[setdiff(diff_covarn, c("SE_mu", "significant"))], tab[setdiff(names(tab), c("xvar", "Df"))])
            pred_df$SE_contr <- tab$ses
          } else {
            pred_df <- plt_df_contr
            pred_df$SE_contr <- stderr_contr
            pred_df$df <- tab$Df
            pred_df$SE <- NULL
            pred_df$SE_mu <- NULL
            pred_df$includes0 <- tab$includes0
            pred_df$significant <- tab$significant
          }
        }
      }
    }else if (length(vars)==2) {	 
      if (is.null(plotord) || all(plotord%in%c("NULL", ""))) plotord <- 1:2
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = fact1))))+         	  
        geom_line(aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = fact1))), linewidth=1) +
        labs(x=covariate, y=response, col=fact1) + 
        facet_wrap(formula(paste("~", fact2)), scales=scales)+		
        theme_bw()
      if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))
      if (ci) plt <- plt +
        geom_ribbon(data=pred_df, aes(ymin = LL, ymax = UL, fill=eval(parse(text = fact1))), alpha = 0.3, stat="identity") +
        labs(fill=fact1)
      if (point) plt <- plt + 
        geom_point(data=mod_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = fact1))), position = position_jitter(width = jitterv, height = jitterv))
    }
  }
  #     facet_grid(formula(paste(Variable1, "~", Variable2)))

  if (deriv > 0) {
    ylab <- paste0(if (deriv == 1) "Slope" else "Curvature", " of ", if (response == "response_y") response0 else response,
                   " (d", if (deriv == 2) "2", "/d", covariate, if (deriv == 2) "2", if (!isLMM(semer)) ", link scale", ")")
    if (is.null(contr)) {
      pred_df$includes0 <- pred_df$LL <= 0 & pred_df$UL >= 0
      pred_df$significant <- !pred_df$includes0
      by <- if (!setequal(modelterm, covariate)) unlist(strsplit(modelterm, ":"))
      sig <- pm_deriv_ranges(pred_df, by, covariate)
      attr(pred_df, "significant") <- sig
      if (prt) {
        bandw <- if (simultaneous) "simultaneous" else "pointwise"
        message(if (is.null(sig)) paste0("The ", bandw, " band (level ", level, ") of the ", tolower(sub(" .*", "", ylab)),
                                          " includes 0 at every value of ", covariate, " shown.") else
          paste0(sub(" \\(.*", "", ylab), " (", bandw, " band, level ", level, "): ",
                 paste0(ifelse(nzchar(sig$level), paste0(sig$level, ": "), ""), format(sig$from, digits = 4), " to ",
                        format(sig$to, digits = 4), " (", sig$direction, ")", collapse = "; "), "."))
      }
      if (inherits(plt, "ggplot")) plt <- plt + geom_hline(yintercept = 0, linetype = 2, colour = "grey40") + labs(y = ylab)
    } else if (inherits(plt, "ggplot")) {
      plt <- plt + labs(y = paste0("The difference of the ", tolower(substr(ylab, 1, 1)), substring(ylab, 2)))
      if (!is.null(plt$labels$title))
        plt <- plt + ggtitle(sub("^Difference between", paste("Difference between the", if (deriv == 1) "slopes" else
          "curvatures", "of"), plt$labels$title))
    }
  }
  
  if (prt) {
    options(warn = -1)	  
    if (length(covariate)==2 && threeD && !setequal(modelterm, covariate)) {
      for (i in names(plt)) print(plt[[i]])
    }else print(plt)
    options(warn = 0) 
  }
  rownames(pred_df) <- NULL
  return(invisible(list(plt=plt, pred_df=pred_df)))
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
  # (pm_sr_uncond(), R/semireg-engine.R)
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
      corr <- pm_sr_uncond(semireg)
      if (!is.null(corr)) V <- as.matrix(V) + corr[idx, idx]
    }
    V
  }
  if (!is.null(semireg$smooths)) {
    # smooths written in the formula: the blocks by name, their design for the
    # grid from the stored mgcv smooth objects (R/semireg-smooth.R); random
    # curves (s(x) | g) are random effects, not part of the population curves
    blk <- intersect(names(semireg$u_lst), names(semireg$smooths$maps))
    if (is.null(sm_term)) sm_term <- blk
    bad <- setdiff(sm_term, blk)
    if (length(bad)) stop("'sm_term' must be among the smooth blocks: ", paste(blk, collapse = ", "), call. = FALSE)
    blocks <- blk[blk %in% sm_term]
    sm_Cov <- cov_of(blocks)
    smooth_basis <- function(nd) do.call(cbind, pm_sm_predict(semireg$smooths, nd, blocks))
    ZN <- pm_refgrid_average(KK$refgrid, smooth_basis)
    D_mx <- cbind(K, ZN)
    beta_u <- c(betaHat, do.call("c", unname(semireg$u_lst[blocks])))
    return(list(semer = semer, KK = KK, K = K, ZN = ZN, mod_df = mod_df, pred_df = pred_df, response = response,
                response0 = response0, D_mx = D_mx, beta_u = beta_u, sm_Cov = sm_Cov, sig.epsHat = sig.epsHat,
                cov_idx = cov_idx(blocks),
                sm_term = blocks))
  }
  # the smoothZ form: terms are the names of the smoothZ list; a term with
  # group = TRUE has a block per level (named term.level)
  blk <- names(semireg$u_lst)
  grp <- sub("\\..*", "", blk)
  terms_all <- intersect(names(as.list(semireg$smoothZ_call))[-1], grp)
  if (is.null(sm_term)) sm_term <- terms_all
  bad <- setdiff(sm_term, c(terms_all, blk))
  if (length(bad)) stop("'sm_term' must be among the smooth terms: ", paste(terms_all, collapse = ", "), call. = FALSE)
  blocks <- blk[grp %in% sm_term | blk %in% sm_term]
  sm_term <- intersect(terms_all, unique(sub("\\..*", "", blocks)))
  sm_Cov <- cov_of(blocks)
  u_lst <- semireg$u_lst[blocks]

  knots_lst <- semireg$knots_lst[sm_term]
  range_lst <- semireg$range_lst[sm_term] 
  type_lst <- semireg$type_lst[sm_term]  
  call_lst <- as.list(semireg$smoothZ_call)[sm_term]
  
  call_lstN <- Z_lstN <- vector("list", length(call_lst))
  names(call_lstN) <- names(Z_lstN) <- sm_term 
  
  for (i in sm_term) {
    xcall_lst  <- as.list(call_lst[[i]])
    if (!is.null(xcall_lst$type) && xcall_lst$type=="smspline") xcall_lst$pred <- TRUE
    if (!is.null(xcall_lst$k)) xcall_lst$k <- NULL
    if (!is.null(range_lst[[i]]) && type_lst[[i]]!="Ztps") xcall_lst$range.x <- range_lst[[i]]
    xcall_lst$intKnots <- knots_lst[[i]]
    call_lstN[[i]] <- as.call(xcall_lst) 
  }  
  
  # Smooth basis for a set of grid rows, in the order of the random effects.
  smooth_basis <- function(nd) {
    for (i in sm_term) {
      Z_lstN[[i]] <- with(nd, eval(call_lstN[[i]]))
    }

    Z_lstNN <- list()
    for (i in sm_term){
      if (is.list(Z_lstN[[i]])) Z_lstNN <- c(Z_lstNN, unlist(Z_lstN[i], recursive=FALSE)) else Z_lstNN <- c(Z_lstNN, Z_lstN[i])
    }

    Z_lstNN <- Z_lstNN[blocks]
    do.call("cbind", Z_lstNN)
  }

  # Evaluate the smooth basis on the full reference grid and average it with
  # the same weights as K, so factors that are averaged over (including 'by'
  # factors of the smooths) are treated identically in both parts.
  ZN <- pm_refgrid_average(KK$refgrid, smooth_basis)
  D_mx <- cbind(K,ZN)
  u_hat <- do.call("c", unname(u_lst))
  beta_u <- c(betaHat, u_hat)
  
  list(semer = semer, KK = KK, K = K, ZN = ZN, mod_df = mod_df, pred_df = pred_df, response = response,
       response0 = response0, D_mx = D_mx, beta_u = beta_u, sm_Cov = sm_Cov, sig.epsHat = sig.epsHat,
                cov_idx = cov_idx(blocks),
       sm_term = sm_term)
}

## Derivatives of the curves of semipred(): the rows of the design D(x) on
## the grid replaced by dD/dx (deriv = 1) or d2D/dx2 (deriv = 2), by central
## differences of the whole design (fixed part and every smooth basis; for
## the O'Sullivan basis they equal ZOSull(drv = ) to 1e-8, tested), the same
## grid rebuilt at x + h and x - h. Per unit of the covariate as given (a
## covariate scaled by semireg(nscale = TRUE): divided by its SD).
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
  s <- attr(mod_df, paste(covariate, "sd", sep = "_"))
  if (!is.null(s)) D <- D/s^deriv
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
