#' Predicted Means of a Linear Model with Covariate Variable(s)
#'
#' @description
#' This function obtains predicted means with graph for a new set of covariate values.
#' @usage
#' covariatemeans(model, modelterm=NULL, covariate, as.is=FALSE, covariateV=NULL, 
#'   data=NULL, level=0.05, Df=NULL, trans=NULL, transOff=0, responsen=NULL, trellis=TRUE, 
#'   plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=FALSE,
#'   cov_reduce=mean, exclude=NULL,
#'   weights=c("equal", "proportional", "outer", "cells", "flat"),
#'   df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
#'   vcov_adjust=TRUE, contr=NULL)
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm}, \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's predicted mean to be calculated. 
#'   The \code{modelterm} must be given exactly as it appears in the printed model, e.g. "A" or "A:B".
#' @param covariate Name (in "quotes") of one the covariate variables in the \code{model}.
#' @param as.is A logic value to specify wheather or not using original covariate values' rage for graph, the default is FALSE.
#' @param covariateV A numeric vector when as.is is FALSE, then covariatemeans will produce the result for \code{covariate} at value of \code{covariateV}.
#' @param data In some cases, you need to provide the data set used in model fitting, especially when you have applied some variable trnasformation in the model.
#' @param level A significant level for calculating confident interval. The default value is 0.05. The limits
#'   are for each line on its own: whether the bands of two groups overlap is not a test of their difference
#'   (the lines are estimated from the same model and are correlated; see \code{letterCI} in
#'   \code{\link{predictmeans}}). Use \code{contr} for the band of the difference between two groups, or compare
#'   groups at chosen covariate values with \code{\link{predictmeans}} (\code{covariate = }) and their slopes
#'   with \code{\link{pm_slopes}}.
#' @param Df A degree of freedom for calculating CI of predicted means (you can manually specified ddf here). For the above models, ddf is obtained from the function automatically.
#' @param trans A function object for calculating the back transformed means, e.g. \code{trans=exp}. When it is
#'   not given, the inverse link of the model is used, or the inverse of a transformation of the response in the
#'   model formula (\code{log(y)}, \code{sqrt(y)}, ...; see \code{\link{predictmeans}}), with a message;
#'   \code{trans = I} keeps the model scale. The back-transformed line is the back-transformed fitted value.
#'   (Before 2.0 half the standard error of the fitted value was added inside \code{exp()} for a log.)
#' @param transOff When you use \code{trans=exp(x+1)}, then \code{transOff=1}, the default is 0.
#' @param responsen Name (in "quotes") of the back transformed response variable in the \code{model}.
#' @param trellis A logical variable. If set to TRUE (default), a trellis plots of predicted means with CI will be drawn.
#' @param plotord A numeric (or character) vector specifying the order of plotting for two or three way interaction (e.g. 
#'   \code{plotord = c(2, 1, 3)} or \code{plotord = c("B", "A", "C")} will put the second variable 'B' in 
#'   \code{modelterm} on the \code{X} axis, the first variable 'A' as the grouping variable, and the third one 
#'   'C' as the panel variable). The defaults are \code{c(1, 2)} and \code{c(1, 2, 3)} for two and three way 
#'   interactions.
#' @param mtitle The main title in the graph.
#' @param ci A logical variable to indicate whether to print confidence interval. The default value is TRUE.
#' @param point A logical variable to indicate whether to print raw data points. The default value is TRUE.
#' @param jitterv A degree of jitter in x and y direction in the graph. The default is zero.
#' @param newwd A logical variable to indicate whether to print graph in a new window. The default value is FALSE.
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
#' @param df_method,vcov_adjust As in \code{\link{predictmeans}}: by default Kenward-Roger df (each point of the curves
#'   has its own) and adjusted standard errors for models fitted by REML (\code{lmer}, \code{lme}, \code{gls}).
#'   predictmeans 1.x used the unadjusted covariance and, for \code{lmer}, one df per term; \code{vcov_adjust = FALSE}
#'   and, for \code{lme} and \code{gls}, \code{df_method = "containment"} give its results.
#' @param contr Two levels of \code{modelterm} (numbers or names, e.g. \code{c(1, 3)} or
#'   \code{c("Golden Rain", "Victory")}; for a term with several factors the combined labels such as
#'   \code{"a:x"}): the difference between the first and the second along the covariate, with its confidence
#'   band, is returned and plotted instead of the lines. Whether the bands of the two lines overlap does not say
#'   whether the lines differ (they are estimated from the same model and are correlated); the band of the
#'   difference does. It is computed from the contrast of the two lines, not from their limits:
#'   \itemize{
#'   \item On the model (link) scale: \eqn{d = (K_a - K_b)\hat\beta} with
#'     \eqn{SE^2 = (K_a - K_b) V (K_a - K_b)'} and limits \eqn{d \pm t\,SE} at the contrast's df (columns
#'     \code{diff}, \code{ses}, \code{Df}, \code{LL}, \code{UL}). This is the scale of the test: the lines
#'     differ at a covariate value when this band leaves out zero. For a log link or a log response the table
#'     also has the \code{ratio} (exp of \code{diff} and of its limits), for a logit link the
#'     \code{odds.ratio}.
#'   \item On the response scale (models with a link, a transformed response or \code{trans}): the difference
#'     of the back-transformed lines, \eqn{h(\eta_a) - h(\eta_b)}, with the delta-method standard error from
#'     the gradient \eqn{h'(\eta_a) K_a - h'(\eta_b) K_b} (columns \code{response.diff},
#'     \code{response.ses}, \code{response.LL}, \code{response.UL}). This is what is plotted, as the lines are
#'     drawn on that scale; \code{trans = I} keeps the model scale. For a generalized linear mixed model it is the
#'     difference for a typical cluster (random effects at zero), as the lines are.
#'   }
#' @return
#' \item{plt}{The plot (a \code{ggplot}).}
#' \item{pltdf}{The data of the plot: the fitted values with their limits, or with \code{contr} the differences
#'   on the model scale and, when back-transformed, on the response scale (see \code{contr}).}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   data(Oats, package="nlme")
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' # library(lme4)
#' # fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#'   covariatemeans(fm, "Variety", covariate="nitro")
#'   covariatemeans(fm, "Variety", covariate="nitro", covariateV=seq(0, 0.6, 0.1))$pltdf
#' @name covariatemeans
NULL

covariatemeans <- function (model, modelterm=NULL, covariate, as.is=FALSE, covariateV=NULL, data=NULL, level=0.05, Df=NULL, trans=NULL, transOff=0, responsen=NULL, trellis=TRUE, plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=FALSE, cov_reduce=mean, exclude=NULL,
                            weights=c("equal", "proportional", "outer", "cells", "flat"),
                            df_method=c("auto", "Kenward-Roger", "Satterthwaite", "containment", "between-within", "KR-expected", "KR-average", "inner-outer", "residual"),
                            vcov_adjust=TRUE, contr=NULL) {
  vcov_adjust_default <- missing(vcov_adjust)
  df_method <- match.arg(df_method)
  if (inherits(model, "aovlist")) model <- pm_aov_carrier(model)     # by strata (R/aovlist.R)

  if (is.null(modelterm) || all(modelterm%in%c("NULL", ""))) {
    modelterm <- covariate
    trellis=FALSE
  }

  if (!is.null(covariateV) && length(covariate) > 1){
      stopifnot(is.matrix(covariateV), dim(covariateV)[2] == length(covariate))
      colnames(covariateV) <- covariate
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
  weights <- match.arg(weights)
  ctr.matrix <- Kmatrix(model, modelterm, covariate, covariateV, data, cov_reduce=cov_reduce, exclude=exclude, weights=weights)
  KK <- ctr.matrix$K
  pltdf <- ctr.matrix$fctnames
  response <- ctr.matrix$response
  preddf <- ctr.matrix$preddf
  inf <- pm_inference(model, df_method, vcov_adjust, vcov_adjust_default)     # R/inference.R
  mp <- inf$mp
  bhat <- mp$coef

  # the difference between two levels along the covariate, with its own band
  if (!is.null(contr)) {
    return(pm_cov_difference(contr, KK[, mp$estimable, drop = FALSE], pltdf, vars, covariate, modelterm, mp, inf,
                             level, Df, model, trans, mtitle, newwd))
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

  if (is.null(Df) || all(Df%in%c("NULL", ""))) {
    Df <- inf$df(KK, modelterm)          # each point its own df (Kenward-Roger, Satterthwaite)
    if (any(Df==0)) {
      stop("You need provide Df for this model!")
    }
  }

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
  if (is.null(trans)) {
    pltdf$Mean <- pltdf$yhat
    pltdf$LL <- pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses
    pltdf$UL <- pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses
  } else {
    # (1.x added half the SE of the fitted value inside exp() for a log:
    # not a bias adjustment, and the line then depended on its own precision)
    pltdf$Mean <- trans(pltdf$yhat)-transOff
    pltdf$LL <- trans(pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses)-transOff
    pltdf$UL <- trans(pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses)-transOff
  }

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
    pltdf <- pltdf[!pltdf$factors%in%ndf0$factors, ]
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
        mdf$bky <- trans(mdf[, response]) - transOff
        response <- paste("Transformed", response)
      }
    } else {       ## Transformed y within modelling
      response <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
      if (length(setdiff(response, names(mdf)))!=0) {
        if (is.null(responsen) || all(responsen%in%c("NULL", "")))  {
          stop("Please provide suitable name for response variable using option 'responsen'!")
        }
        response <- responsen
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
    if (newwd) {
      dev.new()
    }
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
    print(plt)
  }else{
    if (length(vars)==1) {
      if (newwd) {
        dev.new()
      }
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
      print(plt)
    }
    if (length(vars)==2) {
      if (newwd) {
        dev.new()
      }
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
      print(plt)
    }
    if (length(vars)==3) {
      if (newwd) {
        dev.new()
      }
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
      print(plt)
    }
  }

  return(invisible(list(plt=plt, pltdf=pltdf)))
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
                              mtitle, newwd) {
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
  V <- mp$vcov
  out <- data.frame(xvar = x[a], diff = as.numeric(L %*% mp$coef), ses = sqrt(rowSums((L %*% V) * L)))
  out$Df <- if (is.null(Df) || all(Df %in% c("NULL", ""))) inf$df(L, modelterm) else Df
  crit <- stats::qt(1 - level/2, df = out$Df)
  out$LL <- out$diff - crit*out$ses
  out$UL <- out$diff + crit*out$ses

  # the back-transformation: the inverse link, the response transformation
  # in the formula, or 'trans' (trans = I: the model scale only)
  link <- pm_link_trans(model)
  if (is.null(trans)) trans <- link
  if (is.null(trans)) trans <- pm_response_trans(model)$fun
  back <- !is.null(trans) && !identical(trans, I)
  link_name <- if (!is.null(link) && identical(trans, link)) {
    tryCatch(stats::family(model)$link, error = function(e) NA_character_)
  } else NA_character_
  if (back) {
    # a log: the difference is the log of a ratio; a logit: of an odds ratio
    if (identical(trans, exp) || identical(link_name, "log")) {
      out$ratio <- exp(out$diff); out$ratio.LL <- exp(out$LL); out$ratio.UL <- exp(out$UL)
    } else if (identical(trans, stats::plogis) || identical(link_name, "logit")) {
      out$odds.ratio <- exp(out$diff); out$odds.ratio.LL <- exp(out$LL); out$odds.ratio.UL <- exp(out$UL)
    }
    # the difference on the response scale, by the delta method
    eta_a <- as.numeric(Ka %*% mp$coef)
    eta_b <- as.numeric(Kb %*% mp$coef)
    d1 <- pm_back_funs(trans, trans, model, !is.na(link_name))$d1
    G <- d1(eta_a) * Ka - d1(eta_b) * Kb
    out$response.diff <- as.numeric(trans(eta_a)) - as.numeric(trans(eta_b))
    out$response.ses <- sqrt(rowSums((G %*% V) * G))
    out$response.LL <- out$response.diff - crit*out$response.ses
    out$response.UL <- out$response.diff + crit*out$response.ses
  }
  names(out)[1] <- covariate
  if (is.null(mtitle) || mtitle %in% c("NULL", "")) {
    mtitle <- paste0("Difference between '", modelterm, "' ", levs[k[1]], " and ", levs[k[2]], " along ", covariate,
                     " with ", (1 - level)*100, "% CI", if (back) " (response scale)" else "")
  }
  # the plot: the difference of the lines as covariatemeans() draws them
  pd <- if (back) data.frame(x = out[[covariate]], y = out$response.diff, lo = out$response.LL, hi = out$response.UL) else
    data.frame(x = out[[covariate]], y = out$diff, lo = out$LL, hi = out$UL)
  plt <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::labs(title = paste0(pm_wrap_title(mtitle, 11), "\n"), x = covariate,
                  y = paste0("Difference (", levs[k[1]], " - ", levs[k[2]], ")"),
                  subtitle = if (back) pm_wrap_title(paste0(
                    "Delta-method band of the difference of the back-transformed lines",
                    if (inherits(model, c("merMod", "glmmTMB", "lme"))) " (random effects at zero)" else "",
                    "; the test is on the model scale (columns diff, LL, UL)"), 11, rel = 1)) +
    ggplot2::theme_bw()
  if (newwd) grDevices::dev.new()
  print(plt)
  invisible(list(plt = plt, pltdf = out))
}
