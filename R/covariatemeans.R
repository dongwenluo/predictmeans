covariatemeans <- function (model, modelterm=NULL, covariate, as.is=FALSE, covariateV=NULL, data=NULL, level=0.05, Df=NULL, trans=NULL, transOff=0, responsen=NULL, trellis=TRUE, plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=FALSE, cov_reduce=mean, exclude=NULL,
                            weights=c("equal", "proportional", "outer", "cells", "flat")) {

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
  mp <- mymodelparm(model)
  bhat <- mp$coef

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
    if (inherits(model, "lme")) {
      Df <- terms(model$fixDF)[modelterm]
    } else if (inherits(model, "lmerMod")) {
      # Df <- try(median(df_term(model, modelterm, covariate), na.rm=TRUE))
	    Df <- mp$df[modelterm]
      if(inherits(Df, "try-error")) {
        stop("You need provide Df for this model!")
      }
    } else {
      df_custom <- pm_df(model, KK)
      Df <- if (!is.null(df_custom)) rep_len(as.numeric(df_custom), nrow(KK)) else mp$df
    }
    if (any(Df==0)) {
      stop("You need provide Df for this model!")
    }
  }

  if (is.null(trans)) {
    trans <- pm_link_trans(model)   # NULL for an identity link
  }

  Mean <- LL <- UL <- xvar <- factors <- bky <- NULL
  if (is.null(trans)) {
    pltdf$Mean <- pltdf$yhat
    pltdf$LL <- pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses
    pltdf$UL <- pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses
  } else {
    pltdf$Mean <- trans(pltdf$yhat)-transOff
    if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)) {
      pltdf$Mean <- exp(pltdf$yhat+pltdf$ses/2)-transOff
    }
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
        labs(title=paste(mtitle, "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
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
        labs(title=paste(mtitle, "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
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
        labs(title=paste(mtitle, "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
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
        labs(title=paste(mtitle, "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
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
        labs(title=paste(mtitle, "\n", sep=""), x=paste("\n", covariate, sep=""), y=paste(response, "\n"))+
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
