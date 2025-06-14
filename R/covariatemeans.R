#' Predicted Means of a Linear Model with Covariate Variable(s)
#'
#' This function obtains predicted means with graph for a new set of covariate
#' values.
#'
#'
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, and \code{lmer}.
#' @param modelterm Name (in "quotes") for indicating which factor term's
#' predicted mean to be calculated.  The \code{modelterm} must be given exactly
#' as it appears in the printed model, e.g. "A" or "A:B".
#' @param covariate Name (in "quotes") of one the covariate variables in the
#' \code{model}.
#' @param as.is A logic value to specify wheather or not using original
#' covariate values' rage for graph, the default is FALSE.
#' @param covariateV A numeric vector when as.is is FALSE, then covariatemeans
#' will produce the result for \code{covariate} at value of \code{covariateV}.
#' @param data In some cases, you need to provide the data set used in model
#' fitting, especially when you have applied some variable trnasformation in
#' the model.
#' @param level A significant level for calculating confident interval. The
#' default value is 0.05.
#' @param Df A degree of freedom for calculating CI of predicted means (you can
#' manually specified ddf here). For the above models, ddf is obtained from the
#' function automatically.
#' @param trans A function object for calculating the back transformed means,
#' e.g. \code{trans=exp}.
#' @param transOff When you use \code{trans=exp(x+1)}, then \code{transOff=1},
#' the default is 0.
#' @param responsen Name (in "quotes") of the back transformed response
#' variable in the \code{model}.
#' @param trellis A logical variable. If set to TRUE (default), a trellis plots
#' of predicted means with CI will be drawn.
#' @param plotord A numeric (or character) vector specifying the order of
#' plotting for two or three way interaction (e.g.  \code{plotord = c(2, 1, 3)}
#' or \code{plotord = c("B", "A", "C")} will put the second variable 'B' in
#' \code{modelterm} on the \code{X} axis, the first variable 'A' as the
#' grouping variable, and the third one 'C' as the panel variable). The
#' defaults are \code{c(1, 2)} and \code{c(1, 2, 3)} for two and three way
#' interactions.
#' @param mtitle The main title in the graph.
#' @param ci A logical variable to indicate whether to print confidence
#' interval. The default value is TRUE.
#' @param point A logical variable to indicate whether to print raw data
#' points. The default value is TRUE.
#' @param jitterv A degree of jitter in x and y direction in the graph. The
#' default is zero.
#' @param newwd A logical variable to indicate whether to print graph in a new
#' window. The default value is FALSE.
#' @return \item{Predicted Means}{A table of predicted means.}
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
#' library(predictmeans)
#' data(Oats, package="nlme")
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'   covariatemeans(fm, "Variety", covariate="nitro")
#'   covariatemeans(fm, "Variety", covariate="nitro", covariateV=seq(0, 0.6, 0.1))$pltdf
#' @importFrom rlang inherits_any
#' @importFrom stats aggregate
#' @importFrom utils edit
#' @export

covariatemeans <- function(model,
                           modelterm = NULL,
                           covariate,
                           as.is = FALSE,
                           covariateV = NULL,
                           data = NULL,
                           level = 0.05,
                           Df = NULL,
                           trans = NULL,
                           transOff = 0,
                           responsen = NULL,
                           trellis = TRUE,
                           plotord = NULL,
                           mtitle = NULL,
                           ci = TRUE,
                           point = TRUE,
                           jitterv = 0,
                           newwd = FALSE) {
  if (is.null(modelterm) || all(modelterm  %in%  c("NULL", ""))) {
    modelterm <- covariate
    trellis <- FALSE
  }

  if (!is.null(covariateV) && length(covariate) > 1) {
    stopifnot(is.matrix(covariateV),
              dim(covariateV)[2] == length(covariate))
    colnames(covariateV) <- covariate
  }

  vars <- unlist(strsplit(modelterm, "\\:"))
  if (!is.null(plotord) && !unique(plotord  %in%  c("NULL", ""))) {
    if (length(plotord) != length(vars)) {
      stop(paste("plotord must be a character vector of length", length(vars)))
    }
    if (is.character(plotord)) {
      i <- match(plotord, vars)
      if (any(is.na(i))) {
        stop(paste(
          "The entries",
          paste0(plotord[is.na(i)], sep = ", "),
          "do not match any term in ",
          modelterm
        ))
      } else {
        plotord <- i
      }
    }
  }
  ctr.matrix <- Kmatrix(model, modelterm, covariate, covariateV, data)
  KK <- ctr.matrix$K
  pltdf <- ctr.matrix$fctnames
  response <- ctr.matrix$response
  preddf <- ctr.matrix$preddf
  mp <- mymodelparm(model)
  bhat <- mp$coef

  if (!setequal(modelterm, covariate) && as.is) {
    agg_formula <- formula(paste(covariate, "~", paste(vars, collapse = "+")))
    min_df <- aggregate(agg_formula, preddf, min, na.rm = TRUE)
    names(min_df)[ncol(min_df)] <- "min_value"
    max_df <- aggregate(agg_formula, preddf, max, na.rm = TRUE)
    names(max_df)[ncol(max_df)] <- "max_value"
    pltdf <- merge(
      merge(pltdf, min_df, by = vars, sort = FALSE),
      max_df,
      by = vars,
      sort = FALSE
    )
    pltdf <- pltdf[pltdf[[covariate]] >= pltdf[["min_value"]] &
                     pltdf[[covariate]] <= pltdf[["max_value"]], ]
    KK <- KK[as.numeric(rownames(pltdf)), ]
  }

  # We'll work only with the non-NA elements of bhat
  KK <- KK[, mp$estimable, drop = FALSE]
  pltdf$yhat <- as.numeric(KK %*% bhat)
  pltdf$ses <- sqrt(base::diag(KK %*% tcrossprod(mp$vcov, KK)))

  if (is.null(Df) || all(Df %in% c("NULL", ""))) {
    if (inherits(model, "lme")) {
      Df <- terms(model$fixDF)[modelterm]
    } else if (inherits(model, "lmerMod")) {
      # Df <- try(median(df_term(model, modelterm, covariate), na.rm=TRUE))
      # Df <- mp$df[modelterm]
      # if(inherits(Df, "try-error")) stop("You need provide Df for this model!")
      L_term <- get_contrasts_type1(model)[[modelterm]]
      Df <- mean(df_term(model, ctrmatrix = L_term))
      if (is.nan(Df)) {
        stop("You need provide Df for this model term!")
      }
    } else {
      Df <- mp$df
    }
    if (Df == 0) {
      stop("You need provide Df for this model!")
    }
  }

  if (inherits(model, "glm") && is.null(trans)) {
    trans <- model$family$linkinv
  }
  if (inherits(model, "glmerMod") && is.null(trans)) {
    trans <- slot(model, "resp")$family$linkinv
  }
  if (inherits(model, "glmmTMB") && is.null(trans)) {
    trans <- model$modelInfo$family$linkinv
  }

  Mean <- LL <- UL <- xvar <- factors <- bky <- NULL
  if (is.null(trans)) {
    pltdf$Mean <- pltdf$yhat
    pltdf$LL <- pltdf$yhat - qt(1 - level / 2, df = Df) * pltdf$ses
    pltdf$UL <- pltdf$yhat + qt(1 - level / 2, df = Df) * pltdf$ses
  } else {
    pltdf$Mean <- trans(pltdf$yhat) - transOff
    if (identical(trans, make.link("log")$linkinv) || identical(trans, exp)) {
      pltdf$Mean <- exp(pltdf$yhat + pltdf$ses / 2) - transOff
    }
    pltdf$LL <- trans(pltdf$yhat - qt(1 - level / 2, df = Df) * pltdf$ses) -
      transOff
    pltdf$UL <- trans(pltdf$yhat + qt(1 - level / 2, df = Df) * pltdf$ses) -
      transOff
  }

  # pltdf$yhat <- pltdf$ses <- NULL
  pltdf$yhat <- NULL

  if (setequal(modelterm, covariate)) {
    pltdf$factors <- factor(1)
  } else {
    pltdf$factors <- factor(do.call("paste", c(pltdf[, vars, drop = FALSE], sep =
                                                 ":")))
  }
  colnames(pltdf)[colnames(pltdf) == covariate] <- "xvar"

  # delete empty factor combinations
  if (!is.null(data)) {
    mdf <- data
  } else {
    mdf <- model.frame(model)
  }
  if (length(setdiff(response, names(mdf))) != 0) {
    mdf[, response] <- eval(parse(text = response), mdf)
  }
  mdf <- cbind(mdf, preddf[, !names(preddf) %in% names(mdf), drop = FALSE])

  ndf <- data.frame(table(mdf[, vars, drop = FALSE]))
  if (any(ndf$Freq == 0)) {
    ndf0 <- ndf[ndf$Freq == 0, , drop = FALSE]
    ndf0$factors <- factor(do.call("paste", c(ndf0[, vars, drop = FALSE], sep =
                                                ":")))
    pltdf <- pltdf[!pltdf$factors %in% ndf0$factors, ]
  }

  if (is.null(mtitle) || mtitle %in% c("NULL", "")) {
    mtitle <- paste("Fitted and observed relationship with",
                    (1 - level) * 100,
                    "% CI")
  }
  if (is.null(trans)) {
    mdf$bky <- mdf[, response]
  } else {
    if (length(setdiff(response, names(mdf))) == 0) {
      ## Transformed y before modelling
      if (inherits_any(model, c("glm", "glmerMod", "glmmTMB"))) {
        if (inherits(mdf[, response], "factor")) {
          mdf$bky <- as.numeric(mdf[, response]) - 1
        } else if (!is.null(dim(mdf[, response]))) {
          mdf$bky <- mdf[, response][, 1] / rowSums(mdf[, response])
          response <- "Probability"
        } else {
          mdf$bky <- mdf[, response]
        }
        if (any(
          identical(trans, function(x) x, ignore.environment = TRUE),
          identical(trans, I, ignore.environment = TRUE)
        )) {
          if (inherits(model, "glm")) {
            mdf$bky <- ifelse(mdf$bky >= 1 |
                                mdf$bky <= 0,
                              NA,
                              model$family$linkfun(mdf$bky))
          }
          if (inherits(model, "glmerMod")) {
            mdf$bky <- ifelse(
              mdf$bky >= 1 |
                mdf$bky <= 0,
              NA,
              slot(model, "resp")$family$linkfun(mdf$bky)
            )
          }
          if (inherits(model, "glmmTMB")) {
            mdf$bky <- ifelse(mdf$bky >= 1 |
                                mdf$bky <= 0,
                              NA,
                              model$modelInfo$family$linkfun(mdf$bky))
          }
          response <- "Response"
        }
      } else {
        mdf$bky <- trans(mdf[, response])
        response <- paste("Transformed", response)
      }
    } else {
      ## Transformed y within modelling
      response <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
      if (length(setdiff(response, names(mdf))) != 0) {
        if (is.null(responsen) || all(responsen %in% c("NULL", ""))) {
          stop(
            "Please provide suitable name for response variable using option 'responsen'!"
          )
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
    mdf$factors <- do.call("paste", c(mdf[, vars, drop = FALSE], sep = ":"))
  }
  names(mdf)[names(mdf) == covariate] <- "xvar"
  pltdf[, c("Mean", "LL", "UL")] <- lapply(pltdf[, c("Mean", "LL", "UL")], as.numeric)

  if (!trellis) {
    if (newwd) {
      dev.new()
    }
    if (setequal(modelterm, covariate)) {
      plt <- ggplot(pltdf, aes(xvar, Mean)) +
        labs(
          title = paste(mtitle, "\n", sep = ""),
          x = paste("\n", covariate, sep = ""),
          y = paste(response, "\n")
        ) +
        geom_line(linewidth = 0.5) +
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(
          aes(ymin = LL, ymax = UL),
          alpha = 0.2,
          data = pltdf,
          stat = "identity"
        )
      }
      if (point) {
        plt <- plt + geom_point(
          aes(x = xvar, y = bky),
          position = position_jitter(width = jitterv, height = jitterv),
          data = mdf
        )
      }
    } else {
      plt <- ggplot(pltdf, aes(xvar, Mean, colour = factors)) +
        labs(
          title = paste(mtitle, "\n", sep = ""),
          x = paste("\n", covariate, sep = ""),
          y = paste(response, "\n")
        ) +
        geom_line(linewidth = 0.5) +
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(
          aes(
            ymin = LL,
            ymax = UL,
            fill = factors
          ),
          alpha = 0.2,
          data = pltdf,
          stat = "identity"
        )
      }
      if (point) {
        plt <- plt + geom_point(
          aes(x = xvar, y = bky),
          position = position_jitter(width = jitterv, height = jitterv),
          data = mdf
        )
      }
      plt <- plt + guides(col = guide_legend(modelterm), fill = guide_legend(modelterm))
    }
    print(plt)
  } else {
    if (length(vars) == 1) {
      if (newwd) {
        dev.new()
      }
      plt <- ggplot(pltdf, aes(xvar, Mean, colour = factors)) +
        labs(
          title = paste(mtitle, "\n", sep = ""),
          x = paste("\n", covariate, sep = ""),
          y = paste(response, "\n")
        ) +
        geom_line(linewidth = 0.5) +
        facet_wrap(~  factors) +
        theme_bw()
      if (ci) {
        plt <- plt + geom_ribbon(
          aes(
            ymin = LL,
            ymax = UL,
            fill = factors
          ),
          alpha = 0.2,
          data = pltdf,
          stat = "identity"
        )
      }
      if (point) {
        plt <- plt + geom_point(
          aes(x = xvar, y = bky),
          position = position_jitter(width = jitterv, height = jitterv),
          data = mdf
        )
      }
      plt <- plt + guides(col = guide_legend(modelterm), fill = guide_legend(modelterm))
      if (setequal(modelterm, covariate)) {
        plt <- plt + theme(legend.position = "none")
      }
      print(plt)
    }
    if (length(vars) == 2) {
      if (newwd) {
        dev.new()
      }
      if (is.null(plotord) || all(plotord  %in%  c("NULL", ""))) {
        plotord <- 1:2
      }
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      plt <- ggplot(pltdf, aes(xvar, Mean, colour = factor(.data[[fact1]]))) +
        labs(
          title = paste(mtitle, "\n", sep = ""),
          x = paste("\n", covariate, sep = ""),
          y = paste(response, "\n")
        ) +
        geom_line(linewidth = 0.5) +
        facet_grid(as.formula(paste("~", fact2, sep = ""))) +
        theme_bw()
      if (ci) {
        plt <- plt +
          geom_ribbon(
            aes(
              ymin = LL,
              ymax = UL,
              fill = factor(.data[[fact1]])
            ),
            alpha = 0.2,
            data = pltdf,
            stat = "identity"
          )
      }
      if (point) {
        plt <- plt +
          geom_point(
            aes(x = xvar, y = bky),
            position = position_jitter(width = jitterv, height = jitterv),
            data = mdf
          )
      }
      plt <- plt + guides(col = guide_legend(fact1), fill = guide_legend(fact1))
      print(plt)
    }
    if (length(vars) == 3) {
      if (newwd) {
        dev.new()
      }
      if (is.null(plotord) || all(plotord %in% c("NULL", ""))) {
        plotord <- 1:3
      }
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      fact3 <- (vars[plotord])[3]
      plt <- ggplot(pltdf, aes(xvar, Mean, colour = factor(.data[[fact1]]))) +
        labs(
          title = paste(mtitle, "\n", sep = ""),
          x = paste("\n", covariate, sep = ""),
          y = paste(response, "\n")
        ) +
        geom_line(linewidth = 0.5) +
        facet_grid(as.formula(paste(fact2, "~", fact3, sep = ""))) +
        theme_bw()
      if (ci) {
        plt <- plt +
          geom_ribbon(
            aes(
              ymin = LL,
              ymax = UL,
              fill = factor(.data[[fact1]])
            ),
            alpha = 0.2,
            data = pltdf,
            stat = "identity"
          )
      }
      if (point) {
        plt <- plt +
          geom_point(
            aes(x = xvar, y = bky),
            position = position_jitter(width = jitterv, height = jitterv),
            data = mdf
          )
      }
      plt <- plt + guides(col = guide_legend(fact1), fill = guide_legend(fact1))
      print(plt)
    }
  }

  return(invisible(list(plt = plt, pltdf = pltdf)))
}
