#' Diagnostic Plots for a Linear (Mixed) Model
#'
#' This function produces diagnostic plots for linear models including 'aov',
#' 'lm', 'glm', 'gls', 'lme' and 'lmer'.
#'
#'
#' @param model Model object returned by \code{aov}, \code{lm}, \code{glm},
#' \code{gls}, \code{lme}, and \code{lmer}.
#' @param group Name (in "quotes") for indicating the variable used to show
#' grouping in the residual vs predicted plot. If variable is a term in the
#' model, then group will be a name of the variable such as \code{group="A"},
#' otherwise group will be the actual variable such as \code{group=data$A}.
#' @param level An integer 1, 2, etc. used to specify a level of the random
#' effect for plotting. The default value is 1.
#' @param slope A logical variable. If set to TRUE, a Q-Q plot of random slope
#' will be drawn.
#' @param id A logical variable. If set to TRUE, outliers in the residual vs
#' fitted plot can be identified interactively.
#' @param newwd A logical variable to indicate whether to print graph in a new
#' window. The default is FALSE.
#' @param ask logical. If TRUE (and the R session is interactive) the user is
#' asked for input, before a new figure is drawn.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'
#' ## Note that the order of levels of nested random effects is oposite
#' ## between lme and lmer objects.
#'
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' residplot(fm, level=2)    #lme: level=2 for random effect "Block:Variety"
#'
#' \dontrun{
#' library(lme4)
#' fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#' residplot(fm) # lmer: By default level=1 for random effect "Block:Variety"
#' }
#'
#' @importFrom stats lm.influence ppoints predict residuals
#' @export
residplot <- function(model, group="none", level=1, slope=FALSE, id=FALSE, newwd=FALSE, ask=FALSE) {

  if (inherits(model, "aovlist")) {
    model <- aovlist_lmer(model)
  }
  if (newwd) {
    dev.new()
  }
  if (inherits(model, "lm") || inherits(model, "aov")) {
    op <- par(mfrow=c(2, 2), cex=0.7, mar=c(5, 5, 4, 2), mex=0.8, ask=ask)
    plot(model, cex.caption=0.8, which=1:4, col="blue")
    par(op)
  }
  if (inherits(model, "glm")) {  # the code below is from function 'glm.diag.plots' in package 'boot'
    glm.diag <- function (glmfit) {
      w <- if (is.null(glmfit$prior.weights)) {
        rep(1, length(glmfit$residuals))
      } else {
        glmfit$prior.weights
      }
      sd <- switch(family(glmfit)$family[1L], gaussian = sqrt(glmfit$deviance/glmfit$df.residual),
                   Gamma = sqrt(sum(w * (glmfit$y/fitted(glmfit) - 1)^2)/glmfit$df.residual),
                   1)
      dev <- residuals(glmfit, type = "deviance")/sd
      pear <- residuals(glmfit, type = "pearson")/sd
      h <- rep(0, length(w))
      h[w != 0] <- lm.influence(glmfit)$hat
      p <- glmfit$rank
      rp <- pear/sqrt(1 - h)
      rd <- dev/sqrt(1 - h)
      cook <- (h * rp^2)/((1 - h) * p)
      res <- sign(dev) * sqrt(dev^2 + h * rp^2)
      list(res = res, rd = rd, rp = rp, cook = cook, h = h, sd = sd)
    }
    glmdiag <- glm.diag(model)
    subset <- seq_along(glmdiag$h)
    par(mfrow=c(2, 2), cex=0.7, mar=c(5, 5, 4, 2), mex=0.8)
    x1 <- predict(model)
    plot(x1, glmdiag$res, col="blue", xlab = "Linear predictor", ylab = "Residuals")
    pars <- vector(4L, mode = "list")
    pars[[1L]] <- par("usr")
    y2 <- glmdiag$rd
    x2 <- qnorm(ppoints(length(y2)))[rank(y2)]
    plot(x2, y2, col="blue", ylab = "Quantiles of standard normal", xlab = "Ordered deviance residuals")
    abline(0, 1, lty = 2)
    pars[[2L]] <- par("usr")
    hh <- glmdiag$h/(1 - glmdiag$h)
    plot(hh, glmdiag$cook, col="blue", xlab = "h/(1-h)", ylab = "Cook statistic")
    rx <- range(hh)
    ry <- range(glmdiag$cook)
    rank.fit <- model$rank
    nobs <- rank.fit + model$df.residual
    cooky <- 8/(nobs - 2 * rank.fit)
    hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
    if ((cooky >= ry[1L]) && (cooky <= ry[2L])) {
      abline(h = cooky, lty = 2)
    }
    if ((hy >= rx[1L]) && (hy <= rx[2L])) {
      abline(v = hy, lty = 2)
    }
    pars[[3L]] <- par("usr")
    plot(subset, glmdiag$cook, col="blue", type="h", xlab = "Case", ylab = "Cook statistic")
    if ((cooky >= ry[1L]) && (cooky <= ry[2L])) {
      abline(h = cooky, lty = 2)
    }
    xx <- list(x1, x2, hh, subset)
    yy <- list(glmdiag$res, y2, glmdiag$cook, glmdiag$cook)
    pars[[4L]] <- par("usr")
    labels <- names(x1)
    while (id) {
      cat("****************************************************\n")
      cat("Please Input a screen number (1,2,3 or 4)\n")
      cat("0 will terminate the function \n")
      num <- as.numeric(readline())
      if ((length(num) > 0L) && ((num == 1) || (num == 2) ||
                                 (num == 3) || (num == 4))) {
        cat(paste("Interactive Identification for screen",
                  num, "\n"))
        cat("left button = Identify, center button = Exit\n")
        nm <- num + 1
        par(mfg = c(trunc(nm/2), 1 + nm%%2, 2, 2))
        par(usr = pars[[num]])
        identify(xx[[num]], yy[[num]], labels)
      } else {
        id <- FALSE
      }
    }
    par(mfrow = c(1, 1))
  }
  if (inherits(model, "gls")) {
    rsplot.gls(model, group, id, ask)
  }
  if (inherits(model, "lme") || inherits(model, "lmerMod")
      || inherits(model, "merModLmerTest") || inherits(model, "glmerMod")
      || inherits(model, "glmmTMB")) {
    rsplot.lme(model, group, level, slope, id, ask)
  }
}
