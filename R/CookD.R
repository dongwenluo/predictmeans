#' Cook's Distances of a Linear (Mixed) Model
#'
#' @description
#' Cook's distances of the observations (or groups of observations) of a model fitted by \code{aov}, \code{lm},
#'   \code{glm}, \code{gls}, \code{lme}, \code{lmer}, \code{glmer}, \code{glmmTMB} or \code{mgcv::gam}, with a plot.
#' @details
#' Cook's distance of an observation (or of a group of observations, \code{group}) is
#' \eqn{D = (\hat\beta - \hat\beta_{(-i)})' V^{-1} (\hat\beta - \hat\beta_{(-i)}) / p}, where
#' \eqn{\hat\beta_{(-i)}} are the fixed effects of the model refitted without it, \eqn{V} the covariance of
#' \eqn{\hat\beta} and \eqn{p} their number (for a \code{gam}, the effective degrees of freedom, as
#' \code{mgcv}'s \code{cooks.distance()}; the smoothing parameters are re-estimated for each deletion). For \code{lm} and \code{aov} models and single observations this is
#' the exact closed form of \code{stats::cooks.distance()}. Other models are refitted for each deletion: exact also
#' for \code{glm} (\code{stats::cooks.distance()} for a \code{glm} is a one-step approximation), and for
#' mixed models and \code{gls} fitted by maximum likelihood (the REML criterion changes with the fixed effects).
#' A deletion that leaves a coefficient without data gives \code{NA}. Refitting a generalized linear mixed model
#' for every observation takes time; \code{group} deletes whole clusters, usually the more useful check for
#' mixed models.
#' @usage
#' CookD(model, group=NULL, plot=TRUE, id=3, ...)
#' @param model Model object returned by \code{aov} (without \code{Error()}), \code{lm}, \code{glm}, \code{gls},
#'   \code{lme}, \code{lmer}, \code{glmer}, \code{glmmTMB} or \code{gam}.
#' @param group Optional name (in "quotes") of a variable whose levels are deleted one at a time (e.g. the
#'   clusters of a mixed model); the default \code{NULL} deletes single observations.
#' @param plot \code{TRUE} (the default): draw the distances (\pkg{ggplot2}, as the Cook's distance panel of
#'   \code{\link{residplot}}).
#' @param id The number of the largest Cook's distances labelled in the plot (\code{idn} before 2.0, now as in
#'   \code{\link{residplot}}). The default value is 3.
#' @param ... Not used: \code{idn}, \code{newwd} (before 2.0) stop with their 2.0 spelling.
#' @return Invisibly, the Cook's distances (named by observation or group).
#' @seealso \code{\link{residplot}}, \code{\link{pm_check}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#' CookD(fm)
#' CookD(fm, group = "Block")   # delete whole blocks
#' @name CookD
NULL

#' @export
CookD <- function (model, group=NULL, plot=TRUE, id=3, ...) {
  # arguments removed in 2.0 stop with their 2.0 spelling (R/defunct.R)
  pm_defunct("CookD", match.call(), strict = TRUE, formals_of = formals(sys.function()))

  if (inherits(model, "aovlist")) {
    stop("CookD() does not handle aov models with Error(): fit the error strata as random effects (lmer or lme).",
         call. = FALSE)
  }
  if (!inherits(model, c("lm", "gls", "lme", "lmerMod", "glmerMod", "glmmTMB"))) {
    stop("CookD() handles lm, aov, glm, gls, lme, lmer, glmer, glmmTMB and gam models.", call. = FALSE)
  }
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  # mixed models and gls by maximum likelihood (the REML criterion changes
  # with the fixed effects), each deletion refitted
  if (inherits(model, "gls") || inherits(model, "lme")) {
    mdf <- as.data.frame(pm_get_data(model))
    model <- pm_refit(model, data=mdf, method="ML", check=TRUE)
  } else if (inherits(model, "lmerMod")) {
    model <- pm_refit(model, REML=FALSE, check=TRUE)
    mdf <- model.frame(model)
  } else if (inherits(model, "glmerMod")) {
    mdf <- pm_model_rows(model)
  } else if (inherits(model, "glmmTMB")) {
    if (isTRUE(model$modelInfo$REML)) model <- quiet(pm_refit(model, REML = FALSE, check = TRUE))
    mdf <- pm_model_rows(model)
  } else {
    mdf <- as.data.frame(pm_get_data(model))
  }
  mp <- mymodelparm(model, df = 0)
  beta0  <- mp$coef
  vcovb <- mp$vcov
  vb.inv <- solve(vcovb)

  refit_coef <- function(keep) {
    fit <- tryCatch(quiet(pm_refit(model, data = mdf[keep, , drop = FALSE])), error = function(e) NULL)
    b <- if (is.null(fit)) NULL else tryCatch(mymodelparm(fit, df = 0)$coef, error = function(e) NULL)
    # a coefficient that the deletion leaves without data: no distance
    out <- rep(NA_real_, length(beta0))
    if (!is.null(b)) out <- unname(b[names(beta0)])
    out
  }
  if (is.null(group) || group%in%c("NULL", "")) {
    rn <- rownames(mdf)
    if (inherits(model, "lm") && !inherits(model, "glm")) {
      # lm, aov: the exact closed form, (b - b(-i))' X'X (b - b(-i)) / (p s^2)
      CookD <- stats::cooks.distance(model)
      CookD[!is.finite(CookD)] <- NA
      names(CookD) <- rn
      LOObeta <- NULL
    } else {
      LOObeta <- vapply(seq_along(rn), function(i) refit_coef(seq_along(rn) != i), numeric(length(beta0)))
    }
  } else {
    rn <- unique(as.character(mdf[, group]))
    LOObeta <- vapply(rn, function(x) refit_coef(as.character(mdf[, group]) != x), numeric(length(beta0)))
  }

  if (!is.null(LOObeta)) {
    rK <- t(matrix(LOObeta, length(beta0)) - beta0)
    # p: the number of coefficients; for a gam the effective df (as mgcv's cooks.distance())
    p <- if (inherits(model, "gam")) sum(model$edf) else length(beta0)
    CookD <- rowSums((rK %*% vb.inv) * rK)/p
  }
  names(CookD) <- rn
  if (plot) print(pm_cookd_plot(CookD, id, group))
  return(invisible(CookD))
}

## The plot of CookD(): a spike per observation or group, the 'id' largest
## labelled, as the Cook's distance panel of residplot().
pm_cookd_plot <- function(cd, id = 3, group = NULL) {
  grp <- !(is.null(group) || group %in% c("NULL", ""))
  d <- data.frame(case = seq_along(cd), cook = unname(cd), label = names(cd))
  ok <- is.finite(d$cook)
  k <- if (isTRUE(id)) 3L else if (isFALSE(id) || is.null(id)) 0L else as.integer(id)
  top <- if (k > 0 && any(ok)) utils::head(order(d$cook, decreasing = TRUE, na.last = NA), k) else integer(0)
  col <- "steelblue"
  p <- ggplot2::ggplot(d[ok, ], ggplot2::aes(x = .data$case, y = .data$cook)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$case, yend = 0), colour = col) +
    ggplot2::geom_point(data = d[top, , drop = FALSE], colour = col, size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), data = d[top, , drop = FALSE], vjust = -0.6, size = 3,
                       colour = "grey20") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = "Cook's distance", x = if (grp) paste("Level of", group) else "Observation",
                  y = "Cook's distance") +
    ggplot2::theme_bw()
  # few groups: name them on the axis
  if (grp && nrow(d) <= 30) p <- p + ggplot2::scale_x_continuous(breaks = d$case, labels = d$label)
  p
}
