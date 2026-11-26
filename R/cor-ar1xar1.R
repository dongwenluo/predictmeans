## Separable AR1 x AR1 correlation for field trials: an nlme corStruct
## (dev/review-ar1ar1.md, D1-D7). corr(plot i, plot j) = phi_r^|dr| phi_c^|dc|
## (times 1 - nugget off the diagonal), from the integer row and column of
## each plot over the whole field; missing plots need nothing (the matrix of
## the plots present is the sub-matrix). Parameters on the optimiser's scale:
## atanh(phi_r), atanh(phi_c), qlogis(nugget).

#' Separable AR1 x AR1 correlation for field trials
#'
#' A correlation structure for \code{nlme::gls()} and \code{nlme::lme()}: the
#' correlation of two plots is \eqn{\phi_r^{|dr|}\phi_c^{|dc|}}, with \eqn{dr} and
#' \eqn{dc} their distances in rows and columns (Gilmour, Cullis and Verbyla 1997),
#' as Genstat's \code{VSTRUCTURE AR,AR} and ASReml's \code{ar1(row):ar1(col)}.
#'
#' @param value Starting values of \eqn{\phi_r} and \eqn{\phi_c} (in (-1, 1)), and of
#'   the nugget share (in (0, 1)) when \code{nugget = TRUE}.
#' @param form \code{~ row + col} or \code{~ row + col | group}: integer row and
#'   column positions of the plots over the whole field (one AR1 x AR1 per group,
#'   the same \eqn{\phi}).
#' @param nugget An independent plot error besides the AR1 x AR1 part: the
#'   correlation is \eqn{(1 - n)\phi_r^{|dr|}\phi_c^{|dc|}} off the diagonal, n the
#'   share of the plot variance.
#' @param fixed Keep the starting values (not estimated)?
#' @details
#' The row and column are integer positions over the whole field, so plots
#' missing from the grid need nothing: the correlation of the plots present is
#' the sub-matrix. In \code{gls()} the structure covers all plots (or each
#' group). In \code{lme()} a correlation acts only within the groups of the
#' innermost random effect: with \code{random = ~ 1 | rep} the AR1 x AR1 would be
#' within each replicate. For replicate (block) effects plus an AR1 x AR1 over
#' the whole field, as Genstat and ASReml fit it, put the random effects in one
#' group: \code{lme(y ~ variety, random = list(all = pdIdent(~ rep - 1)),
#' correlation = corAR1xAR1(form = ~ row + col | all))} with \code{all = factor(1)}.
#'
#' With \code{nugget = TRUE} the plot variance is split into an AR1 x AR1 part and
#' an independent part, the nugget share \eqn{n}; Genstat's \code{units} term
#' with the AR1 x AR1 term scaled by a ratio \eqn{g} is the same model with
#' \eqn{n = 1/(1 + g)} (with AR1 x AR1 as the residual and units as a random term,
#' \eqn{n = g/(1 + g)}). Each likelihood uses a dense matrix of the plots of a
#' group: fine to about 2000 plots.
#'
#' A smooth field trend is the other common spatial model:
#' \code{semireg(y ~ variety + t2(row, col) + (1 | rowf) + (1 | colf), data = d)}, a
#' smooth surface over the field as SpATS fits, often as good as AR1 x AR1 in
#' practice (on Slate Hall the AR1 x AR1 means are the more precise: average SED
#' 0.590 against 0.668). \code{\link{residvariogram}()} shows whether a fit leaves
#' spatial correlation.
#'
#' Checked against Genstat (REML Guide's Slate Hall trial): \eqn{\phi}, variances,
#' the 25 variety means and all 300 SEDs to 1e-6, the Wald statistic, and the
#' d.d.f. of the variety test with \code{df_method = "KR-average"}.
#' @return An object of class \code{c("corAR1xAR1", "corStruct")} for the
#'   \code{correlation} argument; \code{coef(fit$modelStruct$corStruct,
#'   unconstrained = FALSE)} gives \eqn{\phi_r}, \eqn{\phi_c} (and the nugget share).
#' @references
#' Gilmour, A. R., Cullis, B. R. and Verbyla, A. P. (1997). Accounting for natural and extraneous variation in
#' the analysis of field experiments. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 2,
#' 269-293. \doi{10.2307/1400446}
#' @seealso \code{\link{residvariogram}}, \code{\link{predictmeans}}
#' @examples
#' \donttest{
#' d <- read.csv(system.file("extdata", "slatehall.csv", package = "predictmeans"))
#' d$variety <- factor(d$variety)
#' fit <- nlme::gls(yield ~ variety, data = d, correlation = corAR1xAR1(form = ~ fieldrow + fieldcolumn))
#' coef(fit$modelStruct$corStruct, unconstrained = FALSE)    # Genstat: 0.4586, 0.6838
#' pm <- predictmeans(fit, "variety")
#' head(pm$means)
#' residvariogram(fit, row = "fieldrow", col = "fieldcolumn", data = d, nsim = 19)   # a few seconds
#' }
#' @export
corAR1xAR1 <- function(value = c(0.5, 0.5), form = ~ row + col, nugget = FALSE, fixed = FALSE) {
  if (!inherits(form, "formula")) stop("'form' must be a formula, e.g. ~ row + col.", call. = FALSE)
  np <- 2L + nugget
  if (length(value) == 2L && nugget) value <- c(value, 0.1)
  if (length(value) != np) stop("'value' must have ", np, " numbers.", call. = FALSE)
  if (any(abs(value[1:2]) >= 1)) stop("the starting phi must be in (-1, 1).", call. = FALSE)
  if (nugget && (value[3] <= 0 || value[3] >= 1)) stop("the starting nugget share must be in (0, 1).", call. = FALSE)
  par <- c(atanh(value[1:2]), if (nugget) stats::qlogis(value[3]))
  attr(par, "formula") <- form
  attr(par, "fixed") <- isTRUE(fixed)
  attr(par, "nugget") <- isTRUE(nugget)
  class(par) <- c("corAR1xAR1", "corStruct")
  par
}

## phi_r, phi_c and the nugget share from the optimiser's scale.
pm_ar1x_par <- function(object) {
  v <- as.vector(unclass(object))
  list(phi = tanh(v[1:2]), nugget = if (isTRUE(attr(object, "nugget"))) stats::plogis(v[3]) else 0)
}

## The correlation matrix of the plots at rows r, columns c.
pm_ar1x_mat <- function(rc, p) {
  R <- p$phi[1]^abs(outer(rc[, 1], rc[, 1], "-")) * p$phi[2]^abs(outer(rc[, 2], rc[, 2], "-"))
  if (p$nugget > 0) {
    R <- (1 - p$nugget) * R
    diag(R) <- 1
  }
  R
}

#' @exportS3Method nlme::Initialize
Initialize.corAR1xAR1 <- function(object, data, ...) {
  form <- nlme::getCovariateFormula(object)
  vars <- all.vars(form)
  if (length(vars) != 2) stop("corAR1xAR1(): 'form' must name two variables, the row and the column.", call. = FALSE)
  rc <- as.matrix(data[, vars])
  if (!is.numeric(rc) || any(abs(rc - round(rc)) > 1e-8, na.rm = TRUE)) {
    stop("corAR1xAR1(): the row and column must be integer positions (numbers, not factors).", call. = FALSE)
  }
  object <- utils::getFromNamespace("Initialize.corStruct", "nlme")(object, data)
  grp <- attr(object, "groups")
  # (one matrix without groups, a list by group otherwise, as nlme's classes)
  cov <- if (is.null(grp)) rc else lapply(split(seq_len(nrow(rc)), grp), function(i) rc[i, , drop = FALSE])
  if (max(if (is.list(cov)) vapply(cov, nrow, integer(1)) else nrow(cov)) > 2000) {
    message("corAR1xAR1(): over 2000 plots in a group: each likelihood uses a dense matrix of that size (slow); ",
            "a smooth field trend, semireg(y ~ ... + t2(row, col)), is faster.")
  }
  attr(object, "covariate") <- cov
  attr(object, "factor") <- nlme::corFactor(object)
  attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")
  object
}

#' @exportS3Method nlme::corMatrix
corMatrix.corAR1xAR1 <- function(object, covariate = nlme::getCovariate(object), corr = TRUE, ...) {
  p <- pm_ar1x_par(object)
  one <- function(rc) {
    R <- pm_ar1x_mat(rc, p)
    if (corr) return(R)
    U <- chol(R)
    val <- t(backsolve(U, diag(nrow(R))))
    attr(val, "logDet") <- -sum(log(diag(U)))
    val
  }
  if (!is.list(covariate)) return(one(covariate))
  val <- lapply(covariate, one)
  if (!corr) attr(val, "logDet") <- sum(vapply(val, function(v) attr(v, "logDet"), numeric(1)))
  val
}

#' @exportS3Method nlme::corFactor
corFactor.corAR1xAR1 <- function(object, ...) {
  val <- corMatrix(object, corr = FALSE)
  lD <- attr(val, "logDet")
  if (is.list(val)) val <- unlist(lapply(val, as.vector))
  else val <- as.vector(val)
  attr(val, "logDet") <- lD
  val
}

#' @exportS3Method stats::coef
coef.corAR1xAR1 <- function(object, unconstrained = TRUE, ...) {
  if (unconstrained) {
    if (isTRUE(attr(object, "fixed"))) return(numeric(0))
    v <- as.vector(unclass(object))
    names(v) <- c("atanh(phi_row)", "atanh(phi_col)", if (isTRUE(attr(object, "nugget"))) "logit(nugget)")
    return(v)
  }
  p <- pm_ar1x_par(object)
  c(phi_row = p$phi[1], phi_col = p$phi[2], if (isTRUE(attr(object, "nugget"))) c(nugget = p$nugget))
}

#' @exportS3Method nlme::`coef<-`
`coef<-.corAR1xAR1` <- function(object, ..., value) {
  if (length(value) != length(object)) stop("cannot change the length of the parameter of a corAR1xAR1 object.", call. = FALSE)
  object[] <- value
  attr(object, "factor") <- NULL
  attr(object, "factor") <- nlme::corFactor(object)
  attr(object, "logDet") <- NULL
  attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")
  object
}
