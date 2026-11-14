#' Checks of a Fitted Model
#'
#' @description
#' A table of checks of a model's assumptions and structure, with a plain note on each: the fit itself (singular
#' or not converged), outliers and influential observations, normality of the residuals and of the random effects,
#' equal variances across the levels of each factor, overdispersion (Poisson and binomial models), interactions,
#' whether a factor with numeric levels follows a straight line, and coefficients that cannot be estimated. For any
#' model that \code{\link{predictmeans}} handles.
#' @usage pm_check(model, plot = FALSE)
#' @param model A fitted model, as for \code{\link{residplot}}.
#' @param plot \code{TRUE}: also draw the diagnostic plots of \code{\link{residplot}}.
#' @details
#' The checks inform; they do not decide. Rows marked \code{look = TRUE} are worth a look at the plots of
#' \code{\link{residplot}}. Keep the terms the design calls for (blocks, whole plots, the treatment structure):
#' dropping a term because its test is not significant makes the later tests and intervals too optimistic. Tests of
#' normality and of equal variances reject small departures in large data sets and miss large ones in small data
#' sets; the plots say how large a departure is.
#' \describe{
#' \item{Fit}{singular fits (a variance component estimated as zero) and convergence warnings.}
#' \item{Outliers}{residuals beyond 3 (of the residuals \code{\link{residplot}} uses: standardized, normalized or
#'   randomized quantile), against the number expected for normal residuals; worth a look when that many, or one as
#'   large as the largest (Bonferroni over the observations), has probability below 0.01.}
#' \item{Influence}{(\code{lm}, \code{aov}, \code{glm}) the largest Cook's distance; above 0.5 is worth a look
#'   (\code{\link{CookD}} for mixed models).}
#' \item{Normality}{Shapiro-Wilk tests of the residuals and of each term's random effects (intercepts).}
#' \item{Equal variances}{for each factor: the ratio of the largest to the smallest standard deviation of the
#'   residuals within its levels, and the Brown-Forsythe test (ANOVA of the absolute deviations from the medians).}
#' \item{Overdispersion}{Poisson and binomial (with more than one trial) models: Pearson chi-square / residual df,
#'   about 1 when the variance function is right.}
#' \item{Interaction}{the Type III test of each interaction term (as \code{permmodels(type = "III")}).}
#' \item{Straight line}{for a factor with numeric levels (e.g. nitrogen rates 0, 0.2, 0.4, 0.6) and at least three
#'   of them: the joint test of the polynomial contrasts above the linear one (with the actual spacing of the
#'   levels), from \code{\link{pm_joint}}.}
#' \item{Estimability}{coefficients dropped as aliased (empty cells or confounded factors).}
#' }
#' @return A data frame of class \code{"pm_check"}: \code{check}, \code{result}, \code{p.value}, \code{look} and
#'   \code{note}; printed as a table with the notes below it.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#' pm_check(fm)
#'
#' # counts: overdispersion
#' gm <- glm(breaks ~ wool * tension, family = poisson, data = warpbreaks)
#' pm_check(gm)
#' @name pm_check
NULL

pm_check <- function(model, plot = FALSE) {
  if (inherits(model, "semireg")) model <- model$semer
  out <- list()
  add <- function(check, result, p = NA_real_, note = "", look = FALSE) {
    out[[length(out) + 1]] <<- data.frame(check = check, result = result, p.value = p, look = look, note = note,
                                          stringsAsFactors = FALSE)
  }
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  d <- quiet(pm_resid_data(model))
  r <- d$resid
  ok <- is.finite(r)
  n <- sum(ok)

  ## the fit itself
  fit_issue <- pm_check_fit(model)
  if (length(fit_issue)) add("Fit", fit_issue, note = paste("Standard errors and df treat the estimates as known;",
                                                            "consider simplifying the random effects."), look = TRUE)
  else add("Fit", "no warnings")

  ## outliers
  big <- which(ok & abs(r) > 3)
  expct <- n * 2 * stats::pnorm(-3)
  ord <- big[order(abs(r[big]), decreasing = TRUE)]
  shown <- utils::head(ord, 3)
  # worth a look when so many beyond 3, or the largest so large, is unlikely
  # for normal residuals (p < 0.01; the largest: Bonferroni over n)
  pmax_r <- if (n) 1 - (1 - 2 * stats::pnorm(-max(abs(r[ok]))))^n else 1
  p_count <- stats::ppois(length(big) - 1, expct, lower.tail = FALSE)
  add("Outliers (|residual| > 3)",
      paste0(length(big), " of ", n, " (expected ", format(round(expct, 1), nsmall = 1), ")",
             if (length(shown)) paste0(": ", paste0(d$label[shown], " (", round(r[shown], 1), ")", collapse = ", ")),
             if (length(ord) > 3) ", ..."),
      note = if (length(big)) "Check these observations (data entry, unusual units) before anything else." else "",
      look = length(big) > 0 && (p_count < 0.01 || pmax_r < 0.01))

  ## influence (models with closed-form Cook's distances)
  if (inherits(model, "lm") && !inherits(model, c("aovlist", "gam"))) {
    cd <- quiet(stats::cooks.distance(model))
    i <- which.max(cd)
    add("Influence (Cook's distance)", paste0("largest ", round(cd[i], 3), " (", d$label[i], ")"),
        note = if (cd[i] > 0.5) "An observation changes the estimates noticeably: refit without it to see how much."
        else "", look = cd[i] > 0.5)
  }

  ## normality of the residuals and of the random effects
  if (n >= 3) {
    rs <- r[ok]
    if (n > 5000) rs <- rs[round(seq(1, n, length.out = 5000))]
    sw <- stats::shapiro.test(rs)$p.value
    add(paste0("Normality of residuals", if (identical(attr(d, "rlab"), "Randomized quantile residuals"))
      " (quantile)" else ""), "Shapiro-Wilk", sw,
      note = if (sw < 0.01) "See the Q-Q plot: skewness may call for a transformation or another distribution." else "",
      look = sw < 0.01)
  }
  for (j in seq_len(pm_check_nre(model))) {
    re <- quiet(pm_resid_ranef(model, j, FALSE))
    v <- re$value
    if (length(v) < 3 || length(unique(v)) < 3) next
    sw <- stats::shapiro.test(v)$p.value
    add(paste0("Normality of random effects: ", re$term), paste0("Shapiro-Wilk, ", length(v), " levels"), sw,
        note = if (length(v) < 10) "Few levels: little power." else "", look = sw < 0.01)
  }

  ## equal variances by factor
  frame <- pm_check_frame(model, length(r))
  if (!is.null(frame)) {
    fx <- pm_check_factors(model, frame)
    for (f in fx) {
      g <- factor(frame[[f]])[ok]
      if (nlevels(g) < 2 || any(table(g) < 2)) next
      sds <- tapply(r[ok], g, stats::sd)
      dev <- abs(r[ok] - stats::ave(r[ok], g, FUN = stats::median))
      # (no test when the deviations are equal within every level, e.g. two observations per level)
      bf <- if (stats::var(dev - stats::ave(dev, g)) < 1e-12 * max(1, stats::var(dev))) NA_real_ else
        tryCatch(quiet(stats::anova(stats::lm(dev ~ g))[1, "Pr(>F)"]), error = function(e) NA_real_)
      ratio <- max(sds) / min(sds)
      add(paste0("Equal variances: ", f),
          paste0("SD ratio ", round(ratio, 2), " (", names(which.max(sds)), " / ", names(which.min(sds)), ")"), bf,
          note = if (isTRUE(bf < 0.01)) paste0("The spread differs between levels of ", f, ": consider a variance ",
                                               "model (nlme: weights = varIdent(form = ~ 1 | ", f, ")) or a ",
                                               "transformation.") else "",
          look = isTRUE(bf < 0.01))
    }
  }

  ## overdispersion
  od <- pm_check_overdispersion(model)
  if (!is.null(od)) {
    add("Overdispersion (Pearson chi-square / df)", format(round(od$ratio, 2)), od$p,
        note = if (od$ratio > 1.5 && od$p < 0.01) paste("More variation than the distribution allows: SEs are too",
                                                       "small. Consider a negative binomial model or an",
                                                       "observation-level random effect.") else "",
        look = od$ratio > 1.5 && od$p < 0.01)
  }

  ## interactions and straight lines
  tl <- tryCatch(attr(stats::terms(if (inherits(model, "merMod")) stats::formula(model, fixed.only = TRUE)
                                   else model), "term.labels"), error = function(e) NULL)
  if (inherits(model, "aovlist")) tl <- NULL
  inter <- tl[grepl(":", tl)]
  if (length(inter)) {
    test <- if (inherits(model, c("lmerMod", "lm")) && !inherits(model, "glm")) "F" else
      if (inherits(model, "glm")) "LR" else "Chisq"
    a <- tryCatch(quiet(pm_anova(model, "III", test)), error = function(e) NULL)
    pcol <- if (!is.null(a)) grep("^Pr", colnames(a), value = TRUE)[1]
    for (it in inter) {
      p <- if (!is.null(a) && it %in% rownames(a)) a[it, pcol] else NA_real_
      add(paste0("Interaction: ", it), "Type III test", p,
          note = if (is.na(p)) "" else if (p < 0.05) paste("The effect of one factor depends on the other: compare",
                                                          "the levels of one within each level of the other (atvar",
                                                          "or pm_joint(by = ...)).")
          else "No evidence of interaction; keep the term if the design calls for it.")
    }
  }
  for (f in pm_check_numeric_factors(model, tl, frame)) {
    lev <- levels(factor(frame[[f]]))
    x <- as.numeric(lev)
    P <- stats::poly(x, length(lev) - 1)
    C <- t(P[, -1, drop = FALSE])
    colnames(C) <- lev
    jt <- tryCatch(quiet(pm_joint(model, f, ctrmatrix = C)), error = function(e) NULL)
    if (is.null(jt)) next
    p <- jt$p.value[1]
    add(paste0("Straight line in ", f), paste0("curvature, ", length(lev) - 2, " df"), p,
        note = if (p >= 0.05) paste0("The means of ", f, " are near a straight line: ", f, " could also enter as a ",
                                     "number (a response line); as a factor it compares the levels.")
        else paste0("The means of ", f, " are not on a straight line."))
  }

  ## estimability
  al <- pm_check_aliased(model)
  if (al > 0) add("Estimability", paste(al, "coefficient(s) aliased"),
                  note = "Empty cells or confounded factors: see the 'estimable' column of predictmeans() means.",
                  look = TRUE)

  res <- do.call(rbind, out)
  attr(res, "model") <- class(model)[1]
  attr(res, "n") <- n
  class(res) <- c("pm_check", "data.frame")
  if (plot) residplot(model)
  res
}

## print() of pm_check(): the table, then the notes.
print.pm_check <- function(x, ...) {
  cat("Checks of a ", attr(x, "model"), " model (", attr(x, "n"), " observations)\n\n", sep = "")
  tab <- data.frame(" " = ifelse(x$look, "!", ""), check = x$check, result = x$result,
                    p.value = ifelse(is.na(x$p.value), "", format.pval(x$p.value, digits = 2, eps = 1e-4)),
                    check.names = FALSE)
  print(tab, row.names = FALSE, right = FALSE)
  notes <- which(nzchar(x$note))
  if (length(notes)) {
    cat("\n")
    for (i in notes) cat(strwrap(paste0(if (x$look[i]) "! " else "- ", x$check[i], ": ", x$note[i]), exdent = 4),
                         sep = "\n")
  }
  cat("\nChecks inform, they do not decide; see residplot() for the plots.\n")
  invisible(x)
}

## Warnings about the fit: singular, not converged.
pm_check_fit <- function(model) {
  out <- character(0)
  if (inherits(model, "merMod")) {
    if (isTRUE(pm_singular(model))) out <- c(out, "singular (a variance component is zero)")
    msg <- model@optinfo$conv$lme4$messages
    if (length(msg)) out <- c(out, paste("convergence:", msg[1]))
  } else if (inherits(model, "glmmTMB")) {
    if (isFALSE(model$sdr$pdHess)) out <- c(out, "Hessian not positive definite")
  } else if (inherits(model, "lme")) {
    if (is.character(model$apVar)) out <- c(out, "variance parameters: approximate covariance not positive definite")
  }
  paste(out, collapse = "; ")[length(out) > 0]
}

## Number of random-effect terms.
pm_check_nre <- function(model) {
  re <- tryCatch(if (inherits(model, "lme")) {
    r <- nlme::ranef(model); if (is.data.frame(r)) 1L else length(r)
  } else if (inherits(model, "merMod")) length(lme4::ranef(model))
  else if (inherits(model, "glmmTMB")) length(lme4::ranef(model)[["cond"]]) else 0L, error = function(e) 0L)
  re
}

## The data rows of the residuals (NULL if they cannot be matched).
pm_check_frame <- function(model, n) {
  dat <- tryCatch({
    if (inherits(model, c("lme", "gls"))) {
      d <- pm_get_data(model)
      v <- intersect(all.vars(stats::formula(model)), names(d))
      d[stats::complete.cases(d[v]), , drop = FALSE]
    } else stats::model.frame(model)
  }, error = function(e) NULL)
  if (is.null(dat) || nrow(dat) != n) NULL else dat
}

## The factors of the fixed effects (variables in the frame).
pm_check_factors <- function(model, frame) {
  vars <- tryCatch(all.vars(stats::delete.response(stats::terms(
    if (inherits(model, "merMod")) stats::formula(model, fixed.only = TRUE) else model))), error = function(e) character(0))
  vars <- intersect(vars, names(frame))
  vars[vapply(vars, function(v) is.factor(frame[[v]]) || is.character(frame[[v]]), logical(1))]
}

## Factors that are fixed main-effect terms with at least three numeric levels.
pm_check_numeric_factors <- function(model, tl, frame) {
  if (is.null(frame) || !length(tl)) return(character(0))
  f <- intersect(pm_check_factors(model, frame), tl)
  f[vapply(f, function(v) {
    lev <- levels(factor(frame[[v]]))
    length(lev) >= 3 && !anyNA(suppressWarnings(as.numeric(lev)))
  }, logical(1))]
}

## Pearson chi-square / residual df for Poisson and binomial (trials > 1) models.
pm_check_overdispersion <- function(model) {
  fam <- pm_qresid_family(model)
  if (is.null(fam) || fam == "negbin") return(NULL)
  if (fam == "binomial") {
    w <- tryCatch(stats::weights(model), error = function(e) NULL)
    if (inherits(model, "glmmTMB")) {
      yy <- stats::model.response(stats::model.frame(model))
      w <- if (is.matrix(yy)) rowSums(yy) else NULL
    }
    if (is.null(w) || all(w <= 1)) return(NULL)
  }
  pr <- tryCatch(stats::residuals(model, type = "pearson"), error = function(e) NULL)
  dfr <- tryCatch(stats::df.residual(model), error = function(e) NA)
  if (is.null(pr) || !is.finite(dfr) || dfr <= 0) return(NULL)
  x2 <- sum(pr^2, na.rm = TRUE)
  list(ratio = x2 / dfr, p = stats::pchisq(x2, dfr, lower.tail = FALSE))
}

## Number of aliased (dropped) fixed-effect coefficients.
pm_check_aliased <- function(model) {
  if (inherits(model, "merMod")) {
    dr <- attr(lme4::getME(model, "X"), "col.dropped")
    return(length(dr))
  }
  b <- tryCatch(pm_coef(model), error = function(e) NULL)
  if (is.null(b)) 0L else sum(is.na(b))
}
