#' Checks of a Fitted Model
#'
#' @description
#' A table of checks of a model's assumptions and structure, grouped in sections, with a plain note on each: the fit
#' itself (singular or not converged); the residuals (outliers, influence, normality, a transformation); the
#' variance (equal across factor levels, constant with the fitted values); the distribution of generalized models
#' (over- and underdispersion, excess zeros, separation); the design (balance and empty cells, aliasing,
#' collinearity, few levels of a random effect); and the terms (interactions, straight lines). For any model that
#' \code{\link{predictmeans}} handles.
#' @usage pm_check(model, plot = FALSE)
#' @param model A fitted model, as for \code{\link{residplot}}.
#' @param plot \code{TRUE}: also draw the diagnostic plots of \code{\link{residplot}}.
#' @details
#' The checks inform; they do not decide. Rows marked \code{!} are worth a look at the plots of
#' \code{\link{residplot}}; rows marked \code{i} are information. Keep the terms the design calls for (blocks, whole
#' plots, the treatment structure): dropping a term because its test is not significant makes the later tests and
#' intervals too optimistic. Tests of normality and of equal variances reject small departures in large data sets
#' and miss large ones in small data sets; the plots say how large a departure is.
#' \describe{
#' \item{Fit}{singular fits (a variance component estimated as zero) and convergence warnings.}
#' \item{Outliers}{residuals beyond 3 (of the residuals \code{\link{residplot}} uses: standardized, normalized or
#'   randomized quantile), against the number expected for normal residuals; worth a look when that many, or one as
#'   large as the largest (Bonferroni over the observations), has probability below 0.01.}
#' \item{Influence}{(\code{lm}, \code{aov}, \code{glm}) the largest Cook's distance; above 0.5 is worth a look
#'   (\code{\link{CookD}} for mixed models).}
#' \item{Normality}{Shapiro-Wilk tests of the residuals and of each term's random effects (intercepts).}
#' \item{Transformation of the response}{(\code{lm} and \code{aov} with a positive response that is not
#'   transformed in the formula) the Box-Cox power the residuals favour: the profile likelihood estimate of
#'   \eqn{\lambda} and its 95\% interval (Box and Cox, 1964; as \code{MASS::boxcox()}). Worth a look when the
#'   interval leaves out 1 (no transformation); the note names the usual transformation within the interval.}
#' \item{Equal variances}{for each factor: the ratio of the largest to the smallest standard deviation of the
#'   residuals within its levels, and the Brown-Forsythe test (ANOVA of the absolute deviations from the medians).}
#' \item{Spread vs fitted values}{(gaussian models) Spearman correlation of the absolute residuals with the fitted
#'   values; worth a look when \eqn{|\rho| \ge 0.2} and p < 0.01 (the variance grows with the mean).}
#' \item{Dispersion}{(Poisson, binomial with more than one trial, negative binomial; \code{glm}, \code{glmer},
#'   \code{glmmTMB}) Pearson chi-square / residual df, about 1 when the variance function is right; overdispersion
#'   (ratio > 1.5, p < 0.01) makes SEs too small, underdispersion (ratio < 0.7, p < 0.01) too large. For mixed
#'   models the residual df counts the variance parameters, an approximation. A 0/1 response carries no information
#'   on dispersion; for families with an estimated dispersion (Gamma, quasi) its estimate is shown.}
#' \item{Zero counts}{(Poisson, negative binomial) the observed zeros against those simulated from the fitted model
#'   (200 simulations, fixed seed; the user's random numbers are not changed); worth a look when there are more
#'   (p < 0.01).}
#' \item{Separation}{(binomial) fitted probabilities of 0 or 1, or coefficients beyond 10 on the logit scale or with
#'   standard errors above 50: the estimates and Wald tests of these terms are unreliable.}
#' \item{Balance of cells}{of the crossed fixed factors: balanced, unbalanced or with empty cells; unbalanced data
#'   make Type I, II and III tests differ (see \code{\link{anova.lmerMod}}, section 'Which type of test').}
#' \item{Estimability}{coefficients dropped as aliased (empty cells or confounded factors).}
#' \item{Collinearity}{(models of main effects only) the generalized variance inflation factor of each term,
#'   \eqn{GVIF^{1/(2 df)}} (Fox and Monette, 1992); above \eqn{\sqrt{10}} is worth a look.}
#' \item{Random-effect levels}{a grouping factor with fewer than 5 levels: its variance is poorly estimated.}
#' \item{Interaction}{the Type III test of each interaction term (as \code{permmodels(type = "III")}).}
#' \item{Straight line}{for a factor with numeric levels and at least three of them: the joint test of the
#'   polynomial contrasts above the linear one (with the actual spacing of the levels), from \code{\link{pm_joint}}.}
#' }
#' @return A data frame of class \code{"pm_check"}: \code{section}, \code{check}, \code{result}, \code{p.value},
#'   \code{look}, \code{status} (\code{"!"} worth a look, \code{"i"} information, \code{"OK"}) and \code{note};
#'   printed by section with the notes below. \code{summary()} keeps only the rows worth a look. The rows worth a
#'   look are also kept as advice (\code{\link{advice}}).
#' @references
#' Fox, J. and Monette, G. (1992). Generalized collinearity diagnostics. \emph{Journal of the American Statistical
#' Association}, 87, 178-183.
#' @seealso \code{\link{residplot}}, \code{\link{CookD}}, \code{\link{advice}}, \code{\link{pm_boxcox}}
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
#' summary(pm_check(gm))
#' @name pm_check
NULL

## The order of the sections of pm_check().
pm_check_sections <- c("Fit", "Residuals", "Variance", "Distribution", "Design", "Terms")

#' @export
pm_check <- function(model, plot = FALSE) {
  if (inherits(model, "semireg")) model <- model$semer
  out <- list()
  add <- function(section, check, result, p = NA_real_, note = "", look = FALSE, info = FALSE) {
    out[[length(out) + 1]] <<- data.frame(section = section, check = check, result = result, p.value = p,
                                          look = look, status = if (look) "!" else if (info) "i" else "OK",
                                          note = note, stringsAsFactors = FALSE)
  }
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  d <- quiet(pm_resid_data(model))
  r <- d$resid
  ok <- is.finite(r)
  n <- sum(ok)

  ## the fit itself
  fit_issue <- pm_check_fit(model)
  if (length(fit_issue)) add("Fit", "Fit", fit_issue, note = paste("Standard errors and df treat the estimates as",
                                                                  "known; consider simplifying the random effects."),
                             look = TRUE)
  else add("Fit", "Fit", "no warnings")

  ## outliers
  big <- which(ok & abs(r) > 3)
  expct <- n * 2 * stats::pnorm(-3)
  ord <- big[order(abs(r[big]), decreasing = TRUE)]
  shown <- utils::head(ord, 3)
  # worth a look when so many beyond 3, or the largest so large, is unlikely
  # for normal residuals (p < 0.01; the largest: Bonferroni over n)
  pmax_r <- if (n) 1 - (1 - 2 * stats::pnorm(-max(abs(r[ok]))))^n else 1
  p_count <- stats::ppois(length(big) - 1, expct, lower.tail = FALSE)
  add("Residuals", "Outliers (|residual| > 3)",
      paste0(length(big), " of ", n, " (expected ", format(round(expct, 1), nsmall = 1), ")",
             if (length(shown)) paste0(": ", paste0(d$label[shown], " (", round(r[shown], 1), ")", collapse = ", ")),
             if (length(ord) > 3) ", ..."),
      note = if (length(big)) "Check these observations (data entry, unusual units) before anything else." else "",
      look = length(big) > 0 && (p_count < 0.01 || pmax_r < 0.01))

  ## influence (models with closed-form Cook's distances)
  if (inherits(model, "lm") && !inherits(model, c("aovlist", "gam"))) {
    cd <- quiet(stats::cooks.distance(model))
    i <- which.max(cd)
    add("Residuals", "Influence (Cook's distance)", paste0("largest ", round(cd[i], 3), " (", d$label[i], ")"),
        note = if (cd[i] > 0.5) "An observation changes the estimates noticeably: refit without it to see how much."
        else "", look = cd[i] > 0.5)
  }

  ## normality of the residuals and of the random effects
  if (n >= 3) {
    rs <- r[ok]
    if (n > 5000) rs <- rs[round(seq(1, n, length.out = 5000))]
    sw <- stats::shapiro.test(rs)$p.value
    add("Residuals", paste0("Normality of residuals", if (identical(attr(d, "rlab"), "Randomized quantile residuals"))
      " (quantile)" else ""), "Shapiro-Wilk", sw,
      note = if (sw < 0.01) "See the Q-Q plot: skewness may call for a transformation or another distribution." else "",
      look = sw < 0.01)
  }
  # the Box-Cox power the residuals suggest (lm with a positive, untransformed response)
  bc <- pm_boxcox_profile(model)
  if (!is.null(bc)) {
    worth <- bc$upper < 1 || bc$lower > 1                      # the interval leaves out "no transformation"
    add("Residuals", "Transformation of the response (Box-Cox)",
        paste0("lambda ", format(round(bc$lambda, 2)), " (95% interval ", format(round(bc$lower, 2)), " to ",
               format(round(bc$upper, 2)), if (bc$edge) ", at the edge of -3 to 3" else "", ")"),
        note = if (!worth) "" else if (!is.na(bc$pick)) {
          paste0("The residuals favour a transformation: ", bc$name, " is within the interval",
                 if (bc$pick != 0) paste0(" (pm_boxcox(y, ", bc$pick, ") in the formula)") else "",
                 ". For counts or proportions a generalized linear model is usually better.")
        } else {
          paste0("The residuals favour a transformation: pm_boxcox(y, ", format(round(bc$lambda, 2)),
                 ") in the formula. For counts or proportions a generalized linear model is usually better.")
        }, look = worth)
  }
  for (j in seq_len(pm_check_nre(model))) {
    re <- quiet(pm_resid_ranef(model, j, FALSE))
    v <- re$value
    if (length(v) < 3 || length(unique(v)) < 3) next
    sw <- stats::shapiro.test(v)$p.value
    add("Residuals", paste0("Normality of random effects: ", re$term), paste0("Shapiro-Wilk, ", length(v), " levels"),
        sw, note = if (length(v) < 10) "Few levels: little power." else "", look = sw < 0.01)
  }

  ## equal variances by factor, and with the fitted values
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
      add("Variance", paste0("Equal variances: ", f),
          paste0("SD ratio ", round(ratio, 2), " (", names(which.max(sds)), " / ", names(which.min(sds)), ")"), bf,
          note = if (isTRUE(bf < 0.01)) paste0("The spread differs between levels of ", f, ": consider a variance ",
                                               "model (nlme: weights = varIdent(form = ~ 1 | ", f, ")) or a ",
                                               "transformation.") else "",
          look = isTRUE(bf < 0.01))
    }
  }
  sf <- pm_check_spread_fitted(model, d, ok)
  if (!is.null(sf)) {
    worth <- abs(sf$rho) >= 0.2 && sf$p < 0.01
    add("Variance", "Spread vs fitted values", paste0("Spearman rho ", format(round(sf$rho, 2), nsmall = 2)), sf$p,
        note = if (!worth) "" else paste("The spread of the residuals", if (sf$rho > 0) "grows" else "shrinks",
                                         "with the fitted values: a log or square-root transformation, or a",
                                         "variance function (nlme: weights = varPower()), may fit better."),
        look = worth)
  }

  ## the distribution of generalized models
  od <- pm_check_dispersion(model)
  if (!is.null(od)) {
    if (od$kind == "binary") {
      add("Distribution", "Dispersion", "0/1 response: not identifiable", info = TRUE,
          note = paste("Overdispersion cannot be judged from 0/1 data; with repeated measures on units, a random",
                       "effect for the units allows for it."))
    } else if (od$kind == "estimated") {
      add("Distribution", "Dispersion (estimated)", paste0(od$family, ": ", format(signif(od$ratio, 3))), info = TRUE,
          note = "The dispersion is estimated by the model, so the SEs allow for it.")
    } else {
      over <- od$ratio > 1.5 && od$p_over < 0.01
      under <- od$ratio < 0.7 && od$p_under < 0.01
      add("Distribution", "Dispersion (Pearson chi-square / df)",
          paste0(format(round(od$ratio, 2), nsmall = 2), " on ", format(round(od$df, 1)), " df",
                 if (od$glmm) " (approx.)" else ""),
          if (od$ratio >= 1) od$p_over else od$p_under,
          note = if (over) paste("Overdispersion: more variation than the", od$family, "distribution allows, so SEs",
                                 "are too small and tests too liberal.",
                                 if (od$family == "negbin") "Check for a missing term or excess zeros." else
                                   "Consider a negative binomial model, an observation-level random effect, or a quasi family (glm).")
          else if (under) paste("Underdispersion: less variation than the distribution allows, so SEs are too large",
                                "and tests conservative. Consider a quasi family or a generalized Poisson model.")
          else "",
          look = over || under)
    }
  }
  zi <- pm_check_zeros(model)
  if (!is.null(zi)) {
    worth <- zi$p < 0.01
    add("Distribution", "Zero counts", paste0(zi$obs, " observed, ", format(round(zi$expected, 1)), " expected"), zi$p,
        note = if (!worth) "" else paste("More zeros than the model expects: check overdispersion first (a negative",
                                         "binomial model often accounts for them), then consider a zero-inflated or",
                                         "hurdle model (glmmTMB, ziformula)."),
        look = worth)
  }
  sp <- pm_check_separation(model)
  if (!is.null(sp)) {
    worth <- sp$n01 > 0 || length(sp$big) > 0
    add("Distribution", "Separation (binomial)",
        if (!worth) "none detected" else paste0(sp$n01, " fitted probabilities at 0 or 1",
                                                if (length(sp$big)) paste0("; extreme: ", paste(utils::head(sp$big, 3), collapse = ", "),
                                                                           if (length(sp$big) > 3) ", ...")),
        note = if (!worth) "" else paste("(Quasi-)complete separation: some groups have all or none of the events.",
                                         "Estimates and Wald tests of these terms are unreliable; use likelihood ratio",
                                         "tests, merge levels, or a bias-reduced fit (e.g. brglm2)."),
        look = worth)
  }

  ## the design
  bal <- tryCatch(pm_design_balance(model), error = function(e) NA_character_)
  if (!is.na(bal)) {
    cnt <- as.vector(attr(bal, "counts"))
    add("Design", "Balance of cells",
        switch(bal, balanced = paste0("balanced (", cnt[1], " per cell)"),
               unbalanced = paste0("unbalanced (", min(cnt), " to ", max(cnt), " per cell)"),
               empty = paste0(sum(cnt == 0), " empty of ", length(cnt), " cells")),
        note = switch(bal, balanced = "",
                      unbalanced = paste("Type I, II and III tests differ: Type III tests equally weighted cell means",
                                         "(meaningful with interactions); Type II is more powerful when interactions",
                                         "are negligible; Type I depends on the order of the terms. See",
                                         "?anova.lmerMod."),
                      empty = paste("No type of test answers a main-effect question over all cells; Type III uses",
                                    "the cells with data. Prefer comparisons of the cell means (pm_joint(),",
                                    "predictmeans(pairwise = TRUE)).")),
        info = bal != "balanced")
  }
  al <- pm_check_aliased(model)
  if (al > 0) add("Design", "Estimability", paste(al, "coefficient(s) aliased"),
                  note = "Empty cells or confounded factors: see the 'estimable' column of predictmeans() means.",
                  look = TRUE)
  gv <- tryCatch(pm_check_gvif(model), error = function(e) NULL)
  if (length(gv)) {
    worth <- max(gv) > sqrt(10)
    add("Design", "Collinearity (GVIF^(1/(2 df)))",
        paste0("largest ", format(round(max(gv), 2), nsmall = 2), " (", names(gv)[which.max(gv)], ")"),
        note = if (!worth) "" else paste("Terms are strongly correlated: each one's test, adjusted for the others,",
                                         "has little power and its estimates are unstable. Consider dropping or",
                                         "combining correlated covariates."),
        look = worth)
  }
  ng <- pm_check_re_levels(model)
  for (g in names(ng)[ng < 5]) {
    add("Design", paste0("Random-effect levels: ", g), paste(ng[[g]], "levels"), info = TRUE,
        note = paste("A variance from fewer than 5 levels is poorly estimated (often as zero); tests and intervals",
                     "that use it are approximate. Keep it if the design calls for it."))
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
      add("Terms", paste0("Interaction: ", it), "Type III test", p,
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
    add("Terms", paste0("Straight line in ", f), paste0("curvature, ", length(lev) - 2, " df"), p,
        note = if (p >= 0.05) paste0("The means of ", f, " are near a straight line: ", f, " could also enter as a ",
                                     "number (a response line); as a factor it compares the levels.")
        else paste0("The means of ", f, " are not on a straight line."))
  }

  res <- do.call(rbind, out)
  res <- res[order(match(res$section, pm_check_sections)), , drop = FALSE]
  row.names(res) <- NULL
  attr(res, "model") <- class(model)[1]
  attr(res, "n") <- n
  class(res) <- c("pm_check", "data.frame")
  res <- pm_check_advice(res)
  if (plot) residplot(model)
  res
}

## The notes of the rows worth a look as advice rows (R/advice.R): id
## "check_<check>" (lower case, non-alphanumerics as _), level "caution".
## print.pm_check() keeps its own layout (notes beside their rows).
pm_check_advice <- function(res) {
  k <- which(res$look & nzchar(res$note))
  if (!length(k)) return(res)
  ids <- paste0("check_", gsub("[^a-z0-9]+", "_", tolower(res$check[k])))
  ids <- make.unique(sub("_$", "", ids), sep = "_")
  do.call(pm_add_advice, c(list(res), Map(function(i, t) pm_advice(i, t, "caution"), ids,
                                         paste0(res$check[k], ": ", res$note[k]))))
}

## print() of pm_check(): the rows by section, the notes, a summary line.
#' @exportS3Method print pm_check
print.pm_check <- function(x, ...) {
  if (!all(c("section", "check", "result", "p.value", "look", "status", "note") %in% names(x))) {
    return(print(as.data.frame(unclass(x), stringsAsFactors = FALSE), ...))
  }
  cat("Checks of a ", attr(x, "model"), " model (", attr(x, "n"), " observations)\n", sep = "")
  if (!nrow(x)) {
    cat("\nNo rows.\n")
    return(invisible(x))
  }
  pv <- ifelse(is.na(x$p.value), "", format.pval(x$p.value, digits = 2, eps = 1e-4))
  w_chk <- min(max(nchar(x$check)), 42)
  w_p <- max(7, nchar(pv))
  w_res <- max(20, getOption("width") - w_chk - w_p - 9)
  clip <- function(s, w) ifelse(nchar(s) > w, paste0(substr(s, 1, w - 3), "..."), s)
  cut <- nchar(x$result) > w_res
  for (s in unique(x$section)) {
    cat("\n", s, "\n", sep = "")
    for (k in which(x$section == s)) {
      cat(sprintf("  %-2s %-*s  %-*s %*s", x$status[k], w_chk, clip(x$check[k], w_chk), w_res,
                  clip(x$result[k], w_res), w_p, pv[k]), "\n", sep = "")
    }
  }
  notes <- which(nzchar(x$note) | cut)
  if (length(notes)) {
    cat("\nNotes\n")
    for (i in notes) {
      txt <- paste0(if (x$status[i] == "OK") "-" else x$status[i], " ", x$check[i], ": ", if (cut[i]) paste0(x$result[i], ". "), x$note[i])
      cat(strwrap(txt, width = 0.95 * getOption("width"), exdent = 4), sep = "\n")
    }
  }
  lk <- x$check[x$look]
  cat("\n", if (length(lk)) paste0(length(lk), " worth a look: ", paste(lk, collapse = "; "), ".")
      else "Nothing flagged.", "\n", sep = "")
  cat("(! worth a look, i information, OK no concern.) Checks inform, they do not decide; see residplot().\n")
  invisible(x)
}

## summary() of pm_check(): only the rows worth a look.
#' @exportS3Method summary pm_check
summary.pm_check <- function(object, ...) {
  out <- object[object$look, , drop = FALSE]
  row.names(out) <- NULL
  attributes(out)[c("model", "n", "advice")] <- attributes(object)[c("model", "n", "advice")]
  class(out) <- class(object)
  out
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

## Number of levels of each grouping factor of the random effects (NULL if none).
pm_check_re_levels <- function(model) {
  tryCatch({
    if (inherits(model, "merMod")) lme4::ngrps(model)
    else if (inherits(model, "lme")) {
      g <- model$dims$ngrps
      g[setdiff(names(g), c("X", "y"))]
    } else if (inherits(model, "glmmTMB")) vapply(model$modelInfo$reTrms$cond$flist, nlevels, 1L)
  }, error = function(e) NULL)
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

## The family name of a model ("gaussian" for lm, lme, gls, aov).
pm_check_family <- function(model) {
  if (!inherits(model, c("glm", "glmerMod", "glmmTMB", "gam"))) return("gaussian")
  tryCatch(stats::family(model)$family, error = function(e) "")
}

## Spearman correlation of |residual| with the fitted values (gaussian models).
pm_check_spread_fitted <- function(model, d, ok) {
  if (pm_check_family(model) != "gaussian" || is.null(d$fitted)) return(NULL)
  f <- d$fitted[ok]
  a <- abs(d$resid[ok])
  if (length(unique(f)) < 3 || length(a) < 10) return(NULL)
  ct <- tryCatch(suppressWarnings(stats::cor.test(a, f, method = "spearman", exact = FALSE)), error = function(e) NULL)
  if (is.null(ct) || !is.finite(ct$estimate)) return(NULL)
  list(rho = unname(ct$estimate), p = ct$p.value)
}

## Dispersion of a generalized model (glm, glmer, glmmTMB). kind: "tested"
## (Poisson, binomial with trials, negative binomial: Pearson chi-square /
## residual df, both tails), "binary" (0/1 response), "estimated" (the
## family estimates it), or NULL (gaussian, other classes).
pm_check_dispersion <- function(model) {
  if (!inherits(model, c("glm", "glmerMod", "glmmTMB"))) return(NULL)
  fam <- pm_check_family(model)
  if (fam %in% c("gaussian", "")) return(NULL)
  kind <- pm_qresid_family(model)
  if (is.null(kind)) {
    phi <- tryCatch(if (inherits(model, "glm")) summary(model)$dispersion else stats::sigma(model)^2,
                    error = function(e) NA_real_)
    if (!is.finite(phi)) return(NULL)
    return(list(kind = "estimated", family = fam, ratio = phi))
  }
  if (kind == "binomial") {
    w <- tryCatch(stats::weights(model), error = function(e) NULL)
    if (inherits(model, "glmmTMB")) {
      yy <- stats::model.response(stats::model.frame(model))
      w <- if (is.matrix(yy)) rowSums(yy) else NULL
    }
    if (is.null(w) || all(w <= 1)) return(list(kind = "binary"))
  }
  pr <- tryCatch(stats::residuals(model, type = "pearson"), error = function(e) NULL)
  dfr <- tryCatch(stats::df.residual(model), error = function(e) NA)
  if (is.null(pr) || !is.finite(dfr) || dfr <= 0) return(NULL)
  x2 <- sum(pr^2, na.rm = TRUE)
  list(kind = "tested", family = kind, ratio = x2 / dfr, df = dfr, glmm = !inherits(model, "glm"),
       p_over = stats::pchisq(x2, dfr, lower.tail = FALSE), p_under = stats::pchisq(x2, dfr))
}

## Observed zeros against zeros simulated from a Poisson or negative binomial
## fit (fixed seed; the user's stream is left as it was). p: the share of
## simulations with at least as many zeros (+1 in both counts).
pm_check_zeros <- function(model, nsim = 200) {
  kind <- pm_qresid_family(model)
  if (is.null(kind) || !kind %in% c("poisson", "negbin")) return(NULL)
  y <- tryCatch(stats::model.response(stats::model.frame(model)), error = function(e) NULL)
  if (!is.numeric(y) || is.matrix(y)) return(NULL)
  sims <- tryCatch(suppressWarnings(suppressMessages(pm_with_seed(20261121, stats::simulate(model, nsim = nsim)))),
                   error = function(e) NULL)
  if (is.null(sims)) return(NULL)
  z <- vapply(sims, function(s) sum(as.matrix(s) == 0), numeric(1))
  obs <- sum(y == 0)
  if (obs == 0 && mean(z) < 0.5) return(NULL)                   # no zeros, none expected: nothing to say
  list(obs = obs, expected = mean(z), p = (sum(z >= obs) + 1) / (length(z) + 1))
}

## Signs of separation in a binomial model: fitted probabilities at 0 or 1,
## coefficients beyond 10 (logit scale) or with SE above 50.
pm_check_separation <- function(model) {
  if (!identical(pm_qresid_family(model), "binomial")) return(NULL)
  mu <- tryCatch(as.vector(stats::fitted(model)), error = function(e) NULL)
  if (is.null(mu)) return(NULL)
  b <- tryCatch(pm_coef(model), error = function(e) NULL)
  b <- b[!is.na(b)]
  se <- tryCatch(sqrt(diag(as.matrix(pm_vcov(model)))), error = function(e) NULL)
  big <- if (length(se) == length(b)) names(b)[abs(b) > 10 | se > 50] else character(0)
  list(n01 = sum(mu < 1e-8 | mu > 1 - 1e-8), big = big)
}

## Generalized variance inflation factors, GVIF^(1/(2 df)), of the terms of
## a model of main effects only (Fox and Monette 1992), from the correlation
## of the estimates. NULL with fewer than two terms, an interaction or no
## intercept.
pm_check_gvif <- function(model) {
  if (inherits(model, c("aovlist", "gam"))) return(NULL)
  tl <- attr(stats::terms(if (inherits(model, "merMod")) stats::formula(model, fixed.only = TRUE) else model),
             "term.labels")
  if (length(tl) < 2 || any(grepl(":", tl))) return(NULL)
  b <- pm_coef(model)
  b <- b[!is.na(b)]
  V <- as.matrix(pm_vcov(model))
  if (nrow(V) != length(b)) return(NULL)
  X <- if (inherits(model, "merMod")) lme4::getME(model, "X") else stats::model.matrix(model)
  asg <- attr(X, "assign")[match(names(b), colnames(X))]
  if (anyNA(asg) || !any(asg == 0)) return(NULL)
  keep <- asg != 0
  R <- stats::cov2cor(V[keep, keep, drop = FALSE])
  asg <- asg[keep]
  dR <- det(R)
  terms_in <- unique(asg)
  g <- vapply(terms_in, function(j) {
    i <- asg == j
    gv <- det(R[i, i, drop = FALSE]) * det(R[!i, !i, drop = FALSE]) / dR
    gv^(1 / (2 * sum(i)))
  }, numeric(1))
  stats::setNames(g, tl[terms_in])
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
