## Advice in outputs (rule no. 1, "Explained"): one engine for every function.
## Advice is a data frame (id, level, text) stored as attr(x, "advice");
## printing is one shared block, switched off with
## options(predictmeans.advice = FALSE). Rules are triggered by facts of the
## fit and data only; each id is stable so code and tests can refer to it.

pm_advice_levels <- c("caution", "note")

## One row of advice. 'level': "caution" (the result may mislead as it
## stands) or "note" (a choice the user may want to know about).
pm_advice <- function(id, text, level = c("note", "caution")) {
  level <- match.arg(level)
  stopifnot(is.character(id), length(id) == 1, nzchar(id), is.character(text), length(text) == 1, nzchar(text))
  data.frame(id = id, level = level, text = text, stringsAsFactors = FALSE)
}

## Add rows of advice to an object (rows from pm_advice(); NULL or zero rows
## are ignored). Rows with an id already present replace the old one.
pm_add_advice <- function(x, ...) {
  rows <- Filter(function(r) !is.null(r) && NROW(r) > 0, list(...))
  if (!length(rows)) return(x)
  new <- do.call(rbind, rows)
  old <- attr(x, "advice")
  if (!is.null(old)) new <- rbind(old[!old$id %in% new$id, , drop = FALSE], new)
  new <- new[order(match(new$level, pm_advice_levels)), , drop = FALSE]
  rownames(new) <- NULL
  attr(x, "advice") <- new
  x
}

## The advice of an object: a data frame with columns id, level, text
## (zero rows when there is none).
pm_advice_get <- function(x) {
  a <- attr(x, "advice")
  if (is.null(a)) a <- data.frame(id = character(0), level = character(0), text = character(0),
                                  stringsAsFactors = FALSE)
  a
}

#' Advice on a Result
#'
#' @description
#' Results of predictmeans functions carry short advice when this fit and these data call for a choice or a
#' caution the user may not know about: an unbalanced design, empty cells, a covariate aliased with a factor, a
#' singular fit, an extrapolated mean, few degrees of freedom, many tests, too few permutations, and so on. Advice
#' is printed at the end of the result as an "Advice:" block (\code{!} a caution: the result may mislead as it
#' stands; \code{i} a note: a choice worth knowing). It is never generic, and no test is run that was not asked
#' for. \code{advice()} returns it as a data frame, for code and AI assistants.
#'
#' \code{options(predictmeans.advice = FALSE)} stops the printing; the advice is still stored.
#' @param x A result of \code{predictmeans()}, \code{contrastmeans()}, \code{anova()} of an \code{lmer} or
#'   \code{semireg} fit, \code{pm_check()}, \code{permmodels()}, \code{bootmodels()}, \code{R2_model()},
#'   \code{varcomp()}, \code{se_ranef()}, \code{semireg()}, \code{semipred()}, \code{covariatemeans()} or \code{pm_slopes()}.
#' @return A data frame with columns \code{id} (stable, listed below), \code{level} (\code{"caution"} or
#'   \code{"note"}) and \code{text}; zero rows when there is no advice.
#' @section Advice ids:
#' \tabular{lll}{
#'   \strong{id} \tab \strong{from} \tab \strong{when} \cr
#'   \code{nonestimable_means} \tab predictmeans, covariatemeans \tab some means cannot be estimated (NA) \cr
#'   \code{empty_cells} \tab predictmeans \tab means average over the observed cells only \cr
#'   \code{aliased_covariate} \tab predictmeans \tab a covariate is constant within a factor's levels \cr
#'   \code{large_se} \tab predictmeans \tab a mean's SE is over 3 times the median SE \cr
#'   \code{lsd_varies} \tab predictmeans \tab the average LSD judges some pair differently from its own test \cr
#'   \code{bias_adjust} \tab predictmeans, covariatemeans \tab back-transformed means are a typical cluster's or a median, not the mean \cr
#'   \code{contrast_common_cells} \tab contrastmeans \tab contrasts are estimated on the cells their means share \cr
#'   \code{small_df} \tab predictmeans, contrastmeans \tab degrees of freedom below 5 \cr
#'   \code{many_tests_unadjusted} \tab predictmeans, contrastmeans \tab 6 or more requested tests with \code{adj = "none"} \cr
#'   \code{singular_fit} \tab predictmeans, contrastmeans, R2_model, semireg \tab a variance at zero or a correlation at +-1 \cr
#'   \code{anova_unbalanced_typeI}, \code{_typeII}, \code{_typeIII} \tab anova \tab unbalanced data: what the type tests \cr
#'   \code{anova_nested_in_data} \tab anova \tab a factor nested in another in the data, without the nesting in the model \cr
#'   \code{anova_empty_cells} \tab anova \tab empty cells: no type tests a main effect over all cells \cr
#'   \code{anova_covariate_zero} \tab anova \tab Type III with a covariate in an interaction (tests at covariate = 0) \cr
#'   \code{check_<name>} \tab pm_check \tab a check worth a look (the note of its row) \cr
#'   \code{few_permutations} \tab permmodels, bootmodels \tab fewer than 999 resamples and a p-value at the floor or near 0.05 \cr
#'   \code{r2_latent_scale} \tab R2_model \tab a non-gaussian model: R2 on the link scale \cr
#'   \code{variance_at_zero} \tab varcomp, se_ranef \tab a variance estimated as zero \cr
#'   \code{few_re_levels} \tab varcomp, se_ranef \tab a random effect with fewer than 5 levels \cr
#'   \code{variance_test_boundary} \tab varcomp \tab how to test a variance against zero (\code{test = TRUE}) \cr
#'   \code{ranef_conditional_fixed} \tab se_ranef \tab conditional SEs of a random term that contains a fixed factor \cr
#'   \code{smooth_linear} \tab semireg \tab a smooth estimated as a straight line (variance 0) \cr
#'   \code{extrapolation} \tab semipred, covariatemeans, pm_slopes \tab predictions beyond the range of the data \cr
#'   \code{pointwise_regions} \tab semipred, covariatemeans \tab regions of difference read from pointwise limits \cr
#'   \code{smooth_other_covariate} \tab semipred \tab a smooth in the prediction does not involve \code{covariate} (adds a constant) \cr
#'   \code{nonestimable_slopes} \tab pm_slopes \tab the covariate is aliased with a factor: slopes NA, their differences kept \cr
#'   \code{slopes_equal_by_model} \tab pm_slopes \tab no covariate x term interaction: the slopes are equal by assumption \cr
#'   \code{slopes_at_value} \tab pm_slopes \tab a curved term or a link: the slopes depend on \code{at} \cr
#' }
#' @examples
#' oats <- nlme::Oats
#' oats$nitro <- factor(oats$nitro)
#' fm <- lme4::lmer(log(yield) ~ nitro * Variety + (1 | Block), data = oats)
#' res <- predictmeans(fm, "nitro:Variety", pairwise = TRUE)
#' advice(res)
#' @export
advice <- function(x) pm_advice_get(x)

## Print the advice of an object as one block (nothing when there is none
## or when options(predictmeans.advice = FALSE)).
pm_print_advice <- function(x, width = 0.95 * getOption("width")) {
  txt <- pm_format_advice(pm_advice_get(x), width)
  if (nzchar(txt)) cat(txt)
  invisible(x)
}

## The printed block of advice rows ("" when none or switched off): the one
## layout used by every print method and by anova headings.
pm_format_advice <- function(a, width = 0.95 * getOption("width")) {
  if (is.null(a) || !nrow(a) || isFALSE(getOption("predictmeans.advice", TRUE))) return("")
  mark <- ifelse(a$level == "caution", "! ", "i ")
  lines <- vapply(seq_len(nrow(a)), function(i) {
    paste(strwrap(a$text[i], width = width, initial = mark[i], exdent = 2), collapse = "\n")
  }, "")
  paste0("\nAdvice:\n", paste(lines, collapse = "\n"), "\n")
}

## ---- Rules of predictmeans() and contrastmeans() ----
## Thresholds in one place.
pm_advice_limits <- list(se_ratio = 3, min_df = 5, many_tests = 6)

## Rules shared by means and contrasts: singular fit, small df, many tests
## without adjustment. 'df': the df of the estimates; 'ntests': number of
## tests reported with p-values (0 when none).
pm_advice_common <- function(model, df, ntests, adj, what = "pairs") {
  lim <- pm_advice_limits
  rows <- list()
  if (isTRUE(pm_singular(model))) {
    rows$sing <- pm_advice("singular_fit", paste(
      "The fit is singular (a variance estimated as zero or a correlation as +-1): SEs and df treat it as known.",
      "Consider simplifying the random effects."), "caution")
  }
  d <- df[is.finite(df)]
  if (length(d) && min(d) < lim$min_df) {
    rows$df <- pm_advice("small_df", paste0(
      "Few degrees of freedom (", format(signif(min(d), 3)), "): intervals are wide and p-values depend on the ",
      "df method; report it (df_method) and interpret with care."), "note")
  }
  if (ntests >= lim$many_tests && identical(adj, "none")) {
    rows$many <- pm_advice("many_tests_unadjusted", paste0(
      ntests, " ", what, " tested without adjustment: some will be significant by chance. If they are a family, ",
      "use adj = \"holm\" (any design) or \"tukey\"/\"mvt\" (all pairs)."), "note")
  }
  rows
}

## Advice for a predictmeans() result. 'v': the flat list of the stages
## (pm_flat()); 'parts': pm_tidy_parts(). Returns rows of pm_advice() or NULL.
pm_means_advice <- function(v, parts) {
  means <- parts$means
  s <- parts$settings
  model <- v$model
  km <- v$ctr.matrix
  rows <- list()
  if (any(!means$estimable)) {
    rows$na <- pm_advice("nonestimable_means", paste0(
      sum(!means$estimable), " mean(s) cannot be estimated (NA): the data do not identify them under this model. ",
      "Their estimable differences are kept; see vignette(\"missing-cells\")."), "caution")
  }
  if (!is.null(km$empty_cells) && identical(v$weights, "equal")) {
    aff <- names(which(km$empty_cells$affected))
    rows$empty <- pm_advice("empty_cells", paste0(
      "Some means average over the observed cells only (empty cells): ", paste(utils::head(aff, 5), collapse = ", "),
      if (length(aff) > 5) ", ...", ". They are not comparable with the others as averages; their pairs are ",
      "compared on common cells (column 'over')."), "caution")
  }
  lc <- km$level_covariates
  if (length(lc)) {
    rows$lc <- pm_advice("aliased_covariate", paste0(
      paste(names(lc), collapse = ", "), " is constant within ", paste(unique(vapply(lc, `[[`, "", "factor")),
                                                                       collapse = ", "),
      ": each mean uses each level's own value. For levels without data at every level averaged over, ",
      "weights = \"equal\" extrapolates their regression line, weights = \"cells\" averages the observed cells; ",
      "choose by the question."), "note")
  }
  se <- means$std.error[means$estimable]
  if (length(se) >= 3 && all(is.finite(se))) {
    big <- se > pm_advice_limits$se_ratio * stats::median(se)
    if (any(big)) {
      lab <- pm_advice_labels(means[means$estimable, , drop = FALSE], s$modelterm)[big]
      rows$se <- pm_advice("large_se", paste0(
        "Mean(s) with SE over ", pm_advice_limits$se_ratio, " times the median: ", paste(utils::head(lab, 5),
                                                                                       collapse = ", "),
        if (length(lab) > 5) ", ...", ". Few data or extrapolation (empty cells, covariate far from its data); ",
        "check the counts (n) before comparing them."), "caution")
    }
  }
  bt <- s$back_transform
  if (!is.null(bt) && !isTRUE(is.finite(s$bias_adjust)) && !identical(s$average, "response")) {
    glmm <- grepl("typical cluster", bt$estimate)
    med <- grepl("^median", bt$estimate)
    if (glmm || med) {
      rows$bias <- pm_advice("bias_adjust", paste0(
        "'response' is the ", if (glmm) "response of a typical cluster (random effects at zero)" else
          "median of the original scale (the geometric mean for a log)",
        ", not the mean. If the mean on the original scale is wanted (totals, averages to report), use ",
        "bias_adjust = TRUE", if (pm_advice_slopes(model)) " with 'sigma' (the model has random slopes)", "."), "note")
    }
  }
  # LSDs that differ between the pairs compared: advice only when the one
  # average LSD (the printed one, the plot's bar) would judge some pair
  # differently from its own test
  lsd <- v$LSD
  pr <- parts$pairs
  if (isTRUE(v$compare) && isFALSE(attr(lsd, "equal")) && is.data.frame(pr) && !is.null(pr$lsd)) {
    own <- if (!is.null(pr$adj.conf.high)) (pr$adj.conf.high - pr$adj.conf.low)/2 else pr$lsd
    keep <- is.finite(pr$estimate) & is.finite(own)
    if (!is.null(pr$by)) keep <- keep & !is.na(pr$by)
    avg <- lsd[["Aveg.LSD"]]
    wrong <- sum(keep & ((abs(pr$estimate) > avg) != (abs(pr$estimate) > own)))
    if (wrong > 0) {
      rows$lsd <- pm_advice("lsd_varies", paste0(
        wrong, " of ", sum(keep), " pairs are judged differently by the average LSD than by their own test ",
        "(LSDs ", format(signif(lsd[["Min.LSD"]], 3)), " to ", format(signif(lsd[["Max.LSD"]], 3)),
        "): use the p-values or pairs$lsd."), "note")
    }
  }
  pairs <- parts$pairs
  # only comparisons the user asked for (the pairs table is always there)
  ntests <- if (isTRUE(v$compare) && is.data.frame(pairs) && !is.null(pairs$p.value)) sum(!is.na(pairs$p.value)) else 0
  rows <- c(rows, pm_advice_common(model, means$df, ntests, s$adjust, "pairs"))
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## Advice for a contrastmeans() result.
pm_contrast_advice <- function(model, tab, adj) {
  rows <- list()
  if (!is.null(tab$over) && any(!is.na(tab$over))) {
    rows$over <- pm_advice("contrast_common_cells", paste0(
      "Contrasts with empty cells are estimated on the cells their means share (column 'over'): ",
      paste(utils::head(tab$contrast[!is.na(tab$over)], 5), collapse = ", "), "."), "note")
  }
  rows <- c(rows, pm_advice_common(model, tab$df, if (nrow(tab) > 1) sum(!is.na(tab$p.value)) else 0, adj, "contrasts"))
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## Labels of the rows of a means table (the factor columns of 'modelterm').
pm_advice_labels <- function(means, modelterm) {
  vars <- intersect(strsplit(modelterm, ":")[[1]], names(means))
  do.call(paste, c(lapply(means[vars], as.character), sep = ":"))
}

## Does a mixed model have random slopes (bias_adjust then needs 'sigma')?
pm_advice_slopes <- function(model) {
  isTRUE(tryCatch({
    if (inherits(model, "merMod")) {
      vc <- as.data.frame(lme4::VarCorr(model))
      any(vc$grp != "Residual" & (vc$var1 != "(Intercept)" | !is.na(vc$var2)))
    } else if (inherits(model, "glmmTMB")) {
      any(vapply(glmmTMB::VarCorr(model)$cond, function(m) nrow(as.matrix(m)) != 1, logical(1)))
    } else if (inherits(model, "lme")) {
      any(vapply(nlme::pdMatrix(model$modelStruct$reStruct), nrow, integer(1)) != 1)
    } else FALSE
  }, error = function(e) FALSE))
}

## ---- Rules of permmodels()/bootmodels(), R2_model() and varcomp() ----

## permmodels()/bootmodels(): too few permutations to resolve the p-values
## that matter. A p-value of k/(nperm + 1) has Monte Carlo SE about
## sqrt(p(1 - p)/nperm); flag when a p-value at or near the floor, or one
## within 2 Monte Carlo SEs of 0.05, is reported.
pm_perm_advice <- function(out, nperm, exact = FALSE) {
  if (is.null(nperm) || !isTRUE(nperm > 0)) return(NULL)
  col <- function(tab) {
    if (is.null(tab)) return(numeric(0))
    j <- grep("Perm_p_value|Boot_p_value", colnames(tab))
    if (!length(j)) numeric(0) else stats::na.omit(as.numeric(as.matrix(tab)[, j[1]]))
  }
  p <- c(col(out$COEFFICENTS), col(out$ANOVA))
  if (!length(p)) return(NULL)
  floor_p <- if (isTRUE(exact)) 1/nperm else 1/(nperm + 1)
  mcse <- sqrt(0.05 * 0.95 / nperm)
  rows <- list()
  at_floor <- any(p <= 2 * floor_p)
  near_cut <- any(abs(p - 0.05) < 2 * mcse)
  if (nperm < 999 && (at_floor || near_cut)) {
    rows$n <- pm_advice("few_permutations", paste0(
      "Only ", nperm, " resamples: ", if (at_floor) paste0("p-values at the floor (", format(signif(floor_p, 2)),
                                                           ") only say 'as small as can be shown'")
      else "p-values near 0.05 are uncertain", " (Monte Carlo SE about ", format(signif(mcse, 2)),
      " at p = 0.05). Use nperm = 999 or more (4999 for reporting); 'seed' makes it reproducible."), "caution")
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## R2_model(): what the numbers are and are not.
pm_r2_advice <- function(model, method) {
  model <- pm_fit(model)
  rows <- list()
  fam <- tryCatch(stats::family(model)$family, error = function(e) "gaussian")
  if (!identical(fam, "gaussian") && !inherits(model, c("lmerMod", "lme", "gls"))) {
    rows$latent <- pm_advice("r2_latent_scale", paste(
      "For a non-gaussian model the R2 is on the link (latent) scale: not comparable with",
      "an lm R2 or across families and links; use it to compare models of the same family."), "note")
  }
  if (isTRUE(pm_singular(model))) {
    rows$sing <- pm_advice("singular_fit", paste(
      "The fit is singular: a variance estimated as zero gives no share to its random effect, so the random part",
      "of R2 may be understated. Consider simplifying the random effects."), "caution")
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## varcomp(): estimates at zero (boundary), few levels, tests on the boundary.
pm_varcomp_advice <- function(model, tab, test = FALSE) {
  model <- pm_fit(model)
  rows <- list()
  zero <- isTRUE(pm_singular(model))
  if (!zero && is.data.frame(tab) && "vcov" %in% names(tab)) {
    re <- !grepl("^residual|^cor|^var\\.|Residual", rownames(tab))
    zero <- any(abs(tab$vcov[re]) < 1e-8 * max(1, abs(tab$vcov), na.rm = TRUE), na.rm = TRUE)
  }
  if (zero) {
    rows$sing <- pm_advice("variance_at_zero", paste(
      "A variance is estimated as zero (or a correlation as +-1): its interval and SE are not reliable, and",
      "'no variation' is not shown. Test it (test = TRUE) or simplify the random effects."), "caution")
  }
  ng <- pm_check_re_levels(model)
  few <- ng[!is.na(ng) & ng < 5]
  if (length(few)) {
    rows$few <- pm_advice("few_re_levels", paste0(
      "Few levels for ", paste0(names(few), " (", few, ")", collapse = ", "), ": its variance is poorly estimated ",
      "(often as zero) and its interval wide; consider it as a fixed effect if its levels are of interest."), "note")
  }
  if (!test && !is.null(pm_vc_engine(model)) && length(rows)) {
    rows$test <- pm_advice("variance_test_boundary", paste(
      "To test a variance against zero use varcomp(test = TRUE): it uses the boundary (mixture or bootstrap)",
      "reference, not the naive chi-square, which is conservative."), "note")
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## ---- Rules of semireg() and semipred() ----

## The variances of the smooth blocks of a semireg fit (named by block;
## NULL when they cannot be read).
pm_sr_block_vars <- function(fit) {
  blocks <- names(fit$refit$smoothZ)
  if (!length(blocks)) return(NULL)
  tryCatch({
    semer <- fit$semer
    if (inherits(semer, "glmmTMB")) {
      v <- glmmTMB::VarCorr(semer)$cond
      stats::setNames(vapply(blocks, function(b) as.matrix(v[[b]])[1, 1], numeric(1)), blocks)
    } else {
      v <- as.data.frame(lme4::VarCorr(semer))
      stats::setNames(v$vcov[match(blocks, v$grp)], blocks)
    }
  }, error = function(e) NULL)
}

## semireg(): smooths estimated as straight lines, other variances at zero.
pm_semireg_advice <- function(fit) {
  rows <- list()
  bv <- pm_sr_block_vars(fit)
  flat <- names(bv)[!is.na(bv) & bv < 1e-8 * max(1, abs(bv), na.rm = TRUE)]
  if (length(flat)) {
    rows$flat <- pm_advice("smooth_linear", paste0(
      "Smooth term(s) ", paste(flat, collapse = ", "), " estimated as a straight line (variance 0): the data show ",
      "no curvature there. The fit is still valid; anova(fit) tests each smooth."), "note")
  }
  semer <- fit$semer
  if (isTRUE(pm_singular(semer)) && !length(flat)) {
    rows$sing <- pm_advice("singular_fit", paste(
      "A random-effect variance is estimated as zero (or a correlation as +-1): SEs and df treat it as known.",
      "Consider simplifying the random effects."), "caution")
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## semipred(): predictions outside the data, pointwise bands read as regions.
## 'x': the result; 'fit': the semireg fit.
pm_semipred_advice <- function(x, fit) {
  s <- x$settings
  rows <- list()
  tab <- if (!is.null(x$means)) x$means else if (!is.null(x$slopes)) x$slopes else x$pairs
  cv <- intersect(s$covariate, names(tab))
  out_cv <- character(0)
  for (v in cv) {
    d <- fit$data[[v]]
    if (!is.numeric(d) || !is.numeric(tab[[v]])) next
    r <- range(d, na.rm = TRUE)
    tol <- 1e-8 * max(1, diff(r))
    if (any(tab[[v]] < r[1] - tol | tab[[v]] > r[2] + tol, na.rm = TRUE)) out_cv <- c(out_cv, v)
  }
  if (length(out_cv)) {
    rows$extra <- pm_advice("extrapolation", paste0(
      "Predictions beyond the range of the data for ", paste(out_cv, collapse = ", "), ": a spline there follows ",
      "its straight-line part only and the band understates the uncertainty; restrict 'at' to the data."), "caution")
  }
  sig <- x$significant
  if (!identical(s$band, "simultaneous") && is.data.frame(sig) && nrow(sig) > 0) {
    rows$band <- pm_advice("pointwise_regions", paste(
      "The regions where the", if (!is.null(x$pairs)) "curves differ" else "slope or curvature differs from 0",
      "come from pointwise limits: over a whole curve some will be by chance. band = \"simultaneous\"",
      "controls the error over the curve."), "note")
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## ---- Rules of covariatemeans() and pm_slopes() ----

## Covariates of 'tab' whose values lie outside the range of 'data'.
pm_advice_outside <- function(tab, data, covs) {
  out <- character(0)
  for (v in intersect(covs, names(tab))) {
    d <- data[[v]]
    if (!is.numeric(d) || !is.numeric(tab[[v]])) next
    r <- range(d, na.rm = TRUE)
    tol <- 1e-8 * max(1, diff(r))
    if (any(tab[[v]] < r[1] - tol | tab[[v]] > r[2] + tol, na.rm = TRUE)) out <- c(out, v)
  }
  out
}

## covariatemeans(): lines beyond the data, pointwise regions of difference,
## back-transformed lines that are not means, singular fit.
pm_covariatemeans_advice <- function(x, model) {
  s <- x$settings
  rows <- list()
  tab <- if (!is.null(x$means)) x$means else x$pairs
  dat <- tryCatch(pm_get_data(model), error = function(e) NULL)
  out_cv <- if (is.data.frame(dat) && is.data.frame(tab)) pm_advice_outside(tab, dat, s$covariate) else character(0)
  if (length(out_cv)) {
    rows$extra <- pm_advice("extrapolation", paste0(
      "Lines beyond the range of the data for ", paste(out_cv, collapse = ", "), ": the model's straight line (or ",
      "polynomial) is assumed there and nothing checks it; restrict the covariate values to the data."), "caution")
  }
  if (is.data.frame(x$means) && anyNA(x$means$estimate)) {
    rows$na <- pm_advice("nonestimable_means", paste0(
      sum(is.na(x$means$estimate)), " point(s) of the curves cannot be estimated (NA): the data do not identify ",
      "them under this model."), "caution")
  }
  sig <- x$significant
  if (!is.null(x$pairs) && !identical(s$band, "simultaneous") && is.data.frame(sig) && nrow(sig) > 0) {
    rows$band <- pm_advice("pointwise_regions", paste(
      "The regions where the lines differ come from pointwise limits: over a whole line some will be by chance.",
      "band = \"simultaneous\" controls the error over the line."), "note")
  }
  bt <- s$back_transform
  if (!is.null(bt) && (grepl("typical cluster", bt$estimate) || grepl("^median", bt$estimate))) {
    rows$bias <- pm_advice("bias_adjust", paste0(
      "The back-transformed lines are the ", if (grepl("typical cluster", bt$estimate))
        "response of a typical cluster (random effects at zero)" else "median of the original scale",
      ", not the mean; bias_adjust = TRUE gives the mean", if (pm_advice_slopes(model)) " (with 'sigma': random slopes)",
      "."), "note")
  }
  rows <- c(rows, pm_advice_common(model, numeric(0), 0, "none"))
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## pm_slopes(): slopes that are the same by the model's form, the point
## where slopes are compared, the scale, and the common rules.
pm_slopes_advice <- function(x, model) {
  s <- x$settings
  rows <- list()
  vars <- strsplit(s$modelterm, ":")[[1]]
  # terms containing the covariate, also inside poly(), ns(), I(x^2), ...
  tl <- tryCatch(attr(stats::terms(model), "term.labels"), error = function(e) character(0))
  has_cov <- vapply(tl, function(t) s$covariate %in% all.vars(str2lang(t)), logical(1))
  terms_cov <- tl[has_cov]
  inter <- vapply(terms_cov, function(t) any(vars %in% all.vars(str2lang(t))), logical(1))
  # curved in the covariate: a function of it other than itself (poly, ns, I(x^2), log) or x:x
  poly <- any(vapply(terms_cov, function(t) any(vapply(strsplit(t, ":")[[1]], function(p) {
    p <- trimws(p)
    s$covariate %in% all.vars(str2lang(p)) && p != s$covariate
  }, logical(1))), logical(1)))
  if (isTRUE(any(!x$slopes$estimable))) {
    rows$na <- pm_advice("nonestimable_slopes", paste0(
      sum(!x$slopes$estimable), " slope(s) cannot be estimated (NA): ", s$covariate, " is aliased with a factor of ",
      "the model (e.g. constant within its levels)",
      if (!is.null(x$pairs)) ", so only differences between slopes are estimable; they are in the pairs", "."),
      "caution")
  }
  if (nrow(x$slopes) > 1 && !any(inter)) {
    rows$same <- pm_advice("slopes_equal_by_model", paste0(
      "The model has no ", s$covariate, " x ", s$modelterm, " interaction, so the slopes are equal by assumption; ",
      "add the interaction to estimate and compare different slopes."), "caution")
  }
  nonlin <- poly || !is.null(pm_link_trans(model))
  if (nonlin) {
    rows$at <- pm_advice("slopes_at_value", paste0(
      "The slope changes along ", s$covariate, " (", if (poly) "a curved term" else "a non-identity link (the ",
      "slopes are on the link scale)", "): these are at ", s$covariate, " = ", format(signif(s$at, 4)),
      "; give 'at' for other values."), "note")
  }
  dat <- tryCatch(pm_get_data(model), error = function(e) NULL)
  if (is.data.frame(dat) && is.numeric(dat[[s$covariate]])) {
    r <- range(dat[[s$covariate]], na.rm = TRUE)
    if (any(s$at < r[1] | s$at > r[2])) {
      rows$extra <- pm_advice("extrapolation", paste0(
        "'at' is beyond the range of ", s$covariate, " in the data: the slope there is an extrapolation."), "caution")
    }
  }
  pairs <- x$pairs
  rows <- c(rows, pm_advice_common(model, x$slopes$df, if (is.data.frame(pairs)) nrow(pairs) else 0, s$adjust,
                                   "pairs of slopes"))
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}
