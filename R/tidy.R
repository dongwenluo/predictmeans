## Stable, flat result tables for predictmeans().
##
## predictmeans() keeps its legacy list slots unchanged and adds three slots
## with fixed names, built here from the values it has already computed:
##   means     one row per predicted mean
##   pairs     one row per pairwise comparison
##   settings  the options that produced them
## Values are not rounded; rounding is left to printing.

## 'e': the results of the stages of predictmeans() in one list (pm_flat()).
pm_tidy_parts <- function(e) {
  g <- function(nm) e[[nm]]
  vars <- g("vars")
  K <- g("K")
  KK <- g("KK")
  label <- g("label")
  mt <- g("mt")
  model <- g("model")
  adj <- g("adj")
  permlist <- g("permlist")
  perm <- is.list(permlist)
  atvar <- g("atvar")
  if (!is.null(atvar) && all(atvar %in% c("NULL", ""))) atvar <- NULL
  slevel <- g("slevel")
  trans <- g("trans")
  nK <- nrow(K)

  ## Means
  est <- as.vector(mt$pm)
  se <- as.vector(mt$ses)
  df <- if (perm) {
    rep(NA_real_, nK)
  } else if (!is.null(g("pm_df_raw"))) {
    as.vector(g("pm_df_raw"))
  } else if (!is.null(mt$Df)) {
    rep_len(as.numeric(mt$Df), nK)
  } else {
    rep(NA_real_, nK)
  }
  # with a permlist the tests are of comparisons: the limits of each mean are
  # the model's t limits at conf_level with the model's df ('ci_df' of
  # pm_compare(), pm_perm_ci_df()); 'df' stays NA, as no test uses it
  ci_df <- if (perm) g("ci_df")
  if (perm && is.null(ci_df)) ci_df <- rep(Inf, nK)
  crit <- pm_ci_crit(slevel, if (perm) ci_df else df)
  means <- label
  rownames(means) <- NULL
  for (v in names(means)) if (is.matrix(means[[v]])) means[[v]] <- as.vector(means[[v]])
  # means over empty cells (R/empty-cells.R): the cells each one averages, as
  # 'over' of the pairs (NA: all its cells)
  ec_over <- g("ctr.matrix")$empty_cells$over
  if (any(!is.na(ec_over))) means$over <- unname(ec_over[rownames(K)])
  means$estimate <- est
  means$std.error <- se
  means$df <- df
  means$conf.low <- est - crit*se
  means$conf.high <- est + crit*se
  # comparison intervals (R/comparison.R): overlap answers the pairwise tests
  cmp <- g("cmp_int")
  if (!is.null(cmp)) {
    means$comparison.low <- est - cmp$mult*se
    means$comparison.high <- est + cmp$mult*se
    # none for means over empty cells: their pairs are compared on common
    # cells, which no interval around the mean can show (R/comparison.R)
    off <- rownames(K) %in% cmp$left_out
    means$comparison.low[off] <- means$comparison.high[off] <- NA_real_
  }

  back <- function(x) {
    if (identical(trans, stats::make.link("log")$linkinv) || identical(trans, exp)) {
      exp(x)
    } else {
      as.numeric(trans(x))
    }
  }
  resp <- if (is.function(trans)) back(est) else est
  if (!isTRUE(all.equal(resp, est))) {
    lo <- back(means$conf.low)
    hi <- back(means$conf.high)
    means$response <- resp
    # delta method: SE on the response scale = |d back / d eta| x SE
    means$response.std.error <- abs(pm_back_deriv(trans, est, model, isTRUE(g("link_scale"))))*se
    means$response.conf.low <- pmin(lo, hi)
    means$response.conf.high <- pmax(lo, hi)
  }
  # bias adjustment and averaging on the response scale (R/response-scale.R)
  avg <- if (is.null(g("average"))) "link" else g("average")
  bs <- pm_bias_setup(g("bias_adjust"), g("sigma"), model, isTRUE(g("link_scale")),
                      need = !is.null(means$response))
  ba <- bs$on
  ba_sigma <- NA_real_
  if (ba || identical(avg, "response")) {
    if (is.null(means$response)) {
      warning("'bias_adjust' and 'average = \"response\"' apply to back-transformed means (a non-identity link ",
              "or 'trans'); ignored.", call. = FALSE)
      avg <- "link"
    } else {
      bf <- pm_back_funs(back, trans, model, isTRUE(g("link_scale")))
      if (ba) ba_sigma <- sqrt(bs$s2)
      # the mean of the back-transformed response (R/back-transform.R)
      bias <- pm_bias_funs(bf, bs$s2, bs$method)
      hb <- bias$h
      db <- bias$d
      if (identical(avg, "response")) {
        mp_ <- g("mp")
        ra <- pm_response_average(g("ctr.matrix")$refgrid, match(rownames(K), rownames(KK)), mp_$coef, mp_$vcov,
                                  hb, db)
        means$response <- ra$mu
        means$response.std.error <- ra$se
        means$response.conf.low <- ra$mu - crit*ra$se
        means$response.conf.high <- ra$mu + crit*ra$se
      } else {
        means$response <- hb(est)
        means$response.std.error <- abs(db(est))*se
        lo <- hb(means$conf.low)
        hi <- hb(means$conf.high)
        means$response.conf.low <- pmin(lo, hi)
        means$response.conf.high <- pmax(lo, hi)
      }
    }
  }

  key <- function(df) do.call(paste, c(lapply(df[vars], as.character), sep = "\r"))
  key_m <- key(means)
  mtab <- g("meanTable")
  if (is.list(mtab) && !is.data.frame(mtab)) mtab <- mtab$Table
  if (is.data.frame(mtab) && !is.null(mtab$LetterGrp)) {
    means$group <- as.character(mtab$LetterGrp)[match(key_m, key(mtab))]
  }
  ndf <- g("ndf")
  if (is.data.frame(ndf)) means$n <- ndf$Freq[match(key_m, key(ndf))]
  est_flag <- attr(KK, "estimable")
  means$estimable <- if (is.null(est_flag)) TRUE else unname(est_flag[rownames(K)])

  ## Covariance of the means (model scale), for coef() and vcov()
  vcovm <- g("vcovm")
  V <- NULL
  if (!is.null(vcovm) && ncol(K) == ncol(vcovm)) {
    V <- K %*% tcrossprod(as.matrix(vcovm), K)
    V[!means$estimable, ] <- NA
    V[, !means$estimable] <- NA
    dimnames(V) <- list(rownames(K), rownames(K))
  }

  ## Pairs: every statistic from the pair table of the compare stage
  ## (pm_pair_table(), R/predictmeans-stages.R), as the printed matrices
  pairs <- NULL
  pt <- g("pair_tab")
  if (nK > 1 && is.data.frame(pt)) {
    cm <- pt$estimate
    rnK <- rownames(K)
    pairs <- data.frame(level1 = rnK[pt$j], level2 = rnK[pt$i], stringsAsFactors = FALSE)
    if (!is.null(atvar)) pairs$by <- pt$by
    po <- g("pair_over")
    if (length(po) == nrow(pairs) && any(!is.na(po))) pairs$over <- po   # compared on common cells
    pairs$estimate <- cm
    pairs$std.error <- pt$std.error
    pairs$df <- pt$df
    pairs$statistic <- pt$statistic
    pairs$p.value <- pt$p.value
    pairs$p.adj <- pt$p.adj
    # equivalence tests (TOST); with a permlist they are t-based too
    eqm <- pm_equiv_margins(g("equivalence"))
    if (!is.null(eqm)) {
      tdf <- if (perm) {
        dd <- g("Df_diff")
        if (length(dd)) rep_len(as.numeric(dd), length(cm)) else Inf
      } else pt$df
      pairs$p.equiv <- pm_tost(cm, pt$std.error, tdf, eqm)
    }
    # limits of each pair (t, or the permutation quantile of |t|), and the
    # simultaneous limits that agree with p.adj (R/p-adjust.R)
    crit <- pt$crit
    ci_se <- pt$ci_se              # HC3 SEs with a robust permlist
    pairs$conf.low <- cm - crit*ci_se
    pairs$conf.high <- cm + crit*ci_se
    pairs$lsd <- crit*ci_se
    if (!is.null(pt$crit_adj)) {
      pairs$adj.conf.low <- cm - pt$crit_adj*ci_se
      pairs$adj.conf.high <- cm + pt$crit_adj*ci_se
    }
    pairs$p.method <- if (perm) {
      if (identical(attr(permlist, "resampling"), "parametric bootstrap")) "bootstrap" else "permutation"
    } else if (adj == "tukey") "tukey" else "t"

    # log scale: the difference is the log of the ratio of the means;
    # logit scale: the log of the odds ratio
    kind <- pm_ratio_kind(trans, model, isTRUE(g("link_scale")))
    if (!is.null(kind)) {
      crit_p <- crit
      pre <- if (kind == "log") "ratio" else "odds.ratio"
      pairs[[pre]] <- exp(cm)
      pairs[[paste0(pre, ".conf.low")]] <- exp(cm - crit_p*ci_se)
      pairs[[paste0(pre, ".conf.high")]] <- exp(cm + crit_p*ci_se)
    }
    # class probabilities (R/categorical.R): pairs of two classes compare different things, left out
    cv <- if (inherits(model, "pm_catfit")) model$var
    if (!is.null(cv) && cv %in% names(label)) {
      cl <- as.character(label[[cv]])
      pairs <- pairs[cl[pt$i] == cl[pt$j], , drop = FALSE]
    }
  }

  ## Settings
  # class probabilities (R/categorical.R) are averaged on the probability (response) scale
  if (inherits(model, "pm_catfit") && model$what != "link") avg <- "response"
  df_method <- if (perm) {
    "permutation"
  } else if (length(g("Df")) && isTRUE(g("df_user"))) {
    "user"
  } else {
    g("inf")$label              # R/inference.R
  }
  cov_cols <- setdiff(names(label), vars)
  settings <- list(
    modelterm = g("modelterm"),
    level = slevel,
    adjust = adj,
    atvar = atvar,
    df_method = df_method,
    df_cases = pm_df_cases(g("inf"), user = perm || identical(df_method, "user")),           # inner-outer: what the residual stratum counts (R/df-innerouter.R)
    vcov_adjust = isTRUE(g("vcov_adjust")),
    vcov = if (is.null(g("vcov_label"))) "model" else g("vcov_label"),
    equivalence = pm_equiv_margins(g("equivalence")),
    comparison = cmp,
    back_transform = pm_back_label(means, g("trans_by"), isTRUE(g("link_scale")), model, ba, avg),
    response = pm_response_label_info(model, trans, isTRUE(g("link_scale")), !is.null(means$response),
                                      from_formula = grepl("model formula", paste(g("trans_by"), collapse = ""))),
    bias_adjust = ba_sigma,
    singular = pm_singular(model),
    average = avg,
    weights_by = g("ctr.matrix")$weights_by,     # what "proportional" etc. count: rows or frequency weights
    nperm = if (perm) pm_perm_n(permlist) else NA_integer_,
    covariates = if (length(cov_cols)) unique(label[cov_cols]) else NULL,
    component = if (is.null(pm_zi_component(model))) NA_character_ else g("mean_type"),
    model_class = class(pm_cat_fit(model))[1],     # (the fit, not predictmeans' carrier: R/categorical.R)
    version = as.character(utils::packageVersion("predictmeans"))
  )
  zi <- pm_zi_component(model)
  if (!is.null(zi)) {
    message(if (identical(g("mean_type"), "overall")) paste0(
      "The model has zero inflation: the means are the overall means (1 - pi) mu (mean_type = \"overall\"), ",
      "averaged over the cells of the grid; 'estimate' is their log, the ratios are of these means.") else
      paste(attr(zi, "note"), "mean_type = \"overall\" gives the overall means."))
  }

  # with 'atvar': the means grouped by its levels (atvar first, then the
  # other factors, in modelterm order), the pairs by level in the same order,
  # pairs across levels last
  if (!is.null(atvar)) {
    o <- do.call(order, unname(as.list(means[c(atvar, setdiff(vars, atvar))])))
    means <- means[o, , drop = FALSE]
    rownames(means) <- NULL
    if (!is.null(V)) V <- V[o, o, drop = FALSE]
    if (!is.null(pairs)) {
      by_lev <- unique(do.call(paste, c(lapply(means[atvar], as.character), sep = " : ")))
      po <- order(match(pairs$by, by_lev), match(pairs$level1, rownames(K)[o]), match(pairs$level2, rownames(K)[o]),
                  na.last = TRUE)
      pairs <- pairs[po, , drop = FALSE]
      rownames(pairs) <- NULL
    }
  }

  list(means = means, pairs = pairs, settings = settings, vcov = V)
}

## Derivative of the back-transformation at the means (link scale): the
## family's mu.eta for a non-identity link, exp for a log, else the user's
## 'trans' by central differences (R/deriv.R). The one place for h' (also
## pm_back_funs()).
pm_back_deriv <- function(trans, eta, model, link_scale) {
  if (link_scale) {
    fam <- tryCatch(stats::family(model), error = function(e) NULL)
    link <- if (is.list(fam)) fam$link else NULL
    if (is.character(link) && link %in% c("logit", "probit", "cauchit", "cloglog", "identity", "log", "sqrt",
                                          "1/mu^2", "inverse")) {
      return(stats::make.link(link)$mu.eta(eta))
    }
    if (is.list(fam) && is.function(fam$mu.eta)) return(fam$mu.eta(eta))
  }
  if (identical(trans, exp) || identical(trans, stats::make.link("log")$linkinv)) return(exp(eta))
  d <- pm_fd_fun(function(x) as.numeric(trans(x)))
  vapply(eta, function(e) if (is.na(e)) NA_real_ else d(e), numeric(1))
}

## "log" when differences of the means are log ratios (log link or
## trans = exp), "logit" for log odds ratios, else NULL.
pm_ratio_kind <- function(trans, model, link_scale) {
  if (link_scale) {
    fam <- tryCatch(stats::family(model), error = function(e) NULL)
    link <- if (is.list(fam)) fam$link else ""
    if (identical(link, "log")) return("log")
    if (identical(link, "logit")) return("logit")
    return(NULL)
  }
  if (identical(trans, exp)) "log" else NULL
}

## Add the stable slots and set the class.
pm_finish <- function(outputlist, parts) {
  outputlist$means <- parts$means
  outputlist$pairs <- parts$pairs
  outputlist$settings <- parts$settings
  attr(outputlist, "vcov") <- parts$vcov
  class(outputlist) <- c("predictmeans", "pdmlist")
  outputlist
}

## The means (model scale) and their covariance matrix, named by the mean
## labels ("level1:level2"). With these, multcomp::glht() and
## marginaleffects::hypotheses() can test any contrast of the means.
#' @exportS3Method coef predictmeans
coef.predictmeans <- function(object, ...) {
  stats::setNames(object$means$estimate, rownames(attr(object, "vcov")))
}
#' @exportS3Method vcov predictmeans
vcov.predictmeans <- function(object, ...) {
  V <- attr(object, "vcov")
  if (is.null(V)) stop("This result has no covariance matrix of the means.", call. = FALSE)
  V
}

## Methods for marginaleffects (registered when it is loaded), so that
## marginaleffects::hypotheses() works on the means; its delta method changes
## the coefficients with set_coef().
#' @exportS3Method marginaleffects::get_coef predictmeans
get_coef.predictmeans <- function(model, ...) coef.predictmeans(model)
#' @exportS3Method marginaleffects::set_coef predictmeans
set_coef.predictmeans <- function(model, coefs, ...) {
  model$means$estimate <- unname(coefs)
  model
}
#' @exportS3Method marginaleffects::get_vcov predictmeans
get_vcov.predictmeans <- function(model, ...) vcov.predictmeans(model)

#' @exportS3Method as.data.frame predictmeans
as.data.frame.predictmeans <- function(x, row.names = NULL, optional = FALSE,
                                       what = c("means", "pairs"), ...) {
  pm_result_table(x, if (!missing(what)) what, c("means", "pairs"), row.names)
}

## The tables of a classed result as data frames: tidy() and as.data.frame()
## of predictmeans(), covariatemeans() / semipred() (pm_curves), pm_slopes()
## and contrastmeans(), in one place. 'tables': the result's table names;
## without 'what', the first of them that the result has.
pm_result_table <- function(x, what, tables, row.names = NULL) {
  what <- if (is.null(what)) tables[!vapply(tables, function(t) is.null(x[[t]]), logical(1))][1] else
    match.arg(what, tables)
  out <- if (is.na(what)) NULL else x[[what]]
  if (is.null(out)) stop("There is no '", if (is.na(what)) tables[1] else what, "' table in this result.", call. = FALSE)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' Predicted Means and Pairwise Comparisons as Data Frames
#'
#' @description
#' Return the \code{means} or \code{pairs} table of a \code{\link{predictmeans}} result (\code{glance()}: one row describing it). The tables have fixed column
#' names (see the Value section of \code{\link{predictmeans}}) and unrounded values, which makes them easy to use in
#' further analysis, reports, Shiny apps or AI tools. \code{tidy()} is the generic from package \pkg{generics} (also
#' used by \pkg{broom}); \code{as.data.frame()} needs no other package.
#' \code{coef()} and \code{vcov()} return the means (on the model scale, as in the \code{estimate} column) and
#' their covariance matrix (Kenward-Roger adjusted when the means are), named by the mean labels. With them,
#' \code{multcomp::glht()} or \code{marginaleffects::hypotheses()} can test any contrast of the means, e.g.
#' \code{hypotheses(pm, hypothesis = "n2 - n0 = 0")} (normal-based). Recent versions of \pkg{marginaleffects} need
#' the class listed first: \code{options(marginaleffects_model_classes = "predictmeans")}.
#' @usage
#' \method{tidy}{predictmeans}(x, what = c("means", "pairs"), ...)
#'
#' \method{glance}{predictmeans}(x, ...)
#'
#' \method{as.data.frame}{predictmeans}(x, row.names = NULL, optional = FALSE,
#'               what = c("means", "pairs"), ...)
#'
#' \method{coef}{predictmeans}(object, ...)
#'
#' \method{vcov}{predictmeans}(object, ...)
#' @param x A result of \code{\link{predictmeans}}.
#' @param object A result of \code{\link{predictmeans}}.
#' @param what Which table: \code{"means"} (the default) or \code{"pairs"}.
#' @param row.names Optional row names for the result.
#' @param optional Not used.
#' @param \dots Not used.
#' @return
#' A data frame; for \code{coef()} a named vector and for \code{vcov()} a matrix. \code{glance()} gives one row:
#' \code{modelterm}, \code{n.means}, \code{n.pairs}, \code{nobs} (observations in the cells), \code{df.method},
#' \code{adjust}, \code{level}, \code{vcov}, \code{vcov.adjust}, \code{sed.min}, \code{sed.mean}, \code{sed.max}
#' (standard errors of the differences), \code{lsd.mean}, \code{nperm} and \code{model.class}.
#' @seealso \code{\link{predictmeans}}, \code{\link{autoplot.predictmeans}}, \code{\link{advice}}
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
#'   as.data.frame(pm)
#'   head(as.data.frame(pm, what = "pairs"))
#'   if (requireNamespace("generics", quietly = TRUE)) generics::glance(pm)
#'
#'   # Any contrast of the means from coef() and vcov(): Victory - Golden Rain at nitro 0.6
#'   b <- coef(pm)
#'   L <- (names(b) == "0.6:Victory") - (names(b) == "0.6:Golden Rain")
#'   c(estimate = sum(L*b), SE = sqrt(drop(L %*% vcov(pm) %*% L)))
#' @name tidy.predictmeans
#' @aliases as.data.frame.predictmeans coef.predictmeans vcov.predictmeans glance.predictmeans
NULL

#' @exportS3Method generics::tidy predictmeans
tidy.predictmeans <- function(x, what = c("means", "pairs"), ...) {
  pm_result_table(x, if (!missing(what)) what, c("means", "pairs"))
}

## One row describing a result (the generics/broom glance()).
#' @exportS3Method generics::glance predictmeans
glance.predictmeans <- function(x, ...) {
  s <- x$settings
  m <- x$means
  p <- x$pairs
  rng <- function(v) if (length(v) && any(is.finite(v))) range(v[is.finite(v)]) else c(NA_real_, NA_real_)
  sed <- rng(if (is.null(p)) NULL else p$std.error)
  lsd <- if (is.null(p) || is.null(p$lsd)) NA_real_ else mean(p$lsd[is.finite(p$lsd)])
  num <- function(v) if (is.null(v) || !length(v)) NA_real_ else as.numeric(v)[1]
  chr <- function(v) if (is.null(v) || !length(v)) NA_character_ else paste(v, collapse = ", ")
  data.frame(modelterm = chr(s$modelterm), n.means = nrow(m), n.pairs = if (is.null(p)) 0L else nrow(p),
             nobs = if (is.null(m$n)) NA_real_ else sum(m$n, na.rm = TRUE), df.method = chr(s$df_method),
             adjust = chr(s$adjust), level = num(s$level), vcov = chr(s$vcov), vcov.adjust = isTRUE(s$vcov_adjust),
             sed.min = sed[1], sed.mean = if (is.null(p)) NA_real_ else mean(p$std.error[is.finite(p$std.error)]),
             sed.max = sed[2], lsd.mean = lsd, nperm = num(s$nperm), model.class = chr(s$model_class),
             stringsAsFactors = FALSE)
}
