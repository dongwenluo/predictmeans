## Stable, flat result tables for predictmeans().
##
## predictmeans() keeps its legacy list slots unchanged and adds three slots
## with fixed names, built here from the values it has already computed:
##   means     one row per predicted mean
##   pairs     one row per pairwise comparison
##   settings  the options that produced them
## Values are not rounded; rounding is left to printing.

pm_tidy_parts <- function(e) {
  g <- function(nm) get0(nm, envir = e, inherits = FALSE)
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
  transOff <- g("transOff")
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
  crit <- if (perm) rep(2, nK) else stats::qt(1 - slevel/2, df)
  means <- label
  rownames(means) <- NULL
  for (v in names(means)) if (is.matrix(means[[v]])) means[[v]] <- as.vector(means[[v]])
  means$estimate <- est
  means$std.error <- se
  means$df <- df
  means$conf.low <- est - crit*se
  means$conf.high <- est + crit*se

  back <- function(x) {
    if (identical(trans, stats::make.link("log")$linkinv) || identical(trans, exp)) {
      exp(x) - transOff
    } else {
      as.numeric(trans(x)) - transOff
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

  ## Pairs
  pairs <- NULL
  if (nK > 1 && !is.null(g("cm"))) {
    ij <- do.call(rbind, lapply(2:nK, function(i) cbind(j = seq_len(i - 1), i = i)))
    cm <- as.vector(g("cm"))
    dses <- as.vector(g("dses"))
    pdf <- if (perm) rep(NA_real_, nrow(ij)) else rep_len(as.numeric(g("Df_diff")), nrow(ij))
    stat <- cm/dses
    p <- if (perm) as.vector(g("t.p.values")) else 2*stats::pt(-abs(stat), pdf)

    rnK <- rownames(K)
    pairs <- data.frame(level1 = rnK[ij[, "j"]], level2 = rnK[ij[, "i"]], stringsAsFactors = FALSE)
    by <- rep(NA_character_, nrow(ij))
    if (!is.null(atvar)) {
      A <- label[ij[, "j"], atvar, drop = FALSE]
      B <- label[ij[, "i"], atvar, drop = FALSE]
      same <- Reduce(`&`, lapply(atvar, function(v) as.character(A[[v]]) == as.character(B[[v]])))
      by[same] <- do.call(paste, c(lapply(A[same, , drop = FALSE], as.character), sep = " : "))
      pairs$by <- by
    }
    pairs$estimate <- cm
    pairs$std.error <- dses
    pairs$df <- pdf
    pairs$statistic <- stat
    pairs$p.value <- p

    adjust_group <- function(idx) {
      n <- length(idx)
      if (adj == "tukey") {
        if (n < 2) return(p[idx])
        # number of means k from the k(k-1)/2 pairs among them
        stats::ptukey(sqrt(2)*abs(stat[idx]), (1 + sqrt(1 + 8*n))/2, pdf[idx], lower.tail = FALSE)
      } else if (adj %in% pm_adjust_extra_methods) {
        if (n < 2) return(p[idx])
        pm_adjust_extra(adj, stat[idx], p[idx], pdf[idx], rank = (1 + sqrt(1 + 8*n))/2 - 1,
                        R = if (adj == "mvt") g("pair_cor")[idx, idx, drop = FALSE])
      } else {
        stats::p.adjust(p[idx], adj)
      }
    }
    p_adj <- rep(NA_real_, length(p))
    if ((adj %in% pm_adjust_extra_methods || (perm && adj %in% pm_perm_maxt_methods)) &&
        length(g("p.extra")) == length(p)) {
      p_adj <- g("p.extra")              # computed with the legacy tables
    } else if (is.null(atvar)) {
      p_adj <- adjust_group(seq_along(p))
    } else {
      for (b in unique(stats::na.omit(by))) {
        idx <- which(by == b)
        p_adj[idx] <- adjust_group(idx)
      }
    }
    pairs$p.adj <- p_adj
    # equivalence tests (TOST); with a permlist they are t-based too
    eqm <- pm_equiv_margins(g("equivalence"))
    if (!is.null(eqm)) {
      tdf <- if (perm) {
        dd <- g("Df_diff")
        if (length(dd)) rep_len(as.numeric(dd), length(cm)) else Inf
      } else pdf
      pairs$p.equiv <- pm_tost(cm, dses, tdf, eqm)
    }
    lsd_level <- g("level")
    # critical values: t quantiles, or with a permlist those of the permutation
    # distribution of |t| (simultaneous within each family with a max-|t| adj)
    pT <- g("perm_T")
    crit <- if (perm && !is.null(pT)) {
      maxt_fam <- if (adj %in% pm_perm_maxt_methods) { if (is.null(atvar)) rep("all", length(cm)) else by } else NULL
      pm_perm_crit(pT$T, lsd_level, family = maxt_fam)
    } else if (perm) {
      rep(NA_real_, length(cm))
    } else stats::qt(1 - lsd_level/2, pdf)
    # a robust permlist (permmodels(robust = TRUE)) studentises with HC3 SEs
    ci_se <- if (!is.null(pT$se)) pT$se else dses
    pairs$conf.low <- cm - crit*ci_se
    pairs$conf.high <- cm + crit*ci_se
    pairs$lsd <- crit*ci_se
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
  }

  ## Settings
  df_method <- if (perm) {
    "permutation"
  } else if (length(g("Df")) && isTRUE(g("df_user"))) {
    "user"
  } else if (!is.null(g("df_type"))) {
    # "KR-SAS" is Kenward-Roger with the observed information
    if (identical(g("df_method"), "KR-SAS") && identical(g("df_type"), "Kenward-Roger")) "KR-SAS" else g("df_type")
  } else if (inherits(model, "lme")) {
    "containment (nlme)"
  } else {
    "residual"
  }
  cov_cols <- setdiff(names(label), vars)
  settings <- list(
    modelterm = g("modelterm"),
    level = slevel,
    adjust = adj,
    atvar = atvar,
    df_method = df_method,
    vcov_adjust = isTRUE(g("vcov_adjust")),
    vcov = if (is.null(g("vcov_label"))) "model" else g("vcov_label"),
    equivalence = pm_equiv_margins(g("equivalence")),
    nperm = if (perm) pm_perm_n(permlist) else NA_integer_,
    covariates = if (length(cov_cols)) unique(label[cov_cols]) else NULL,
    model_class = class(model)[1],
    version = as.character(utils::packageVersion("predictmeans"))
  )

  list(means = means, pairs = pairs, settings = settings, vcov = V)
}

## Derivative of the back-transformation at the means (link scale): the
## family's mu.eta for a non-identity link, else the user's 'trans'
## (numerically).
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
  if (identical(trans, exp)) return(exp(eta))
  vapply(eta, function(e) {
    if (is.na(e)) NA_real_ else numDeriv::grad(function(x) as.numeric(trans(x)), e)
  }, numeric(1))
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
coef.predictmeans <- function(object, ...) {
  stats::setNames(object$means$estimate, rownames(attr(object, "vcov")))
}
vcov.predictmeans <- function(object, ...) {
  V <- attr(object, "vcov")
  if (is.null(V)) stop("This result has no covariance matrix of the means.", call. = FALSE)
  V
}

## Methods for marginaleffects (registered when it is loaded), so that
## marginaleffects::hypotheses() works on the means; its delta method changes
## the coefficients with set_coef().
get_coef.predictmeans <- function(model, ...) coef.predictmeans(model)
set_coef.predictmeans <- function(model, coefs, ...) {
  model$means$estimate <- unname(coefs)
  model
}
get_vcov.predictmeans <- function(model, ...) vcov.predictmeans(model)

as.data.frame.predictmeans <- function(x, row.names = NULL, optional = FALSE,
                                       what = c("means", "pairs"), ...) {
  what <- match.arg(what)
  out <- x[[what]]
  if (is.null(out)) stop("There are no ", what, " in this result.", call. = FALSE)
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' Predicted Means and Pairwise Comparisons as Data Frames
#'
#' @description
#' Return the \code{means} or \code{pairs} table of a \code{\link{predictmeans}} result. The tables have fixed column
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
#' A data frame; for \code{coef()} a named vector and for \code{vcov()} a matrix.
#' @author
#' Dongwen Luo, Siva Ganesh and John Koolaard
#' @examples
#'   library(predictmeans)
#'   Oats <- as.data.frame(nlme::Oats)
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#'   pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey", plot = FALSE)
#'   as.data.frame(pm)
#'   head(as.data.frame(pm, what = "pairs"))
#'
#'   # Any contrast of the means from coef() and vcov(): Victory - Golden Rain at nitro 0.6
#'   b <- coef(pm)
#'   L <- (names(b) == "0.6:Victory") - (names(b) == "0.6:Golden Rain")
#'   c(estimate = sum(L*b), SE = sqrt(drop(L %*% vcov(pm) %*% L)))
#' @name tidy.predictmeans
#' @aliases as.data.frame.predictmeans coef.predictmeans vcov.predictmeans
NULL

tidy.predictmeans <- function(x, what = c("means", "pairs"), ...) {
  what <- match.arg(what)
  out <- x[[what]]
  if (is.null(out)) stop("There are no ", what, " in this result.", call. = FALSE)
  out
}
