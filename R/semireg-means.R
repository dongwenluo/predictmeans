## predictmeans() for semireg() fits: means of a factor term including the
## smooth terms (fixed part plus the smooth's random effects, as semipred()),
## with covariates at their means or at the values in 'covariate' (a named
## numeric vector), and pairwise comparisons, with the df of pm_sr_ddf()
## (gaussian fits: Satterthwaite for the combinations of fixed and smooth
## coefficients; generalized fits: inner-outer df of the design's strata, z
## without strata), as semipred(). Returns the stable tables (means, pairs,
## settings).
pm_semireg_means <- function(semireg, modelterm, covariate = NULL, level = 0.05, adj = "none", atvar = NULL,
                             unconditional = getOption("predictmeans.semi_unconditional"), vcov_adjust = TRUE) {
  adj <- pm_match_adj(adj)
  mod_df <- semireg$data
  cov_names <- covV <- NULL
  if (!is.null(covariate)) {
    if (!is.numeric(covariate) || is.null(names(covariate))) {
      stop("For semireg models 'covariate' must be a named numeric vector, e.g. c(construction.date = 1980).",
           call. = FALSE)
    }
    cov_names <- names(covariate)
    covV <- covariate
    covV <- if (length(covV) == 1) unname(covV) else matrix(covV, nrow = 1, dimnames = list(NULL, cov_names))
  }
  des <- pm_semi_design(semireg, modelterm, cov_names, covV, NULL, mod_df, unconditional = unconditional)
  vars <- unlist(strsplit(modelterm, ":"))
  lab <- des$pred_df[vars]
  rownames(lab) <- NULL
  levs <- do.call(paste, c(lapply(lab, as.character), sep = ":"))
  D <- as.matrix(des$D_mx)
  C <- as.matrix(des$sig.epsHat^2 * des$sm_Cov)
  est <- as.numeric(D %*% des$beta_u)
  V <- D %*% C %*% t(D)
  dimnames(V) <- list(levs, levs)
  se <- sqrt(diag(V))
  n <- length(est)
  CM <- if (n > 1) pm_pair_matrix(n)
  # df of the means and of their differences (one computation for both)
  dfs <- if (!is.null(des$cov_idx)) pm_sr_ddf(semireg, rbind(D, if (n > 1) CM %*% D), des$cov_idx, modelterm) else
    structure(rep(Inf, n + if (n > 1) nrow(CM) else 0), method = "normal (semireg)")
  df_label <- attr(dfs, "method")
  df_cases <- attr(dfs, "df_cases")
  dfs <- as.numeric(dfs)
  dfs[is.na(dfs)] <- Inf
  # pairs that are contrasts of the fixed effects only (the smooth parts cancel, as between levels of a factor
  # added to the smooth): Kenward-Roger adjusted SEs and df of the gaussian REML fit, as predictmeans() of an
  # lmer fit (decision of 2026-11-26); means, which use the smooth's random coefficients, keep Satterthwaite
  kr <- if (isTRUE(vcov_adjust) && n > 1 && !is.null(des$cov_idx)) pm_sr_kr(semireg, CM %*% D, des$cov_idx)
  if (!is.null(kr)) {
    dfs[-seq_len(n)] <- kr$df
    df_label <- paste0("Kenward-Roger (pairs, fixed effects; adjusted SEs); ", df_label, " (means)")
  }
  df_m <- dfs[seq_len(n)]
  tq <- stats::qt(1 - level/2, df_m)
  means <- cbind(lab, estimate = est, std.error = se, df = df_m, conf.low = est - tq*se, conf.high = est + tq*se,
                 estimable = rep(TRUE, n))           # (the column of predictmeans()' means; summary() reads it)
  rownames(means) <- NULL

  semer <- des$semer
  fam <- tryCatch(stats::family(semer), error = function(e) NULL)
  if (is.list(fam) && !identical(fam$link, "identity")) {
    means$response <- fam$linkinv(est)
    means$response.std.error <- abs(fam$mu.eta(est))*se
    means$response.conf.low <- fam$linkinv(means$conf.low)
    means$response.conf.high <- fam$linkinv(means$conf.high)
  }

  pairs <- NULL
  if (n > 1) {
    pp <- pm_pair_index(n)
    df_p <- dfs[-seq_len(n)]
    pe <- drop(CM %*% est)
    pV <- if (is.null(kr)) CM %*% V %*% t(CM) else kr$V
    pse <- sqrt(diag(pV))
    stat <- pe/pse
    p <- 2*stats::pt(-abs(stat), df_p)
    pairs <- data.frame(level1 = levs[pp$j], level2 = levs[pp$i], stringsAsFactors = FALSE)
    by <- rep(NA_character_, length(p))
    if (!is.null(atvar)) {
      A <- lab[pp$j, atvar, drop = FALSE]
      B <- lab[pp$i, atvar, drop = FALSE]
      same <- Reduce(`&`, lapply(atvar, function(v) as.character(A[[v]]) == as.character(B[[v]])))
      by[same] <- do.call(paste, c(lapply(A[same, , drop = FALSE], as.character), sep = " : "))
      pairs$by <- by
    }
    pairs$estimate <- pe
    pairs$std.error <- pse
    pairs$df <- df_p
    pairs$statistic <- stat
    pairs$p.value <- p
    # all the pairs, or the pairs within each level of atvar (R/p-adjust.R)
    fams <- pm_families(if (is.null(atvar)) rep("all", length(p)) else by, est = !is.na(stat),
                        ij = cbind(pp$j, pp$i), L = CM %*% D)
    R <- if (adj == "mvt") stats::cov2cor(pV)
    p_adj <- pm_adjust_p(adj, stat, p, df_p, fams, R = R)
    pairs$p.adj <- p_adj
    # limits of each pair, and the simultaneous limits that agree with p.adj
    # (the columns of predictmeans()$pairs; R/p-adjust.R)
    ct <- stats::qt(1 - level/2, df_p)
    pairs$conf.low <- pe - ct*pse
    pairs$conf.high <- pe + ct*pse
    pairs$lsd <- ct*pse
    crit_adj <- pm_adjust_crit(adj, level, df_p, fams, R = R)
    if (!is.null(crit_adj)) {
      pairs$adj.conf.low <- pe - crit_adj*pse
      pairs$adj.conf.high <- pe + crit_adj*pse
    }
    pairs$p.method <- if (adj == "tukey") "tukey" else if (any(is.finite(dfs))) "t" else "z"
    # a log link: ratios of the means; a logit link: odds ratios
    kind <- if (is.list(fam)) switch(fam$link, log = "ratio", logit = "odds.ratio", NULL)
    if (!is.null(kind)) {
      pairs[[kind]] <- exp(pe)
      pairs[[paste0(kind, ".conf.low")]] <- exp(pairs$conf.low)
      pairs[[paste0(kind, ".conf.high")]] <- exp(pairs$conf.high)
    }
    if (is.null(atvar)) {             # letters among all means
      P <- matrix(1, n, n, dimnames = list(levs, levs))
      P[cbind(pp$j, pp$i)] <- P[cbind(pp$i, pp$j)] <- p_adj
      means$group <- unname(pm_letter_groups(P, rep(TRUE, n), level))
    }
  }
  zi <- pm_zi_component(semireg)
  if (!is.null(zi)) message(attr(zi, "note"))
  settings <- list(modelterm = modelterm, level = level, adjust = adj, atvar = atvar, df_method = df_label,
                   df_cases = if (is.null(df_cases)) NA_character_ else df_cases,
                   vcov_adjust = !is.null(kr), nperm = NA_integer_, covariates = covariate, weights_by = NA_character_,
                   component = if (is.null(pm_zi_component(semireg))) NA_character_ else "conditional", model_class = "semireg",
                   version = as.character(utils::packageVersion("predictmeans")))
  out <- structure(list(means = means, pairs = pairs, settings = settings), vcov = V,
                   class = c("pm_semimeans", "predictmeans"))
  # (no singular-fit line: a smooth at a straight line is a variance at 0, not random effects to simplify)
  # (ntests 0: the pairs are always built here, not asked for)
  rows <- pm_advice_common(semireg$semer, c(means$df, if (is.data.frame(pairs)) pairs$df), 0, adj, "pairs",
                           settings$df_cases, singular = FALSE)
  if (length(rows)) out <- pm_add_advice(out, do.call(rbind, unname(rows)))
  out
}

#' @exportS3Method print pm_semimeans
print.pm_semimeans <- function(x, digits = 4, ...) {
  cat("Predicted means of", sQuote(x$settings$modelterm), "from a semireg model (smooth terms included)\n\n")
  print(x$means, digits = digits, row.names = FALSE)
  if (!is.null(x$pairs)) {
    cat("\nPairwise comparisons (", if (any(is.finite(x$pairs$df))) paste0("t tests, ", x$settings$df_method, " df")
        else "normal-based",
        "; p.adj: ", x$settings$adjust, ")\n\n", sep = "")
    print(x$pairs, digits = digits, row.names = FALSE)
  }
  pm_print_advice(x)
  invisible(x)
}

## Kenward-Roger for the pairs of semireg means (gaussian lme4 engine): when
## every row of L (columns: the coefficients 'idx' of the fit, fixed effects
## first) uses fixed effects only, list(V: their KR-adjusted covariance, df:
## KR df), from the REML mixed model of the fit (pm_vcov_kr(), R/kr.R); NULL
## otherwise (a row with smooth random coefficients, a glmmTMB or ML fit, or
## a failed adjustment: then Satterthwaite, with a message for the failure).
pm_sr_kr <- function(semireg, L, idx) {
  semer <- semireg$semer
  if (!inherits(semer, "lmerMod") || !isTRUE(lme4::getME(semer, "is_REML"))) return(NULL)
  p <- length(lme4::fixef(semer))
  L <- as.matrix(L)
  fx <- idx <= p
  if (any(abs(L[, !fx, drop = FALSE]) > 1e-10)) return(NULL)
  Lf <- matrix(0, nrow(L), p)
  Lf[, idx[fx]] <- L[, fx, drop = FALSE]
  Vadj <- pm_try_vcov_kr(semer)
  if (is.null(Vadj)) return(NULL)
  V0 <- as.matrix(stats::vcov(semer))
  list(V = Lf %*% as.matrix(Vadj) %*% t(Lf), df = as.numeric(pm_kr_ddf(Lf, V0, Vadj)))
}
