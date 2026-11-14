## predictmeans() for semireg() fits: means of a factor term including the
## smooth terms (fixed part plus the smooth's random effects, as semipred()),
## with covariates at their means or at the values in 'covariate' (a named
## numeric vector), and pairwise comparisons. Tests are normal-based, as in
## semipred(). Returns the stable tables (means, pairs, settings).
pm_semireg_means <- function(semireg, modelterm, covariate = NULL, level = 0.05, adj = "none", atvar = NULL) {
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
    for (i in cov_names) {     # to the centred / scaled values semireg() used
      m <- attr(mod_df, paste(i, "mean", sep = "_"))
      s <- attr(mod_df, paste(i, "sd", sep = "_"))
      if (!is.null(m)) covV[i] <- covV[i] - m
      if (!is.null(s)) covV[i] <- covV[i]/s
    }
    covV <- if (length(covV) == 1) unname(covV) else matrix(covV, nrow = 1, dimnames = list(NULL, cov_names))
  }
  des <- pm_semi_design(semireg, modelterm, cov_names, covV, NULL, mod_df)
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
  z <- stats::qnorm(1 - level/2)
  means <- cbind(lab, estimate = est, std.error = se, df = Inf, conf.low = est - z*se, conf.high = est + z*se)
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
  n <- length(est)
  if (n > 1) {
    pp <- pm_pair_index(n)
    CM <- pm_pair_matrix(n)
    pe <- drop(CM %*% est)
    pV <- CM %*% V %*% t(CM)
    pse <- sqrt(diag(pV))
    stat <- pe/pse
    p <- 2*stats::pnorm(-abs(stat))
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
    pairs$df <- Inf
    pairs$statistic <- stat
    pairs$p.value <- p
    # all the pairs, or the pairs within each level of atvar (R/p-adjust.R)
    fams <- pm_families(if (is.null(atvar)) rep("all", length(p)) else by, est = !is.na(stat),
                        ij = cbind(pp$j, pp$i), L = CM %*% D)
    p_adj <- pm_adjust_p(adj, stat, p, Inf, fams, R = if (adj == "mvt") stats::cov2cor(pV))
    pairs$p.adj <- p_adj
    pairs$p.method <- if (adj == "tukey") "tukey" else "z"
    if (is.null(atvar)) {             # letters among all means
      P <- matrix(1, n, n, dimnames = list(levs, levs))
      P[cbind(pp$j, pp$i)] <- P[cbind(pp$i, pp$j)] <- p_adj
      means$group <- unname(pm_letter_groups(P, rep(TRUE, n), level))
    }
  }
  settings <- list(modelterm = modelterm, level = level, adjust = adj, atvar = atvar, df_method = "normal (semireg)",
                   vcov_adjust = FALSE, nperm = NA_integer_, covariates = covariate, model_class = "semireg",
                   version = as.character(utils::packageVersion("predictmeans")))
  structure(list(means = means, pairs = pairs, settings = settings), vcov = V,
            class = c("pm_semimeans", "predictmeans"))
}

print.pm_semimeans <- function(x, digits = 4, ...) {
  cat("Predicted means of", sQuote(x$settings$modelterm), "from a semireg model (smooth terms included)\n\n")
  print(x$means, digits = digits, row.names = FALSE)
  if (!is.null(x$pairs)) {
    cat("\nPairwise comparisons (normal-based; p.adj: ", x$settings$adjust, ")\n\n", sep = "")
    print(x$pairs, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
