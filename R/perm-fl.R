## Freedman-Lane permutation tests (permmodels(method = "freedman-lane")) for
## every model class, from one linear model:
##   lm, aov      the (weighted) least-squares model itself;
##   lmer, lme,   the GLS model with the covariance of the data V fixed at
##   gls          the fit: X and y decorrelated with the order-invariant
##                square root of V (R/whiten.R), so least squares on them
##                gives the GLS estimates and the model's Wald F and t;
##   glm, glmer,  the working (link-scale) model of pm_working(), decorrelated
##   glmmTMB      the same way (V = A A' + diag(working variances)).
## For each test (a term of the ANOVA table, or one coefficient) the reduced
## model without it (for type I and II tests also without the terms that
## contain it) gives fitted values and residuals; the residuals are
## permuted, added back, and the statistic computed again. The statistics
## come from QR projections of whole blocks of permuted responses, so no model
## is refitted. The same permutations are used for every test. The result
## has 'perm', a permutation object with which predictmeans() and
## contrastmeans() test pairs and contrasts in the same way (pm_perm_fl_t()).
## dev/permutation-plan.md, steps 3 and 6-8.

## The linear model for the permutations: X (non-aliased columns), z, the
## term each column belongs to (assign), the terms, and 'ord' when the rows
## were reordered (nlme models, by group).
pm_fl_design <- function(model) {
  if (inherits(model, "lm") && !inherits(model, "glm")) {
    mf <- stats::model.frame(model)
    y <- stats::model.response(mf)
    if (!is.null(dim(y))) stop("method = \"freedman-lane\" needs a single response.", call. = FALSE)
    X <- stats::model.matrix(model)
    asgn <- attr(X, "assign")
    b <- stats::coef(model)
    ok <- !is.na(b)
    X <- X[, ok, drop = FALSE]
    asgn <- asgn[ok]
    off <- stats::model.offset(mf)
    z <- if (is.null(off)) y else y - off
    w <- stats::weights(model)
    if (!is.null(w)) {
      keep <- w > 0
      X <- X[keep, , drop = FALSE] * sqrt(w[keep])
      z <- z[keep] * sqrt(w[keep])
    }
    return(list(X = X, z = as.vector(z), asgn = asgn, Terms = stats::terms(model), ord = NULL, kind = "lm"))
  }
  if (inherits(model, c("glm", "glmerMod", "glmmTMB"))) {
    wk <- pm_working(model)
    X0 <- if (inherits(model, "glm")) stats::model.matrix(model) else lme4::getME(model, "X")
    asgn <- attr(X0, "assign")[match(colnames(wk$X), colnames(X0))]
    # rows with zero prior weight carry no information
    keep <- is.finite(wk$resvar) & wk$resvar > 0
    Xk <- wk$X[keep, , drop = FALSE]
    zk <- wk$z[keep]
    if (ncol(wk$A)) {
      W <- pm_whiten(wk$A[keep, , drop = FALSE], wk$resvar[keep])
      X <- W$solve(Xk)
      z <- W$solve(zk)
    } else {
      X <- Xk / sqrt(wk$resvar[keep])
      z <- zk / sqrt(wk$resvar[keep])
    }
    colnames(X) <- colnames(wk$X)
    return(list(X = X, z = as.vector(z), asgn = asgn, Terms = stats::terms(model), ord = NULL, kind = "working"))
  }
  if (inherits(model, "lmerMod")) {
    X0 <- lme4::getME(model, "X")
    resp <- model@resp
    W <- pm_whiten(lme4::getME(model, "Z") %*% lme4::getME(model, "Lambda"), 1 / resp$weights)
    X <- W$solve(as.matrix(X0))
    colnames(X) <- colnames(X0)
    return(list(X = X, z = as.vector(W$solve(resp$y - resp$offset)), asgn = attr(X0, "assign"),
                Terms = stats::terms(model), ord = NULL, kind = "gls"))
  }
  if (inherits(model, c("lme", "gls"))) {
    d <- as.data.frame(pm_get_data(model))
    Terms <- stats::terms(model)
    X0 <- stats::model.matrix(Terms, d)
    b <- pm_fixef(model)
    X <- X0[, names(b), drop = FALSE]
    asgn <- attr(X0, "assign")[match(names(b), colnames(X0))]
    y <- as.numeric(nlme::getResponse(model))
    if (nrow(X) != length(y)) stop("Please remove missing values from the data and refit the model.", call. = FALSE)
    V_list <- build_Sigma_mats(model, sigma_scale = FALSE)
    o <- get_order(A = V_list, B = cbind(X, y))
    ord <- o[[2]]
    W <- pm_whiten_blocks(V_list)
    Xw <- W$solve(o[[1]][, colnames(X), drop = FALSE])
    colnames(Xw) <- colnames(X)
    return(list(X = Xw, z = as.vector(W$solve(o[[1]][, "y"])), asgn = asgn, Terms = Terms, ord = ord, kind = "gls"))
  }
  stop("method = \"freedman-lane\" is available for lm, aov, glm, lme, gls, lmer, glmer and glmmTMB models.",
       call. = FALSE)
}

pm_perm_fl <- function(model, nperm, type, exact, seed, prt, perms = NULL) {
  des <- pm_fl_design(model)
  X <- des$X
  z <- des$z
  asgn <- des$asgn
  n <- nrow(X)
  df_res <- n - ncol(X)
  if (df_res < 1) stop("residual df = 0", call. = FALSE)
  if (!is.null(perms) && !is.null(des$ord)) {
    # the rows were reordered by group: map the user's permutations
    perms <- pm_perm_check_index(perms, n)
    perms <- matrix(order(des$ord)[perms[des$ord, , drop = FALSE]], nrow = n)
  }

  Terms <- des$Terms
  labs <- attr(Terms, "term.labels")
  fac <- attr(Terms, "factors")
  all_cols <- seq_len(ncol(X))
  cols_of <- function(terms) which(asgn %in% terms)
  qr_cache <- list()
  get_qr <- function(cols) {
    key <- paste(cols, collapse = ",")
    if (is.null(qr_cache[[key]])) qr_cache[[key]] <<- qr(X[, cols, drop = FALSE])
    qr_cache[[key]]
  }
  rss <- function(cols, Y) {
    if (!length(cols)) return(colSums(Y^2))
    colSums(qr.resid(get_qr(cols), Y)^2)
  }
  rank_of <- function(cols) if (length(cols)) get_qr(cols)$rank else 0L

  ## one test: columns of the reduced model (red) and the two nested models
  ## whose difference in RSS is the test's sum of squares (A within B)
  tests <- list()
  for (j in seq_along(labs)) {
    if (!any(asgn == j)) next
    rel <- pm_relatives(labs[j], labs, fac)
    red <- setdiff(all_cols, cols_of(c(j, rel)))
    AB <- switch(type,
      I = , "1" = list(cols_of(0:(j - 1)), cols_of(0:j)),
      II = , "2" = list(red, setdiff(all_cols, cols_of(rel))),
      list(setdiff(all_cols, cols_of(j)), all_cols))
    if (type %in% c("III", "3")) red <- AB[[1]]
    tests[[labs[j]]] <- list(red = red, A = AB[[1]], B = AB[[2]])
  }
  for (k in all_cols) {
    if (colnames(X)[k] == "(Intercept)") next
    tests[[paste0("coef:", colnames(X)[k])]] <- list(red = setdiff(all_cols, k), A = setdiff(all_cols, k), B = all_cols)
  }
  stat <- function(tst, Y) {
    dfk <- rank_of(tst$B) - rank_of(tst$A)
    rB <- rss(tst$B, Y)
    rF <- if (length(tst$B) == length(all_cols)) rB else rss(all_cols, Y)
    ((rss(tst$A, Y) - rB) / dfk) / (rF / df_res)
  }
  for (k in names(tests)) {
    r <- tests[[k]]$red
    f0 <- if (length(r)) qr.fitted(get_qr(r), z) else rep(0, n)
    tests[[k]]$f0 <- f0
    tests[[k]]$e0 <- z - f0
    tests[[k]]$value <- stat(tests[[k]], matrix(z))
    tests[[k]]$df <- rank_of(tests[[k]]$B) - rank_of(tests[[k]]$A)
    tests[[k]]$obs <- round(tests[[k]]$value, 6)
    tests[[k]]$count <- 0
  }

  perms <- pm_perm_source(n, seed, perms)
  chunk <- max(1L, min(nperm, floor(2e6 / n)))
  done <- 0L
  while (done < nperm) {
    m <- min(chunk, nperm - done)
    idx <- perms(m)
    for (k in names(tests)) {
      tk <- tests[[k]]
      Y <- tk$f0 + matrix(tk$e0[idx], n)
      tests[[k]]$count <- tk$count + sum(round(stat(tk, Y), 6) >= tk$obs)
    }
    done <- done + m
  }
  pval <- vapply(tests, function(tk) (tk$count + !exact) / (nperm + !exact), numeric(1))

  tabs <- if (des$kind == "lm") pm_fl_tables_lm(model, type, pval) else pm_fl_tables(X, z, df_res, tests, pval, type, des$kind)
  COEF <- tabs$COEF
  ANOVA <- tabs$ANOVA
  if (prt) {
    if (!is.null(COEF)) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEF)
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "Freedman-Lane permutations.\n")
    }
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep = ""))
    print(ANOVA)
    cat("\nNote: Perm_p_value of F test is obtained using", sQuote(nperm), "Freedman-Lane permutations.\n\n")
  }
  list(permlist = NULL, COEFFICENTS = COEF, ANOVA = ANOVA,
       perm = pm_perm_fl_object(X, z, nperm, exact, attr(perms, "start"), colnames(X), attr(perms, "matrix")))
}

## Tables for lm and aov: the model's own coefficient and ANOVA tables with
## the permutation p-values added (the statistics are the same).
pm_fl_tables_lm <- function(model, type, pval) {
  COEF <- NULL
  if (!inherits(model, "aov")) {        # as with method = "raw"
    ctab <- stats::coef(summary(model))
    COEF <- round(cbind(ctab, Perm_p_value = unname(pval[paste0("coef:", rownames(ctab))])), 4)
  }
  aT <- if (type %in% c("I", "1")) stats::anova(model) else pm_anova(model, type = type, test.statistic = "F")
  ANOVA <- round(cbind(aT, Perm_p_value = unname(pval[rownames(aT)])), 4)
  attr(ANOVA, "heading") <- attr(aT, "heading")
  ANOVA[is.na(ANOVA)] <- ""
  list(COEF = COEF, ANOVA = ANOVA)
}

## Tables for the other classes, from the decorrelated linear model: the
## GLS (or working-scale) estimates, standard errors and t values, and for
## each term its F statistic, which is the one permuted. For lmer, lme and
## gls these equal the model's Wald t and F tests (unadjusted covariance;
## type I: sequential in the GLS metric, as anova() of lme4 and nlme); for
## GLMs and GLMMs they are Wald tests on the working scale with the
## dispersion estimated from the working residuals.
pm_fl_tables <- function(X, z, df_res, tests, pval, type, kind) {
  q <- qr(X)
  b <- qr.coef(q, z)
  s2 <- sum(qr.resid(q, z)^2) / df_res
  se <- sqrt(diag(chol2inv(qr.R(q))))[order(q$pivot)] * sqrt(s2)
  tv <- b / se
  p_coef <- unname(pval[paste0("coef:", colnames(X))])
  COEF <- round(cbind(Estimate = b, "Std. Error" = se, "t value" = tv, Perm_p_value = p_coef), 4)
  rownames(COEF) <- colnames(X)
  terms <- names(tests)[!startsWith(names(tests), "coef:")]
  ANOVA <- data.frame(Df = vapply(tests[terms], `[[`, numeric(1), "df"),
                      "F value" = round(vapply(tests[terms], `[[`, numeric(1), "value"), 4),
                      Perm_p_value = round(unname(pval[terms]), 4), check.names = FALSE, row.names = terms)
  tl <- c(I = "I", "1" = "I", II = "II", "2" = "II", III = "III", "3" = "III")[[type]]
  attr(ANOVA, "heading") <- paste0("Type ", tl, " F tests ",
                                   if (kind == "working") "on the working (link) scale, " else "",
                                   "with the covariance of the data fixed at the fitted model")
  list(COEF = COEF, ANOVA = ANOVA)
}
