## Heteroscedasticity-robust Freedman-Lane tests for lm, aov and glm models
## (permmodels(robust = TRUE); plan step 12).
##
## With unequal error variances that the model ignores, permuting residuals
## is not exact, and permutation tests of unstudentised (or classically
## studentised) statistics can reject far too often (19-23% at the 5% level in
## dev/sim-perm-robust.R). Studentising with a heteroscedasticity-consistent
## (sandwich) standard error makes the Freedman-Lane permutation test
## asymptotically valid (DiCiccio and Romano, 2017, JASA 112, 1211-1220), and
## flipping the signs of the reduced model's residuals instead of permuting
## them (a wild bootstrap with Rademacher weights; Davidson and Flachaire 2008;
## Winkler et al. 2014) keeps each residual with its own unit, and so its
## variance: in the simulation it was closer to 5% in every case with unequal
## variances (e.g. 5.8% against 7.3%; skewed errors 9.4% against 15.9%) and
## 5.5% with equal variances. So robust = TRUE flips signs. The statistics are
## HC3 Wald statistics: |t| = |L b| / sqrt(L S L') for one row, F = (L b)'
## (L S L')^(-1) (L b) / q for q rows, S = (X'X)^-1 X' diag(r^2 / (1 - h)^2) X
## (X'X)^-1, computed for blocks of resampled responses at once (no refits).
## Each test is a hypothesis matrix L (type I: rows of the Cholesky factor of
## X'X, as lmerTest; type II: car's hypotheses in the least-squares metric;
## type III: the term's columns; a coefficient or a pair: one row), and its
## Freedman-Lane reduced model is the null space of L. glm: the same on the
## working scale (the IRLS weights fixed at the fit).

## The parts of the least-squares fit used for every statistic.
pm_hc_parts <- function(X) {
  qrX <- qr(X)
  XtXi <- solve(crossprod(X))
  h <- rowSums((X %*% XtXi) * X)
  list(X = X, qr = qrX, XtXi = XtXi, h = pmin(h, 1 - 1e-8), df_res = nrow(X) - qrX$rank)
}

## Squared HC3 residuals of each column of Y (n x B).
pm_hc_w <- function(P, Y) qr.resid(P$qr, Y)^2 / (1 - P$h)^2

## Robust |t| of each row of L for each column of Y (k x B).
pm_hc_absT <- function(P, L, Y) {
  A <- L %*% P$XtXi %*% t(P$X)
  abs(L %*% qr.coef(P$qr, Y)) / sqrt((A^2) %*% pm_hc_w(P, Y))
}

## Robust Wald F of all rows of L (independent rows) for each column of Y.
pm_hc_F <- function(P, L, Y) {
  q <- nrow(L)
  A <- L %*% P$XtXi %*% t(P$X)
  est <- L %*% qr.coef(P$qr, Y)
  W <- pm_hc_w(P, Y)
  if (q == 1) return(as.vector(est^2 / ((A^2) %*% W)))
  vapply(seq_len(ncol(Y)), function(b) {
    S <- A %*% (W[, b] * t(A))
    as.numeric(crossprod(est[, b], solve(S, est[, b]))) / q
  }, numeric(1))
}

## The model with L b = 0: fitted values and residuals of z.
pm_hc_reduced <- function(X, z, L) {
  q <- qr(t(L))
  N <- qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
  f0 <- if (ncol(N)) qr.fitted(qr(X %*% N), z) else rep(0, length(z))
  list(f0 = f0, e0 = z - f0)
}

pm_perm_fl_robust <- function(model, nperm, type, exact, seed, prt, perms = NULL) {
  if (!inherits(model, "lm")) {
    stop("robust = TRUE is for lm, aov and glm models.", call. = FALSE)
  }
  if (!is.null(perms)) {
    stop("'perms' cannot be used with robust = TRUE (the residuals' signs are flipped, not permuted).",
         call. = FALSE)
  }
  exact <- FALSE
  des <- pm_fl_design(model)
  X <- des$X
  z <- des$z
  asgn <- des$asgn
  n <- nrow(X)
  p <- ncol(X)
  P <- pm_hc_parts(X)
  if (P$df_res < 1) stop("residual df = 0", call. = FALSE)
  labs <- attr(des$Terms, "term.labels")
  fac <- attr(des$Terms, "factors")
  I.p <- diag(p)
  L1 <- if (type %in% c("I", "1")) pm_type1_L(X, labs, asgn) else NULL
  tests <- list()
  for (j in seq_along(labs)) {
    cols <- which(asgn == j)
    if (!length(cols)) next
    L <- switch(type,
      I = , "1" = L1[[labs[j]]],
      II = , "2" = pm_type2_hyp(labs[j], labs, asgn, fac, I.p, rep(TRUE, p), P$XtXi),
      I.p[cols, , drop = FALSE])
    tests[[labs[j]]] <- list(L = pm_indep_rows(L), F = TRUE)
  }
  for (k in seq_len(p)) {
    if (colnames(X)[k] == "(Intercept)") next
    tests[[paste0("coef:", colnames(X)[k])]] <- list(L = I.p[k, , drop = FALSE], F = FALSE)
  }
  stat <- function(tk, Y) if (tk$F) pm_hc_F(P, tk$L, Y) else as.vector(pm_hc_absT(P, tk$L, Y))
  for (k in names(tests)) {
    tk <- tests[[k]]
    R0 <- pm_hc_reduced(X, z, tk$L)
    tk$f0 <- R0$f0
    tk$e0 <- R0$e0
    tk$value <- stat(tk, matrix(z))
    tk$df <- nrow(tk$L)
    tests[[k]] <- tk
  }
  # the permuted statistics (kept for permstats() and permcombine())
  Tm <- matrix(NA_real_, length(tests), nperm, dimnames = list(names(tests), NULL))
  src <- pm_perm_source(n, seed, flip = TRUE)
  chunk <- max(1L, min(nperm, floor(2e6 / n)))
  done <- 0L
  while (done < nperm) {
    m <- min(chunk, nperm - done)
    sg <- src(m)
    for (k in names(tests)) {
      tk <- tests[[k]]
      Tm[k, done + seq_len(m)] <- stat(tk, tk$f0 + tk$e0 * sg)
    }
    done <- done + m
  }
  obs <- vapply(tests, `[[`, numeric(1), "value")
  pval <- stats::setNames(pm_perm_p(obs, Tm, exact), names(tests))

  # tables: least-squares (working-scale) estimates with HC3 standard errors
  b <- qr.coef(P$qr, z)
  w0 <- as.vector(pm_hc_w(P, matrix(z)))
  A <- P$XtXi %*% t(X)
  se <- sqrt(as.vector((A^2) %*% w0))
  COEF <- if (inherits(model, "aov")) NULL else {
    out <- round(cbind(Estimate = b, "Robust SE" = se, "t value" = b / se,
                       Perm_p_value = unname(pval[paste0("coef:", colnames(X))])), 4)
    rownames(out) <- colnames(X)
    out
  }
  terms <- names(tests)[!startsWith(names(tests), "coef:")]
  ANOVA <- data.frame(Df = vapply(tests[terms], `[[`, numeric(1), "df"),
                      "F value" = round(vapply(tests[terms], `[[`, numeric(1), "value"), 4),
                      Perm_p_value = round(unname(pval[terms]), 4), check.names = FALSE, row.names = terms)
  tl <- c(I = "I", "1" = "I", II = "II", "2" = "II", III = "III", "3" = "III")[[type]]
  attr(ANOVA, "heading") <- paste0("Type ", tl, " robust (HC3) Wald F tests",
                                   if (des$kind == "working") " on the working (link) scale" else "",
                                   " (Freedman-Lane residuals, sign flips)")
  if (prt) {
    if (!is.null(COEF)) {
      cat("\nCoefficients of (fixed) effects (robust standard errors):\n")
      print(COEF)
      cat("\nNote: Perm_p_value of robust t test is obtained using", sQuote(nperm),
          "sign flips of the Freedman-Lane residuals (wild bootstrap).\n")
    }
    cat("\nANOVA:\n")
    cat(attr(ANOVA, "heading"), "\n", sep = "")
    print(ANOVA)
    cat("\nNote: Perm_p_value of robust F test is obtained using", sQuote(nperm),
        "sign flips of the Freedman-Lane residuals (wild bootstrap).\n\n")
  }
  perm <- pm_perm_fl_object(X, z, nperm, exact, attr(src, "start"), colnames(X))
  perm$robust <- TRUE
  perm$flip <- TRUE
  perm$stats <- pm_perm_statset(model, obs, Tm, pval, vapply(tests, `[[`, numeric(1), "df"),
                                ifelse(startsWith(names(tests), "coef:"), "t", "F"), lapply(tests, `[[`, "L"),
                                P$df_res, n)
  list(permlist = NULL, COEFFICENTS = COEF, ANOVA = ANOVA, perm = perm)
}

## Pairs and contrasts from a robust permutation object: robust |t| of each
## row of L with its own reduced model (as pm_perm_fl_t()), or the joint
## distribution and robust Wald F with one reduced model for the family (as
## pm_perm_fl_T(); also the robust standard errors of the rows).
pm_perm_fl_t_robust <- function(perm, L) {
  X <- perm$X
  z <- perm$z
  n <- nrow(X)
  P <- pm_hc_parts(X)
  red <- lapply(seq_len(nrow(L)), function(j) {
    R0 <- pm_hc_reduced(X, z, L[j, , drop = FALSE])
    R0$obs <- round(as.vector(pm_hc_absT(P, L[j, , drop = FALSE], matrix(z))), 6)
    R0
  })
  count <- valid <- numeric(nrow(L))
  src <- pm_perm_source(n, state = perm$start, flip = TRUE)
  chunk <- max(1L, min(perm$nperm, floor(2e6 / n)))
  done <- 0L
  while (done < perm$nperm) {
    m <- min(chunk, perm$nperm - done)
    sg <- src(m)
    for (j in seq_len(nrow(L))) {
      tj <- round(as.vector(pm_hc_absT(P, L[j, , drop = FALSE], red[[j]]$f0 + red[[j]]$e0 * sg)), 6)
      ok <- !is.na(tj)
      count[j] <- count[j] + sum(tj[ok] >= red[[j]]$obs)
      valid[j] <- valid[j] + sum(ok)
    }
    done <- done + m
  }
  (count + !perm$exact) / (valid + !perm$exact)
}

pm_perm_fl_T_robust <- function(perm, L) {
  X <- perm$X
  z <- perm$z
  n <- nrow(X)
  P <- pm_hc_parts(X)
  Lr <- pm_indep_rows(L)
  R0 <- pm_hc_reduced(X, z, L)
  src <- pm_perm_source(n, state = perm$start, flip = TRUE)
  chunk <- max(1L, min(perm$nperm, floor(2e6 / n)))
  done <- 0L
  Tl <- Fl <- list()
  while (done < perm$nperm) {
    m <- min(chunk, perm$nperm - done)
    Y <- R0$f0 + R0$e0 * src(m)
    Tl[[length(Tl) + 1]] <- pm_hc_absT(P, L, Y)
    Fl[[length(Fl) + 1]] <- pm_hc_F(P, Lr, Y)
    done <- done + m
  }
  # se: robust (HC3) standard errors of the rows, for the confidence limits
  A <- L %*% P$XtXi %*% t(X)
  list(obs = as.vector(pm_hc_absT(P, L, matrix(z))), T = do.call(cbind, Tl),
       Fobs = pm_hc_F(P, Lr, matrix(z)), F = unlist(Fl),
       se = sqrt(as.vector((A^2) %*% pm_hc_w(P, matrix(z)))))
}
