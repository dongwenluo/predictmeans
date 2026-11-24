## Type II and III tests for the model terms, as car::Anova() computes them.
##
## predictmeans uses these in permmodels(); they are ported from package car
## (Fox & Weisberg, GPL >= 2) for the model classes predictmeans supports, so
## car is not needed. The results equal car::Anova() (see the tests). Classes
## other than lm, glm, lme, gls and merMod get Wald chi-square tests from
## pm_coef() and pm_vcov(), so this also works through the model adapter.
## Kenward-Roger F tests (lmer, test_statistic = "F") use the observed
## information, as everywhere in predictmeans (R/inference.R);
## kr_information = "expected" gives car's (pbkrtest's) values. 'type3_map':
## pm_type3_map() of a fit with the same design, which permutation refits
## (only the response changes) compute once.

## The type of the tests, the one place it is decided: "I"/"II"/"III",
## 1/2/3 or "1"/"2"/"3" (any case) give "I", "II" or "III". 'allowed' limits
## the types a caller supports. Callers pass type[1] (the default vector's first).
pm_anova_type <- function(type, allowed = c("I", "II", "III")) {
  if (length(type) != 1 || is.na(type)) stop("'type' must be one of I, II, III (or 1, 2, 3).", call. = FALSE)
  t <- toupper(trimws(as.character(type)))
  t <- c("1" = "I", "2" = "II", "3" = "III", I = "I", II = "II", III = "III")[t]
  if (is.na(t)) stop("'type' must be one of I, II, III (or 1, 2, 3), not '", type, "'.", call. = FALSE)
  if (!t %in% allowed) {
    stop("Type ", t, " tests are not available here; use type ", paste(allowed, collapse = " or "), ".", call. = FALSE)
  }
  unname(t)
}

pm_anova <- function(model, type = c("II", "III", "2", "3"),
                     test_statistic = c("Chisq", "F", "LR", "Wald"), kr_information = "observed",
                     type3_map = NULL, ddf = c("Kenward-Roger", "Satterthwaite", "KR-average")) {
  ddf <- match.arg(ddf)
  type <- pm_anova_type(type[1], c("II", "III"))
  test_statistic <- match.arg(test_statistic)

  if (inherits(model, "glm")) {
    if (!test_statistic %in% c("LR", "Wald", "F")) test_statistic <- "LR"
    return(switch(test_statistic,
                  LR = if (type == "II") pm_anova_glm_II(model, "LR") else pm_anova_glm_III(model, "LR", type3_map),
                  F = if (type == "II") pm_anova_glm_II(model, "F") else pm_anova_glm_III(model, "F", type3_map),
                  Wald = pm_anova_wald(model, type, layout = "default", type3_map = type3_map)))
  }
  if (inherits(model, "lm")) return(pm_anova_lm(model, type, type3_map))
  if (inherits(model, "lme")) return(pm_anova_wald(model, type, layout = "lme", type3_map = type3_map))
  if (inherits(model, "merMod")) {
    test <- if (test_statistic == "F") "F" else "Chisq"
    return(pm_anova_wald(model, type, layout = "mer", test = test, kr_information = kr_information,
                         type3_map = type3_map, ddf = ddf))
  }
  pm_anova_wald(model, type, layout = "default", type3_map = type3_map)
}

## car::relatives(): the other terms that contain 'term'.
pm_relatives <- function(term, names, factors) {
  is_relative <- function(term1, term2) all(!(factors[, term1] & (!factors[, term2])))
  if (length(names) == 1) return(NULL)
  which_term <- which(term == names)
  (seq_along(names))[-which_term][vapply(names[-which_term], function(t2) is_relative(term, t2), logical(1))]
}

## car::ConjComp(): basis of the part of span(Z) conjugate to span(X) under 'ip'.
pm_conj_comp <- function(X, Z = diag(nrow(X)), ip = diag(nrow(X))) {
  xq <- qr(t(Z) %*% ip %*% X)
  if (xq$rank == 0) return(Z)
  Z %*% qr.Q(xq, complete = TRUE)[, -(1:xq$rank), drop = FALSE]
}

## Type II hypothesis matrix for one term (car's hyp.term()).
pm_type2_hyp <- function(term, names, assign, fac, I.p, not_aliased, V) {
  which_term <- which(term == names)
  subs_term <- which(assign == which_term)
  rel <- pm_relatives(term, names, fac)
  subs_rel <- unlist(lapply(rel, function(r) which(assign == r)))
  h1 <- I.p[subs_rel, , drop = FALSE][, not_aliased, drop = FALSE]
  h2 <- I.p[c(subs_rel, subs_term), , drop = FALSE][, not_aliased, drop = FALSE]
  H <- if (nrow(h1) == 0) h2 else t(pm_conj_comp(t(h1), t(h2), V))
  H[!apply(H, 1, function(x) all(x == 0)), , drop = FALSE]
}

pm_wald_chisq <- function(H, b, V) pm_wald_quad(H, b, V)

pm_response_name <- function(model) deparse(attr(stats::terms(model), "variables")[[2]])

pm_anova_result <- function(df, row_names, col_names, heading) {
  names(df) <- col_names
  row.names(df) <- row_names
  class(df) <- c("anova", "data.frame")
  attr(df, "heading") <- heading
  df
}

## Wald tests: lme (Chisq, Df, Pr), mer (Chisq or KR F), others (Df, Chisq, Pr).
pm_anova_wald <- function(model, type, layout = c("default", "lme", "mer"), test = "Chisq",
                          kr_information = "observed", type3_map = NULL, ddf = "Kenward-Roger") {
  layout <- match.arg(layout)
  tt <- stats::terms(model)
  fac <- attr(tt, "factors")
  intercept <- attr(tt, "intercept") == 1
  b_all <- pm_coef(model)
  not_aliased <- !is.na(b_all)
  b <- b_all[not_aliased]
  map <- NULL
  V <- if (layout == "default" && inherits(model, "glm")) {
    as.matrix(stats::vcov(model))
  } else {
    pm_vcov(model)
  }
  V <- as.matrix(V)
  if (nrow(V) == length(b_all) && length(b) < length(b_all)) V <- V[not_aliased, not_aliased, drop = FALSE]
  PhiA <- NULL
  satt <- layout == "mer" && test == "F" && ddf == "Satterthwaite"
  kr_avg <- layout == "mer" && test == "F" && ddf == "KR-average"
  if (kr_avg) {
    # Genstat: Wald/q from the unadjusted covariance, df Kenward-Roger with the average information
    PhiA <- pm_vcov_kr(model, "average")
  } else if (layout == "mer" && test == "F" && !satt) {
    PhiA <- pm_vcov_kr(model, kr_information)   # "expected": as car::Anova() (pbkrtest)
    V <- as.matrix(PhiA)
    attributes(V) <- attributes(V)[c("dim", "dimnames")]
  }
  # Satterthwaite: the unadjusted covariance, df by lmerTest's rule (pm_inference())
  inf_s <- if (satt) pm_inference(model, df_method = "Satterthwaite", vcov_adjust = FALSE)
  ftest <- function(H, term) {
    if (satt) {
      r <- inf_s$ftest(H, term)
      return(c(stat = r[[3]], ddf = r[[2]]))
    }
    kr <- pm_kr_ftest(model, H, PhiA)
    if (kr_avg) return(c(stat = pm_waldF(H, b, V), ddf = kr[["ddf"]]))
    c(stat = kr[["stat"]], ddf = kr[["ddf"]])
  }
  mm <- stats::model.matrix(model)
  assign <- attr(mm, "assign")
  if (layout == "lme") {
    valid <- colnames(mm) %in% names(stats::coef(model))
    assign[!not_aliased[seq_along(assign)]] <- NA
    assign <- assign[valid]
  } else {
    assign[!not_aliased] <- NA
  }
  p <- length(b_all)
  I.p <- diag(p)
  term_names <- labels(tt)
  names_all <- if (intercept) c("(Intercept)", term_names) else term_names

  if (type == "II") {
    names <- term_names
    stat <- df <- res_df <- pval <- rep(0, length(names))
    for (i in seq_along(names)) {
      H <- pm_type2_hyp(names[i], names, assign, fac, I.p, not_aliased, V)
      if (nrow(H) == 0) {
        stat[i] <- NA; df[i] <- 0
        next
      }
      if (test == "F") {
        kr <- ftest(H, names[i])
        stat[i] <- kr[["stat"]]; res_df[i] <- kr[["ddf"]]
      } else {
        stat[i] <- pm_wald_chisq(H, b, V)
      }
      df[i] <- nrow(H)
    }
  } else {
    map <- if (is.null(type3_map)) pm_type3_map(model, names(b)) else type3_map
    off <- intercept && !isTRUE(map$cells)               # empty cells: no intercept test
    names <- if (off) names_all else term_names
    stat <- df <- res_df <- pval <- rep(0, length(names))
    for (term in seq_along(names)) {
      H <- pm_type3_hyp(map, term - off)       # sum-to-zero coding, whatever the model's
      if (nrow(H) == 0) {
        stat[term] <- NA; df[term] <- 0
        next
      }
      if (test == "F") {
        kr <- ftest(H, names[term])
        stat[term] <- kr[["stat"]]; res_df[term] <- kr[["ddf"]]
      } else {
        stat[term] <- pm_wald_chisq(H, b, V)
      }
      df[term] <- nrow(H)
    }
  }
  pval <- if (test == "F") stats::pf(stat, df, res_df, lower.tail = FALSE) else stats::pchisq(stat, df, lower.tail = FALSE)
  resp <- paste("Response:", pm_response_name(model))

  if (layout == "mer") {
    kind <- if (test != "F") "Wald chisquare tests" else if (ddf == "KR-average")
      "Wald F tests (unadjusted F) with Kenward-Roger df, average information" else paste0("Wald F tests with ", ddf, " df")
    heading <- c(paste0("Analysis of Deviance Table (Type ", type, " ", kind, ")\n"), resp, pm_cells_note(map))
    if (test == "F") {
      return(pm_anova_result(data.frame(stat, df, res_df, pval), names, c("F", "Df", "Df.res", "Pr(>F)"), heading))
    }
    return(pm_anova_result(data.frame(stat, df, pval), names, c("Chisq", "Df", "Pr(>Chisq)"), heading))
  }
  heading <- c(paste0("Analysis of Deviance Table (Type ", type, " tests)\n"), resp, pm_cells_note(map))
  if (layout == "lme") {
    return(pm_anova_result(data.frame(stat, df, pval), names, c("Chisq", "Df", "Pr(>Chisq)"), heading))
  }
  keep <- df > 0
  pm_anova_result(data.frame(df[keep], stat[keep], pval[keep]), names[keep], c("Df", "Chisq", "Pr(>Chisq)"), heading)
}

## lm: F tests with sums of squares (car's Anova.II.lm / Anova.III.lm).
pm_anova_lm <- function(model, type, type3_map = NULL) {
  if (stats::df.residual(model) == 0) stop("residual df = 0")
  if (stats::deviance(model) < sqrt(.Machine$double.eps)) stop("residual sum of squares is 0 (within rounding error)")
  if (any(is.na(stats::coef(model))) && type == "II") {
    # aliased coefficients: sums of squares by model comparison, as car does
    out <- pm_anova_glm_II(model, "F")
    attr(out, "heading") <- c("Anova Table (Type II tests)\n", attr(out, "heading")[2])
    return(out)
  }
  tt <- stats::terms(model)
  fac <- attr(tt, "factors")
  intercept <- attr(tt, "intercept") == 1
  b_all <- stats::coef(model)
  not_aliased <- !is.na(b_all)
  b <- b_all[not_aliased]
  V <- as.matrix(stats::vcov(model))
  if (nrow(V) == length(b_all)) V <- V[not_aliased, not_aliased, drop = FALSE]
  s2 <- summary(model)$sigma^2
  map <- NULL
  if (type == "III") map <- if (is.null(type3_map)) pm_type3_map(model, names(b)) else type3_map
  assign <- model$assign
  assign[!not_aliased] <- NA
  I.p <- diag(length(b_all))
  error_df <- stats::df.residual(model)
  error_SS <- stats::deviance(model)
  term_names <- labels(tt)
  if (type == "II") {
    names <- term_names
    SS <- df <- rep(0, length(names))
    for (i in seq_along(names)) {
      H <- pm_type2_hyp(names[i], names, assign, fac, I.p, not_aliased, V)
      if (nrow(H) == 0) { SS[i] <- NA; df[i] <- 0; next }
      SS[i] <- s2 * pm_wald_chisq(H, b, V)
      df[i] <- nrow(H)
    }
    error_SS <- s2 * error_df
  } else {
    off <- intercept && !isTRUE(map$cells)               # empty cells: no intercept test
    names <- if (off) c("(Intercept)", term_names) else term_names
    SS <- df <- rep(0, length(names))
    for (term in seq_along(names)) {
      H <- pm_type3_hyp(map, term - off)       # sum-to-zero coding, whatever the model's
      if (nrow(H) == 0) { SS[term] <- NA; df[term] <- 0; next }
      SS[term] <- s2 * pm_wald_chisq(H, b, V)
      df[term] <- nrow(H)
    }
  }
  f <- (SS/df)/(error_SS/error_df)
  pval <- stats::pf(f, df, error_df, lower.tail = FALSE)
  res <- data.frame(c(SS, error_SS), c(df, error_df), c(f, NA), c(pval, NA))
  pm_anova_result(res, c(names, "Residuals"), c("Sum Sq", "Df", "F value", "Pr(>F)"),
                  c(paste0("Anova Table (Type ", type, " tests)\n"), paste("Response:", pm_response_name(model)),
                    pm_cells_note(map)))
}

## Pieces of a glm (or lm treated as a gaussian glm) needed for refits.
pm_glm_parts <- function(model) {
  X <- stats::model.matrix(model)
  y <- model$y
  if (is.null(y)) y <- stats::model.response(stats::model.frame(model), "numeric")
  wt <- model$prior.weights
  if (is.null(wt)) wt <- stats::weights(model)
  if (is.null(wt)) wt <- rep(1, length(y))
  fam <- if (inherits(model, "glm")) model$family else stats::gaussian()
  ctrl <- if (!is.null(model$control)) model$control else stats::glm.control()
  list(X = X, y = y, wt = wt, offset = model$offset, family = fam, control = ctrl, asgn = attr(X, "assign"))
}

## glm type II: LR or F tests by refitting without the term (car's
## Anova.II.LR.glm / Anova.II.F.glm; F uses the Pearson error estimate).
pm_anova_glm_II <- function(model, test = c("LR", "F")) {
  test <- match.arg(test)
  tt <- stats::terms(model)
  fac <- attr(tt, "factors")
  names <- labels(tt)
  g <- pm_glm_parts(model)
  which_nms <- function(name) which(g$asgn == which(names == name))
  fit <- function(exclude) {
    stats::glm.fit(g$X[, -exclude, drop = FALSE], g$y, g$wt, offset = g$offset,
                   family = g$family, control = g$control)
  }
  n <- length(names)
  stat <- df <- rep(0, n)
  if (test == "LR") {
    dispersion <- summary(model, corr = FALSE)$dispersion
  } else {
    df_res <- stats::df.residual(model)
    error_SS <- sum(stats::residuals(model, "pearson")^2, na.rm = TRUE)
    dispersion <- error_SS/df_res
  }
  for (term in seq_len(n)) {
    rels <- names[pm_relatives(names[term], names, fac)]
    mod1 <- fit(as.vector(unlist(lapply(c(names[term], rels), which_nms))))
    dev1 <- stats::deviance(mod1)
    if (length(rels) == 0) {
      dev2 <- stats::deviance(model)
      dfres2 <- stats::df.residual(model)
    } else {
      mod2 <- fit(as.vector(unlist(lapply(rels, which_nms))))
      dev2 <- stats::deviance(mod2)
      dfres2 <- stats::df.residual(mod2)
    }
    df[term] <- stats::df.residual(mod1) - dfres2
    stat[term] <- if (df[term] == 0) NA else (dev1 - dev2)
  }
  resp <- paste("Response:", pm_response_name(model))
  if (test == "LR") {
    LR <- stat/dispersion
    p <- stats::pchisq(LR, df, lower.tail = FALSE)
    return(pm_anova_result(data.frame(LR, df, p), names, c("LR Chisq", "Df", "Pr(>Chisq)"),
                           c("Analysis of Deviance Table (Type II tests)\n", resp)))
  }
  f <- stat/(dispersion * df)
  p <- stats::pf(f, df, df_res, lower.tail = FALSE)
  res <- data.frame(c(stat, error_SS), c(df, df_res), c(f, NA), c(p, NA))
  pm_anova_result(res, c(names, "Residuals"), c("Sum Sq", "Df", "F value", "Pr(>F)"),
                  c("Analysis of Deviance Table (Type II tests)\n", resp,
                    "Error estimate based on Pearson residuals \n"))
}

## glm type III: LR or F tests from drop1() (car's Anova.III.LR.glm / Anova.III.F.glm).
## The term's columns are dropped from the sum-to-zero design (pm_type3_map()),
## so the tests do not depend on the model's contrasts.
pm_anova_glm_III <- function(model, test = c("LR", "F"), type3_map = NULL) {
  test <- match.arg(test)
  if (any(is.na(stats::coef(model)))) {
    stop("Type III likelihood ratio tests need a complete design; with empty cells (aliased coefficients) use ",
         "test_statistic = \"Wald\" (the hypotheses on the cells with data) or type = \"II\".", call. = FALSE)
  }
  Source <- labels(stats::terms(model))
  n <- length(Source)
  stat <- df <- rep(0, n)
  g <- pm_glm_parts(model)
  map <- if (is.null(type3_map)) pm_type3_map(model, names(stats::coef(model))) else type3_map
  for (term in seq_len(n)) {
    drop <- map$assign == term
    f1 <- stats::glm.fit(map$Xs[, !drop, drop = FALSE], g$y, g$wt, offset = g$offset, family = g$family,
                         control = g$control)
    df[term] <- sum(drop)
    stat[term] <- stats::deviance(f1)
  }
  resp <- paste("Response:", pm_response_name(model))
  if (test == "LR") {
    dispersion <- summary(model, corr = FALSE)$dispersion
    LR <- ifelse(df == 0, NA, stat/dispersion - stats::deviance(model)/dispersion)
    p <- stats::pchisq(LR, df, lower.tail = FALSE)
    return(pm_anova_result(data.frame(LR, df, p), Source, c("LR Chisq", "Df", "Pr(>Chisq)"),
                           c("Analysis of Deviance Table (Type III tests)\n", resp)))
  }
  df_res <- stats::df.residual(model)
  error_SS <- sum(stats::residuals(model, "pearson")^2, na.rm = TRUE)
  dispersion <- error_SS/df_res
  SS <- stat - stats::deviance(model)
  f <- ifelse(df == 0, NA, (SS/df)/dispersion)
  p <- stats::pf(f, df, df_res, lower.tail = FALSE)
  res <- data.frame(c(SS, error_SS), c(df, df_res), c(f, NA), c(p, NA))
  pm_anova_result(res, c(Source, "Residuals"), c("Sum Sq", "Df", "F values", "Pr(>F)"),
                  c("Analysis of Deviance Table (Type III tests)\n", resp,
                    "Error estimate based on Pearson residuals \n"))
}

## Type I (sequential) hypotheses of lmerTest: for each term, the rows of the
## Doolittle factor of X'X in the term's columns, i.e. the rows of the
## Cholesky factor R of X'X scaled to a unit diagonal (the scaling does not
## change a Wald F, but the Satterthwaite df of several rows depend on it).
## Returns a list of L matrices named by term.
pm_type1_L <- function(X, labels, asgn = attr(X, "assign")) {
  R <- chol(crossprod(as.matrix(X)))
  R <- R / diag(R)
  out <- lapply(seq_along(labels), function(j) R[asgn == j, , drop = FALSE])
  names(out) <- labels
  out[vapply(out, nrow, 1L) > 0]
}

## The tables lmerTest prints for an lmer fit, from predictmeans' own
## Satterthwaite df (R/satterthwaite.R; equal to lmerTest's): the
## coefficients (Estimate, Std. Error, df, t value, Pr(>|t|)) and the type I
## ANOVA (Sum Sq, Mean Sq, NumDF, DenDF, F value, Pr(>F); sums of squares as
## lmerTest's, F x sigma^2 x NumDF). No lmerTest object is needed.
pm_lmer_tables <- function(model) {
  b <- lme4::fixef(model)
  V <- as.matrix(stats::vcov(model))
  parts <- pm_satt_parts(model)
  se <- sqrt(diag(V))
  dfc <- pm_satt_ddf(diag(length(b)), parts)
  tv <- b/se
  coefT <- cbind(Estimate = b, "Std. Error" = se, df = dfc, "t value" = tv, "Pr(>|t|)" = 2*stats::pt(-abs(tv), dfc))
  rownames(coefT) <- names(b)
  Ls <- pm_type1_L(lme4::getME(model, "X"), attr(stats::terms(model), "term.labels"))
  Fv <- pm_type1_F(Ls, b, V)[, 1]
  q <- vapply(Ls, nrow, 1L)
  ddf <- vapply(Ls, function(L) pm_satt_Fddf(L, parts), numeric(1))
  ms <- Fv * stats::sigma(model)^2
  aT <- data.frame(ms * q, ms, q, ddf, Fv, stats::pf(Fv, q, ddf, lower.tail = FALSE), row.names = names(Ls))
  names(aT) <- c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "Pr(>F)")
  attr(aT, "heading") <- "Type I Analysis of Variance Table with Satterthwaite's method"
  class(aT) <- c("anova", "data.frame")
  list(coef = coefT, anova = aT)
}

## Type I (sequential) F tests of an lmer fit. basis = "design" (default):
## hypotheses built from X'X (pm_type1_L()), as SAS PROC MIXED HTYPE=1 and
## lmerTest; they depend on the design only. basis = "fitted": from
## X'V^-1X, the GLS fit (the rows of the Cholesky factor of vcov^-1 that go
## with each term), as Genstat's "sequentially adding terms": each term
## after the earlier ones in the mixed model's own fit; the parts add up to
## the whole Wald statistic, and the hypotheses move with the estimated
## variances. The two agree for the first and the last term
## (dev/sas/type1-exam8.2.sas). ddf: "Kenward-Roger" (adjusted covariance and
## scaled F, = SAS DDFM=KR), "Satterthwaite" (= SAS DDFM=SATTERTHWAITE),
## "KR-average" (Wald/q unadjusted with Kenward-Roger df from the average
## information, as Genstat). A failed adjustment falls back to Satterthwaite.
pm_lmer_type1 <- function(model, ddf = "Kenward-Roger", basis = "design") {
  X <- lme4::getME(model, "X")
  labs <- attr(stats::terms(model), "term.labels")
  b <- lme4::fixef(model)
  V <- as.matrix(stats::vcov(model))
  Ls <- if (basis == "fitted") pm_type1_L_fitted(V, attr(X, "assign"), labs) else pm_type1_L(X, labs)
  PhiA <- switch(ddf, "Kenward-Roger" = pm_try_vcov_kr(model),
                 "KR-average" = pm_try_vcov_kr(model, "average"), NULL)
  if (is.null(PhiA)) ddf <- "Satterthwaite"
  parts <- if (ddf == "Satterthwaite") pm_satt_parts(model)
  ft <- vapply(Ls, function(L) {
    if (ddf == "Satterthwaite") return(c(pm_waldF(L, b, V), nrow(L), pm_satt_Fddf(L, parts)))
    k <- pm_kr_ftest(model, L, PhiA)
    c(if (ddf == "KR-average") pm_waldF(L, b, V) else k[["stat"]], k[["ndf"]], k[["ddf"]])
  }, numeric(3))
  Fv <- ft[1, ]; q <- ft[2, ]; dd <- ft[3, ]
  ms <- Fv * stats::sigma(model)^2
  aT <- data.frame(ms * q, ms, q, dd, Fv, stats::pf(Fv, q, dd, lower.tail = FALSE), row.names = names(Ls))
  names(aT) <- c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "Pr(>F)")
  method <- switch(ddf, "Kenward-Roger" = "Kenward-Roger's method",
                   "KR-average" = "Kenward-Roger df (average information), unadjusted F",
                   "Satterthwaite's method")
  attr(aT, "heading") <- paste0("Type I Analysis of Variance Table with ", method,
                                if (basis == "fitted") "\nSequential hypotheses from the fitted model (X'V^-1 X), as Genstat")
  class(aT) <- c("anova", "data.frame")
  aT
}

## Type I hypotheses from the GLS fit: with R the upper Cholesky factor of
## vcov^-1 (= X'V^-1X / sigma^2), (R b)[term k] are the sequential parts, so
## the rows of R for term k are its hypothesis. 'assign': term of each
## coefficient (0 = intercept).
pm_type1_L_fitted <- function(V, assign, labs) {
  R <- chol(solve(V))
  dimnames(R) <- list(colnames(V), colnames(V))
  k <- sort(unique(assign[assign > 0]))
  stats::setNames(lapply(k, function(j) R[assign == j, , drop = FALSE]), labs[k])
}

## Wald F of each hypothesis L b = 0 (q = rows of L), with covariance V.
## Returns a one-column matrix "F value" with a row per hypothesis.
pm_type1_F <- function(Ls, b, V) {
  V <- as.matrix(V)
  Fv <- vapply(Ls, function(L) pm_waldF(L, b, V), numeric(1))
  matrix(Fv, ncol = 1, dimnames = list(names(Ls), "F value"))
}

## Type III hypotheses that do not depend on the model's contrasts.
##
## A term's Type III hypothesis is that its coefficients are zero when the
## factors are coded with sum-to-zero contrasts (SAS's Type 3 for complete
## data; car::Anova() gives it only for models fitted with such contrasts:
## with R's default treatment contrasts it tests each main effect at the
## first level of the other factors). Both codings span the same fitted
## values, so the sum-to-zero coefficients are a linear function of the
## model's own, beta_sum = B beta with B = (X_sum' X_sum)^-1 X_sum' X, and the
## hypothesis of term j is B[term j, ] beta = 0: no refit.
## Returns list(B, assign (term of each row of B; 0 = intercept), Xs (the
## sum-to-zero design, rows as the model frame)).
pm_type3_map <- function(model, coef_names) {
  tt <- if (inherits(model, "merMod")) stats::terms(model, fixed.only = TRUE) else stats::terms(model)
  tt <- stats::delete.response(tt)
  # an offset is not part of any hypothesis, and glmmTMB's frame does not hold
  # its variables (only "offset(o)"): leave it out of the rebuilt design
  if (!is.null(attr(tt, "offset"))) {
    tl <- attr(tt, "term.labels")
    tt <- stats::terms(stats::reformulate(if (length(tl)) tl else "1", intercept = attr(tt, "intercept") == 1))
  }
  dat <- if (inherits(model, c("lme", "gls"))) pm_get_data(model) else stats::model.frame(model)
  mf <- stats::model.frame(tt, dat, na.action = stats::na.omit)
  facs <- names(mf)[vapply(mf, function(v) is.factor(v) || is.character(v), logical(1))]
  ctr <- pm_contrasts(model)
  if (inherits(model, "merMod")) ctr <- attr(lme4::getME(model, "X"), "contrasts")
  ctr <- ctr[intersect(names(ctr), facs)]
  X <- stats::model.matrix(tt, mf, contrasts.arg = if (length(ctr)) ctr)
  if (!all(coef_names %in% colnames(X))) stop("cannot rebuild the model matrix for Type III tests", call. = FALSE)
  X <- X[, coef_names, drop = FALSE]
  Xs <- stats::model.matrix(tt, mf, contrasts.arg = if (length(facs)) stats::setNames(rep(list("contr.sum"), length(facs)), facs))
  qs <- qr(Xs)
  # empty cells (aliased coefficients): no sum-to-zero map; the hypotheses
  # on the cells with data instead
  if (qs$rank < ncol(Xs) || qr(X)$rank < ncol(Xs)) return(pm_type3_map_cells(model, coef_names))
  B <- qr.coef(qs, X)
  if (max(abs(Xs %*% B - X)) > 1e-8 * max(1, abs(X))) stop("cannot map the model's coefficients to sum-to-zero contrasts")
  dimnames(B) <- list(colnames(Xs), coef_names)
  list(B = B, assign = attr(Xs, "assign"), Xs = Xs)
}

## The Type III hypothesis matrix of term j (0 = intercept), from pm_type3_map().
pm_type3_hyp <- function(map, j) {
  if (isTRUE(map$cells)) return(if (j == 0) matrix(0, 0, length(map$coef)) else map$hyp[[j]])
  H <- map$B[map$assign == j, , drop = FALSE]
  H[apply(abs(H), 1, max) > 1e-12, , drop = FALSE]
}

## Type III contrasts of a nested term: one whose factors are not all in
## lower-order terms of the model (contcompf:standard without standard, or
## A + A:B), as SAS's B(A). The factors without a lower-order term are nested
## in the others: the hypothesis is equal cells within each combination of
## the outer factors, on the estimable cells (combinations absent by design
## are not empty cells). NULL for a term that is not nested (crossed rule).
pm_type3_nested_C <- function(term, labs, fct, est) {
  vars <- unlist(strsplit(term, ":"))
  if (length(vars) < 2) return(NULL)
  lower <- labs[vapply(strsplit(labs, ":"), function(s) length(s) < length(vars) && all(s %in% vars), logical(1))]
  outer <- intersect(vars, unlist(strsplit(lower, ":")))
  if (length(outer) == length(vars)) return(NULL)
  key <- if (length(outer)) do.call(paste, c(lapply(fct[outer], as.character), sep = "\r")) else rep("", length(est))
  rows <- lapply(unique(key[est]), function(g) {
    i <- which(est & key == g)
    if (length(i) < 2) return(NULL)
    m <- matrix(0, length(i) - 1, length(est))
    for (k in seq_len(length(i) - 1)) m[k, i[c(k, k + 1)]] <- c(1, -1)
    m
  })
  do.call(rbind, c(list(matrix(0, 0, length(est))), rows))
}

## The heading line of a Type III table on the cells with data.
pm_cells_note <- function(map) {
  if (!isTRUE(map$cells) || isFALSE(map$empty)) return(character(0))
  "\nEmpty cells: each term is tested on the cells with data (see ?anova.lmerMod)"
}

## Type III hypotheses when cells are empty (aliased coefficients), where the
## sum-to-zero map does not exist. Each term's hypothesis is written on its
## means (Kmatrix(), equal weights), as for complete data:
##  * its contrasts: for each factor of the term the differences of
##    consecutive levels, multiplied across the factors (main effects: the
##    differences among the term's means; interactions: the interaction
##    contrasts of the cells);
##  * those of them that use no empty cell of the term (an interaction with
##    an empty cell loses one df per independent empty cell, as SAS Type III
##    and IV);
##  * means that average over empty cells of other factors: compared on the
##    cells all the means have (pm_common_rows(), R/empty-cells.R), the
##    hypotheses of pm_joint().
## A term that is not made of factors, or whose means have no cell in common,
## is not tested (NA). Returns list(cells = TRUE, hyp (per term, in the
## columns 'coef_names'), coef, over (the common cells, per term)).
pm_type3_map_cells <- function(model, coef_names) {
  tt <- if (inherits(model, "merMod")) stats::terms(model, fixed.only = TRUE) else stats::terms(model)
  labs <- labels(tt)
  dat <- tryCatch(pm_get_data(model), error = function(e) stats::model.frame(model))
  over <- stats::setNames(rep(NA_character_, length(labs)), labs)
  hyp <- lapply(labs, function(term) {
    none <- matrix(0, 0, length(coef_names), dimnames = list(NULL, coef_names))
    vars <- unlist(strsplit(term, ":"))
    if (!all(vars %in% names(dat)) || !all(vapply(dat[vars], function(v) is.factor(v) || is.character(v), logical(1)))) {
      return(none)
    }
    km <- suppressWarnings(suppressMessages(Kmatrix(model, term)))
    K <- km$K
    if (!all(coef_names %in% colnames(K))) return(none)
    K <- K[, coef_names, drop = FALSE]
    est <- attr(km$K, "estimable")
    if (is.null(est)) est <- rep(TRUE, nrow(K))
    est <- est & !apply(K, 1, function(r) any(is.na(r)))
    K[!est, ] <- 0
    fct <- km$fctnames[vars]
    C <- pm_type3_nested_C(term, labs, fct, est)
    if (is.null(C)) {
      # the term's contrasts: differences of consecutive levels, multiplied
      D <- lapply(vars, function(v) {
        lv <- levels(factor(fct[[v]]))
        m <- matrix(0, length(lv) - 1, length(lv))
        if (length(lv) > 1) for (k in seq_len(length(lv) - 1)) m[k, c(k, k + 1)] <- c(1, -1)
        list(m = m, idx = match(as.character(fct[[v]]), lv))
      })
      combos <- expand.grid(lapply(D, function(d) seq_len(nrow(d$m))))
      if (!nrow(combos)) return(none)
      C <- t(vapply(seq_len(nrow(combos)), function(r) {
        Reduce(`*`, lapply(seq_along(D), function(f) D[[f]]$m[combos[r, f], D[[f]]$idx]))
      }, numeric(nrow(K))))
      if (nrow(K) == 1) C <- matrix(C, ncol = 1)
      # only the contrasts that use no empty cell of the term
      if (any(!est)) {
        W <- MASS::Null(C[, !est, drop = FALSE])
        if (!ncol(W)) return(none)
        C <- t(W) %*% C
      }
    } else if (!nrow(C)) {
      return(none)
    }
    # means over empty cells of other factors: on the cells all of them have
    cr <- pm_common_rows(km$empty_cells, C, C %*% K, coef_names, together = TRUE)
    if (any(cr$na)) {
      over[term] <<- "none"
      return(none)
    }
    if (any(cr$sub)) over[term] <<- cr$over[which(cr$sub)[1]]
    H <- cr$L
    q <- qr(t(H))
    H[q$pivot[seq_len(q$rank)], , drop = FALSE]
  })
  list(cells = TRUE, hyp = hyp, coef = coef_names, over = over,
       empty = identical(as.vector(pm_design_balance(model)), "empty"))
}

#' ANOVA Tables for lmer Fits with Kenward-Roger or Satterthwaite df
#'
#' @description
#' \code{anova()} of a single \code{lmer} fit (and of a gaussian \code{\link{semireg}} fit, \code{fit$semer}) gives
#' F tests with denominator degrees of freedom, computed by the package's own code: Type III tests with
#' Kenward-Roger df by default, as SAS PROC MIXED \code{DDFM=KR}. This replaces the tables of \pkg{lmerTest}, which
#' predictmeans attached before 2.0. With several models, \code{anova(m1, m2)} is lme4's likelihood ratio comparison,
#' unchanged.
#' @details
#' Type III hypotheses do not depend on the contrasts the model was fitted with (sum-to-zero coding, as SAS Type 3).
#' With empty cells (aliased coefficients) each term is tested on the cells with data: an interaction by its
#' contrasts that use no empty cell (one df fewer per empty cell), a main effect whose means average over empty
#' cells on the cells its means share, as \code{\link{pm_joint}}; the table says so, and the intercept is not tested;
#' Type II are those of \code{car::Anova()}; Type I (sequential) are by default those of SAS (\code{HTYPE=1}) and
#' lmerTest, built from \eqn{X'X}, so they depend on the design only. Genstat's "sequentially adding terms" builds
#' them from \eqn{X'V^{-1}X} (\code{type1_basis = "fitted"}): each term after the earlier ones in the mixed model's
#' own fit, with parts that add up to the whole Wald statistic, but hypotheses that move with the estimated
#' variances. With unbalanced data the two differ for the terms between the first and the last.
#' Kenward-Roger F tests use the adjusted covariance of the fixed effects (\code{pbkrtest}'s, observed information);
#' a variance component estimated as zero is left out of the adjustment, as in SAS. \code{ddf = "lme4"} gives lme4's
#' own table (sums of squares and F values without df or p-values).
#'
#' With \pkg{lmerTest} attached, \code{lmer()} is lmerTest's and its fits are of class \code{lmerModLmerTest}: their
#' \code{anova()} and \code{summary()} are lmerTest's (Satterthwaite by default), and this method applies to fits of
#' \code{lme4::lmer()}. Everything else in predictmeans works the same with either class.
#' @param object An \code{lmer} fit (class \code{lmerMod}).
#' @param ... Other fitted models: then lme4's comparison of the models.
#' @param type \code{"III"} (default), \code{"II"} or \code{"I"}; also \code{3}, \code{2}, \code{1}.
#' @param ddf \code{"Kenward-Roger"} (default; adjusted covariance and scaled F, as SAS \code{DDFM=KR}),
#'   \code{"Satterthwaite"}, \code{"KR-average"} (F = Wald/numerator df from the unadjusted covariance, with
#'   Kenward-Roger df from the average information, as Genstat) or \code{"lme4"}.
#' @param type1_basis For \code{type = "I"}: \code{"design"} (default) builds the sequential hypotheses from
#'   \eqn{X'X}, as SAS and lmerTest; \code{"fitted"} from \eqn{X'V^{-1}X}, as Genstat's "sequentially adding
#'   terms". \code{type = "I", type1_basis = "fitted", ddf = "KR-average"} gives Genstat's table, and
#'   \code{type = "II", ddf = "KR-average"} its "dropping individual terms".
#' @param advice Logical: add a note on the choice of type when the data are unbalanced or have empty cells
#'   (default \code{TRUE}).
#' @section Which type of test:
#' With balanced data (equal numbers in all cells) Types I, II and III give the same tests. With unbalanced data:
#' \describe{
#' \item{Type I}{sequential: each term adjusted for the terms before it, so the order in the formula matters. Use
#'   it when the order is meaningful, e.g. a block or a covariate to be removed first.}
#' \item{Type II}{each term adjusted for all terms that do not contain it (the principle of marginality). The most
#'   powerful tests of main effects when the interactions are negligible; with an interaction present the main
#'   effects are not well defined.}
#' \item{Type III}{each term adjusted for all others, testing equally weighted cell means; meaningful with an
#'   interaction present, and here independent of the model's contrasts. A covariate in an interaction makes the
#'   other terms be tested at covariate = 0: centre it.}
#' }
#' With empty cells no type answers the question of a main effect over all cells: Type III uses the cells with
#' data (and says so). State the comparisons of interest on the cell means instead, with \code{\link{pm_joint}} or
#' \code{\link{predictmeans}}.
#' @return An \code{anova} table: for each term the numerator df, the denominator df (column heading names the
#'   method), the F value and its p-value. With \code{advice = TRUE} a note on the type is kept in
#'   \code{attr(, "advice")} (read it with \code{\link{advice}}).
#' @seealso \code{\link{anova.semireg}}, \code{\link{pm_joint}}, \code{\link{df_term}}, \code{\link{permmodels}}
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = Oats)
#' anova(fm)
#' anova(fm, ddf = "Satterthwaite")
#' anova(fm, type = "I")
#' @exportS3Method anova lmerMod
anova.lmerMod <- function(object, ..., type = c("III", "II", "I"),
                          ddf = c("Kenward-Roger", "Satterthwaite", "KR-average", "lme4"),
                          type1_basis = c("design", "fitted"), advice = TRUE) {
  # several models: lme4's comparison, by NextMethod() so that lme4 sees the
  # original call (model names), also when lmerTest's anova() hands on to here
  if (length(list(...))) return(NextMethod())
  type <- pm_anova_type(type[1])
  ddf <- match.arg(ddf)
  type1_basis <- match.arg(type1_basis)
  if (ddf == "lme4") return(utils::getS3method("anova", "merMod", envir = asNamespace("lme4"))(object))
  out <- if (type == "I") pm_anova_type1(object, ddf, type1_basis) else pm_anova(object, type = type, test_statistic = "F", ddf = ddf)
  if (isTRUE(advice)) out <- pm_anova_add_advice(out, object, type)
  out
}

## Type I (sequential) table of a mixed model: lmer, lmerTest's Satterthwaite
## F; other models, sequential Wald chi-square tests on the hypotheses of
## pm_type1_L() (no refits).
pm_anova_type1 <- function(model, ddf = "Kenward-Roger", basis = "design") {
  if (inherits(model, "lmerMod")) return(pm_lmer_type1(model, ddf, basis))
  X <- if (inherits(model, "merMod")) lme4::getME(model, "X") else stats::model.matrix(model)
  tl <- attr(stats::terms(if (inherits(model, "merMod")) stats::formula(model, fixed.only = TRUE) else model),
             "term.labels")
  b <- pm_coef(model)
  b <- b[!is.na(b)]
  V <- as.matrix(pm_vcov(model))
  Ls <- pm_type1_L(X[, names(b), drop = FALSE], tl, attr(X, "assign")[match(names(b), colnames(X))])
  q <- vapply(Ls, nrow, 1L)
  chi <- pm_type1_F(Ls, b, V)[, 1] * q
  pm_anova_result(data.frame(chi, q, stats::pchisq(chi, q, lower.tail = FALSE)), names(Ls),
                  c("Chisq", "Df", "Pr(>Chisq)"),
                  c("Analysis of Deviance Table (Type I Wald chisquare tests)\n",
                    paste("Response:", pm_response_name(model))))
}

## Balance of the fixed factors' cells: "balanced" (all cells observed,
## equal counts), "unbalanced" or "empty" (a cell of the crossed fixed
## factors has no data). NA when there are fewer than two factors in the
## fixed terms (one factor: types I, II and III agree for it anyway).
## Attribute "covariate": a numeric variable is in an interaction term.
pm_design_balance <- function(model) {
  tt <- tryCatch(stats::terms(if (inherits(model, "merMod")) stats::formula(model, fixed.only = TRUE) else model),
                 error = function(e) NULL)
  if (is.null(tt)) return(NA_character_)
  dat <- tryCatch(if (inherits(model, c("lme", "gls"))) pm_get_data(model) else stats::model.frame(model),
                  error = function(e) NULL)
  if (is.null(dat)) return(NA_character_)
  vars <- intersect(all.vars(stats::delete.response(tt)), names(dat))
  isfac <- vapply(dat[vars], function(v) is.factor(v) || is.character(v), logical(1))
  facs <- vars[isfac]
  inter <- grep(":", attr(tt, "term.labels"), value = TRUE)
  cov_in <- any(vapply(vars[!isfac], function(v) any(vapply(strsplit(inter, ":"), function(s) v %in% s,
                                                            logical(1))), logical(1)))
  out <- if (length(facs) < 2) NA_character_ else {
    n <- table(lapply(dat[facs], factor))
    # empty: a cell of an interaction term without data; combinations of
    # factors that only have main effects together, or of nested factors
    # (absent by design), are not cells
    nested <- pm_nested_pairs(dat, facs)
    empty <- any(vapply(pm_term_sets(tt), function(v) {
      v <- intersect(v, facs)
      if (length(v) < 2) return(FALSE)
      g <- as.data.frame(table(lapply(dat[v], factor)), stringsAsFactors = FALSE)
      any(g$Freq == 0 & !pm_absent_by_design(g, nested))
    }, logical(1)))
    if (empty) "empty" else if (length(unique(as.vector(n))) == 1) "balanced" else "unbalanced"
  }
  structure(out, covariate = cov_in, counts = if (length(facs) >= 2) n)
}

## Plain-language advice on the type of tests, for the design at hand; one
## place for the wording (anova(), pm_check()). Rows of pm_advice() (R/advice.R).
pm_anova_advice_rows <- function(model, type) {
  bal <- pm_design_balance(model)
  rows <- list()
  if (isTRUE(attr(bal, "covariate")) && type == "III") {
    rows$cov <- pm_advice("anova_covariate_zero", paste(
      "A numeric covariate is in an interaction: Type III tests of the other terms are at",
      "covariate = 0; centre the covariate if 0 is not a meaningful value."), "caution")
  }
  if (!is.na(bal) && bal != "balanced") {
    rows$type <- pm_advice(paste0("anova_unbalanced_type", type), switch(type,
      I = paste("Unbalanced data: Type I tests are sequential, each term adjusted only for the terms before it, so",
                "the order of the terms in the formula changes the results. Use them only when that order is meaningful",
                "(e.g. a covariate or block first); otherwise Type II or III."),
      II = paste("Unbalanced data: Type II tests each term adjusted for all terms that do not contain it. They are the",
                 "more powerful tests of main effects when the interactions are negligible; when an interaction is",
                 "present the main effects are hard to interpret, so look at the interaction first."),
      III = paste("Unbalanced data: Type III tests compare equally weighted cell means (contrasts as sum-to-zero),",
                  "meaningful with an interaction present. With no interaction, Type II tests of main effects are",
                  "more powerful.")), if (type == "I") "caution" else "note")
    if (bal == "empty") {
      rows$empty <- pm_advice("anova_empty_cells", paste(
        "Empty cells: no type of test answers the question of a main effect over all cells.",
        "Type III tests use the cells with data; better, state the comparisons of interest",
        "on the cell means (pm_joint(), predictmeans(pairwise = TRUE))."), "caution")
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, unname(rows))
}

## The texts of the advice (character; kept for pm_check() and tests).
pm_anova_advice <- function(model, type) {
  r <- pm_anova_advice_rows(model, type)
  if (is.null(r)) character(0) else r$text
}

## Store the advice on an anova table and, unless
## options(predictmeans.advice = FALSE), append it to the table's heading
## (print.anova() prints the heading).
pm_anova_add_advice <- function(tab, model, type) {
  rows <- tryCatch(pm_anova_advice_rows(model, type), error = function(e) NULL)
  if (is.null(rows)) return(tab)
  tab <- pm_add_advice(tab, rows)
  txt <- pm_format_advice(pm_advice_get(tab))
  if (nzchar(txt)) attr(tab, "heading") <- c(attr(tab, "heading"), txt)
  tab
}

## lme4's own table for one model (sequential sums of squares and F values,
## no df), whatever anova.lmerMod() gives: the permutation and bootstrap code
## compare these F values. Other models: stats::anova().
pm_anova_seq <- function(model, ...) {
  if (inherits(model, "merMod")) return(utils::getS3method("anova", "merMod", envir = asNamespace("lme4"))(model, ...))
  stats::anova(model, ...)
}
