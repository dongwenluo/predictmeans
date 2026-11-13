## Type II and III tests for the model terms, as car::Anova() computes them.
##
## predictmeans uses these in permmodels(); they are ported from package car
## (Fox & Weisberg, GPL >= 2) for the model classes predictmeans supports, so
## car is not needed. The results equal car::Anova() (see the tests). Classes
## other than lm, glm, lme, gls and merMod get Wald chi-square tests from
## pm_coef() and pm_vcov(), so this also works through the model adapter.
## Kenward-Roger F tests (lmer, test.statistic = "F") use the observed
## information, as everywhere in predictmeans (R/inference.R);
## kr_information = "expected" gives car's (pbkrtest's) values. 'type3_map':
## pm_type3_map() of a fit with the same design, which permutation refits
## (only the response changes) compute once.

pm_anova <- function(model, type = c("II", "III", "2", "3"),
                     test.statistic = c("Chisq", "F", "LR", "Wald"), kr_information = "observed",
                     type3_map = NULL) {
  type <- as.character(type)
  type <- match.arg(type)
  type <- c("2" = "II", "3" = "III", II = "II", III = "III")[[type]]
  test.statistic <- match.arg(test.statistic)

  if (inherits(model, "glm")) {
    if (!test.statistic %in% c("LR", "Wald", "F")) test.statistic <- "LR"
    return(switch(test.statistic,
                  LR = if (type == "II") pm_anova_glm_II(model, "LR") else pm_anova_glm_III(model, "LR", type3_map),
                  F = if (type == "II") pm_anova_glm_II(model, "F") else pm_anova_glm_III(model, "F", type3_map),
                  Wald = pm_anova_wald(model, type, layout = "default", type3_map = type3_map)))
  }
  if (inherits(model, "lm")) return(pm_anova_lm(model, type, type3_map))
  if (inherits(model, "lme")) return(pm_anova_wald(model, type, layout = "lme", type3_map = type3_map))
  if (inherits(model, "merMod")) {
    test <- if (test.statistic == "F") "F" else "Chisq"
    return(pm_anova_wald(model, type, layout = "mer", test = test, kr_information = kr_information,
                         type3_map = type3_map))
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

pm_wald_chisq <- function(H, b, V) {
  d <- H %*% b
  as.vector(t(d) %*% solve(H %*% V %*% t(H)) %*% d)
}

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
                          kr_information = "observed", type3_map = NULL) {
  layout <- match.arg(layout)
  tt <- stats::terms(model)
  fac <- attr(tt, "factors")
  intercept <- attr(tt, "intercept") == 1
  b_all <- pm_coef(model)
  not_aliased <- !is.na(b_all)
  if (type == "III" && !all(not_aliased)) stop("there are aliased coefficients in the model")
  b <- b_all[not_aliased]
  V <- if (layout == "default" && inherits(model, "glm")) {
    as.matrix(stats::vcov(model))
  } else {
    pm_vcov(model)
  }
  PhiA <- NULL
  if (layout == "mer" && test == "F") {
    PhiA <- pm_vcov_kr(model, kr_information)   # "expected": as car::Anova() (pbkrtest)
    V <- as.matrix(PhiA)
    attributes(V) <- attributes(V)[c("dim", "dimnames")]
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
        kr <- pm_kr_ftest(model, H, PhiA)
        stat[i] <- kr[["stat"]]; res_df[i] <- kr[["ddf"]]
      } else {
        stat[i] <- pm_wald_chisq(H, b, V)
      }
      df[i] <- nrow(H)
    }
  } else {
    names <- names_all
    stat <- df <- res_df <- pval <- rep(0, length(names))
    map <- if (is.null(type3_map)) pm_type3_map(model, names(b)) else type3_map
    for (term in seq_along(names)) {
      H <- pm_type3_hyp(map, term - intercept)       # sum-to-zero coding, whatever the model's
      if (nrow(H) == 0) {
        stat[term] <- NA; df[term] <- 0
        next
      }
      if (test == "F") {
        kr <- pm_kr_ftest(model, H, PhiA)
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
    kind <- if (test == "F") "Wald F tests with Kenward-Roger df" else "Wald chisquare tests"
    heading <- c(paste0("Analysis of Deviance Table (Type ", type, " ", kind, ")\n"), resp)
    if (test == "F") {
      return(pm_anova_result(data.frame(stat, df, res_df, pval), names, c("F", "Df", "Df.res", "Pr(>F)"), heading))
    }
    return(pm_anova_result(data.frame(stat, df, pval), names, c("Chisq", "Df", "Pr(>Chisq)"), heading))
  }
  heading <- c(paste0("Analysis of Deviance Table (Type ", type, " tests)\n"), resp)
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
  if (type == "III" && !all(not_aliased)) stop("there are aliased coefficients in the model")
  b <- b_all[not_aliased]
  V <- as.matrix(stats::vcov(model))
  s2 <- summary(model)$sigma^2
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
    names <- if (intercept) c("(Intercept)", term_names) else term_names
    SS <- df <- rep(0, length(names))
    for (term in seq_along(names)) {
      H <- pm_type3_hyp(map, term - intercept)       # sum-to-zero coding, whatever the model's
      if (nrow(H) == 0) { SS[term] <- NA; df[term] <- 0; next }
      SS[term] <- s2 * pm_wald_chisq(H, b, V)
      df[term] <- nrow(H)
    }
  }
  f <- (SS/df)/(error_SS/error_df)
  pval <- stats::pf(f, df, error_df, lower.tail = FALSE)
  res <- data.frame(c(SS, error_SS), c(df, error_df), c(f, NA), c(pval, NA))
  pm_anova_result(res, c(names, "Residuals"), c("Sum Sq", "Df", "F value", "Pr(>F)"),
                  c(paste0("Anova Table (Type ", type, " tests)\n"), paste("Response:", pm_response_name(model))))
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
  if (any(is.na(stats::coef(model)))) stop("there are aliased coefficients in the model")
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
## Cholesky factor R of X'X (Doolittle's t(L) is R with rows scaled, which
## does not change a Wald test). Returns a list of L matrices named by term.
pm_type1_L <- function(X, labels, asgn = attr(X, "assign")) {
  R <- chol(crossprod(as.matrix(X)))
  out <- lapply(seq_along(labels), function(j) R[asgn == j, , drop = FALSE])
  names(out) <- labels
  out[vapply(out, nrow, 1L) > 0]
}

## Wald F of each hypothesis L b = 0 (q = rows of L), with covariance V.
## Returns a one-column matrix "F value" with a row per hypothesis.
pm_type1_F <- function(Ls, b, V) {
  V <- as.matrix(V)
  Fv <- vapply(Ls, function(L) {
    d <- L %*% b
    as.numeric(crossprod(d, solve(L %*% V %*% t(L), d))) / nrow(L)
  }, numeric(1))
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
  if (qs$rank < ncol(Xs) || qr(X)$rank < ncol(Xs)) stop("there are aliased coefficients in the model")
  B <- qr.coef(qs, X)
  if (max(abs(Xs %*% B - X)) > 1e-8 * max(1, abs(X))) stop("cannot map the model's coefficients to sum-to-zero contrasts")
  dimnames(B) <- list(colnames(Xs), coef_names)
  list(B = B, assign = attr(Xs, "assign"), Xs = Xs)
}

## The Type III hypothesis matrix of term j (0 = intercept), from pm_type3_map().
pm_type3_hyp <- function(map, j) {
  H <- map$B[map$assign == j, , drop = FALSE]
  H[apply(abs(H), 1, max) > 1e-12, , drop = FALSE]
}
