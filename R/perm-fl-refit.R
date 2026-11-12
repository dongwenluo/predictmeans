## Freedman-Lane permutations for linear mixed models (lmer, lme, gls) with
## the variance parameters re-estimated for each permuted data set
## (permmodels(method = "freedman-lane"), refit = TRUE, the default).
##
## With V fixed at the fit (refit = FALSE, R/perm-fl.R) the permuted
## statistics behave as if V were known, while the observed one has V
## estimated; for comparisons whose error has few degrees of freedom (whole
## plots of a split plot with few blocks) that rejected too often (study 1b:
## 6.2% for whole-plot pairs with 6 blocks). Here each test's reduced model is
## fitted as a mixed model, its residuals are decorrelated with its own V,
## permuted, correlated again and added to its fitted values, and the full
## model is refitted to each permuted data set (Freedman-Lane with Lee-Braun
## whitening; plan section 3.3: 4.6-5.3%).
##
## One set of refits per term: the term's F statistic, and the t statistics
## of the term's coefficients (each under the term's null hypothesis). Pairs
## and contrasts in predictmeans()/contrastmeans(): one set of refits for the
## family, from the model with all its comparisons zero (pm_perm_fl_refit_t()).

## What the refits need, for lmer and nlme models. Rows are used in the
## model's own order ('ord' reorders them by group for nlme, as the
## whitening of nlme models needs; identity for lmer).
pm_lmm_parts <- function(model) {
  ctl <- lme4::lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  if (inherits(model, "lmerMod")) {
    resp <- model@resp
    X <- as.matrix(lme4::getME(model, "X"))
    y <- resp$y - resp$offset
    wts <- resp$weights
    Wof <- function(m) pm_whiten(lme4::getME(m, "Z") %*% lme4::getME(m, "Lambda"), 1 / wts)
    bars <- vapply(pm_findbars(stats::formula(model)),
                   function(b) paste0("(", paste(deparse(b), collapse = ""), ")"), "")
    d <- pm_model_rows(model)
    fit_reduced <- function(XN) {
      dd <- d
      dd$.pm_y <- y
      dd$.pm_X <- XN
      .pm_wts <- wts
      f <- stats::as.formula(paste(".pm_y ~ 0 + .pm_X +", paste(bars, collapse = " + ")))
      environment(f) <- environment()
      m0 <- quiet(lme4::lmer(f, data = dd, weights = .pm_wts, REML = lme4::isREML(model), control = ctl))
      list(f0 = as.vector(XN %*% lme4::fixef(m0)), W = Wof(m0))
    }
    return(list(X = X, y = y, asgn = attr(lme4::getME(model, "X"), "assign"), Terms = stats::terms(model),
                ord = seq_along(y), Wfit = Wof(model), Wof = Wof, fit_reduced = fit_reduced,
                refit_full = function(ys) quiet(lme4::refit(model, ys + resp$offset)),
                coefvcov = function(m) list(b = lme4::fixef(m), V = as.matrix(stats::vcov(m))),
                packages = "lme4"))
  }
  # lme and gls
  d <- as.data.frame(pm_get_data(model))
  Terms <- stats::terms(model)
  X0 <- stats::model.matrix(Terms, d)
  b <- pm_fixef(model)
  X <- X0[, names(b), drop = FALSE]
  y <- as.numeric(nlme::getResponse(model))
  yname <- all.vars(stats::formula(model)[[2]])
  if (length(yname) != 1 || !identical(deparse(stats::formula(model)[[2]]), yname) || nrow(d) != length(y)) {
    stop("Freedman-Lane refits of lme/gls models need a response that is one variable and data without missing ",
         "values; use refit = FALSE.", call. = FALSE)
  }
  V_list <- build_Sigma_mats(model, sigma_scale = FALSE)
  ord <- get_order(A = V_list, B = cbind(X, y))[[2]]
  Wof <- function(m) pm_whiten_blocks(build_Sigma_mats(m, sigma_scale = FALSE))
  fa <- if (inherits(model, "lme")) "fixed" else "model"
  fit_reduced <- function(XN) {
    dd <- d
    dd$.pm_y <- y
    dd$.pm_X <- XN
    # (gls names its formula 'model', so the call is changed here, not by pm_refit())
    cl <- stats::getCall(model)
    cl[[fa]] <- .pm_y ~ 0 + .pm_X
    cl$data <- dd
    m0 <- quiet(eval(cl, pm_fit_env(model)))
    list(f0 = as.vector(XN %*% pm_fixef(m0)), W = Wof(m0))
  }
  refit_full <- function(ys) {
    dd <- d
    dd[[yname]] <- ys
    quiet(pm_refit(model, data = dd))
  }
  list(X = X, y = y, asgn = attr(X0, "assign")[match(names(b), colnames(X0))], Terms = Terms, ord = ord,
       Wfit = pm_whiten_blocks(V_list), Wof = Wof, fit_reduced = fit_reduced, refit_full = refit_full,
       coefvcov = function(m) list(b = pm_fixef(m), V = as.matrix(stats::vcov(m))), packages = "nlme")
}

## F statistic of a test (columns of the nested models A within B) by least
## squares on decorrelated data.
pm_fl_F <- function(Xw, zw, A, B) {
  fit <- function(cols) {
    if (!length(cols)) return(list(rss = sum(zw^2), rank = 0L))
    q <- qr(Xw[, cols, drop = FALSE])
    list(rss = sum(qr.resid(q, zw)^2), rank = q$rank)
  }
  full <- fit(seq_len(ncol(Xw)))
  rA <- fit(A)
  rB <- if (length(B) == ncol(Xw)) full else fit(B)
  ((rA$rss - rB$rss) / (rB$rank - rA$rank)) / (full$rss / (nrow(Xw) - full$rank))
}

pm_perm_fl_refit <- function(model, nperm, type, exact, seed, prt, perms, ncore, only = NULL) {
  # GLMMs: the same with the working linear mixed model (pm_working_lmm())
  glmm <- inherits(model, c("glmerMod", "glmmTMB"))
  if (glmm) model <- pm_working_lmm(model)
  P <- pm_lmm_parts(model)
  ord <- P$ord
  X <- P$X[ord, , drop = FALSE]                  # internal (group) order
  y <- P$y[ord]
  n <- nrow(X)
  asgn <- P$asgn
  if (!is.null(perms)) {
    perms <- pm_perm_check_index(perms, n)
    perms <- matrix(order(ord)[perms[ord, , drop = FALSE]], nrow = n)
  }
  src <- pm_perm_source(n, seed, perms)
  idx <- src(nperm)
  labs <- attr(P$Terms, "term.labels")
  fac <- attr(P$Terms, "factors")
  all_cols <- seq_len(ncol(X))
  cols_of <- function(terms) which(asgn %in% terms)
  Xw <- P$Wfit$solve(X)
  zw <- P$Wfit$solve(y)
  cv0 <- P$coefvcov(model)
  t_obs <- abs(cv0$b) / sqrt(diag(cv0$V))
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
    tests[[labs[j]]] <- list(red = red, A = AB[[1]], B = AB[[2]], cols = cols_of(j))
  }
  back <- order(ord)                              # internal order -> model order
  pval <- numeric(0)
  failures <- 0
  # ('only': test these terms only; used by the simulation scripts)
  if (!is.null(only)) tests <- tests[intersect(names(tests), only)]
  for (k in names(tests)) {
    tk <- tests[[k]]
    R0 <- P$fit_reduced(P$X[, tk$red, drop = FALSE])
    f0 <- R0$f0[ord]
    e0 <- R0$W$solve(y - f0)
    tk$value <- pm_fl_F(Xw, zw, tk$A, tk$B)
    tk$df <- qr(Xw[, tk$B, drop = FALSE])$rank - (if (length(tk$A)) qr(Xw[, tk$A, drop = FALSE])$rank else 0)
    one <- function(i) {
      ys <- f0 + R0$W$mult(e0[idx[, i]])
      m <- P$refit_full(ys[back])
      Ws <- P$Wof(m)
      cv <- P$coefvcov(m)
      c(pm_fl_F(Ws$solve(X), Ws$solve(ys), tk$A, tk$B), abs(cv$b[tk$cols]) / sqrt(diag(cv$V)[tk$cols]))
    }
    sims <- pm_perm_drop(pm_perm_map(seq_len(nperm), one, ncore, P$packages))
    failures <- max(failures, sims$failures)
    S <- matrix(unlist(sims$res), ncol = length(sims$res))
    pval[k] <- pm_perm_p(tk$value, S[1, ], exact)
    pval[paste0("coef:", colnames(X)[tk$cols])] <- pm_perm_p(t_obs[tk$cols], S[-1, , drop = FALSE], exact)
    tests[[k]] <- tk
  }
  tabs <- pm_fl_tables(Xw, zw, n - qr(Xw)$rank, tests, pval, type, if (glmm) "refit-working" else "refit")
  COEF <- tabs$COEF
  ANOVA <- tabs$ANOVA
  failnote <- pm_perm_failnote(failures)
  if (prt) {
    cat("\nCoefficients of (fixed) effects:\n")
    print(COEF)
    cat("\nNote: Perm_p_value of t test is obtained using ", sQuote(nperm), " Freedman-Lane permutations",
        failnote, ".\n", sep = "")
    cat("\nANOVA:\n")
    cat(attr(ANOVA, "heading"), "\n", sep = "")
    print(ANOVA)
    cat("\nNote: Perm_p_value of F test is obtained using ", sQuote(nperm), " Freedman-Lane permutations",
        failnote, ".\n\n", sep = "")
  }
  perm <- structure(list(scheme = "freedman-lane-refit", nperm = as.integer(nperm), exact = exact,
                         start = attr(src, "start"), idx = attr(src, "matrix"), coef_names = colnames(X),
                         wmodel = if (glmm) model),
                    class = "pm_perm")
  out <- list(permlist = NULL, COEFFICENTS = COEF, ANOVA = ANOVA, perm = perm)
  if (failures > 0) attr(out, "failures") <- failures
  out
}

## Freedman-Lane p-values, with refits, of the contrasts in the rows of L:
## the reduced model has all of them zero (one set of refits for the
## family); each row's |t| is compared with its |t| in the refits.
pm_perm_fl_refit_t <- function(perm, L, model, ncore = 1) {
  if (!identical(pm_perm_coef_names(model), perm$coef_names)) {
    stop("'permlist' was made for a different model.", call. = FALSE)
  }
  if (!is.null(perm$wmodel)) model <- perm$wmodel      # a GLMM: its working linear mixed model
  P <- pm_lmm_parts(model)
  L <- matrix(as.numeric(as.matrix(L)), ncol = ncol(as.matrix(L)))
  if (ncol(L) != ncol(P$X)) stop("'permlist' was made for a different model.", call. = FALSE)
  ord <- P$ord
  n <- length(ord)
  q <- qr(t(L))
  N <- qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
  R0 <- P$fit_reduced(P$X %*% N)
  y <- P$y[ord]
  X <- P$X[ord, , drop = FALSE]
  f0 <- R0$f0[ord]
  e0 <- R0$W$solve(y - f0)
  idx <- pm_perm_source(n, perms = perm$idx, state = perm$start)(perm$nperm)
  back <- order(ord)
  tstat <- function(cv) abs(as.vector(L %*% cv$b)) / sqrt(rowSums((L %*% cv$V) * L))
  obs <- tstat(P$coefvcov(model))
  one <- function(i) tstat(P$coefvcov(P$refit_full((f0 + R0$W$mult(e0[idx[, i]]))[back])))
  sims <- pm_perm_drop(pm_perm_map(seq_len(perm$nperm), one, ncore, P$packages))
  pm_perm_p(obs, matrix(unlist(sims$res), nrow = nrow(L)), perm$exact)
}

## The working linear mixed model of a GLMM (glmer, glmmTMB): the working
## response of pm_working() fitted by lmer with the GLMM's fixed and random
## terms and weights 1/(working variance), by REML. Freedman-Lane with refits
## of this model is the GLMM version of pm_perm_fl_refit().
pm_working_lmm <- function(model) {
  wk <- pm_working(model)
  n <- length(wk$z)
  if (any(vapply(pm_re_terms(model), function(tm) tm$nl >= n, logical(1)))) {
    stop("A random effect with one level per observation cannot be used on the working scale; use refit = FALSE.",
         call. = FALSE)
  }
  d <- pm_model_rows(model)
  d$.pm_z <- wk$z
  .pm_wts <- 1 / wk$resvar
  f <- pm_working_formula(stats::formula(model))
  environment(f) <- environment()
  ctl <- lme4::lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")
  suppressWarnings(suppressMessages(lme4::lmer(f, data = d, weights = .pm_wts, REML = TRUE, control = ctl)))
}
