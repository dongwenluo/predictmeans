## Tests of variance components, one per row of the varcomp() table
## (varcomp(test = TRUE); dev/predictmeans-2.0-plan.md). Each row's
## hypothesis is a model one would fit, refitted:
##   a term with one variance ((1 | g), a smooth block)     the term dropped              kind "drop"
##   the variance of a slope in (1 + x | g)                  (1 | g)                       "reduce"
##   the intercept variance of (1 + x | g), or a variance    the whole term dropped        "drop_term"
##     of a structured term (ar1(), cs(), ...)
##   a covariance of (1 + x | g)                             (1 + x || g), diag() in glmmTMB "zero_cov"
##   the correlation of ar1(), ou(), cs()                    homdiag() / diag()            "zero_cor"
##   the residual variance                                   no test                       "none"
## Step 1 (this file): the hypotheses, the refits and the likelihood ratio
## statistics. 'kind' and 'npar' (the number of covariance parameters set
## to zero) say which reference distribution a p-value needs (boundary
## mixtures, or the exact restricted LRT for one variance of a linear mixed
## model).

## The model behind a fit: "semireg" (semireg()/semireg_tmb(), refitted by
## pm_sr_refit()), "mer" (lmer/glmer), "glmmTMB"; NULL for other classes.
pm_vc_engine <- function(model) {
  if (inherits(model, "semireg")) return("semireg")
  if (inherits(model, "merMod")) return("mer")
  if (inherits(model, "glmmTMB")) return("glmmTMB")
  NULL
}

pm_vc_dp <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = " ")

## The random terms of a model formula as written: text, left-hand side,
## grouping text, components (column names of the term's design; NULL when
## they can't be computed) and structure ("us" for a plain bar).
pm_vc_formula_terms <- function(f, mf) {
  sp <- pm_special_terms(f)
  bars <- pm_plain_bars(f, sp)
  one <- function(e, special) {
    bar <- if (special) e[[2]] else e
    lhs <- bar[[2]]
    comps <- tryCatch(colnames(stats::model.matrix(stats::as.formula(paste("~", pm_vc_dp(lhs))), mf)),
                      error = function(err) NULL)
    list(text = if (special) pm_vc_dp(e) else paste0("(", pm_vc_dp(bar), ")"), lhs = lhs, group = pm_vc_dp(bar[[3]]),
         comps = comps, struct = if (special) as.character(e[[1]]) else "us")
  }
  c(lapply(bars, one, special = FALSE), lapply(sp, one, special = TRUE))
}

## The fitted random terms, each with its varcomp() group name, structure
## code and the formula term it comes from (matched by grouping and
## components). For a semireg fit the smooth blocks are marked.
pm_vc_terms <- function(model) {
  eng <- pm_vc_engine(model)
  fit <- if (eng == "semireg") model$semer else model
  if (eng == "semireg") {
    f <- model$refit$formula
    blocks <- names(model$refit$smoothZ)
    mf <- model$refit$data
  } else {
    f <- stats::formula(fit)
    blocks <- character(0)
    mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  }
  ft <- pm_vc_formula_terms(f, mf)
  re <- pm_re_terms(fit)
  tmb <- inherits(fit, "glmmTMB")
  vc_names <- if (tmb) names(glmmTMB::VarCorr(fit)$cond) else names(lme4::getME(fit, "cnms"))
  codes <- if (tmb) vapply(glmmTMB::VarCorr(fit)$cond, function(v) names(attr(v, "blockCode")), "") else
    rep("us", length(re))
  used <- logical(length(ft))
  terms <- vector("list", length(re))
  for (j in seq_along(re)) {
    tm <- re[[j]]
    cand <- which(!used & vapply(ft, function(x) identical(x$group, tm$group), logical(1)))
    if (length(cand) > 1) {
      exact <- cand[vapply(ft[cand], function(x) identical(x$comps, tm$comps), logical(1))]
      size <- cand[vapply(ft[cand], function(x) length(x$comps) == length(tm$comps), logical(1))]
      cand <- if (length(exact)) exact else if (length(size)) size else cand
    }
    hit <- if (length(cand)) cand[1] else NA_integer_
    if (!is.na(hit)) used[hit] <- TRUE
    npar <- if (tmb) as.numeric(fit$modelInfo$reStruc$condReStruc[[j]]$blockNumTheta) else
      length(tm$comps)*(length(tm$comps) + 1)/2
    terms[[j]] <- list(group = vc_names[j], comps = tm$comps, k = length(tm$comps), code = codes[j], npar = npar,
                       smooth = tm$group %in% blocks && length(tm$comps) == 1,
                       term = if (!is.na(hit)) ft[[hit]], index = hit)
  }
  # glmmTMB with dispformula = ~0: a structured term is the residual covariance
  resid_struct <- tmb && identical(pm_vc_dp(fit$modelInfo$allForm$dispformula), "~0")
  list(terms = terms, formula_terms = ft, formula = f, engine = eng, fit = fit, blocks = blocks,
       resid_struct = resid_struct)
}

## One hypothesis per varcomp() row: list(row, label, kind, npar, change,
## block); 'change' maps the formula's random-term texts to those of the
## reduced model, 'block' is a smooth block to drop (semireg).
pm_vc_hypotheses <- function(model) {
  info <- pm_vc_terms(model)
  tmb <- inherits(info$fit, "glmmTMB")
  out <- list()
  cur <- NA_integer_                      # the random term being worked on
  add <- function(row, label, kind, npar = NA, change = NULL, block = NULL) {
    out[[length(out) + 1L]] <<- list(row = row, label = label, kind = kind, npar = npar, change = change,
                                      block = block, term = cur)
  }
  for (cur in seq_along(info$terms)) {
    tm <- info$terms[[cur]]
    g <- tm$group
    k <- tm$k
    ft <- tm$term
    if (tm$smooth) {
      add(paste(g, tm$comps, sep = "."), paste("smooth", g, "dropped"), "drop", 1, block = g)
      next
    }
    if (is.null(ft)) {
      for (cn in tm$comps) add(paste(g, cn, sep = "."), "not available (term not found in the formula)", "none")
      next
    }
    # (made by factories: functions made in the loop itself would all see the
    # last term's index)
    drop_it <- pm_vc_drop_fun(tm$index)
    make <- function(new) pm_vc_replace_fun(tm$index, new)
    with_struct <- function(s) paste0(s, "(", pm_vc_dp(ft$lhs), " | ", ft$group, ")")
    lab_drop <- paste("drop", ft$text)
    code <- tm$code
    if (k == 1) {                           # one variance: the term dropped
      add(paste(g, tm$comps, sep = "."), lab_drop, "drop", 1, drop_it)
      next
    }
    if (code == "us") {                     # (1 + x | g): k components, unstructured
      tl <- tryCatch(stats::terms(stats::as.formula(paste("~", pm_vc_dp(ft$lhs)))), error = function(e) NULL)
      labs <- if (is.null(tl)) character(0) else attr(tl, "term.labels")
      int <- if (is.null(tl)) 1L else attr(tl, "intercept")
      for (cn in tm$comps) {
        row <- paste(g, cn, sep = ".")
        rest <- setdiff(labs, cn)
        if (cn == "(Intercept)" || (int == 0 && cn %in% labs && !length(rest))) {
          add(row, lab_drop, "drop_term", k*(k + 1)/2, drop_it)
        } else if (cn %in% labs) {
          new <- if (int == 1 && !length(rest)) paste0("(1 | ", ft$group, ")") else
            paste0("(", paste(c(if (int == 0) "0", rest), collapse = " + "), " | ", ft$group, ")")
          add(row, paste(ft$text, "->", new), "reduce", k, make(new))
        } else {
          add(row, "not available (a level of a factor in the term)", "none")
        }
      }
      nocov <- if (tmb) with_struct("diag") else paste0("(", pm_vc_dp(ft$lhs), " || ", ft$group, ")")
      for (b in seq_len(k - 1)) for (a in (b + 1):k) {
        add(paste(g, tm$comps[a], tm$comps[b], sep = "."), paste(ft$text, "->", nocov), "zero_cov", k*(k - 1)/2,
            make(nocov))
      }
      next
    }
    # a structured glmmTMB term: its variances by dropping the whole term;
    # the correlation of ar1(), ou(), cs() set to 0
    zero <- switch(code, ar1 = , ou = , homcs = with_struct("homdiag"), hetar1 = , cs = with_struct("diag"), NULL)
    if (code %in% c("ar1", "hetar1", "cs", "homcs", "toep", "diag")) {
      # rows as pm_tmb_varcomp() names them; with dispformula = ~0 the term is
      # the residual covariance, and its variances are not tested (as the
      # residual variance)
      vrows <- if (code %in% c("ar1", "homcs")) paste0(g, ".var") else paste(g, tm$comps, sep = ".")
      for (r in vrows) {
        if (info$resid_struct) add(r, "no test (the residual covariance: dispformula = ~0)", "none")
        else add(r, lab_drop, "drop_term", tm$npar, drop_it)
      }
      if (!is.null(zero)) add(paste0("cor.", g), paste(ft$text, "->", zero), "zero_cor", 1, make(zero))
      if (code == "toep") for (l in seq_len(k - 1)) add(paste0("cor.", g, ".lag", l), "not available for toep()", "none")
    } else {
      # the full matrix (ou(), exp(), ...): variances and covariances
      E <- pm_vc_entries(stats::setNames(list(matrix(1, k, k, dimnames = list(tm$comps, tm$comps))), g))
      for (r in seq_along(E)) {
        if (attr(E, "type")[r] == "var") {
          if (info$resid_struct) add(names(E)[r], "no test (the residual covariance: dispformula = ~0)", "none")
          else add(names(E)[r], lab_drop, "drop_term", tm$npar, drop_it)
        }
        else if (!is.null(zero)) add(names(E)[r], paste(ft$text, "->", zero), "zero_cor", 1, make(zero))
        else add(names(E)[r], paste("not available for", paste0(code, "()")), "none")
      }
    }
  }
  cur <- NA_integer_
  fam <- tryCatch(stats::family(info$fit)$family, error = function(e) "gaussian")
  if (identical(fam, "gaussian")) add("residual", "no test (the residual variance)", "none")
  structure(out, info = info)
}

pm_vc_drop_fun <- function(i) {
  force(i)
  function(txt) txt[-i]
}
pm_vc_replace_fun <- function(i, new) {
  force(i)
  force(new)
  function(txt) {
    txt[i] <- new
    txt
  }
}

## Refit with the random terms 'txt' (texts as in the formula; for a semireg
## fit without the smooth placeholders, and 'drop_block' a smooth block to
## leave out). Returns the fit: merMod, glmmTMB, or lm/glm when no random
## term is left (an lm wrapped to give the REML log-likelihood of a REML fit).
pm_vc_refit <- function(model, info, txt, drop_block = NULL) {
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  if (info$engine == "semireg") return(pm_sr_refit(model, drop = drop_block, random = txt))
  fit <- info$fit
  f <- info$formula
  fixed <- pm_nobars(f)
  cl <- stats::getCall(fit)
  env <- pm_fit_env(fit)
  if (!length(txt) && inherits(fit, "merMod")) {
    keep <- names(cl) %in% c("", "data", "weights", "subset", "offset", "na.action", "family", "contrasts")
    cl <- cl[keep]
    cl$formula <- fixed
    if (inherits(fit, "glmerMod")) {
      cl[[1]] <- quote(stats::glm)
      return(quiet(eval(cl, env)))
    }
    cl[[1]] <- quote(stats::lm)
    return(structure(list(fit = quiet(eval(cl, env))), class = "pm_lm_reml",
                     REML = isTRUE(lme4::getME(fit, "is_REML"))))
  }
  f_new <- stats::as.formula(paste(pm_vc_dp(f[[2]]), "~", paste(c(pm_vc_dp(fixed[[3]]), txt), collapse = " + ")),
                             env = environment(f))
  if (methods::is(fit, "semerMod")) return(quiet(stats::update(fit, formula = f_new)))
  cl$formula <- f_new
  quiet(eval(cl, env))
}

## The log-likelihood of a fit (glmmTMB: from its objective, which is there
## also when glmmTMB flags the Hessian and logLik() gives NA).
pm_vc_loglik <- function(m) {
  if (inherits(m, "glmmTMB")) return(-m$fit$objective)
  as.numeric(stats::logLik(m))
}

## The tests of the rows of varcomp(): hypotheses, refits (one per distinct
## reduced model) and likelihood ratio statistics. A data frame with one row
## per varcomp() row: row, hypothesis, kind, npar, logLik (full),
## logLik.null (reduced model), LRT.
pm_vc_test <- function(model, rows = NULL, method = c("auto", "mixture"), nsim = 10000L, seed = 2026L) {
  method <- match.arg(method)
  if (is.null(pm_vc_engine(model))) {
    stop("Tests of variance components are available for lmer, glmer, glmmTMB and semireg fits.", call. = FALSE)
  }
  hyp <- pm_vc_hypotheses(model)
  info <- attr(hyp, "info")
  if (!is.null(rows)) hyp <- hyp[vapply(hyp, function(h) h$row %in% rows, logical(1))]
  txt0 <- vapply(info$formula_terms, `[[`, "", "text")
  ph <- paste0("(1 | ", info$blocks, ")")       # smooth placeholders of a semireg fit
  ll_full <- pm_vc_loglik(info$fit)
  cache <- list()
  rows_out <- lapply(hyp, function(h) {
    ll0 <- NA_real_
    if (!identical(h$kind, "none")) {
      txt <- if (is.null(h$change)) txt0 else h$change(txt0)
      txt <- txt[!txt %in% ph]
      key <- paste(c("model:", h$block, sort(txt)), collapse = " ; ")     # (never "": an empty name can't be looked up)
      if (is.null(cache[[key]])) {
        red <- tryCatch(pm_vc_refit(model, info, txt, drop_block = h$block), error = function(e) NULL)
        cache[[key]] <<- if (is.null(red)) NA_real_ else pm_vc_loglik(red)
      }
      ll0 <- cache[[key]]
    }
    data.frame(row = h$row, hypothesis = h$label, kind = h$kind, npar = h$npar, term = h$term, logLik = ll_full,
               logLik.null = ll0, LRT = if (is.na(ll0)) NA_real_ else max(0, 2*(ll_full - ll0)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows_out)
  rownames(out) <- NULL
  pm_vc_pvalues(out, info, method = method, nsim = nsim, seed = seed)
}

## ---- step 2: p-values ----
##
## Each row's reference distribution (Self and Liang 1987; Stram and Lee
## 1994; Crainiceanu and Ruppert 2004):
##   drop       one variance: when it is the only random effect of a gaussian
##              lme4 REML fit, the exact null distribution of the restricted
##              LRT (simulated); otherwise 0.5 chi2(0) + 0.5 chi2(1). With
##              other random effects the exact distribution is not exact:
##              taken in the model with only the tested term (Greven et al.
##              2008) or with the other variances at their null estimates, it
##              rejected 7% in a crossed design and 0.5-4% for deviation curves
##              of s(x, f1*f2), where the mixture gave 5% and 0-2%; alone it
##              rejected 4.3%, the mixture 2-3% (dev/sim-vc-test.R)
##   reduce     a slope of a k-component term: 0.5 chi2(k-1) + 0.5 chi2(k)
##   drop_term  a whole multi-parameter term: chi2(npar), conservative (SAS
##              also reports a plain chi-square here)
##   zero_cov,  covariances / a correlation 0 (not on the boundary):
##   zero_cor   chi2(npar)
## Adds columns p.value and test (the distribution used).
pm_vc_pvalues <- function(tab, info, method = "auto", nsim = 10000L, seed = 2026L) {
  fit <- info$fit
  re <- if (inherits(fit, "lmerMod")) pm_re_terms(fit)
  exact_ok <- identical(method, "auto") && inherits(fit, "lmerMod") && isTRUE(lme4::getME(fit, "is_REML")) &&
    length(re) == 1
  pchi <- function(x, df) if (df == 0) as.numeric(x <= 0) else stats::pchisq(x, df, lower.tail = FALSE)
  mix <- function(x, d1, d2) 0.5*pchi(x, d1) + 0.5*pchi(x, d2)
  tab$p.value <- NA_real_
  tab$test <- NA_character_
  tab$note <- ""
  # where the reference distribution may be inaccurate (too liberal or too
  # conservative): few levels of the tested term; binary data in small
  # clusters (Laplace approximation)
  re_all <- tryCatch(pm_re_terms(fit), error = function(e) NULL)
  nobs <- tryCatch(stats::nobs(fit), error = function(e) NA_real_)
  fam <- tryCatch(stats::family(fit), error = function(e) NULL)
  binary <- !is.null(fam) && fam$family == "binomial" &&
    tryCatch(all(stats::model.response(stats::model.frame(fit)) %in% c(0, 1)), error = function(e) FALSE)
  for (i in seq_len(nrow(tab))) {
    x <- tab$LRT[i]
    k <- tab$kind[i]
    np <- tab$npar[i]
    if (is.na(x) || k == "none") next
    j <- tab$term[i]
    smooth <- !is.na(j) && isTRUE(info$terms[[j]]$smooth)
    nl <- if (!is.na(j) && !is.null(re_all)) re_all[[j]]$nl else NA
    notes <- c(if (!smooth && !is.na(nl) && nl < 6) sprintf("%d levels: asymptotic reference", nl),
               if (binary && !smooth && !is.na(nl) && nobs/nl < 5) "binary data, small clusters (Laplace)")
    if (k == "drop" && exact_ok && !is.na(tab$term[i])) {
      Zt <- lme4::getME(fit, "Z")[, re[[tab$term[i]]]$cols, drop = FALSE]
      null <- pm_rlrt_null(lme4::getME(fit, "X"), Zt, w = stats::weights(fit), nsim = nsim, seed = seed)
      tab$p.value[i] <- (sum(null >= x - 1e-8) + 1)/(length(null) + 1)     # never 0
      tab$test[i] <- paste0("exact restricted LRT (", nsim, " simulations)")
      notes <- character(0)                         # exact whatever the number of levels
    } else if (k == "drop") {
      tab$p.value[i] <- mix(x, 0, 1)
      tab$test[i] <- "0.5 chi2(0) + 0.5 chi2(1)"
    } else if (k == "reduce") {
      tab$p.value[i] <- mix(x, np - 1, np)
      tab$test[i] <- sprintf("0.5 chi2(%d) + 0.5 chi2(%d)", np - 1, np)
    } else if (k == "drop_term") {
      if (is.na(np)) next
      tab$p.value[i] <- pchi(x, np)
      tab$test[i] <- sprintf("chi2(%d), conservative", np)
    } else if (k %in% c("zero_cov", "zero_cor")) {
      tab$p.value[i] <- pchi(x, np)
      tab$test[i] <- sprintf("chi2(%d)", np)
    }
    tab$note[i] <- paste(notes, collapse = "; ")
  }
  tab
}

## The exact null distribution of the restricted likelihood ratio statistic
## for one variance component in a linear mixed model y = X b + Z u + e,
## u ~ N(0, lambda sigma^2 I) (Crainiceanu and Ruppert 2004, Theorem 1):
##   RLRT = sup_lambda (n - p) log(1 + N(lambda)/D(lambda)) - sum_s log(1 + lambda mu_s)
##   N = sum_s lambda mu_s/(1 + lambda mu_s) w_s^2,
##   D = sum_s w_s^2/(1 + lambda mu_s) + sum_{s > K} w_s^2,
## mu_s the eigenvalues of Z' P0 Z (P0 = I - X (X'X)^- X'), w_s iid N(0, 1).
## The supremum over a grid of lambda (log scale, by the range of mu), as
## RLRsim. Prior weights w: X and Z scaled by sqrt(w). nsim draws.
pm_rlrt_null <- function(X, Z, w = NULL, nsim = 10000L, seed = 2026L, ngrid = 200L) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  if (!is.null(w) && length(w) == nrow(X)) {
    X <- X*sqrt(w)
    Z <- Z*sqrt(w)
  }
  qx <- qr(X)
  n <- nrow(X)
  p <- qx$rank
  R <- qr.resid(qx, Z)
  mu <- eigen(crossprod(R), symmetric = TRUE, only.values = TRUE)$values
  mu <- mu[mu > 1e-10*max(mu)]
  K <- length(mu)
  if (!K) return(rep(0, nsim))
  lam <- c(0, exp(seq(log(1e-5/max(mu)), log(1e5/min(mu)), length.out = ngrid)))
  # keep the caller's random numbers as they were
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  W <- matrix(stats::rchisq(nsim*K, 1), nsim, K)
  rest <- if (n - p - K > 0) stats::rchisq(nsim, n - p - K) else rep(0, nsim)
  best <- rep(0, nsim)
  for (l in lam[-1]) {
    a <- l*mu/(1 + l*mu)
    Nn <- drop(W %*% a)
    Dn <- drop(W %*% (1/(1 + l*mu))) + rest
    val <- (n - p)*log1p(Nn/Dn) - sum(log1p(l*mu))
    best <- pmax(best, val)
  }
  best
}
