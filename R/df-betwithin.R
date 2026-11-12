## Between-within denominator df (SAS PROC MIXED DDFM=BETWITHIN) for lmer,
## lme and gls models.
##
## SAS's rule (Schluchter and Elashoff 1990): the residual df N - rank(X) are
## split into between-subject and within-subject parts. An effect that does
## not change within any subject gets the between-subject df, the number of
## subjects minus the rank of the between-subject columns of X (the intercept
## included). An effect that changes within a subject gets within-subject df;
## when two or more such effects contain factors, the within-subject df are
## split into subject-by-effect interactions: the df of an effect is the rank
## its subject-by-effect columns add to X, the subject and the interactions
## of the subject with the effect's lower-order parts. Only the variables
## that change within subjects form the interaction, so Group x Time with
## Group constant within subjects shares Subject x Time with Time. Ranks make
## missing cells count. One exception, as in SAS: with no random effects and
## an unstructured covariance (gls with corSymm), every effect gets the
## between-subject df.
##
## Checked against SAS (dev/sas): split plot (subject Block) Variety 10,
## Nitro 15, Variety x Nitro 30 (29 with one plot missing); repeated measures
## with an unstructured covariance: 10 for every effect.
##
## The subject is the coarsest grouping factor in which all the others are
## nested (lmer), the outermost grouping level (lme), or the grouping of the
## correlation structure (gls).

## Subject factor and the data rows used in the fit.
pm_bw_subject <- function(model) {
  if (inherits(model, "merMod")) {
    fl <- lme4::getME(model, "flist")
    nested_in <- function(g, f) all(tapply(as.integer(f), g, function(v) length(unique(v))) <= 1, na.rm = TRUE)
    ok <- vapply(seq_along(fl), function(i) {
      all(vapply(seq_along(fl), function(j) i == j || nested_in(fl[[j]], fl[[i]]), logical(1)))
    }, logical(1))
    if (!any(ok)) {
      stop("between-within df need one subject factor that holds all the random effects (nested grouping); ",
           "use df_method = \"containment\", \"Satterthwaite\" or \"Kenward-Roger\" for crossed random effects.",
           call. = FALSE)
    }
    i <- which(ok)[which.min(vapply(fl[ok], nlevels, integer(1)))]
    return(list(subject = factor(fl[[i]]), name = names(fl)[i], data = stats::model.frame(model)))
  }
  gf <- tryCatch(nlme::getGroupsFormula(model, asList = TRUE), error = function(e) NULL)
  dat <- pm_get_data(model)          # the data used in the fit (subset and missing values removed)
  if (is.null(gf) || !length(gf)) return(list(subject = NULL, name = NULL, data = dat))
  f <- gf[[1]]
  s <- interaction(lapply(all.vars(f), function(v) dat[[v]]), drop = TRUE)
  list(subject = s, name = paste(all.vars(f), collapse = ":"), data = dat)
}

## Between-within df of every fixed term (named, "(Intercept)" included),
## with attributes "residual", "between", "within" and "subject".
pm_betwithin_df <- function(model) {
  if (!inherits(model, c("merMod", "lme", "gls"))) {
    stop("between-within df are for lmer, lme and gls models.", call. = FALSE)
  }
  sb <- pm_bw_subject(model)
  tt <- if (inherits(model, "merMod")) stats::terms(model) else stats::terms(stats::formula(model))
  dat <- sb$data
  X <- stats::model.matrix(tt, dat)
  if (nrow(X) != length(stats::fitted(model))) stop("cannot match the model's data to its fit.", call. = FALSE)
  rk <- function(M) qr(M, tol = 1e-9)$rank
  asgn <- attr(X, "assign")
  r_X <- rk(X)
  N <- nrow(X)
  res_df <- N - r_X
  labs <- attr(tt, "term.labels")
  all_labs <- c(if (attr(tt, "intercept") == 1) "(Intercept)", labs)
  if (is.null(sb$subject)) {             # no subject: the residual df for every term
    out <- stats::setNames(rep(res_df, length(all_labs)), all_labs)
    attr(out, "residual") <- res_df
    return(out)
  }
  S <- sb$subject
  fac <- attr(tt, "factors")
  term_vars <- lapply(labs, function(l) rownames(fac)[fac[, l] > 0])
  # variables that change within some subject
  varies <- function(v) {
    x <- dat[[v]]
    if (is.null(x)) x <- eval(str2lang(v), dat)
    any(tapply(seq_along(S), S, function(i) length(unique(x[i])) > 1))
  }
  vv <- unique(unlist(term_vars))
  within_var <- stats::setNames(vapply(vv, varies, logical(1)), vv)
  within_term <- vapply(term_vars, function(v) any(within_var[v]), logical(1))
  # between-subject df: subjects minus the rank of the between-subject columns
  bcols <- asgn %in% c(0, which(!within_term))
  n_subj <- nlevels(droplevels(S))
  between_df <- n_subj - rk(X[, bcols, drop = FALSE])
  within_df <- res_df - between_df
  df <- stats::setNames(rep(between_df, length(labs)), labs)
  unstructured <- inherits(model, "gls") && inherits(model$modelStruct$corStruct, "corSymm")
  if (!unstructured && any(within_term)) {
    is_fac <- vapply(vv, function(v) {
      x <- dat[[v]]
      is.factor(x) || is.character(x) || is.logical(x)
    }, logical(1))
    wv <- lapply(term_vars, function(v) v[within_var[v]])            # the within part of each term
    has_fac <- vapply(seq_along(labs), function(j) within_term[j] && any(is_fac[term_vars[[j]]]), logical(1))
    df[within_term] <- within_df
    if (sum(has_fac) >= 2) {
      dat$.pm_subject <- S
      inter <- function(v) {
        f <- stats::as.formula(paste("~", paste(c(".pm_subject", sprintf("`%s`", v)), collapse = ":"), "- 1"))
        stats::model.matrix(f, dat)
      }
      sets <- unique(wv[has_fac])
      comp <- vapply(sets, function(w) {
        lower <- Filter(function(u) length(u) < length(w) && all(u %in% w), sets)
        lower <- unique(c(lower, unlist(lapply(seq_len(length(w) - 1), function(k) utils::combn(w, k, simplify = FALSE)),
                                        recursive = FALSE)))
        base <- cbind(X, stats::model.matrix(~ .pm_subject - 1, dat),
                      do.call(cbind, c(list(matrix(0, N, 0)), lapply(lower, inter))))
        rk(cbind(base, inter(w))) - rk(base)
      }, numeric(1))
      key <- vapply(sets, paste, character(1), collapse = ":")
      df[has_fac] <- comp[match(vapply(wv[has_fac], paste, character(1), collapse = ":"), key)]
    }
  }
  out <- c(if (attr(tt, "intercept") == 1) c("(Intercept)" = between_df), df)
  attr(out, "residual") <- res_df
  attr(out, "between") <- between_df
  attr(out, "within") <- within_df
  attr(out, "subject") <- sb$name
  out
}
