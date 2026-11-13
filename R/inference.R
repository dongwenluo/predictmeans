## The inference choices of a call, made in one place: which covariance of the
## coefficients is used, which denominator df, and how they are reported.
## predictmeans(), contrastmeans(), pm_joint(), pm_slopes() and
## covariatemeans() all use pm_inference(), so a default or a rule is changed
## here only.
##
## The rules:
##  * df method (pm_df_method()): "auto" is Kenward-Roger for REML lmer, lme
##    and gls fits, Satterthwaite for ML fits (from the ML log-likelihood); "containment" for lme and gls is nlme's df (one df per term for
##    lme, the residual df for gls), as predictmeans 1.x; GLMMs (glmer,
##    glmmTMB) inner-outer df ("residual": the residual df, as 1.x); glm
##    with a known dispersion (Poisson, binomial, negative binomial) z tests
##    ("residual": t with the residual df, as 1.x); other models use the
##    residual df (or a pm_df() method, else normal-based).
##  * Kenward-Roger uses the observed information of the variance parameters
##    (as SAS DDFM=KR): in the simulations of dev/sim-df/ it was never worse
##    than the expected information (pbkrtest, emmeans; "KR-expected") and
##    held 5% where that gave 8.5% with missing values (Kenward and
##    Molenberghs 1998: with missing data only the observed information is
##    valid). "KR-average" uses the average information, as Genstat.
##  * covariance: the Kenward-Roger adjusted covariance when the df are
##    Kenward-Roger and vcov_adjust is TRUE; not with "containment" or
##    "between-within" unless vcov_adjust = TRUE is given (as SAS); not with a
##    permlist (the refits use the unadjusted covariance); the user's 'vcov'
##    (R/vcov.R) replaces the model's covariance.
##  * aov() models with Error(): analysed by strata (R/aovlist.R): stratum
##    mean squares and Satterthwaite df across strata ("multistratum").
##  * df of a row of L: per row with Kenward-Roger and Satterthwaite (and a
##    pm_df() method); with inner-outer df (GLMMs) the smallest df of the
##    coefficients the row uses (R/df-innerouter.R); with one df per term (containment, between-within,
##    nlme's lme rule) the df of the term the rows belong to ('terms'; the
##    smallest if several, e.g. the terms with the covariate for slopes); else
##    the residual df.

## Arguments: df_method as in predictmeans(); vcov_adjust and
## vcov_adjust_default (TRUE when the user left vcov_adjust at its default:
## then it is switched off quietly where it does not apply); vcov, the user's
## covariance; permlist; term_df: also compute the per-term df of lmer models
## (mp$df, for the legacy tables of predictmeans() and covariatemeans()).
##
## Value, a list of class "pm_inference":
##   mp           coefficients, the covariance used, the model's df, estimable flags (mymodelparm())
##   df_type      "Kenward-Roger", "Satterthwaite", "containment", "between-within", or NULL (the model's own rule)
##   vcov_adj     Kenward-Roger adjusted covariance, or the Satterthwaite ingredients ("pm_satt")
##   vcov_adjust  TRUE when the adjusted covariance is used
##   vcov_label   label of the user's covariance, or NULL
##   per_row      TRUE when the df vary by row (Kenward-Roger, Satterthwaite, pm_df())
##   label        how the df are reported (settings$df_method, attr(pm_joint(), "df_method"))
##   df(L, terms)     df of each row of L (coefficient space)
##   ftest(L, terms)  F test of L beta = 0: c(df1, df2, F, p)
pm_inference <- function(model, df_method = "auto", vcov_adjust = TRUE, vcov_adjust_default = TRUE,
                         vcov = NULL, permlist = NULL, term_df = FALSE) {
  ## covariance: which one
  if (df_method %in% c("containment", "between-within") && vcov_adjust_default) vcov_adjust <- FALSE
  vcov_adjust <- pm_perm_vcov_adjust(model, permlist, vcov_adjust, vcov_adjust_default)
  vcov <- pm_vcov_with_perm(vcov, permlist)
  if (!is.null(vcov)) vcov_adjust <- FALSE          # the user's covariance replaces the adjusted one

  ## df method and its ingredients (computed once, shared by all df)
  lmer <- inherits(model, "lmerMod")
  nlme <- inherits(model, c("lme", "gls"))
  aovl <- inherits(model, "pm_aovlist")               # aov() with Error(): by strata (R/aovlist.R)
  if (aovl) {
    if (df_method != "auto") warning("'df_method' is not used for aov models with Error(): they are analysed by ",
                                     "strata.", call. = FALSE)
    if (!is.null(vcov)) warning("'vcov' is not used for aov models with Error().", call. = FALSE)
    vcov <- NULL
    sa <- pm_aov_strata(model)
    dfm <- list(df_type = "Satterthwaite", vcov_adj = sa$parts)
  } else dfm <- pm_df_method(model, df_method)
  df_type <- dfm$df_type
  vcov_adj <- dfm$vcov_adj
  if (any(df_type %in% c("containment", "between-within")) && vcov_adjust && lmer && lme4::getME(model, "is_REML"))
    vcov_adj <- tryCatch(pm_vcov_kr(model), error = function(e) NULL)

  ## coefficients and covariance
  mp <- if (aovl) sa$mp else if (lmer && !term_df) mymodelparm(model, df = 0) else
    mymodelparm(model, df_type = df_type, vcov_adj = vcov_adj)
  adjusted <- FALSE
  if (vcov_adjust) {
    if (is.null(vcov_adj) || inherits(vcov_adj, "pm_satt")) {
      # the default applies to REML fits with Kenward-Roger df only; say so only when asked for
      if (!vcov_adjust_default) warning("'vcov_adjust = TRUE' needs an lmer model fitted by REML; ignored.")
    } else {
      mp$vcov <- as.matrix(vcov_adj)
      adjusted <- TRUE
    }
  }
  V <- pm_vcov_resolve(model, vcov, names(mp$coef))
  if (!is.null(V)) mp$vcov <- V

  ## df
  one_df <- if (identical(df_type, "containment")) {
    pm_containment_df(model)
  } else if (identical(df_type, "between-within")) {
    pm_betwithin_df(model)
  } else if (is.null(df_type) && inherits(model, "lme")) {
    stats::terms(model$fixDF)                          # nlme's df of each term
  }
  io <- dfm$io                                         # inner-outer df of each coefficient
  per_row_kr_satt <- any(df_type %in% c("Kenward-Roger", "Satterthwaite")) && (lmer || nlme || aovl)
  resid_df <- if (identical(df_type, "asymptotic")) Inf else if (is.null(one_df) && !per_row_kr_satt && is.null(io)) {
    d <- mp$df
    if (is.numeric(d) && length(d) == 1 && is.finite(d) && d > 0) d else Inf    # normal-based when unknown
  }
  custom <- is.null(one_df) && !per_row_kr_satt && is.null(io) && !identical(df_type, "asymptotic") &&
    !is.null(tryCatch(pm_df(model, matrix(0, 1, length(mp$coef))), error = function(e) NULL))

  df <- function(L, terms) {
    L <- rbind(L)
    n <- nrow(L)
    if (!is.null(one_df)) {
      d <- one_df[intersect(terms, names(one_df))]
      if (!length(d)) stop("'", paste(terms, collapse = "', '"), "' is not a term of the model")
      return(rep(min(unname(d)), n))
    }
    if (!is.null(io)) {
      # contrasts: the smallest df of the coefficients a row uses; means (rows
      # with the intercept): the df of their term (R/df-innerouter.R)
      d <- io[names(mp$coef)]
      int <- attr(io, "intercept")[match(names(mp$coef), names(io))]
      tdf <- attr(io, "terms")
      tdf <- tdf[intersect(terms, names(tdf))]
      return(apply(L, 1, function(l) {
        nz <- abs(l) > 1e-10
        if (!any(nz)) return(NA_real_)
        if (any(nz & int) && length(tdf) && all(is.finite(tdf))) return(min(tdf))
        min(d[nz])
      }))
    }
    if (per_row_kr_satt) {
      if (nlme || aovl) return(as.numeric(pm_nlme_df(L, model, df_type, vcov_adj)))
      return(as.numeric(df_term(model, ctrmatrix = L, type = df_type, vcov_adj = vcov_adj)))
    }
    if (custom) return(rep_len(as.numeric(pm_df(model, L)), n))
    rep(resid_df, n)
  }

  ftest <- function(L, terms) {
    r <- qr(t(L))
    q <- r$rank
    if (q == 0) return(c(0, NA, NA, NA))
    L <- L[r$pivot[seq_len(q)], , drop = FALSE]
    beta <- mp$coef
    b <- L %*% beta
    wald <- function(Vb) as.numeric(t(b) %*% solve(L %*% Vb %*% t(L), b))/q
    if (identical(df_type, "Kenward-Roger") && !is.null(attr(vcov_adj, "P"))) {   # lmer, lme, gls
      kr <- pm_kr_ftest(model, L, PhiA = vcov_adj)
      ddf <- unname(kr["ddf"])
      Fv <- unname(kr["stat"])
    } else if (identical(df_type, "Satterthwaite") && inherits(vcov_adj, "pm_satt")) {
      # lmerTest's rule for several contrasts (contestMD): split L V L' into
      # its eigen-directions, each with a Satterthwaite df nu_m, and combine
      # E = sum(nu_m / (nu_m - 2)), df2 = 2E / (E - q)
      eg <- eigen(L %*% vcov_adj$vcov_beta %*% t(L), symmetric = TRUE)
      PtL <- crossprod(eg$vectors, L)
      Fv <- sum(drop(PtL %*% beta)^2/eg$values)/q
      nu <- pm_satt_ddf(PtL, vcov_adj)
      E <- sum((nu/(nu - 2))[nu > 2])
      ddf <- if (E > q) 2*E/(E - q) else NA_real_
    } else {
      ddf <- min(df(L, terms))
      Fv <- wald(mp$vcov)
    }
    if (!is.null(V)) Fv <- wald(mp$vcov)    # the given covariance, with the df of the model's own test
    c(q, ddf, Fv, stats::pf(Fv, q, ddf, lower.tail = FALSE))
  }

  label <- if (aovl) {
    "multistratum"                                  # stratum mean squares, Satterthwaite across strata
  } else if (identical(df_type, "asymptotic")) {
    "asymptotic (z)"                                # glm with a known dispersion
  } else if (df_method %in% c("KR-expected", "KR-average") && identical(df_type, "Kenward-Roger")) {
    df_method                                       # Kenward-Roger with the expected or average information
  } else if (!is.null(df_type)) {
    df_type
  } else if (inherits(model, "lme")) {
    "containment (nlme)"
  } else if (custom) {
    "model (pm_df)"
  } else {
    "residual"
  }
  structure(list(mp = mp, df_type = df_type, vcov_adj = vcov_adj, vcov_adjust = adjusted,
                 vcov_label = attr(V, "label"), per_row = per_row_kr_satt || custom || !is.null(io), label = label,
                 df = df, ftest = ftest), class = "pm_inference")
}

## Resolve the df method: "auto" is Kenward-Roger for REML fits (lmer, lme,
## gls), Satterthwaite for ML fits;
## lme and gls also take "Satterthwaite", "KR-expected", "between-within" and
## "containment" (nlme's df); other models keep their own df. Kenward-Roger uses the observed
## information, "KR-expected" the expected one, "KR-average" the average one
## (Genstat). Returns list(df_type, vcov_adj).
pm_df_method <- function(model, df_method) {
  kr_info <- switch(df_method, "KR-expected" = "expected", "KR-average" = "average", "observed")
  glmm <- inherits(model, c("glmerMod", "glmmTMB"))
  if (identical(df_method, "inner-outer") || (glmm && df_method %in% c("auto", "Kenward-Roger", "KR-expected",
                                                                         "KR-average", "Satterthwaite"))) {
    # df by grouping level (R/df-innerouter.R): the default for GLMMs, which have
    # no Kenward-Roger or Satterthwaite df here
    if (glmm && !df_method %in% c("auto", "inner-outer")) {
      warning("'df_method = \"", df_method, "\"' is not available for GLMMs; using \"inner-outer\".", call. = FALSE)
    }
    if (inherits(model, c("merMod", "glmmTMB", "lme"))) {
      io <- tryCatch(pm_inner_outer_df(model), error = function(e) e)
      if (!inherits(io, "error")) return(list(df_type = "inner-outer", vcov_adj = NULL, io = io))
      if (df_method == "inner-outer") warning("inner-outer df not available (", conditionMessage(io), "); ignored.",
                                              call. = FALSE)
      else message("inner-outer df not available (", conditionMessage(io), "); using the residual df.")
      if (glmm) return(list(df_type = NULL, vcov_adj = NULL))
    } else {
      warning("'df_method = \"inner-outer\"' needs random effects (lmer, glmer, glmmTMB, lme); ignored.", call. = FALSE)
      return(list(df_type = NULL, vcov_adj = NULL))
    }
  }
  if (identical(df_method, "residual")) {
    if (inherits(model, c("lmerMod", "lme"))) {
      warning("'df_method = \"residual\"' is for models without random effects and GLMMs; use \"containment\". ",
              "Ignored.", call. = FALSE)
      df_method <- "auto"
    } else return(list(df_type = NULL, vcov_adj = NULL))       # GLMMs as predictmeans 1.x; gls, lm, glm
  }
  if (inherits(model, "glm") && pm_fixed_dispersion(model)) {
    # Poisson, binomial, negative binomial: the dispersion is known, so the
    # Wald tests are z (chi-square) tests, as emmeans and SAS GLIMMIX. In
    # dev/sim-df-glm/ t with the residual df rejected 0.3-4.5% of true nulls
    # at 5%, z 3.6-5.7%. "residual" (above) keeps the t tests of 1.x.
    if (df_method != "auto") {
      warning("'df_method = \"", df_method, "\"' does not apply to glm models; ignored.", call. = FALSE)
    }
    return(list(df_type = "asymptotic", vcov_adj = NULL))
  }
  if (glmm && df_method %in% c("containment", "between-within")) {
    if (inherits(model, "glmerMod")) return(list(df_type = df_method, vcov_adj = NULL))
    warning("'df_method = \"", df_method, "\"' is not available for glmmTMB; using the residual df.", call. = FALSE)
    return(list(df_type = NULL, vcov_adj = NULL))
  }
  if (identical(df_method, "between-within")) {
    if (inherits(model, c("lmerMod", "lme", "gls"))) return(list(df_type = "between-within", vcov_adj = NULL))
    warning("'df_method = \"between-within\"' applies to lmer, lme and gls models; ignored.")
    return(list(df_type = NULL, vcov_adj = NULL))
  }
  if (inherits(model, c("lme", "gls"))) {
    # nlme models (R/satt-nlme.R). "auto": Kenward-Roger for REML fits (in the
    # simulations of dev/sim-df/ nlme's and the residual df reached 16-27% for
    # F tests of unstructured repeated measures); "containment": nlme's df
    # (lme) or the residual df (gls), as predictmeans 1.x
    if (df_method == "containment") return(list(df_type = NULL, vcov_adj = NULL))
    ml <- !identical(model$method, "REML")
    if (ml && df_method %in% c("Kenward-Roger", "KR-expected", "KR-average")) {
      warning("Kenward-Roger df need a REML fit: using Satterthwaite instead")
      df_method <- "Satterthwaite"
    }
    if (df_method == "auto" && ml) {
      # as for ML lmer fits: Satterthwaite from the ML log-likelihood
      message("The model was fitted by ML: using Satterthwaite df.")
      sp <- tryCatch(pm_satt_parts_nlme(model), error = function(e) e)
      if (!inherits(sp, "error")) return(list(df_type = "Satterthwaite", vcov_adj = sp))
      message("Satterthwaite df are not available for this model (", conditionMessage(sp), "); using ",
              if (inherits(model, "lme")) "nlme's df." else "the residual df.")
      return(list(df_type = NULL, vcov_adj = NULL))
    }
    if (df_method == "auto") {
      kr <- tryCatch(pm_vcov_kr_nlme(model, "observed"), error = function(e) e)
      if (!inherits(kr, "error")) return(list(df_type = "Kenward-Roger", vcov_adj = kr))
      message("Kenward-Roger df are not available for this model (", conditionMessage(kr), "); using ",
              if (inherits(model, "lme")) "nlme's df." else "the residual df.")
      return(list(df_type = NULL, vcov_adj = NULL))
    }
    if (df_method == "Satterthwaite") return(list(df_type = "Satterthwaite", vcov_adj = pm_satt_parts_nlme(model)))
    if (df_method %in% c("Kenward-Roger", "KR-expected", "KR-average")) {
      return(list(df_type = "Kenward-Roger", vcov_adj = pm_vcov_kr_nlme(model, kr_info)))
    }
  }
  if (!inherits(model, "lmerMod")) {
    if (df_method != "auto") {
      warning("'df_method' applies to lmer models only (\"between-within\" and \"Satterthwaite\" also to lme and ",
              "gls); ignored.")
    }
    return(list(df_type = NULL, vcov_adj = NULL))
  }
  reml <- lme4::getME(model, "is_REML")
  df_type <- switch(df_method, auto = if (reml) "Kenward-Roger" else "Satterthwaite", "KR-expected" = , "KR-average" = "Kenward-Roger",
                    df_method)
  if (df_type == "Kenward-Roger" && !reml) {
    warning("Kenward-Roger df need a REML fit: using Satterthwaite instead")
    df_type <- "Satterthwaite"
  }
  if (df_method == "auto" && !reml) message("The model was fitted by ML: using Satterthwaite df.")
  vcov_adj <- switch(df_type,
    "Kenward-Roger" = tryCatch(pm_vcov_kr(model, kr_info), error = function(e) NULL),
    "Satterthwaite" = pm_satt_parts(model),
    NULL)
  list(df_type = df_type, vcov_adj = vcov_adj)
}

## The model terms that contain variable 'v' (e.g. the covariate of slopes).
pm_terms_with <- function(model, v) {
  tl <- tryCatch(attr(stats::terms(model), "term.labels"), error = function(e) NULL)
  if (is.null(tl)) tl <- tryCatch(attr(stats::terms(stats::formula(model)), "term.labels"), error = function(e) NULL)
  tl[vapply(strsplit(tl, ":"), function(x) v %in% x, logical(1))]
}

## TRUE for a glm whose dispersion is known (fixed at 1), not estimated:
## Poisson, binomial and negative binomial (MASS::glm.nb) families.
pm_fixed_dispersion <- function(model) {
  fam <- tryCatch(stats::family(model)$family, error = function(e) "")
  fam %in% c("poisson", "binomial") || grepl("^Negative Binomial", fam)
}
