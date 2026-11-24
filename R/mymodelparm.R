## Coefficients, their covariance, the model's df and which coefficients are
## estimable: the one reader of a fit's fixed effects for every front end
## (through pm_inference(), R/inference.R). Coefficients and covariance come
## only from the adapter generics pm_coef() and pm_vcov() (R/adapter.R), so a
## new model class needs no code here. 'estimable' has one entry per
## coefficient of pm_coef() (FALSE for aliased ones, NA there); 'coef' and
## 'vcov' hold the estimable ones. df: the model's residual df, or for lmer
## one df per term of the type I table (pm_model_df()).
pm_model_parm <- function(model, df = NULL, ...) {
  if (inherits(model, "aovlist")) stop("This function does not support objects of class ", sQuote("aovlist"))
  ocoef <- tryCatch(pm_coef(model), error = function(e) NULL)
  if (!is.numeric(ocoef)) stop("no ", sQuote("coef"), " method for ", sQuote("model"), " found!")
  sigma <- tryCatch(as.matrix(pm_vcov(model)), error = function(e) NULL)
  if (is.null(sigma)) stop("no ", sQuote("vcov"), " method for ", sQuote("model"), " found!")
  # gam: every basis coefficient is estimated (penalised); the Bayesian covariance Vp
  estimable <- if (inherits(model, "gam")) rep(TRUE, length(ocoef)) else !is.na(ocoef)
  beta <- ocoef[estimable]
  if (any(!estimable) && nrow(sigma) == length(ocoef)) sigma <- sigma[estimable, estimable, drop = FALSE]
  if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable)) {
    stop("could not extract coefficients and covariance matrix from ", sQuote("model"))
  }
  if (is.null(df)) {
    df <- pm_model_df(model, ...)
  } else if (any(df < 0)) {
    stop(sQuote("df"), " is not positive")
  }
  structure(list(coef = beta, vcov = sigma, df = df, estimable = estimable), class = "mymodelparm")
}

## The model's df for pm_model_parm(): residual df for lm/glm/aov/gls, the
## residual df of the likelihood for glmer/glmmTMB, mgcv's residual df for gam,
## and for lmer one df per term of the type I table (df_type: Kenward-Roger
## for REML fits, Satterthwaite for ML, or containment / between-within).
## Other classes: df.residual(), else Inf (normal-based; see ?pm_coef).
pm_model_df <- function(model, df_type = NULL, vcov_adj = NULL, ...) {
  if (inherits(model, "gam")) return(model$df.residual)
  if (inherits(model, "lm")) {          # aov, lm, glm
    class(model) <- "lm"
    return(summary(model)$df[2])
  }
  if (inherits(model, "gls")) return(model$dims[["N"]] - model$dims[["p"]])
  if (inherits(model, "lmerMod")) {
    if (is.null(df_type)) df_type <- if (lme4::getME(model, "is_REML")) "Kenward-Roger" else "Satterthwaite"
    L_list <- get_contrasts_type1(model)
    if (identical(df_type, "containment")) return(pm_containment_df(model)[names(L_list)])
    if (identical(df_type, "between-within")) return(pm_betwithin_df(model)[names(L_list)])
    return(sapply(L_list, function(L) mean(pm_df_term(model, ctrmatrix = L, type = df_type, vcov_adj = vcov_adj))))
  }
  if (inherits(model, c("glmerMod", "glmmTMB"))) return(summary(model)$AICtab["df.resid"])
  if (inherits(model, "lme")) return(0)
  pm_residual_df(model)
}

## The columns of a K (or L) matrix that go with the model's estimable
## coefficients, in their order: the one place where the columns of a design
## are matched to mp$coef. By name when the coefficients are named (any
## class: a model may list aliased coefficients as NA, as lm() does, or leave
## them out of its design, as lme4 does, or both, as glmmTMB does); by
## position only when there are no names to go by and the sizes agree.
## Anything else stops: a silent mismatch would give wrong means.
pm_K_cols <- function(K, mp) {
  b <- names(mp$coef)
  cn <- colnames(K)
  if (length(b) && length(cn)) {
    if (all(b %in% cn)) return(K[, b, drop = FALSE])
  } else if (ncol(K) == length(mp$estimable)) {
    return(K[, mp$estimable, drop = FALSE])
  }
  miss <- if (length(b) && length(cn)) setdiff(b, cn) else character(0)
  stop("predictmeans cannot match the model's coefficients (", length(mp$coef), " estimable of ",
       length(mp$estimable), ") to the columns of its design (", ncol(K), ")",
       if (length(miss)) paste0("; not in the design: ", paste(utils::head(miss, 5), collapse = ", "),
                                if (length(miss) > 5) ", ..."),
       ". For a new model class, pm_coef() and pm_basis() must give the same column names (see ?pm_coef).", call. = FALSE)
}
