## What bootmodels() needs for each model class: the model fitted without its
## fixed effects (only an intercept, any offset, and the same random effects
## and covariance structure: the global null hypothesis), responses simulated
## from it, and the full model refitted to a simulated response.
##   lmer, glmer    simulate() and lme4::refit()
##   lm, aov, glm   simulate() (gaussian, binomial, Poisson, Gamma, negative
##                  binomial), refitted with the response replaced
##   glmmTMB        simulate(), refitted with the response replaced
##   lme, gls       y = X b0 + L z, L L' the null fit's covariance of the data
##                  (random effects, correlation and variance structures;
##                  pm_nlme_reml(), checked against its vcov())
## Returns list(null_fit, simulate = function(nsim) -> list of responses,
## refit = function(y) -> fitted model, packages).
pm_boot_parts <- function(model) {
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  ff <- stats::formula(model)
  lhs <- paste(deparse(ff[[2]]), collapse = "")
  # random terms as written (a glmmTMB ar1() term stays ar1(), not unstructured)
  bars <- if (inherits(model, c("merMod", "glmmTMB"))) pm_random_terms(ff) else character(0)
  # offsets in the formula stay in the null model (an offset argument stays in the call)
  tt <- stats::terms(stats::delete.response(stats::terms(pm_nobars(ff))))
  offs <- if (is.null(attr(tt, "offset"))) character(0) else
    vapply(attr(tt, "variables")[-1][attr(tt, "offset")], function(e) paste(deparse(e), collapse = ""), "")
  null_f <- stats::as.formula(paste(lhs, "~", paste(c("1", offs, bars), collapse = " + ")), env = environment(ff))
  new_lhs <- function(f) {
    f[[2]] <- quote(.pm_y)
    f
  }

  if (inherits(model, "merMod")) {
    null_fit <- quiet(pm_refit(model, formula = null_f))
    return(list(null_fit = null_fit, simulate = function(nsim) stats::simulate(null_fit, nsim = nsim),
                refit = function(y) quiet(pm_mer_refit(model, y)), packages = "lme4"))
  }
  if (inherits(model, c("lme", "gls"))) {
    fa <- if (inherits(model, "lme")) "fixed" else "model"
    d <- as.data.frame(pm_get_data(model))
    # (gls names its formula 'model', which pm_refit() takes for the fit: change the call here)
    refit_f <- function(f, data = NULL) {
      cl <- stats::getCall(model)
      cl[[fa]] <- f
      if (!is.null(data)) cl$data <- data
      quiet(eval(cl, pm_fit_env(model)))
    }
    null_fit <- refit_f(null_f)
    P <- pm_nlme_reml(null_fit)
    V <- P$Sig(P$th0)
    R <- Matrix::chol(V)                                  # R'R = V
    mu <- as.numeric(P$X %*% pm_fixef(null_fit))
    n <- length(mu)
    f_full <- new_lhs(stats::formula(model))
    return(list(null_fit = null_fit,
                simulate = function(nsim) {
                  Z <- matrix(stats::rnorm(n * nsim), n)
                  Y <- mu + as.matrix(Matrix::crossprod(R, Z))
                  lapply(seq_len(nsim), function(i) Y[, i])
                },
                refit = function(y) {
                  dd <- d
                  dd$.pm_y <- y
                  refit_f(f_full, dd)
                },
                packages = "nlme"))
  }
  if (inherits(model, c("glmmTMB", "lm"))) {
    if (inherits(model, "glm")) {
      fam <- stats::family(model)$family
      if (grepl("^quasi", fam)) {
        stop("bootmodels() simulates from the model: a quasi family has no distribution to simulate from. ",
             "Use permmodels() instead.", call. = FALSE)
      }
    }
    d <- as.data.frame(if (inherits(model, "glmmTMB")) pm_model_rows(model) else pm_get_data(model))
    null_fit <- quiet(pm_refit(model, formula = null_f))
    f_full <- new_lhs(stats::formula(model))
    return(list(null_fit = null_fit, simulate = function(nsim) as.list(stats::simulate(null_fit, nsim = nsim)),
                refit = function(y) {
                  dd <- d
                  dd$.pm_y <- y
                  quiet(pm_refit(model, data = dd, formula = f_full))
                },
                packages = if (inherits(model, "glmmTMB")) "glmmTMB" else "stats"))
  }
  stop("bootmodels() handles lm, aov, glm, gls, lme, lmer, glmer and glmmTMB models.", call. = FALSE)
}

## The ANOVA table whose statistic bootmodels() resamples, and that column.
##   Type I: sequential tests (F for lm, lme, gls, merMod; for glm the
##   deviance with a known dispersion, F otherwise); glmmTMB has none (Type II).
##   Types II and III: pm_anova() (Wald chi-square; lm F; glm likelihood ratio).
pm_boot_anova <- function(m, type) {
  if (type %in% c("I", "1")) {
    if (inherits(m, "glm")) {
      fam <- stats::family(m)$family
      known <- fam %in% c("poisson", "binomial") || grepl("^Negative Binomial", fam)
      # (MASS warns that a negative binomial's theta is not re-estimated for
      # the sequential tests: the same holds for the observed and the
      # simulated statistics, which is what the bootstrap compares)
      aT <- suppressWarnings(stats::anova(m, test = if (known) "Chisq" else "F"))
      return(aT[rownames(aT) != "NULL", , drop = FALSE])
    }
    return(pm_anova_seq(m))
  }
  pm_anova(m, type = type, test_statistic = "Chisq")
}

pm_boot_stat_col <- function(aT) {
  cn <- colnames(aT)
  intersect(c("F value", "F-value", "F", "LR Chisq", "Chisq", "Deviance"), cn)[1]
}

## The coefficient table of a model (estimate, SE, statistic, p-value).
pm_coef_table <- function(model) {
  if (inherits(model, "glmmTMB")) return(summary(model)$coefficients$cond)
  if (inherits(model, c("lme", "gls"))) return(summary(model)$tTable)
  if (inherits(model, "lm") && !inherits(model, "glm")) return(stats::coef(stats::summary.lm(model)))
  stats::coef(summary(model))
}
