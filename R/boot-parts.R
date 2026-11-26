## What bootmodels() needs for each model class: the model fitted without its
## fixed effects (only an intercept, any offset, and the same random effects
## and covariance structure: the global null hypothesis), responses simulated
## from it, and the full model refitted to a simulated response.
##   lmer, glmer    simulate() and lme4::refit()
##   lm, aov, glm   simulate() (gaussian, binomial, Poisson, Gamma, negative
##                  binomial), refitted with the response replaced
##   glmmTMB        simulate(), refitted with the response replaced
##   clmm           the cut points and random intercepts of the null fit:
##                  u drawn per level, the class from F(theta_j - u) (latent
##                  scale; pm_boot_fit.clmm()), refitted with the response replaced
##   lme, gls       y = X b0 + L z, L L' the null fit's covariance of the data
##                  (random effects, correlation and variance structures;
##                  pm_nlme_reml(), checked against its vcov())
## Returns list(null_fit, simulate = function(nsim) -> list of responses,
## refit = function(y) -> fitted model, packages). null = FALSE: simulate
## from the model itself (residvariogram()); null_fit is then the model.
pm_boot_parts <- function(model, null = TRUE) {
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  ff <- stats::formula(model)
  lhs <- paste(deparse(ff[[2]]), collapse = "")
  # random terms as written (a glmmTMB ar1() term stays ar1(), not unstructured)
  bars <- if (inherits(model, c("merMod", "glmmTMB", "clmm"))) pm_random_terms(ff) else character(0)
  # offsets in the formula stay in the null model (an offset argument stays in the call)
  tt <- stats::terms(stats::delete.response(stats::terms(pm_nobars(ff))))
  offs <- if (is.null(attr(tt, "offset"))) character(0) else
    vapply(attr(tt, "variables")[-1][attr(tt, "offset")], function(e) paste(deparse(e), collapse = ""), "")
  null_f <- if (null) stats::as.formula(paste(lhs, "~", paste(c("1", offs, bars), collapse = " + ")), env = environment(ff))
  new_lhs <- function(f) {
    f[[2]] <- quote(.pm_y)
    f
  }

  pm_boot_fit(model, null_f, new_lhs, quiet)
}

## The class-specific part of pm_boot_parts(): the null fit, the simulation
## of responses from it, and the refit of the full model to one response.
pm_boot_fit <- function(model, null_f, new_lhs, quiet) UseMethod("pm_boot_fit")

#' @exportS3Method NULL
pm_boot_fit.default <- function(model, null_f, new_lhs, quiet) {
  stop("simulation from the model handles lm, aov, glm, gls, lme, lmer, glmer and glmmTMB fits.", call. = FALSE)
}

#' @exportS3Method NULL
pm_boot_fit.merMod <- function(model, null_f, new_lhs, quiet) {
  null_fit <- if (is.null(null_f)) model else quiet(pm_refit(model, formula = null_f))
  return(list(null_fit = null_fit, simulate = function(nsim) stats::simulate(null_fit, nsim = nsim),
              refit = function(y) quiet(pm_mer_refit(model, y)), packages = "lme4"))
}

#' @exportS3Method NULL
pm_boot_fit.lme <- function(model, null_f, new_lhs, quiet) {
  fa <- if (inherits(model, "lme")) "fixed" else "model"
  d <- as.data.frame(pm_get_data(model))
  # (gls names its formula 'model', which pm_refit() takes for the fit: change the call here)
  refit_f <- function(f, data = NULL) {
    cl <- stats::getCall(model)
    cl[[fa]] <- f
    if (!is.null(data)) cl$data <- data
    quiet(eval(cl, pm_fit_env(model)))
  }
  null_fit <- if (is.null(null_f)) model else refit_f(null_f)
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

#' @exportS3Method NULL
pm_boot_fit.gls <- pm_boot_fit.lme

#' @exportS3Method NULL
pm_boot_fit.lm <- function(model, null_f, new_lhs, quiet) {
  if (inherits(model, "glm")) {
    fam <- stats::family(model)$family
    if (grepl("^quasi", fam)) {
      stop("a quasi family has no distribution to simulate from (for tests use permmodels()).", call. = FALSE)
    }
  }
  d <- as.data.frame(if (inherits(model, "glmmTMB")) pm_model_rows(model) else pm_get_data(model))
  null_fit <- if (is.null(null_f)) model else quiet(pm_refit(model, formula = null_f))
  f_full <- new_lhs(stats::formula(model))
  return(list(null_fit = null_fit, simulate = function(nsim) as.list(stats::simulate(null_fit, nsim = nsim)),
              refit = function(y) {
                dd <- d
                dd$.pm_y <- y
                quiet(pm_refit(model, data = dd, formula = f_full))
              },
              packages = if (inherits(model, "glmmTMB")) "glmmTMB" else "stats"))
}

#' @exportS3Method NULL
pm_boot_fit.glmmTMB <- pm_boot_fit.lm

## clmm (B4 G3): random intercepts only, one row per observation (case weights
## would need a multinomial draw per row), no offset; the null fit's cut points
## theta and SDs; refits by clmm() with the response replaced.
#' @exportS3Method NULL
pm_boot_fit.clmm <- function(model, null_f, new_lhs, quiet) {
  if (!pm_ord_intercepts_only(model)) {
    stop("bootmodels() simulates clmm fits with random intercepts; this one has random slopes.", call. = FALSE)
  }
  w <- pm_freq_w(model)
  if (!is.null(w) && any(w != 1)) {
    stop("bootmodels() needs a clmm fit with one row per observation (no case weights).", call. = FALSE)
  }
  if (!is.null(model$offset) && any(model$offset != 0)) stop("bootmodels() does not handle offsets in clmm fits.",
                                                             call. = FALSE)
  d <- pm_clmm_rows(model)
  null_fit <- if (is.null(null_f)) model else quiet(pm_refit(model, formula = null_f, data = d))
  gf <- lapply(null_fit$gfList, factor)
  sd <- vapply(null_fit$ST, function(m) as.matrix(m)[1, 1], numeric(1))
  Fd <- pm_ord_link_fd(null_fit$link)
  th <- as.vector(null_fit$Theta)
  lev <- null_fit$y.levels
  n <- length(gf[[1]])
  f_full <- new_lhs(stats::formula(model))
  list(null_fit = null_fit,
       simulate = function(nsim) lapply(seq_len(nsim), function(i) {
         eta <- Reduce(`+`, Map(function(g, s) stats::rnorm(nlevels(g), 0, s)[g], gf, sd))
         cum <- Fd[[1]](outer(-eta, th, "+"))                       # P(Y <= j), rows x (J - 1)
         factor(lev[1 + rowSums(stats::runif(n) > cum)], levels = lev, ordered = TRUE)
       }),
       refit = function(y) {
         dd <- d
         dd$.pm_y <- y
         f <- quiet(pm_refit(model, data = dd, formula = f_full))
         if (!pm_clmm_ok(f)) stop("the refit did not converge", call. = FALSE)       # (a failure, left out)
         f
       },
       packages = "ordinal")
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
  if (inherits(model, "clmm")) {                 # the latent coefficients (pm_model_parm(), R/ordinal.R)
    mp <- pm_model_parm(model, df = 0)
    se <- sqrt(diag(mp$vcov))
    z <- mp$coef / se
    return(cbind(Estimate = mp$coef, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))))
  }
  if (inherits(model, c("lme", "gls"))) return(summary(model)$tTable)
  if (inherits(model, "lm") && !inherits(model, "glm")) return(stats::coef(stats::summary.lm(model)))
  stats::coef(summary(model))
}
