#' Heritability on a Line-Mean (Entry-Mean) Basis
#'
#' The share of the variance of genotype means that is genetic, from a mixed model with random genotype
#' effects: the standard form for balanced trials and the mixed-model forms of Cullis et al. (2006), Piepho and
#' Möhring (2007) and Oakey et al. (2006) for unbalanced trials, incomplete blocks or correlated errors.
#'
#' @details
#' With \eqn{\sigma^2_g} the genotypic variance:
#' \describe{
#'   \item{\code{"standard"}}{\eqn{\sigma^2_g / (\sigma^2_g + \sum_t \sigma^2_t / n_t + \sigma^2_e / n)}, over
#'   the random terms \eqn{t} that contain the genotype (genotype x environment, ...), with \eqn{n_t} the number of
#'   their levels per genotype and \eqn{n} the number of plots per genotype (harmonic means over genotypes; give
#'   \code{n_env}, \code{n_rep} to set them). With one environment and \eqn{r} replicates,
#'   \eqn{\sigma^2_g/(\sigma^2_g + \sigma^2_e/r)}; with \eqn{e} environments,
#'   \eqn{\sigma^2_g/(\sigma^2_g + \sigma^2_{ge}/e + \sigma^2_e/(re))}. Exact for balanced data only.}
#'   \item{\code{"cullis"}}{\eqn{1 - \bar v_{BLUP}/(2\sigma^2_g)}, \eqn{\bar v_{BLUP}} the mean variance of a
#'   difference of two genotype BLUPs (prediction error variance, from the inverse of the mixed model
#'   equations, as \code{se_ranef(type = "prediction")}).}
#'   \item{\code{"piepho"}}{\eqn{\sigma^2_g / (\sigma^2_g + \bar v_{BLUE}/2)}, \eqn{\bar v_{BLUE}} the mean
#'   variance of a difference of two genotype means from the same model refitted with genotype as a fixed
#'   effect (other random terms kept), from their unadjusted covariance as in Piepho and Möhring (2007).
#'   On the oats alpha design of John and Williams (1995), Example 1 of Piepho and Möhring (2007): 0.803 as
#'   published (Cullis 0.809). (\code{inti::H2cal(emmeans = TRUE)} and \pkg{heritable} take emmeans's
#'   Kenward-Roger adjusted covariance: 0.797 there.)}
#'   \item{\code{"oakey"}}{the mean of the non-zero eigenvalues of \eqn{I - C_{gg}/\sigma^2_g}, \eqn{C_{gg}}
#'   the prediction error covariance of the genotype BLUPs (the zero ones are directions confounded with fixed
#'   effects, as the mean of all genotypes with the intercept). With independent genotypes it equals
#'   \code{"cullis"}; it differs with a relationship matrix.}
#'   \item{\code{"family"}}{narrow-sense heritability of a single tree (plant, animal) from families:
#'   \eqn{(\sigma^2_f/\rho) / (\sigma^2_f + \sum_t \sigma^2_t + \sigma^2_e)}, \eqn{\sigma^2_f} the family variance
#'   (\code{genotype} is the family factor), \eqn{\rho} the \code{relationship} among members of a family (1/4
#'   for half-sibs, 1/2 for full-sibs; often taken as 1/3 for open-pollinated progeny, which are partly full-sibs
#'   or selfs), and \eqn{\sigma^2_t} the random terms with the family in them (family x block, i.e. plot; family x
#'   site). The other design terms (blocks, rows, columns) are left out of the phenotypic variance. Can exceed
#'   1 when \eqn{\rho} is too small for the material. Family-mean heritability is the line-mean forms above
#'   with the family factor as \code{genotype}.}
#'   \item{\code{"repeatability"}}{\eqn{\sigma^2_s / (\sigma^2_s + \sum_t \sigma^2_t + \sigma^2_e)} with
#'   \code{genotype} the subject that has repeated records: the share of the variance of a single record due to
#'   the subject.}
#' }
#' The four line-mean forms agree for balanced data with independent genotypes; with imbalance, incomplete blocks or
#' spatial errors they differ, and Schmidt et al. (2019) recommend the mixed-model forms. Narrow-sense
#' heritability from a pedigree or genomic relationship matrix needs a model with that matrix, which lme4
#' cannot fit: \pkg{heritable} computes these forms (and the per-genotype \eqn{H^2_\Delta} of Schmidt et
#' al. 2019) from asreml fits with a relationship matrix. Its standard, Cullis and Oakey values for lmer fits
#' are those of this function (tested); this function adds standard errors, the family and repeatability
#' forms, and advice.
#' @param model An \code{lmer} fit with a random genotype term \code{(1 | genotype)}.
#' @param genotype Name of the genotype (line, entry, family) factor.
#' @param method \code{"cullis"} (default), \code{"piepho"}, \code{"standard"}, \code{"oakey"}, or \code{"all"}.
#' @param n_env,n_rep For \code{"standard"}: the number of environments (levels of the genotype x environment
#'   terms) and of replicates per environment; by default from the data.
#' @param relationship For \code{"family"}: the coefficient of relationship among members of a family.
#' @param nboot For \code{"cullis"}, \code{"piepho"} and \code{"oakey"}: parametric bootstrap draws for their
#'   standard error and percentile limits (each draw refits the model; 0, the default, gives none). The other
#'   forms get theirs by the delta method from the covariance of the variance parameters, as \code{varcomp()}.
#' @param ncore Cores for the bootstrap; results are the same for any number.
#' @param seed Seed for the bootstrap draws.
#' @param conf_level Confidence level of the limits (logit-normal for the delta method, so they stay in (0, 1)).
#' @return A data frame of class \code{"pm_heritability"}: \code{method}, \code{estimate}, \code{std.error},
#'   \code{conf.low}, \code{conf.high}; \code{attr(, "components")} holds the variance components and the mean
#'   variances of differences used.
#' @references
#' Cullis, B. R., Smith, A. B. and Coombes, N. E. (2006). On the design of early generation variety trials
#' with correlated data. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 11, 381-393.
#'
#' Oakey, H., Verbyla, A., Pitchford, W., Cullis, B. and Kuchel, H. (2006). Joint modeling of additive and
#' non-additive genetic line effects in single field trials. \emph{Theoretical and Applied Genetics}, 113,
#' 809-819.
#'
#' Piepho, H.-P. and Möhring, J. (2007). Computing heritability and selection response from unbalanced plant
#' breeding trials. \emph{Genetics}, 177, 1881-1888.
#'
#' Schmidt, P., Hartung, J., Bennewitz, J. and Piepho, H.-P. (2019). Heritability in plant breeding on a
#' genotype-difference basis. \emph{Genetics}, 212, 991-1008.
#' @seealso \code{\link{varcomp}}, \code{\link{se_ranef}}
#' @examples
#' library(predictmeans)
#' Oats$nitro <- factor(Oats$nitro)
#' fm <- lmer(yield ~ nitro + (1 | Block) + (1 | Variety), data = Oats)
#' heritability(fm, "Variety", method = "all")
#' \donttest{
#' heritability(fm, "Variety", nboot = 199, seed = 1)   # SE of the Cullis form by bootstrap
#' }
#' @export
heritability <- function(model, genotype,
                         method = c("cullis", "piepho", "standard", "oakey", "family", "repeatability", "all"),
                         relationship = NULL, n_env = NULL, n_rep = NULL, nboot = 0, ncore = 1, seed = NULL,
                         conf_level = 0.95) {
  method <- match.arg(method)
  methods <- if (method == "all") c("standard", "cullis", "piepho", "oakey") else method
  if (method == "family" && !(is.numeric(relationship) && length(relationship) == 1 && relationship > 0 &&
                              relationship <= 1)) {
    stop("method = \"family\" needs 'relationship', the coefficient of relationship among members of a family ",
         "(1/4 half-sibs, 1/2 full-sibs; often 1/3 for open-pollinated progeny).", call. = FALSE)
  }
  alpha <- pm_conf_alpha(conf_level)
  model <- pm_fit(model)
  if (!inherits(model, "lmerMod")) {
    stop("heritability() handles lmer fits (gaussian, random genotype effects); a ", class(model)[1],
         " model is not yet covered.", call. = FALSE)
  }
  comp <- pm_h2_components(model, genotype, n_env, n_rep)
  est <- pm_h2_point(model, comp, methods, relationship)
  # SE and limits: the delta method on the variance parameters for the forms
  # that are functions of the components; a parametric bootstrap for the
  # forms that also depend on the prediction errors (nboot > 0)
  vfun <- pm_h2_vfun(comp, relationship)
  boot_m <- intersect(methods, c("cullis", "piepho", "oakey"))
  bt <- if (length(boot_m) && nboot > 0) pm_h2_boot(model, comp, boot_m, nboot, ncore, seed)
  rows <- lapply(methods, function(m) {
    if (m %in% names(vfun)) {
      se <- pm_h2_delta(model, vfun[[m]])
      lim <- pm_h2_limits(est[[m]], se, alpha)
    } else if (!is.null(bt)) {
      draws <- bt[, m]
      se <- stats::sd(draws, na.rm = TRUE)
      lim <- unname(stats::quantile(draws, c(alpha/2, 1 - alpha/2), na.rm = TRUE))
    } else {
      se <- NA_real_
      lim <- c(NA_real_, NA_real_)
    }
    data.frame(method = m, estimate = est[[m]], std.error = se, conf.low = lim[1], conf.high = lim[2])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- structure(out, class = c("pm_heritability", "data.frame"), conf_level = conf_level,
                   components = comp[c("s2g", "terms", "s2e", "n_terms", "n_plots")], relationship = relationship,
                   genotype = genotype,
                   nboot = if (!is.null(bt)) sum(stats::complete.cases(bt)) else 0)
  pm_add_advice(out, tryCatch(pm_heritability_advice(out, comp, nboot), error = function(e) NULL))
}

## The estimates of the requested forms from one fit (also each bootstrap refit).
pm_h2_point <- function(model, comp, methods, relationship = NULL, y = NULL) {
  vf <- pm_h2_vfun(comp, relationship)
  v <- pm_h2_vc(comp)
  vapply(methods, function(m) switch(m,
    cullis = 1 - pm_h2_vdiff(comp$C)/(2*comp$s2g),
    piepho = comp$s2g / (comp$s2g + pm_h2_vblue(model, comp$genotype, y)/2),
    oakey = pm_h2_oakey(comp$C, comp$s2g),
    vf[[m]](v)), numeric(1))
}

## The variance components of the fit, named as pm_vc_mer_fun() names them.
pm_h2_vc <- function(comp) {
  stats::setNames(c(comp$s2g, comp$terms, comp$s2e),
                  c(paste0(comp$genotype, ".(Intercept)"), pm_h2_tnames(comp$with_g), "residual"))
}

## The forms that are functions of the variance components alone (v, named as
## pm_vc_mer_fun()): standard (line mean), family (single tree, narrow sense:
## the family variance is relationship x additive variance) and repeatability
## (single record). The terms with the genotype in them (genotype x
## environment, family x block = plot) are in the variance of a single record;
## the other design terms (blocks, rows, columns) are not.
pm_h2_vfun <- function(comp, relationship = NULL) {
  g <- paste0(comp$genotype, ".(Intercept)")
  tt <- pm_h2_tnames(comp$with_g)
  single <- function(v) v[[g]] + sum(v[tt]) + v[["residual"]]
  list(standard = function(v) v[[g]] / (v[[g]] + sum(v[tt] / comp$n_terms) + v[["residual"]] / comp$n_plots),
       family = function(v) (v[[g]] / relationship) / single(v),
       repeatability = function(v) v[[g]] / single(v))
}

## The pieces of every form, from one fit: the genotypic variance, the other
## random terms that contain the genotype with their numbers of levels per
## genotype, the residual, the plots per genotype, and the prediction error
## covariance of the genotype BLUPs (pm_ranef_cov(), R/ranef-cov.R).
pm_h2_components <- function(model, genotype, n_env = NULL, n_rep = NULL) {
  vc <- lme4::VarCorr(model)
  if (!genotype %in% names(vc) || !identical(colnames(vc[[genotype]]), "(Intercept)")) {
    stop("heritability() needs a random intercept for '", genotype, "', (1 | ", genotype, "), in the model.",
         call. = FALSE)
  }
  s2g <- vc[[genotype]][1, 1]
  if (s2g <= 0) stop("The genotypic variance is estimated as 0: the heritability is 0.", call. = FALSE)
  if (stats::family(model)$family != "gaussian") stop("heritability() needs a gaussian lmer fit.", call. = FALSE)
  grp <- strsplit(names(vc), ":", fixed = TRUE)
  with_g <- names(vc)[vapply(grp, function(g) genotype %in% g && length(g) > 1, logical(1))]
  dat <- stats::model.frame(model)
  per_g <- function(vars) {             # levels of 'vars' per genotype (harmonic mean)
    k <- tapply(do.call(paste, c(dat[vars], sep = "\r")), dat[[genotype]], function(x) length(unique(x)))
    k <- k[!is.na(k)]
    length(k) / sum(1/k)
  }
  n_terms <- vapply(with_g, function(t) {
    other <- setdiff(strsplit(t, ":", fixed = TRUE)[[1]], genotype)
    if (!is.null(n_env)) n_env else per_g(other)
  }, numeric(1))
  kp <- table(dat[[genotype]]); kp <- kp[kp > 0]
  n_plots <- if (!is.null(n_env) && !is.null(n_rep)) n_env * n_rep else length(kp) / sum(1/kp)
  s2t <- vapply(with_g, function(t) vc[[t]][1, 1], numeric(1))
  s2e <- stats::sigma(model)^2
  rc <- pm_ranef_cov(model)
  i <- which(rc$tab$group == genotype & rc$tab$term == "(Intercept)")
  p <- length(rc$beta)
  list(s2g = s2g, terms = s2t, s2e = s2e, n_terms = n_terms, n_plots = n_plots, n_env_arg = n_env, n_rep_arg = n_rep,
       unequal = length(unique(as.vector(kp))) > 1, incomplete = pm_h2_incomplete(dat, names(vc), genotype),
       C = rc$V[p + i, p + i, drop = FALSE],
       with_g = with_g, genotype = genotype)
}

## Incomplete blocks: a random term without the genotype whose levels hold
## fewer genotypes than the trial (an alpha or lattice design). The standard
## form then ignores the block adjustment of the genotype means.
pm_h2_incomplete <- function(dat, groups, genotype) {
  ng <- length(unique(dat[[genotype]]))
  any(vapply(groups, function(g) {
    v <- strsplit(g, ":", fixed = TRUE)[[1]]
    if (genotype %in% v || !all(v %in% names(dat))) return(FALSE)
    k <- tapply(dat[[genotype]], do.call(paste, c(dat[v], sep = "\r")), function(x) length(unique(x)))
    any(k < ng)
  }, logical(1)))
}

## Mean variance of a difference between the m entries with covariance C:
## mean over pairs of C_ii + C_jj - 2 C_ij = 2 (m tr(C) - sum(C)) / (m (m - 1)).
pm_h2_vdiff <- function(C) {
  m <- nrow(C)
  2 * (m * sum(diag(C)) - sum(C)) / (m * (m - 1))
}

## vBLUE: the same model with genotype fixed (its other random terms kept), the
## genotype means (Kmatrix()) and their unadjusted covariance; lm() when no
## random term is left.
pm_h2_vblue <- function(model, genotype, y = NULL) {
  f <- stats::formula(model)
  bars <- pm_findbars(f)
  keep <- bars[vapply(bars, function(b) !identical(deparse(b), paste("1 |", genotype)), logical(1))]
  fixed_terms <- attr(stats::terms(pm_nobars(f)), "term.labels")
  rhs <- c(fixed_terms, genotype, vapply(keep, function(b) paste0("(", deparse(b), ")"), ""))
  # the response as a column (also a transformed one, log(y), or a bootstrap draw)
  newf <- stats::as.formula(paste(".pm_y ~", paste(rhs, collapse = " + ")), env = environment(f))
  dat <- stats::model.frame(model)
  dat$.pm_y <- if (is.null(y)) stats::model.response(dat) else y
  dat[[genotype]] <- factor(dat[[genotype]])
  fit <- if (length(keep)) {
    suppressMessages(lme4::lmer(newf, data = dat, REML = lme4::getME(model, "is_REML")))
  } else {
    stats::lm(newf, data = dat)
  }
  km <- suppressMessages(suppressWarnings(Kmatrix(fit, genotype)))
  mp <- pm_model_parm(fit)
  K <- pm_K_cols(km$K, mp)
  pm_h2_vdiff(K %*% mp$vcov %*% t(K))
}

## Oakey et al. (2006): the mean of the non-zero eigenvalues of I - G^-1 C
## (G = s2g I); the zero ones are directions confounded with fixed effects
## (the mean of the genotypes, with the intercept).
pm_h2_oakey <- function(C, s2g) {
  ev <- eigen(diag(nrow(C)) - as.matrix(C)/s2g, symmetric = TRUE, only.values = TRUE)$values
  mean(ev[ev > 1e-8])
}

## Names of random intercept terms as pm_vc_mer_fun() gives them (none for none:
## paste0() would return ".(Intercept)" for an empty vector).
pm_h2_tnames <- function(t) if (length(t)) paste0(t, ".(Intercept)") else character(0)

## SE of a form that is a function of the variance components, by the delta
## method on the variance parameters (theta, sigma) with their covariance from
## the observed information (pm_lmer_varpar(), as varcomp()).
pm_h2_delta <- function(model, f) {
  fun <- pm_vc_mer_fun(model, residual = TRUE)
  vp <- pm_lmer_varpar(model)
  J <- numDeriv::jacobian(function(x) f(fun(x)), vp$varpar)
  sqrt(drop(J %*% vp$vcov_varpar %*% t(J)))
}

## Parametric bootstrap of the forms that depend on prediction errors: draws
## from the fitted model, each refitted (pm_mer_refit()) and its forms
## recomputed. Seeds per draw are made here, so the result is the same for any
## ncore (pm_perm_map(), R/perm-engine.R). Returns a matrix nboot x methods.
pm_h2_boot <- function(model, comp, methods, nboot, ncore = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, nboot)
  ys <- lapply(seeds, function(s) { set.seed(s); stats::simulate(model, nsim = 1)[[1]] })
  one <- function(b) {
    fit <- suppressMessages(suppressWarnings(pm_mer_refit(model, ys[[b]])))
    cb <- tryCatch(pm_h2_components(fit, comp$genotype, comp$n_env_arg, comp$n_rep_arg), error = function(e) NULL)
    if (is.null(cb)) return(stats::setNames(rep(NA_real_, length(methods)), methods))
    pm_h2_point(fit, cb, methods, y = ys[[b]])
  }
  out <- pm_perm_map(seq_len(nboot), one, ncore, packages = "predictmeans")
  out <- lapply(out, function(r) if (is.null(r)) stats::setNames(rep(NA_real_, length(methods)), methods) else r)
  do.call(rbind, out)
}

## Limits on the logit scale, so they stay in (0, 1).
pm_h2_limits <- function(est, se, alpha) {
  if (is.na(se) || est <= 0 || est >= 1) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - alpha/2)
  l <- stats::qlogis(est) + c(-1, 1) * z * se / (est * (1 - est))
  stats::plogis(l)
}

#' @exportS3Method print pm_heritability
print.pm_heritability <- function(x, digits = 4, ...) {
  cmp <- attr(x, "components")
  m <- x$method
  kind <- c(if (any(m %in% c("standard", "cullis", "piepho", "oakey"))) "heritability on a line-mean basis",
            if ("family" %in% m) paste0("narrow-sense heritability of a single plant (relationship ",
                                        format(attr(x, "relationship"), digits = 3), ")"),
            if ("repeatability" %in% m) "repeatability of a single record")
  cat(toupper(substring(kind[1], 1, 1)), substring(paste(kind, collapse = "; "), 2), " (",
      format(100 * attr(x, "conf_level")), "% limits)\n", sep = "")
  cat("Variance of ", attr(x, "genotype"), " ", format(cmp$s2g, digits = digits), "; residual ",
      format(cmp$s2e, digits = digits),
      if (length(cmp$terms)) paste0("; ", paste(names(cmp$terms), format(cmp$terms, digits = digits), sep = " ",
                                                 collapse = ", ")), "\n\n", sep = "")
  y <- x
  class(y) <- "data.frame"
  attr(y, "components") <- attr(y, "conf_level") <- NULL
  print(y, digits = digits, row.names = FALSE)
  nb <- attr(x, "nboot")
  if (isTRUE(nb > 0)) cat("\nSE and limits of the mixed-model forms from", nb, "parametric bootstrap draws (percentile).\n")
  pm_print_advice(x)
  invisible(x)
}
