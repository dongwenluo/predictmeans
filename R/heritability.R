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
#'   (\code{inti::H2cal(emmeans = TRUE)} takes emmeans's Kenward-Roger adjusted covariance, a slightly smaller
#'   value with unbalanced data.)}
#'   \item{\code{"oakey"}}{the mean of the non-zero eigenvalues of \eqn{I - C_{gg}/\sigma^2_g}, \eqn{C_{gg}}
#'   the prediction error covariance of the genotype BLUPs (the zero ones are directions confounded with fixed
#'   effects, as the mean of all genotypes with the intercept). With independent genotypes it equals
#'   \code{"cullis"}; it differs with a relationship matrix.}
#' }
#' The four agree for balanced data with independent genotypes; with imbalance, incomplete blocks or
#' spatial errors they differ, and Schmidt et al. (2019) recommend the mixed-model forms. Narrow-sense
#' heritability from a pedigree or genomic relationship matrix needs a model with that matrix, which lme4
#' cannot fit.
#' @param model An \code{lmer} fit with a random genotype term \code{(1 | genotype)}.
#' @param genotype Name of the genotype (line, entry, family) factor.
#' @param method \code{"cullis"} (default), \code{"piepho"}, \code{"standard"}, \code{"oakey"}, or \code{"all"}.
#' @param n_env,n_rep For \code{"standard"}: the number of environments (levels of the genotype x environment
#'   terms) and of replicates per environment; by default from the data.
#' @param conf_level Confidence level of the limits.
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
#' @export
heritability <- function(model, genotype, method = c("cullis", "piepho", "standard", "oakey", "all"),
                         n_env = NULL, n_rep = NULL, conf_level = 0.95) {
  method <- match.arg(method)
  methods <- if (method == "all") c("standard", "cullis", "piepho", "oakey") else method
  alpha <- pm_conf_alpha(conf_level)
  model <- pm_fit(model)
  if (!inherits(model, "lmerMod")) {
    stop("heritability() handles lmer fits (gaussian, random genotype effects); a ", class(model)[1],
         " model is not yet covered.", call. = FALSE)
  }
  comp <- pm_h2_components(model, genotype, n_env, n_rep)
  rows <- lapply(methods, function(m) {
    est <- switch(m,
      standard = comp$s2g / (comp$s2g + comp$denom_standard),
      cullis = 1 - pm_h2_vdiff(comp$C)/(2*comp$s2g),
      piepho = comp$s2g / (comp$s2g + pm_h2_vblue(model, genotype)/2),
      oakey = pm_h2_oakey(comp$C, comp$s2g))
    se <- if (m == "standard") pm_h2_se_standard(model, comp) else NA_real_
    lim <- pm_h2_limits(est, se, alpha)
    data.frame(method = m, estimate = est, std.error = se, conf.low = lim[1], conf.high = lim[2])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out, class = c("pm_heritability", "data.frame"), conf_level = conf_level,
            components = comp[c("s2g", "terms", "s2e", "n_terms", "n_plots")])
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
  n_plots <- if (!is.null(n_env) && !is.null(n_rep)) n_env * n_rep else {
    k <- table(dat[[genotype]]); k <- k[k > 0]; length(k) / sum(1/k)
  }
  s2t <- vapply(with_g, function(t) vc[[t]][1, 1], numeric(1))
  s2e <- stats::sigma(model)^2
  rc <- pm_ranef_cov(model)
  i <- which(rc$tab$group == genotype & rc$tab$term == "(Intercept)")
  p <- length(rc$beta)
  list(s2g = s2g, terms = s2t, s2e = s2e, n_terms = n_terms, n_plots = n_plots,
       denom_standard = sum(s2t / n_terms) + s2e / n_plots, C = rc$V[p + i, p + i, drop = FALSE],
       with_g = with_g, genotype = genotype)
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
pm_h2_vblue <- function(model, genotype) {
  f <- stats::formula(model)
  bars <- pm_findbars(f)
  keep <- bars[vapply(bars, function(b) !identical(deparse(b), paste("1 |", genotype)), logical(1))]
  fixed_terms <- attr(stats::terms(pm_nobars(f)), "term.labels")
  rhs <- c(fixed_terms, genotype, vapply(keep, function(b) paste0("(", deparse(b), ")"), ""))
  newf <- stats::as.formula(paste(deparse(f[[2]]), "~", paste(rhs, collapse = " + ")), env = environment(f))
  dat <- stats::model.frame(model)
  dat[[genotype]] <- factor(dat[[genotype]])
  fit <- if (length(keep)) {
    lme4::lmer(newf, data = dat, REML = lme4::getME(model, "is_REML"))
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

## SE of the standard form by the delta method on the variance parameters
## (theta, sigma) and their covariance from the observed information
## (pm_lmer_varpar(), as varcomp()).
pm_h2_se_standard <- function(model, comp) {
  fun <- pm_vc_mer_fun(model, residual = TRUE)
  vp <- pm_lmer_varpar(model)
  g <- paste0(comp$genotype, ".(Intercept)")
  tt <- paste0(comp$with_g, ".(Intercept)")
  h2 <- function(x) {
    v <- fun(x)
    v[g] / (v[g] + sum(v[tt] / comp$n_terms) + v["residual"] / comp$n_plots)
  }
  J <- numDeriv::jacobian(h2, vp$varpar)
  sqrt(drop(J %*% vp$vcov_varpar %*% t(J)))
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
  cat("Heritability on a line-mean basis (", format(100 * attr(x, "conf_level")), "% limits)\n", sep = "")
  cat("Genotypic variance ", format(cmp$s2g, digits = digits), "; residual ", format(cmp$s2e, digits = digits),
      if (length(cmp$terms)) paste0("; ", paste(names(cmp$terms), format(cmp$terms, digits = digits), sep = " ",
                                                 collapse = ", ")), "\n\n", sep = "")
  y <- x
  class(y) <- "data.frame"
  attr(y, "components") <- attr(y, "conf_level") <- NULL
  print(y, digits = digits, row.names = FALSE)
  invisible(x)
}
